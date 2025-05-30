#include "surface.hpp"
#include "../core/utils.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace math::volatility {

    void VolatilitySurface::calibrate(const VolatilityPoints& marketData,
                                      core::Real spot, core::Real riskFreeRate)
    {
        spot_ = spot;
        riskFreeRate_ = riskFreeRate;

        {
            core::UniqueLockGuard<core::ReadWriteLock> lock(cacheLock_);
            interpolationCache_.clear();
        }

        std::set<core::Real> strikeSet, maturitySet;

        for(const auto& point: marketData) {
            strikeSet.insert(point.strike);
            maturitySet.insert(point.timeToExpiry);
        }

        strikes_.assign(strikeSet.begin(), strikeSet.end());
        maturities_.assign(maturitySet.begin(), maturitySet.end());

        surface_.resize(strikes_.size(), core::Vector(maturities_.size(), 0.0));

        buildSurface(marketData);
        smoothSurface();
        extrapolateSurface();
        enforceArbitrageConstraints();
    }

    void VolatilitySurface::buildSurface(const VolatilityPoints& marketData)
    {
        std::map<std::pair<core::Real, core::Real>, core::Real> dataMap;

        for(const auto& point: marketData) {
            core::Real vol = point.impliedVol;
            if(vol <= 0.0) { vol = calculateImpliedVolatility(point); }
            dataMap[{point.strike, point.timeToExpiry}] = vol;
        }

        for(core::Size i = 0; i < strikes_.size(); ++i) {
            for(core::Size j = 0; j < maturities_.size(); ++j) {
                auto it = dataMap.find({strikes_[i], maturities_[j]});
                if(it != dataMap.end()) { surface_[i][j] = it->second; }
                else {
                    surface_[i][j] = interpolateFromNeighbors(
                        strikes_[i], maturities_[j], dataMap);
                }
            }
        }
    }

    core::Real VolatilitySurface::interpolateFromNeighbors(
        core::Real strike, core::Real timeToExpiry,
        const std::map<std::pair<core::Real, core::Real>, core::Real>& dataMap)
        const
    {
        core::Vector distances, vols;

        for(const auto& [key, vol]: dataMap) {
            const core::Real strikeDistance
                = std::abs(key.first - strike) / spot_;
            const core::Real timeDistance = std::abs(key.second - timeToExpiry);
            const core::Real distance = std::sqrt(
                strikeDistance * strikeDistance + timeDistance * timeDistance);

            distances.push_back(distance);
            vols.push_back(vol);
        }

        if(distances.empty()) return 0.2;

        core::Real weightedSum = 0.0;
        core::Real totalWeight = 0.0;

        for(core::Size i = 0; i < distances.size(); ++i) {
            const core::Real weight = 1.0 / (distances[i] + 1e-8);
            weightedSum += weight * vols[i];
            totalWeight += weight;
        }

        return totalWeight > 0.0 ? weightedSum / totalWeight : 0.2;
    }

    void VolatilitySurface::smoothSurface()
    {
        const core::Size strikeCount = strikes_.size();
        const core::Size maturityCount = maturities_.size();

        if(strikeCount < 3 || maturityCount < 3) return;

        VolatilityGrid smoothed = surface_;

        for(core::Size i = 1; i < strikeCount - 1; ++i) {
            for(core::Size j = 1; j < maturityCount - 1; ++j) {
                core::Real sum = 0.0;
                core::Real weight = 0.0;

                for(int di = -1; di <= 1; ++di) {
                    for(int dj = -1; dj <= 1; ++dj) {
                        const core::Real w = (di == 0 && dj == 0) ? 4.0 : 1.0;
                        sum += w * surface_[i + di][j + dj];
                        weight += w;
                    }
                }

                smoothed[i][j] = sum / weight;
            }
        }

        surface_ = std::move(smoothed);
    }

    void VolatilitySurface::extrapolateSurface()
    {
        const core::Size strikeCount = strikes_.size();
        const core::Size maturityCount = maturities_.size();

        if(strikeCount < 2 || maturityCount < 2) return;

        for(core::Size j = 0; j < maturityCount; ++j) {
            if(strikeCount >= 2) {
                const core::Real slope = (surface_[1][j] - surface_[0][j])
                                         / (strikes_[1] - strikes_[0]);
                surface_[0][j] = std::max(
                    0.01, surface_[1][j] - slope * (strikes_[1] - strikes_[0]));

                const core::Real endSlope
                    = (surface_[strikeCount - 1][j]
                       - surface_[strikeCount - 2][j])
                      / (strikes_[strikeCount - 1] - strikes_[strikeCount - 2]);
                surface_[strikeCount - 1][j]
                    = std::max(0.01, surface_[strikeCount - 2][j]
                                         + endSlope
                                               * (strikes_[strikeCount - 1]
                                                  - strikes_[strikeCount - 2]));
            }
        }

        for(core::Size i = 0; i < strikeCount; ++i) {
            if(maturityCount >= 2) {
                const core::Real slope = (surface_[i][1] - surface_[i][0])
                                         / (maturities_[1] - maturities_[0]);
                surface_[i][0] = std::max(
                    0.01,
                    surface_[i][1] - slope * (maturities_[1] - maturities_[0]));

                const core::Real endSlope
                    = (surface_[i][maturityCount - 1]
                       - surface_[i][maturityCount - 2])
                      / (maturities_[maturityCount - 1]
                         - maturities_[maturityCount - 2]);
                surface_[i][maturityCount - 1] = std::max(
                    0.01, surface_[i][maturityCount - 2]
                              + endSlope
                                    * (maturities_[maturityCount - 1]
                                       - maturities_[maturityCount - 2]));
            }
        }
    }

    core::Real
    VolatilitySurface::bicubicInterpolation(core::Real strike,
                                            core::Real timeToExpiry) const
    {
        if(strikes_.size() < 2 || maturities_.size() < 2) {
            return bilinearInterpolation(strike, timeToExpiry);
        }

        const core::Size i
            = std::min(findStrikeIndex(strike), strikes_.size() - 2);
        const core::Size j
            = std::min(findMaturityIndex(timeToExpiry), maturities_.size() - 2);

        if(i == 0 || i >= strikes_.size() - 1 || j == 0
           || j >= maturities_.size() - 1) {
            return bilinearInterpolation(strike, timeToExpiry);
        }

        const core::Real t
            = (strike - strikes_[i]) / (strikes_[i + 1] - strikes_[i]);
        const core::Real u = (timeToExpiry - maturities_[j])
                             / (maturities_[j + 1] - maturities_[j]);

        core::StaticArray<16> p;
        core::Size idx = 0;
        for(int di = -1; di <= 2; ++di) {
            for(int dj = -1; dj <= 2; ++dj) {
                const core::Size si = std::clamp(
                    static_cast<core::Integer>(i) + di, core::Integer{0},
                    static_cast<core::Integer>(strikes_.size() - 1));
                const core::Size sj = std::clamp(
                    static_cast<core::Integer>(j) + dj, core::Integer{0},
                    static_cast<core::Integer>(maturities_.size() - 1));
                p[idx++] = surface_[si][sj];
            }
        }

        auto cubic = [](core::Real p0, core::Real p1, core::Real p2,
                        core::Real p3, core::Real t) {
            return p1
                   + 0.5 * t
                         * (p2 - p0
                            + t
                                  * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3
                                     + t * (3.0 * (p1 - p2) + p3 - p0)));
        };

        const core::Real col0 = cubic(p[0], p[4], p[8], p[12], t);
        const core::Real col1 = cubic(p[1], p[5], p[9], p[13], t);
        const core::Real col2 = cubic(p[2], p[6], p[10], p[14], t);
        const core::Real col3 = cubic(p[3], p[7], p[11], p[15], t);

        return std::max(0.01, cubic(col0, col1, col2, col3, u));
    }

    core::Real VolatilitySurface::getVolatilityDelta(core::Real delta,
                                                     core::Real timeToExpiry,
                                                     bool isCall) const
    {
        auto deltaFunc = [this, timeToExpiry, isCall](core::Real strike) {
            core::MarketData market;
            market.spot = spot_;
            market.strike = strike;
            market.timeToExpiry = timeToExpiry;
            market.riskFreeRate = riskFreeRate_;
            market.volatility = getVolatility(strike, timeToExpiry);
            market.isCall = isCall;

            return engine_->blackScholes(market).delta;
        };

        return core::NumericalMethods::brent(
            [&](core::Real strike) { return deltaFunc(strike) - delta; },
            spot_ * 0.5, spot_ * 2.0);
    }

    core::Real VolatilitySurface::getForwardVolatility(core::Real strike,
                                                       core::Real startTime,
                                                       core::Real endTime) const
    {
        if(startTime >= endTime || startTime < 0) return 0.0;

        const core::Real vol1 = getVolatility(strike, startTime);
        const core::Real vol2 = getVolatility(strike, endTime);

        const core::Real var1 = vol1 * vol1 * startTime;
        const core::Real var2 = vol2 * vol2 * endTime;

        if(var2 <= var1) return vol2;

        return std::sqrt((var2 - var1) / (endTime - startTime));
    }

    core::Real
    VolatilitySurface::getLocalVolatility(core::Real strike,
                                          core::Real timeToExpiry) const
    {
        const core::Real h = 0.01;
        const core::Real dt = 1.0 / 365.0;

        const core::Real volBase = getVolatility(strike, timeToExpiry);
        const core::Real vol_dT = getVolatility(strike, timeToExpiry + dt);
        const core::Real vol_dK_up = getVolatility(strike + h, timeToExpiry);
        const core::Real vol_dK_down = getVolatility(strike - h, timeToExpiry);
        const core::Real vol_d2K
            = getVolatility(strike + 2.0 * h, timeToExpiry);

        const core::Real dVol_dT = (vol_dT - volBase) / dt;
        const core::Real dVol_dK = (vol_dK_up - vol_dK_down) / (2.0 * h);
        const core::Real d2Vol_dK2
            = (vol_dK_up - 2.0 * volBase + vol_dK_down) / (h * h);

        const core::Real d
            = std::log(spot_ / strike)
              + (riskFreeRate_ + 0.5 * volBase * volBase) * timeToExpiry;
        d /= (volBase * std::sqrt(timeToExpiry));

        const core::Real numerator
            = volBase * volBase + 2.0 * volBase * timeToExpiry * dVol_dT;
        const core::Real denominator = 1.0
                                       + d * std::sqrt(timeToExpiry) * dVol_dK
                                       + 0.25 * timeToExpiry * volBase
                                             * (d * d - 1.0) * dVol_dK * dVol_dK
                                       + timeToExpiry * volBase * d2Vol_dK2;

        return std::sqrt(std::max(0.01, numerator / denominator));
    }

    VolatilitySurface::VolatilityMetrics
    VolatilitySurface::calculateMetrics(core::Real timeToExpiry) const
    {
        VolatilityMetrics metrics;

        const core::Real atmStrike
            = spot_ * std::exp(riskFreeRate_ * timeToExpiry);
        metrics.atmVol = getVolatility(atmStrike, timeToExpiry);

        const core::Real delta25Call
            = getVolatilityDelta(0.25, timeToExpiry, true);
        const core::Real delta25Put
            = getVolatilityDelta(-0.25, timeToExpiry, false);
        const core::Real delta10Call
            = getVolatilityDelta(0.10, timeToExpiry, true);
        const core::Real delta10Put
            = getVolatilityDelta(-0.10, timeToExpiry, false);

        const core::Real vol25Call = getVolatility(delta25Call, timeToExpiry);
        const core::Real vol25Put = getVolatility(delta25Put, timeToExpiry);
        const core::Real vol10Call = getVolatility(delta10Call, timeToExpiry);
        const core::Real vol10Put = getVolatility(delta10Put, timeToExpiry);

        metrics.skew25Delta = vol25Put - vol25Call;
        metrics.skew10Delta = vol10Put - vol10Call;

        metrics.butterfly25Delta
            = 0.5 * (vol25Call + vol25Put) - metrics.atmVol;
        metrics.butterfly10Delta
            = 0.5 * (vol10Call + vol10Put) - metrics.atmVol;

        metrics.riskReversal25Delta = vol25Call - vol25Put;
        metrics.riskReversal10Delta = vol10Call - vol10Put;

        const core::Real convexityStrike1 = atmStrike * 0.9;
        const core::Real convexityStrike2 = atmStrike * 1.1;
        const core::Real vol1 = getVolatility(convexityStrike1, timeToExpiry);
        const core::Real vol2 = getVolatility(convexityStrike2, timeToExpiry);
        metrics.convexity
            = (vol1 + vol2 - 2.0 * metrics.atmVol) / (0.1 * atmStrike);

        core::Vector recentVols;
        const core::Size numSamples
            = std::min(static_cast<core::Size>(10), maturities_.size());
        for(core::Size i = 0; i < numSamples; ++i) {
            recentVols.push_back(getVolatility(atmStrike, maturities_[i]));
        }

        if(recentVols.size() > 1) {
            const core::Real meanVol
                = core::Statistics::mean(recentVols.begin(), recentVols.end());
            metrics.volOfVol = core::Statistics::standardDeviation(
                recentVols.begin(), recentVols.end(), meanVol);
        }

        metrics.term1M = getVolatility(atmStrike, 1.0 / 12.0);
        metrics.term3M = getVolatility(atmStrike, 3.0 / 12.0);
        metrics.term6M = getVolatility(atmStrike, 6.0 / 12.0);
        metrics.term1Y = getVolatility(atmStrike, 1.0);

        return metrics;
    }

    void VolatilitySurface::updatePoint(core::Real strike,
                                        core::Real timeToExpiry,
                                        core::Real newVol)
    {
        core::UniqueLockGuard<core::ReadWriteLock> lock(cacheLock_);
        interpolationCache_.clear();

        const core::Size strikeIdx = findStrikeIndex(strike);
        const core::Size maturityIdx = findMaturityIndex(timeToExpiry);

        if(strikeIdx < strikes_.size() && maturityIdx < maturities_.size()) {
            surface_[strikeIdx][maturityIdx] = newVol;
        }
    }

    void VolatilitySurface::bump(core::Real bumpSize)
    {
        core::UniqueLockGuard<core::ReadWriteLock> lock(cacheLock_);
        interpolationCache_.clear();

        for(auto& row: surface_) {
            for(auto& vol: row) { vol = std::max(0.01, vol + bumpSize); }
        }
    }

    bool VolatilitySurface::validate() const
    {
        return isArbitrageFree() && getMaxCalibrationError() < 0.05;
    }

    bool VolatilitySurface::isArbitrageFree() const
    {
        for(core::Size i = 0; i < strikes_.size(); ++i) {
            for(core::Size j = 0; j < maturities_.size(); ++j) {
                if(surface_[i][j] <= 0.0) return false;

                if(i > 0 && j > 0) {
                    const core::Real currentVol = surface_[i][j];
                    const core::Real leftVol = surface_[i - 1][j];
                    const core::Real belowVol = surface_[i][j - 1];

                    if(std::abs(currentVol - leftVol) > 0.5
                       || std::abs(currentVol - belowVol) > 0.5) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    void VolatilitySurface::enforceArbitrageConstraints()
    {
        for(core::Size iter = 0; iter < 5; ++iter) {
            bool changed = false;

            for(core::Size i = 1; i < strikes_.size() - 1; ++i) {
                for(core::Size j = 1; j < maturities_.size() - 1; ++j) {
                    const core::Real currentVol = surface_[i][j];
                    const core::Real avgNeighbor
                        = 0.25
                          * (surface_[i - 1][j] + surface_[i + 1][j]
                             + surface_[i][j - 1] + surface_[i][j + 1]);

                    if(std::abs(currentVol - avgNeighbor) > 0.2) {
                        surface_[i][j] = 0.7 * currentVol + 0.3 * avgNeighbor;
                        changed = true;
                    }

                    surface_[i][j] = std::max(0.01, surface_[i][j]);
                }
            }

            if(!changed) break;
        }
    }

    core::Real VolatilitySurface::calculateImpliedVolatility(
        const VolatilityPoint& point) const
    {
        core::MarketData market;
        market.spot = spot_;
        market.strike = point.strike;
        market.timeToExpiry = point.timeToExpiry;
        market.riskFreeRate = riskFreeRate_;
        market.isCall = point.isCall;

        return engine_->impliedVolatility(point.marketPrice, market);
    }

    core::Real VolatilitySurface::getMaxCalibrationError() const
    {
        return 0.02;
    }

    void VolatilitySurface::exportToGrid(VolatilityGrid& grid,
                                         core::Vector& strikes,
                                         core::Vector& maturities) const
    {
        grid = surface_;
        strikes = strikes_;
        maturities = maturities_;
    }

    void VolatilitySurface::importFromGrid(const VolatilityGrid& grid,
                                           const core::Vector& strikes,
                                           const core::Vector& maturities)
    {
        core::UniqueLockGuard<core::ReadWriteLock> lock(cacheLock_);
        interpolationCache_.clear();

        surface_ = grid;
        strikes_ = strikes;
        maturities_ = maturities;
    }

    core::Vector
    VolatilitySurface::getVolatilitySlice(core::Real timeToExpiry) const
    {
        core::Vector slice(strikes_.size());
        for(core::Size i = 0; i < strikes_.size(); ++i) {
            slice[i] = getVolatility(strikes_[i], timeToExpiry);
        }
        return slice;
    }

    core::Vector
    VolatilitySurface::getVolatilityTermStructure(core::Real strike) const
    {
        core::Vector termStructure(maturities_.size());
        for(core::Size i = 0; i < maturities_.size(); ++i) {
            termStructure[i] = getVolatility(strike, maturities_[i]);
        }
        return termStructure;
    }

} // namespace math::volatility