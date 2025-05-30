#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include "../core/threading.hpp"
#include <map>
#include <functional>

namespace math::volatility {

struct VolatilityPoint {
    core::Real strike;
    core::Real timeToExpiry;
    core::Real impliedVol;
    core::Real marketPrice;
    core::Real bid;
    core::Real ask;
    core::Real volume;
    core::Real openInterest;
    bool isCall;
};

using VolatilityGrid = std::vector<std::vector<core::Real>>;
using VolatilityPoints = std::vector<VolatilityPoint>;

class VolatilitySurface {
private:
    std::shared_ptr<core::MathEngine> engine_;
    VolatilityGrid surface_;
    core::Vector strikes_;
    core::Vector maturities_;
    core::Real spot_;
    core::Real riskFreeRate_;
    
    mutable std::map<std::pair<core::Real, core::Real>, core::Real> interpolationCache_;
    mutable core::ReadWriteLock cacheLock_;
    
public:
    explicit VolatilitySurface(std::shared_ptr<core::MathEngine> engine);
    
    void calibrate(const VolatilityPoints& marketData, core::Real spot, core::Real riskFreeRate);
    
    core::Real getVolatility(core::Real strike, core::Real timeToExpiry) const;
    core::Real getVolatilityMoneyness(core::Real moneyness, core::Real timeToExpiry) const;
    core::Real getVolatilityDelta(core::Real delta, core::Real timeToExpiry, bool isCall = true) const;
    
    core::Real getForwardVolatility(core::Real strike, core::Real startTime, core::Real endTime) const;
    core::Real getLocalVolatility(core::Real strike, core::Real timeToExpiry) const;
    
    core::Vector getVolatilitySlice(core::Real timeToExpiry) const;
    core::Vector getVolatilityTermStructure(core::Real strike) const;
    
    struct VolatilityMetrics {
        core::Real atmVol;
        core::Real skew25Delta;
        core::Real skew10Delta;
        core::Real convexity;
        core::Real butterfly25Delta;
        core::Real butterfly10Delta;
        core::Real riskReversal25Delta;
        core::Real riskReversal10Delta;
        core::Real volOfVol;
        core::Real term1M;
        core::Real term3M;
        core::Real term6M;
        core::Real term1Y;
    };
    
    VolatilityMetrics calculateMetrics(core::Real timeToExpiry) const;
    
    void updatePoint(core::Real strike, core::Real timeToExpiry, core::Real newVol);
    void bump(core::Real bumpSize);
    void bumpStrike(core::Real strike, core::Real bumpSize);
    void bumpTenor(core::Real timeToExpiry, core::Real bumpSize);
    
    bool validate() const;
    core::Real getMaxCalibrationError() const;
    
    void exportToGrid(VolatilityGrid& grid, core::Vector& strikes, core::Vector& maturities) const;
    void importFromGrid(const VolatilityGrid& grid, const core::Vector& strikes, 
                       const core::Vector& maturities);
    
private:
    void buildSurface(const VolatilityPoints& marketData);
    void smoothSurface();
    void extrapolateSurface();
    
    core::Real bilinearInterpolation(core::Real strike, core::Real timeToExpiry) const;
    core::Real bicubicInterpolation(core::Real strike, core::Real timeToExpiry) const;
    core::Real splineInterpolation(core::Real strike, core::Real timeToExpiry) const;
    
    core::Size findStrikeIndex(core::Real strike) const;
    core::Size findMaturityIndex(core::Real timeToExpiry) const;
    
    bool isArbitrageFree() const;
    void enforceArbitrageConstraints();
    
    core::Real calculateImpliedVolatility(const VolatilityPoint& point) const;
};

class SVIModel {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit SVIModel(std::shared_ptr<core::MathEngine> engine);
    
    struct SVIParams {
        core::Real a;
        core::Real b;
        core::Real rho;
        core::Real m;
        core::Real sigma;
    };
    
    core::Real volatility(core::Real logMoneyness, const SVIParams& params) const;
    core::Real totalVariance(core::Real logMoneyness, const SVIParams& params) const;
    
    SVIParams calibrate(const core::Vector& logMoneyness, const core::Vector& totalVariance,
                       const core::Vector& weights = {}) const;
    
    bool validateParams(const SVIParams& params) const;
    SVIParams enforceConstraints(const SVIParams& params) const;
    
    core::Real calculateRMSE(const core::Vector& logMoneyness, const core::Vector& totalVariance,
                            const SVIParams& params) const;
    
private:
    core::Real objectiveFunction(const core::Vector& params, const core::Vector& logMoneyness,
                                const core::Vector& totalVariance, const core::Vector& weights) const;
    
    core::Vector paramsToVector(const SVIParams& params) const;
    SVIParams vectorToParams(const core::Vector& vec) const;
};

class SSVIModel {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit SSVIModel(std::shared_ptr<core::MathEngine> engine);
    
    struct SSVIParams {
        core::Real rho;
        core::Real eta;
        core::Real lambda;
    };
    
    struct TimeSliceParams {
        core::Real atmTotalVariance;
        core::Real atmSkew;
        core::Real atmCurvature;
    };
    
    core::Real volatility(core::Real logMoneyness, core::Real timeToExpiry,
                         const SSVIParams& globalParams, 
                         const std::map<core::Real, TimeSliceParams>& timeSlices) const;
    
    SSVIParams calibrateGlobal(const std::map<core::Real, core::Vector>& logMoneynessGrid,
                              const std::map<core::Real, core::Vector>& totalVarianceGrid) const;
    
    TimeSliceParams calibrateTimeSlice(core::Real timeToExpiry,
                                      const core::Vector& logMoneyness,
                                      const core::Vector& totalVariance,
                                      const SSVIParams& globalParams) const;
    
private:
    core::Real phi(core::Real eta) const;
    core::Real psi(core::Real eta, core::Real lambda, core::Real rho) const;
};

class HestonModel {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit HestonModel(std::shared_ptr<core::MathEngine> engine);
    
    struct HestonParams {
        core::Real kappa;
        core::Real theta;
        core::Real sigma;
        core::Real rho;
        core::Real v0;
    };
    
    core::Real characteristicFunction(const std::complex<core::Real>& u, core::Real tau,
                                     const HestonParams& params) const;
    
    core::Real optionPrice(core::Real spot, core::Real strike, core::Real timeToExpiry,
                          core::Real riskFreeRate, const HestonParams& params,
                          bool isCall = true) const;
    
    core::Real impliedVolatility(core::Real spot, core::Real strike, core::Real timeToExpiry,
                                core::Real riskFreeRate, const HestonParams& params,
                                bool isCall = true) const;
    
    HestonParams calibrate(const VolatilityPoints& marketData, core::Real spot,
                          core::Real riskFreeRate) const;
    
    bool validateParams(const HestonParams& params) const;
    HestonParams enforceConstraints(const HestonParams& params) const;
    
    core::Vector simulatePath(const HestonParams& params, core::Real spot,
                             core::Real timeToExpiry, core::Integer timeSteps,
                             std::mt19937& rng) const;
    
private:
    std::complex<core::Real> complexLog(const std::complex<core::Real>& z) const;
    core::Real integrand(core::Real u, core::Real spot, core::Real strike,
                        core::Real timeToExpiry, const HestonParams& params,
                        bool isCall, bool isFirstIntegral) const;
    
    core::Real trapezoidalIntegration(const std::function<core::Real(core::Real)>& func,
                                     core::Real a, core::Real b, core::Size n = 1000) const;
};

class LocalVolatilityModel {
private:
    std::shared_ptr<core::MathEngine> engine_;
    std::shared_ptr<VolatilitySurface> impliedSurface_;
    
public:
    explicit LocalVolatilityModel(std::shared_ptr<core::MathEngine> engine,
                                 std::shared_ptr<VolatilitySurface> impliedSurface);
    
    core::Real localVolatility(core::Real spot, core::Real timeToExpiry) const;
    
    void calibrateFromImpliedSurface();
    
    core::Real optionPrice(core::Real spot, core::Real strike, core::Real timeToExpiry,
                          core::Real riskFreeRate, bool isCall = true) const;
    
    core::Vector simulatePath(core::Real spot, core::Real timeToExpiry,
                             core::Integer timeSteps, std::mt19937& rng) const;
    
private:
    core::Real dupireFormula(core::Real strike, core::Real timeToExpiry) const;
    core::Real numericalDerivative(const std::function<core::Real(core::Real)>& func,
                                  core::Real x, core::Real h = 1e-6) const;
};

inline VolatilitySurface::VolatilitySurface(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)), spot_(0.0), riskFreeRate_(0.0) {}

inline core::Real VolatilitySurface::getVolatility(core::Real strike, core::Real timeToExpiry) const {
    const auto key = std::make_pair(strike, timeToExpiry);
    
    {
        core::SharedLockGuard<core::ReadWriteLock> lock(cacheLock_);
        auto it = interpolationCache_.find(key);
        if (it != interpolationCache_.end()) {
            return it->second;
        }
    }
    
    const core::Real vol = bicubicInterpolation(strike, timeToExpiry);
    
    {
        core::UniqueLockGuard<core::ReadWriteLock> lock(cacheLock_);
        interpolationCache_[key] = vol;
    }
    
    return vol;
}

inline core::Real VolatilitySurface::getVolatilityMoneyness(core::Real moneyness, core::Real timeToExpiry) const {
    const core::Real strike = spot_ * moneyness;
    return getVolatility(strike, timeToExpiry);
}

inline core::Real VolatilitySurface::bilinearInterpolation(core::Real strike, core::Real timeToExpiry) const {
    const core::Size strikeIdx = findStrikeIndex(strike);
    const core::Size maturityIdx = findMaturityIndex(timeToExpiry);
    
    if (strikeIdx >= strikes_.size() - 1 || maturityIdx >= maturities_.size() - 1) {
        return surface_[std::min(strikeIdx, strikes_.size() - 1)]
                      [std::min(maturityIdx, maturities_.size() - 1)];
    }
    
    const core::Real x1 = strikes_[strikeIdx];
    const core::Real x2 = strikes_[strikeIdx + 1];
    const core::Real y1 = maturities_[maturityIdx];
    const core::Real y2 = maturities_[maturityIdx + 1];
    
    const core::Real q11 = surface_[strikeIdx][maturityIdx];
    const core::Real q12 = surface_[strikeIdx][maturityIdx + 1];
    const core::Real q21 = surface_[strikeIdx + 1][maturityIdx];
    const core::Real q22 = surface_[strikeIdx + 1][maturityIdx + 1];
    
    const core::Real wx = (strike - x1) / (x2 - x1);
    const core::Real wy = (timeToExpiry - y1) / (y2 - y1);
    
    return q11 * (1.0 - wx) * (1.0 - wy) +
           q21 * wx * (1.0 - wy) +
           q12 * (1.0 - wx) * wy +
           q22 * wx * wy;
}

inline core::Size VolatilitySurface::findStrikeIndex(core::Real strike) const {
    auto it = std::lower_bound(strikes_.begin(), strikes_.end(), strike);
    return std::distance(strikes_.begin(), it);
}

inline core::Size VolatilitySurface::findMaturityIndex(core::Real timeToExpiry) const {
    auto it = std::lower_bound(maturities_.begin(), maturities_.end(), timeToExpiry);
    return std::distance(maturities_.begin(), it);
}

inline SVIModel::SVIModel(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

inline core::Real SVIModel::volatility(core::Real logMoneyness, const SVIParams& params) const {
    return std::sqrt(totalVariance(logMoneyness, params));
}

inline core::Real SVIModel::totalVariance(core::Real logMoneyness, const SVIParams& params) const {
    const core::Real k = logMoneyness;
    return params.a + params.b * (params.rho * (k - params.m) + 
                                 std::sqrt((k - params.m) * (k - params.m) + params.sigma * params.sigma));
}

inline bool SVIModel::validateParams(const SVIParams& params) const {
    return params.a >= 0.0 && params.b >= 0.0 && 
           params.rho >= -1.0 && params.rho <= 1.0 && 
           params.sigma > 0.0 &&
           params.a + params.b * params.sigma * std::sqrt(1.0 - params.rho * params.rho) >= 0.0;
}

inline HestonModel::HestonModel(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

inline bool HestonModel::validateParams(const HestonParams& params) const {
    return params.kappa > 0.0 && params.theta > 0.0 && params.sigma > 0.0 &&
           params.rho >= -1.0 && params.rho <= 1.0 && params.v0 > 0.0 &&
           2.0 * params.kappa * params.theta >= params.sigma * params.sigma;
}

}