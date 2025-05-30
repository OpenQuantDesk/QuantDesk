#include "calibration.hpp"
#include "../core/utils.hpp"
#include <algorithm>
#include <numeric>

namespace math::volatility {

VolatilityCalibrator::VolatilityCalibrator(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

VolatilityCalibrator::CalibrationResult VolatilityCalibrator::calibrateImpliedSurface(
    const VolatilityPoints& marketData, core::Real spot, core::Real riskFreeRate,
    const CalibrationParams& params) const {
    
    CalibrationResult result;
    result.success = false;
    
    if (marketData.empty()) {
        result.errorMessage = "No market data provided";
        return result;
    }
    
    std::set<core::Real> strikeSet, maturitySet;
    for (const auto& point : marketData) {
        strikeSet.insert(point.strike);
        maturitySet.insert(point.timeToExpiry);
    }
    
    core::Vector strikes(strikeSet.begin(), strikeSet.end());
    core::Vector maturities(maturitySet.begin(), maturitySet.end());
    
    result.strikes = strikes;
    result.maturities = maturities;
    result.impliedVols.resize(strikes.size(), core::Vector(maturities.size(), 0.0));
    
    for (const auto& point : marketData) {
        auto strikeIt = std::find(strikes.begin(), strikes.end(), point.strike);
        auto maturityIt = std::find(maturities.begin(), maturities.end(), point.timeToExpiry);
        
        if (strikeIt != strikes.end() && maturityIt != maturities.end()) {
            core::Size i = std::distance(strikes.begin(), strikeIt);
            core::Size j = std::distance(maturities.begin(), maturityIt);
            
            if (point.impliedVol > 0.0) {
                result.impliedVols[i][j] = point.impliedVol;
            } else {
                core::MarketData market;
                market.spot = spot;
                market.strike = point.strike;
                market.timeToExpiry = point.timeToExpiry;
                market.riskFreeRate = riskFreeRate;
                market.isCall = point.isCall;
                
                result.impliedVols[i][j] = engine_->impliedVolatility(point.marketPrice, market);
            }
        }
    }
    
    interpolateMissingPoints(result.impliedVols, strikes, maturities);
    
    if (params.enforceArbitrage) {
        enforceArbitrageConstraints(result.impliedVols, strikes, maturities, spot, riskFreeRate);
    }
    
    if (params.smoothing > 0.0) {
        applySmoothening(result.impliedVols, params.smoothing);
    }
    
    result.calibrationError = calculateCalibrationError(marketData, result.impliedVols, 
                                                       strikes, maturities, spot, riskFreeRate);
    result.success = result.calibrationError < params.tolerance;
    
    if (!result.success) {
        result.errorMessage = "Calibration error exceeds tolerance: " + 
                             std::to_string(result.calibrationError);
    }
    
    return result;
}

VolatilityCalibrator::ModelCalibrationResult VolatilityCalibrator::calibrateHeston(
    const VolatilityPoints& marketData, core::Real spot, core::Real riskFreeRate,
    const CalibrationParams& params) const {
    
    ModelCalibrationResult result;
    result.success = false;
    
    HestonModel::HestonParams hestonParams;
    hestonParams.kappa = 2.0;
    hestonParams.theta = 0.04;
    hestonParams.sigma = 0.3;
    hestonParams.rho = -0.5;
    hestonParams.v0 = 0.04;
    
    const core::Integer maxIterations = params.maxIterations;
    const core::Real tolerance = params.tolerance;
    core::Real bestError = std::numeric_limits<core::Real>::max();
    HestonModel::HestonParams bestParams = hestonParams;
    
    HestonModel hestonModel(engine_);
    
    for (core::Integer iter = 0; iter < maxIterations; ++iter) {
        core::Real totalError = 0.0;
        core::Size validPoints = 0;
        
        for (const auto& point : marketData) {
            try {
                const core::Real modelPrice = hestonModel.optionPrice(spot, point.strike, 
                                                                    point.timeToExpiry, riskFreeRate,
                                                                    hestonParams, point.isCall);
                const core::Real error = std::abs(modelPrice - point.marketPrice) / point.marketPrice;
                totalError += error * error;
                ++validPoints;
            } catch (...) {
                continue;
            }
        }
        
        if (validPoints == 0) break;
        
        const core::Real rmse = std::sqrt(totalError / validPoints);
        
        if (rmse < bestError) {
            bestError = rmse;
            bestParams = hestonParams;
        }
        
        if (rmse < tolerance) {
            result.success = true;
            break;
        }
        
        optimizeHestonParams(hestonParams, marketData, spot, riskFreeRate, params.stepSize);
        hestonParams = hestonModel.enforceConstraints(hestonParams);
    }
    
    result.hestonParams = bestParams;
    result.calibrationError = bestError;
    result.success = result.success || (bestError < tolerance * 2.0);
    
    if (!result.success) {
        result.errorMessage = "Heston calibration failed. RMSE: " + std::to_string(bestError);
    }
    
    return result;
}

VolatilityCalibrator::ModelCalibrationResult VolatilityCalibrator::calibrateSABR(
    const VolatilityPoints& marketData, core::Real forward, core::Real timeToExpiry,
    const CalibrationParams& params) const {
    
    ModelCalibrationResult result;
    result.success = false;
    
    core::Vector strikes, impliedVols;
    for (const auto& point : marketData) {
        if (std::abs(point.timeToExpiry - timeToExpiry) < 1e-6) {
            strikes.push_back(point.strike);
            impliedVols.push_back(point.impliedVol > 0.0 ? point.impliedVol : 0.2);
        }
    }
    
    if (strikes.size() < 4) {
        result.errorMessage = "Insufficient data points for SABR calibration";
        return result;
    }
    
    StochasticVolatilityModel sabrModel(engine_);
    
    try {
        result.sabrParams = sabrModel.calibrateSABR(strikes, impliedVols, forward, timeToExpiry);
        
        core::Real totalError = 0.0;
        for (core::Size i = 0; i < strikes.size(); ++i) {
            const core::Real modelVol = sabrModel.sabr(forward, strikes[i], timeToExpiry, result.sabrParams);
            const core::Real error = std::abs(modelVol - impliedVols[i]);
            totalError += error * error;
        }
        
        result.calibrationError = std::sqrt(totalError / strikes.size());
        result.success = result.calibrationError < params.tolerance;
        
        if (!result.success) {
            result.errorMessage = "SABR calibration error exceeds tolerance: " + 
                                 std::to_string(result.calibrationError);
        }
    } catch (const std::exception& e) {
        result.errorMessage = "SABR calibration failed: " + std::string(e.what());
    }
    
    return result;
}

void VolatilityCalibrator::interpolateMissingPoints(core::Matrix& impliedVols,
                                                   const core::Vector& strikes,
                                                   const core::Vector& maturities) const {
    const core::Size numStrikes = strikes.size();
    const core::Size numMaturities = maturities.size();
    
    for (core::Size i = 0; i < numStrikes; ++i) {
        for (core::Size j = 0; j < numMaturities; ++j) {
            if (impliedVols[i][j] <= 0.0) {
                core::Real interpolatedVol = 0.0;
                core::Real totalWeight = 0.0;
                
                for (core::Size ii = 0; ii < numStrikes; ++ii) {
                    for (core::Size jj = 0; jj < numMaturities; ++jj) {
                        if (impliedVols[ii][jj] > 0.0) {
                            const core::Real strikeDistance = std::abs(strikes[i] - strikes[ii]) / strikes[i];
                            const core::Real timeDistance = std::abs(maturities[j] - maturities[jj]);
                            const core::Real distance = std::sqrt(strikeDistance * strikeDistance + 
                                                                timeDistance * timeDistance);
                            
                            if (distance > 0.0) {
                                const core::Real weight = 1.0 / (distance + 1e-8);
                                interpolatedVol += weight * impliedVols[ii][jj];
                                totalWeight += weight;
                            }
                        }
                    }
                }
                
                if (totalWeight > 0.0) {
                    impliedVols[i][j] = interpolatedVol / totalWeight;
                } else {
                    impliedVols[i][j] = 0.2;
                }
            }
        }
    }
}

void VolatilityCalibrator::enforceArbitrageConstraints(core::Matrix& impliedVols,
                                                      const core::Vector& strikes,
                                                      const core::Vector& maturities,
                                                      core::Real spot, core::Real riskFreeRate) const {
    const core::Size numStrikes = strikes.size();
    const core::Size numMaturities = maturities.size();
    
    for (core::Integer iter = 0; iter < 5; ++iter) {
        bool changed = false;
        
        for (core::Size i = 1; i < numStrikes - 1; ++i) {
            for (core::Size j = 1; j < numMaturities - 1; ++j) {
                const core::Real currentVol = impliedVols[i][j];
                const core::Real neighborAvg = 0.25 * (impliedVols[i-1][j] + impliedVols[i+1][j] +
                                                      impliedVols[i][j-1] + impliedVols[i][j+1]);
                
                if (std::abs(currentVol - neighborAvg) > 0.2) {
                    impliedVols[i][j] = 0.7 * currentVol + 0.3 * neighborAvg;
                    changed = true;
                }
                
                impliedVols[i][j] = std::max(0.01, impliedVols[i][j]);
            }
        }
        
        for (core::Size j = 0; j < numMaturities; ++j) {
            const core::Real timeToExpiry = maturities[j];
            const core::Real forward = spot * std::exp(riskFreeRate * timeToExpiry);
            
            for (core::Size i = 1; i < numStrikes - 1; ++i) {
                const core::Real strike = strikes[i];
                const core::Real vol = impliedVols[i][j];
                
                const core::Real localVar = vol * vol * timeToExpiry;
                const core::Real d1 = (std::log(forward / strike) + 0.5 * localVar) / (vol * std::sqrt(timeToExpiry));
                const core::Real d2 = d1 - vol * std::sqrt(timeToExpiry);
                
                if (std::abs(d1) > 5.0 || std::abs(d2) > 5.0) {
                    const core::Real newVol = std::abs(std::log(forward / strike)) / std::sqrt(timeToExpiry) / 3.0;
                    impliedVols[i][j] = std::clamp(newVol, 0.05, 2.0);
                    changed = true;
                }
            }
        }
        
        if (!changed) break;
    }
}

void VolatilityCalibrator::applySmoothening(core::Matrix& impliedVols, core::Real factor) const {
    const core::Size numStrikes = impliedVols.size();
    if (numStrikes < 3) return;
    
    const core::Size numMaturities = impliedVols[0].size();
    if (numMaturities < 3) return;
    
    core::Matrix smoothed = impliedVols;
    
    for (core::Size i = 1; i < numStrikes - 1; ++i) {
        for (core::Size j = 1; j < numMaturities - 1; ++j) {
            core::Real sum = 0.0;
            core::Real weight = 0.0;
            
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    const core::Real w = (di == 0 && dj == 0) ? 4.0 : 1.0;
                    sum += w * impliedVols[i + di][j + dj];
                    weight += w;
                }
            }
            
            const core::Real smoothedVol = sum / weight;
            smoothed[i][j] = (1.0 - factor) * impliedVols[i][j] + factor * smoothedVol;
        }
    }
    
    impliedVols = std::move(smoothed);
}

core::Real VolatilityCalibrator::calculateCalibrationError(const VolatilityPoints& marketData,
                                                          const core::Matrix& impliedVols,
                                                          const core::Vector& strikes,
                                                          const core::Vector& maturities,
                                                          core::Real spot, core::Real riskFreeRate) const {
    core::Real totalError = 0.0;
    core::Size validPoints = 0;
    
    for (const auto& point : marketData) {
        auto strikeIt = std::find(strikes.begin(), strikes.end(), point.strike);
        auto maturityIt = std::find(maturities.begin(), maturities.end(), point.timeToExpiry);
        
        if (strikeIt != strikes.end() && maturityIt != maturities.end()) {
            core::Size i = std::distance(strikes.begin(), strikeIt);
            core::Size j = std::distance(maturities.begin(), maturityIt);
            
            core::MarketData market;
            market.spot = spot;
            market.strike = point.strike;
            market.timeToExpiry = point.timeToExpiry;
            market.riskFreeRate = riskFreeRate;
            market.volatility = impliedVols[i][j];
            market.isCall = point.isCall;
            
            const core::Real theoreticalPrice = engine_->blackScholes(market).price;
            const core::Real error = std::abs(theoreticalPrice - point.marketPrice) / point.marketPrice;
            totalError += error * error;
            ++validPoints;
        }
    }
    
    return validPoints > 0 ? std::sqrt(totalError / validPoints) : 1.0;
}

void VolatilityCalibrator::optimizeHestonParams(HestonModel::HestonParams& params,
                                               const VolatilityPoints& marketData,
                                               core::Real spot, core::Real riskFreeRate,
                                               core::Real stepSize) const {
    HestonModel hestonModel(engine_);
    
    const core::Vector paramNames = {"kappa", "theta", "sigma", "rho", "v0"};
    const core::Vector stepSizes = {stepSize, stepSize * 0.1, stepSize, stepSize * 0.1, stepSize * 0.1};
    
    auto calculateError = [&](const HestonModel::HestonParams& testParams) -> core::Real {
        core::Real totalError = 0.0;
        core::Size validPoints = 0;
        
        for (const auto& point : marketData) {
            try {
                const core::Real modelPrice = hestonModel.optionPrice(spot, point.strike,
                                                                    point.timeToExpiry, riskFreeRate,
                                                                    testParams, point.isCall);
                const core::Real error = std::abs(modelPrice - point.marketPrice) / point.marketPrice;
                totalError += error * error;
                ++validPoints;
            } catch (...) {
                return std::numeric_limits<core::Real>::max();
            }
        }
        
        return validPoints > 0 ? std::sqrt(totalError / validPoints) : std::numeric_limits<core::Real>::max();
    };
    
    const core::Real baseError = calculateError(params);
    
    HestonModel::HestonParams testParams = params;
    
    testParams.kappa += stepSizes[0];
    testParams = hestonModel.enforceConstraints(testParams);
    if (calculateError(testParams) < baseError) {
        params.kappa = testParams.kappa;
    } else {
        testParams.kappa = params.kappa - stepSizes[0];
        testParams = hestonModel.enforceConstraints(testParams);
        if (calculateError(testParams) < baseError) {
            params.kappa = testParams.kappa;
        }
    }
    
    testParams = params;
    testParams.theta += stepSizes[1];
    testParams = hestonModel.enforceConstraints(testParams);
    if (calculateError(testParams) < baseError) {
        params.theta = testParams.theta;
    } else {
        testParams.theta = params.theta - stepSizes[1];
        testParams = hestonModel.enforceConstraints(testParams);
        if (calculateError(testParams) < baseError) {
            params.theta = testParams.theta;
        }
    }
    
    testParams = params;
    testParams.sigma += stepSizes[2];
    testParams = hestonModel.enforceConstraints(testParams);
    if (calculateError(testParams) < baseError) {
        params.sigma = testParams.sigma;
    } else {
        testParams.sigma = params.sigma - stepSizes[2];
        testParams = hestonModel.enforceConstraints(testParams);
        if (calculateError(testParams) < baseError) {
            params.sigma = testParams.sigma;
        }
    }
    
    testParams = params;
    testParams.rho += stepSizes[3];
    testParams = hestonModel.enforceConstraints(testParams);
    if (calculateError(testParams) < baseError) {
        params.rho = testParams.rho;
    } else {
        testParams.rho = params.rho - stepSizes[3];
        testParams = hestonModel.enforceConstraints(testParams);
        if (calculateError(testParams) < baseError) {
            params.rho = testParams.rho;
        }
    }
    
    testParams = params;
    testParams.v0 += stepSizes[4];
    testParams = hestonModel.enforceConstraints(testParams);
    if (calculateError(testParams) < baseError) {
        params.v0 = testParams.v0;
    } else {
        testParams.v0 = params.v0 - stepSizes[4];
        testParams = hestonModel.enforceConstraints(testParams);
        if (calculateError(testParams) < baseError) {
            params.v0 = testParams.v0;
        }
    }
}

VolatilityCalibrator::ValidationResult VolatilityCalibrator::validateCalibration(
    const VolatilityPoints& marketData, const CalibrationResult& calibrationResult,
    core::Real spot, core::Real riskFreeRate) const {
    
    ValidationResult result;
    result.isValid = true;
    
    const auto& strikes = calibrationResult.strikes;
    const auto& maturities = calibrationResult.maturities;
    const auto& impliedVols = calibrationResult.impliedVols;
    
    for (core::Size i = 0; i < strikes.size(); ++i) {
        for (core::Size j = 0; j < maturities.size(); ++j) {
            if (impliedVols[i][j] <= 0.0 || impliedVols[i][j] > 5.0) {
                result.isValid = false;
                result.issues.push_back("Invalid volatility value at strike " + 
                                      std::to_string(strikes[i]) + ", maturity " + 
                                      std::to_string(maturities[j]));
            }
        }
    }
    
    for (core::Size j = 0; j < maturities.size(); ++j) {
        const core::Real timeToExpiry = maturities[j];
        for (core::Size i = 1; i < strikes.size() - 1; ++i) {
            const core::Real strike = strikes[i];
            const core::Real vol = impliedVols[i][j];
            
            core::MarketData market;
            market.spot = spot;
            market.strike = strike;
            market.timeToExpiry = timeToExpiry;
            market.riskFreeRate = riskFreeRate;
            market.volatility = vol;
            market.isCall = true;
            
            const core::Greeks greeks = engine_->blackScholes(market);
            
            if (greeks.gamma < 0.0) {
                result.isValid = false;
                result.issues.push_back("Negative gamma detected at strike " + 
                                      std::to_string(strike) + ", maturity " + 
                                      std::to_string(timeToExpiry));
            }
        }
    }
    
    for (const auto& point : marketData) {
        auto strikeIt = std::find(strikes.begin(), strikes.end(), point.strike);
        auto maturityIt = std::find(maturities.begin(), maturities.end(), point.timeToExpiry);
        
        if (strikeIt != strikes.end() && maturityIt != maturities.end()) {
            core::Size i = std::distance(strikes.begin(), strikeIt);
            core::Size j = std::distance(maturities.begin(), maturityIt);
            
            core::MarketData market;
            market.spot = spot;
            market.strike = point.strike;
            market.timeToExpiry = point.timeToExpiry;
            market.riskFreeRate = riskFreeRate;
            market.volatility = impliedVols[i][j];
            market.isCall = point.isCall;
            
            const core::Real theoreticalPrice = engine_->blackScholes(market).price;
            const core::Real relativeError = std::abs(theoreticalPrice - point.marketPrice) / point.marketPrice;
            
            result.priceErrors.push_back(relativeError);
            
            if (relativeError > 0.05) {
                result.largeErrorCount++;
            }
        }
    }
    
    if (!result.priceErrors.empty()) {
        result.averageError = std::accumulate(result.priceErrors.begin(), result.priceErrors.end(), 0.0) 
                             / result.priceErrors.size();
        result.maxError = *std::max_element(result.priceErrors.begin(), result.priceErrors.end());
    }
    
    if (result.averageError > 0.02) {
        result.isValid = false;
        result.issues.push_back("Average pricing error exceeds 2%: " + std::to_string(result.averageError * 100.0) + "%");
    }
    
    if (result.maxError > 0.10) {
        result.isValid = false;
        result.issues.push_back("Maximum pricing error exceeds 10%: " + std::to_string(result.maxError * 100.0) + "%");
    }
    
    return result;
}

}