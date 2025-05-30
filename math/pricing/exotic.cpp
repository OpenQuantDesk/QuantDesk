#include "exotic.hpp"
#include "../core/utils.hpp"
#include <complex>

namespace math::pricing {

RainbowEngine::RainbowEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real RainbowEngine::maxOption(const core::Vector& spots, const core::Vector& strikes,
                                   core::Real timeToExpiry, core::Real riskFreeRate,
                                   const core::Matrix& correlations, const core::Vector& vols,
                                   const core::Vector& dividends) const {
    if (spots.size() != 2) {
        return monteCarloRainbow(spots, strikes, timeToExpiry, riskFreeRate, 
                               correlations, vols, dividends, true, true);
    }
    
    const core::Real S1 = spots[0];
    const core::Real S2 = spots[1];
    const core::Real K = strikes.empty() ? 0.0 : strikes[0];
    const core::Real T = timeToExpiry;
    const core::Real r = riskFreeRate;
    const core::Real v1 = vols[0];
    const core::Real v2 = vols[1];
    const core::Real rho = correlations[0][1];
    const core::Real q1 = dividends.empty() ? 0.0 : dividends[0];
    const core::Real q2 = dividends.size() < 2 ? 0.0 : dividends[1];
    
    const core::Real sqrtT = std::sqrt(T);
    const core::Real discount = std::exp(-r * T);
    
    const core::Real d1 = (std::log(S1 / S2) + (q2 - q1 + 0.5 * (v1 * v1 - 2.0 * rho * v1 * v2 + v2 * v2)) * T) /
                          (std::sqrt(v1 * v1 - 2.0 * rho * v1 * v2 + v2 * v2) * sqrtT);
    
    const core::Real y1 = (std::log(S1 / K) + (r - q1 + 0.5 * v1 * v1) * T) / (v1 * sqrtT);
    const core::Real y2 = (std::log(S2 / K) + (r - q2 + 0.5 * v2 * v2) * T) / (v2 * sqrtT);
    
    const core::Real rho1 = (rho * v2 - v1) / std::sqrt(v1 * v1 - 2.0 * rho * v1 * v2 + v2 * v2);
    const core::Real rho2 = (rho * v1 - v2) / std::sqrt(v1 * v1 - 2.0 * rho * v1 * v2 + v2 * v2);
    
    const core::Real M1 = bivariateCDF(d1, y1, rho1);
    const core::Real M2 = bivariateCDF(-d1, y2, rho2);
    const core::Real M3 = bivariateCDF(d1, y1 - v1 * sqrtT, rho1);
    const core::Real M4 = bivariateCDF(-d1, y2 - v2 * sqrtT, rho2);
    
    const core::Real term1 = S1 * std::exp(-q1 * T) * M1;
    const core::Real term2 = S2 * std::exp(-q2 * T) * M2;
    const core::Real term3 = K * discount * (M3 + M4);
    
    return term1 + term2 - term3;
}

core::Real RainbowEngine::minOption(const core::Vector& spots, const core::Vector& strikes,
                                   core::Real timeToExpiry, core::Real riskFreeRate,
                                   const core::Matrix& correlations, const core::Vector& vols,
                                   const core::Vector& dividends) const {
    if (spots.size() != 2) {
        return monteCarloRainbow(spots, strikes, timeToExpiry, riskFreeRate,
                               correlations, vols, dividends, true, false);
    }
    
    const core::Real maxPrice = maxOption(spots, strikes, timeToExpiry, riskFreeRate,
                                        correlations, vols, dividends);
    
    core::Real vanillaSum = 0.0;
    for (core::Size i = 0; i < spots.size(); ++i) {
        core::MarketData market;
        market.spot = spots[i];
        market.strike = strikes.empty() ? 0.0 : strikes[0];
        market.timeToExpiry = timeToExpiry;
        market.riskFreeRate = riskFreeRate;
        market.volatility = vols[i];
        market.dividendYield = i < dividends.size() ? dividends[i] : 0.0;
        market.isCall = true;
        
        vanillaSum += engine_->blackScholes(market).price;
    }
    
    return vanillaSum - maxPrice;
}

core::Real RainbowEngine::spreadOption(core::Real S1, core::Real S2, core::Real K,
                                      core::Real timeToExpiry, core::Real riskFreeRate,
                                      core::Real vol1, core::Real vol2, core::Real correlation,
                                      core::Real dividend1, core::Real dividend2) const {
    const core::Real T = timeToExpiry;
    const core::Real r = riskFreeRate;
    const core::Real q1 = dividend1;
    const core::Real q2 = dividend2;
    const core::Real rho = correlation;
    
    const core::Real F1 = S1 * std::exp((r - q1) * T);
    const core::Real F2 = S2 * std::exp((r - q2) * T);
    
    const core::Real sigma = std::sqrt(vol1 * vol1 + vol2 * vol2 - 2.0 * rho * vol1 * vol2);
    const core::Real sqrtT = std::sqrt(T);
    
    const core::Real d1 = (std::log((F1 - F2) / K) + 0.5 * sigma * sigma * T) / (sigma * sqrtT);
    const core::Real d2 = d1 - sigma * sqrtT;
    
    const core::Real discount = std::exp(-r * T);
    
    return discount * ((F1 - F2) * engine_->normalCDF(d1) - K * engine_->normalCDF(d2));
}

core::Real RainbowEngine::monteCarloRainbow(const core::Vector& spots, const core::Vector& strikes,
                                           core::Real timeToExpiry, core::Real riskFreeRate,
                                           const core::Matrix& correlations, const core::Vector& vols,
                                           const core::Vector& dividends, bool isCall, bool isMax,
                                           core::Integer numSimulations) const {
    if (spots.empty() || vols.empty()) return 0.0;
    
    const core::Size numAssets = spots.size();
    const core::Real T = timeToExpiry;
    const core::Real r = riskFreeRate;
    const core::Real dt = T;
    
    core::Matrix cholesky = choleskyDecomposition(correlations);
    
    auto& rng = engine_->getRNG();
    std::normal_distribution<core::Real> normal(0.0, 1.0);
    
    core::Real totalPayoff = 0.0;
    
    for (core::Integer sim = 0; sim < numSimulations; ++sim) {
        core::Vector finalPrices(numAssets);
        core::Vector randoms(numAssets);
        
        for (core::Size i = 0; i < numAssets; ++i) {
            randoms[i] = normal(rng);
        }
        
        core::Vector correlatedRandoms(numAssets, 0.0);
        for (core::Size i = 0; i < numAssets; ++i) {
            for (core::Size j = 0; j <= i && j < cholesky[i].size(); ++j) {
                correlatedRandoms[i] += cholesky[i][j] * randoms[j];
            }
        }
        
        for (core::Size i = 0; i < numAssets; ++i) {
            const core::Real q = i < dividends.size() ? dividends[i] : 0.0;
            const core::Real drift = (r - q - 0.5 * vols[i] * vols[i]) * dt;
            const core::Real diffusion = vols[i] * std::sqrt(dt) * correlatedRandoms[i];
            
            finalPrices[i] = spots[i] * std::exp(drift + diffusion);
        }
        
        core::Real payoff = 0.0;
        
        if (isMax) {
            const core::Real maxPrice = *std::max_element(finalPrices.begin(), finalPrices.end());
            const core::Real strike = strikes.empty() ? 0.0 : strikes[0];
            
            if (isCall) {
                payoff = std::max(0.0, maxPrice - strike);
            } else {
                payoff = std::max(0.0, strike - maxPrice);
            }
        } else {
            const core::Real minPrice = *std::min_element(finalPrices.begin(), finalPrices.end());
            const core::Real strike = strikes.empty() ? 0.0 : strikes[0];
            
            if (isCall) {
                payoff = std::max(0.0, minPrice - strike);
            } else {
                payoff = std::max(0.0, strike - minPrice);
            }
        }
        
        totalPayoff += payoff;
    }
    
    const core::Real avgPayoff = totalPayoff / numSimulations;
    return avgPayoff * std::exp(-r * T);
}

BermudanEngine::BermudanEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real BermudanEngine::price(const core::MarketData& market, const core::Vector& exerciseTimes,
                                const ExerciseParams& params) const {
    if (exerciseTimes.empty()) {
        return engine_->blackScholes(market).price;
    }
    
    const core::Size numSteps = params.timeSteps;
    const core::Size numPaths = params.numPaths;
    const core::Real dt = market.timeToExpiry / numSteps;
    
    const core::Real drift = (market.riskFreeRate - market.dividendYield - 
                             0.5 * market.volatility * market.volatility) * dt;
    const core::Real diffusion = market.volatility * std::sqrt(dt);
    
    auto& rng = engine_->getRNG();
    std::normal_distribution<core::Real> normal(0.0, 1.0);
    
    core::Matrix spotPaths(numPaths, core::Vector(numSteps + 1));
    core::Matrix payoffMatrix(numPaths, core::Vector(exerciseTimes.size()));
    
    for (core::Size path = 0; path < numPaths; ++path) {
        spotPaths[path][0] = market.spot;
        
        for (core::Size step = 1; step <= numSteps; ++step) {
            const core::Real z = normal(rng);
            spotPaths[path][step] = spotPaths[path][step - 1] * std::exp(drift + diffusion * z);
        }
        
        for (core::Size ex = 0; ex < exerciseTimes.size(); ++ex) {
            const core::Size timeIndex = static_cast<core::Size>(exerciseTimes[ex] / dt);
            if (timeIndex <= numSteps) {
                const core::Real spot = spotPaths[path][timeIndex];
                
                if (market.isCall) {
                    payoffMatrix[path][ex] = std::max(0.0, spot - market.strike);
                } else {
                    payoffMatrix[path][ex] = std::max(0.0, market.strike - spot);
                }
            }
        }
    }
    
    core::Vector optimalExercise(numPaths);
    for (core::Size path = 0; path < numPaths; ++path) {
        core::Real maxPayoff = 0.0;
        core::Size bestExercise = exerciseTimes.size() - 1;
        
        for (core::Size ex = 0; ex < exerciseTimes.size(); ++ex) {
            const core::Real discountedPayoff = payoffMatrix[path][ex] * 
                                               std::exp(-market.riskFreeRate * exerciseTimes[ex]);
            if (discountedPayoff > maxPayoff) {
                maxPayoff = discountedPayoff;
                bestExercise = ex;
            }
        }
        
        optimalExercise[path] = maxPayoff;
    }
    
    return core::Statistics::mean(optimalExercise.begin(), optimalExercise.end());
}

core::Real BermudanEngine::lsmPrice(const core::MarketData& market, const core::Vector& exerciseTimes,
                                   const ExerciseParams& params) const {
    const core::Size numSteps = params.timeSteps;
    const core::Size numPaths = params.numPaths;
    const core::Real dt = market.timeToExpiry / numSteps;
    
    const core::Real drift = (market.riskFreeRate - market.dividendYield - 
                             0.5 * market.volatility * market.volatility) * dt;
    const core::Real diffusion = market.volatility * std::sqrt(dt);
    
    auto& rng = engine_->getRNG();
    std::normal_distribution<core::Real> normal(0.0, 1.0);
    
    core::Matrix spotPaths(numPaths, core::Vector(numSteps + 1));
    core::Vector cashFlows(numPaths, 0.0);
    
    for (core::Size path = 0; path < numPaths; ++path) {
        spotPaths[path][0] = market.spot;
        
        for (core::Size step = 1; step <= numSteps; ++step) {
            const core::Real z = normal(rng);
            spotPaths[path][step] = spotPaths[path][step - 1] * std::exp(drift + diffusion * z);
        }
    }
    
    for (core::Integer ex = static_cast<core::Integer>(exerciseTimes.size()) - 1; ex >= 0; --ex) {
        const core::Real exerciseTime = exerciseTimes[ex];
        const core::Size timeIndex = static_cast<core::Size>(exerciseTime / dt);
        
        if (timeIndex > numSteps) continue;
        
        core::Vector exerciseValues(numPaths);
        core::Vector continuationValues(numPaths);
        
        for (core::Size path = 0; path < numPaths; ++path) {
            const core::Real spot = spotPaths[path][timeIndex];
            
            if (market.isCall) {
                exerciseValues[path] = std::max(0.0, spot - market.strike);
            } else {
                exerciseValues[path] = std::max(0.0, market.strike - spot);
            }
            
            if (ex == static_cast<core::Integer>(exerciseTimes.size()) - 1) {
                continuationValues[path] = exerciseValues[path];
            } else {
                continuationValues[path] = cashFlows[path] * 
                                         std::exp(-market.riskFreeRate * 
                                                (market.timeToExpiry - exerciseTime));
            }
        }
        
        for (core::Size path = 0; path < numPaths; ++path) {
            if (exerciseValues[path] > continuationValues[path]) {
                cashFlows[path] = exerciseValues[path] * 
                                std::exp(-market.riskFreeRate * exerciseTime);
            }
        }
    }
    
    return core::Statistics::mean(cashFlows.begin(), cashFlows.end());
}

QuantoEngine::QuantoEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real QuantoEngine::quantoOption(const core::MarketData& market, core::Real fxVol,
                                     core::Real correlation, core::Real foreignRate) const {
    const core::Real adjustedRate = market.riskFreeRate - correlation * market.volatility * fxVol;
    
    core::MarketData quantoMarket = market;
    quantoMarket.riskFreeRate = adjustedRate;
    
    return engine_->blackScholes(quantoMarket).price;
}

core::Greeks QuantoEngine::quantoGreeks(const core::MarketData& market, core::Real fxVol,
                                       core::Real correlation, core::Real foreignRate) const {
    const core::Real adjustedRate = market.riskFreeRate - correlation * market.volatility * fxVol;
    
    core::MarketData quantoMarket = market;
    quantoMarket.riskFreeRate = adjustedRate;
    
    return engine_->blackScholes(quantoMarket);
}

CompositeEngine::CompositeEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real CompositeEngine::powerOption(const core::MarketData& market, core::Real power) const {
    if (std::abs(power - 1.0) < core::EPSILON) {
        return engine_->blackScholes(market).price;
    }
    
    const core::Real S = market.spot;
    const core::Real K = market.strike;
    const core::Real T = market.timeToExpiry;
    const core::Real r = market.riskFreeRate;
    const core::Real v = market.volatility;
    const core::Real q = market.dividendYield;
    
    const core::Real adjustedVol = v * std::abs(power);
    const core::Real adjustedRate = r * power + 0.5 * v * v * power * (power - 1.0);
    const core::Real adjustedStrike = std::pow(K, 1.0 / power);
    
    core::MarketData adjustedMarket = market;
    adjustedMarket.volatility = adjustedVol;
    adjustedMarket.riskFreeRate = adjustedRate;
    adjustedMarket.strike = adjustedStrike;
    adjustedMarket.spot = std::pow(S, power);
    
    return engine_->blackScholes(adjustedMarket).price;
}

core::Real CompositeEngine::logOption(const core::MarketData& market) const {
    const core::Real S = market.spot;
    const core::Real K = market.strike;
    const core::Real T = market.timeToExpiry;
    const core::Real r = market.riskFreeRate;
    const core::Real v = market.volatility;
    const core::Real q = market.dividendYield;
    
    const core::Real logS = std::log(S);
    const core::Real logK = std::log(K);
    
    const core::Real mu = r - q - 0.5 * v * v;
    const core::Real expectedLogS = logS + mu * T;
    const core::Real varianceLogS = v * v * T;
    
    const core::Real d1 = (expectedLogS - logK) / std::sqrt(varianceLogS);
    const core::Real d2 = d1 - std::sqrt(varianceLogS);
    
    const core::Real discount = std::exp(-r * T);
    
    if (market.isCall) {
        return discount * (expectedLogS * engine_->normalCDF(d1) - logK * engine_->normalCDF(d2));
    } else {
        return discount * (logK * engine_->normalCDF(-d2) - expectedLogS * engine_->normalCDF(-d1));
    }
}

core::Matrix RainbowEngine::choleskyDecomposition(const core::Matrix& correlation) const {
    const core::Size n = correlation.size();
    core::Matrix cholesky(n, core::Vector(n, 0.0));
    
    for (core::Size i = 0; i < n; ++i) {
        for (core::Size j = 0; j <= i; ++j) {
            if (i == j) {
                core::Real sum = 0.0;
                for (core::Size k = 0; k < j; ++k) {
                    sum += cholesky[j][k] * cholesky[j][k];
                }
                cholesky[j][j] = std::sqrt(correlation[j][j] - sum);
            } else {
                core::Real sum = 0.0;
                for (core::Size k = 0; k < j; ++k) {
                    sum += cholesky[i][k] * cholesky[j][k];
                }
                cholesky[i][j] = (correlation[i][j] - sum) / cholesky[j][j];
            }
        }
    }
    
    return cholesky;
}

core::Real RainbowEngine::bivariateCDF(core::Real x, core::Real y, core::Real rho) const {
    if (std::abs(rho) < core::EPSILON) {
        return engine_->normalCDF(x) * engine_->normalCDF(y);
    }
    
    const core::Integer n = 25;
    core::Real sum = 0.0;
    
    for (core::Integer i = 0; i < n; ++i) {
        const core::Real xi = -5.0 + 10.0 * i / (n - 1);
        const core::Real weight = (i == 0 || i == n - 1) ? 1.0 : ((i % 2 == 0) ? 2.0 : 4.0);
        
        if (xi <= x) {
            const core::Real arg = (y - rho * xi) / std::sqrt(1.0 - rho * rho);
            sum += weight * engine_->normalPDF(xi) * engine_->normalCDF(arg);
        }
    }
    
    return sum * 10.0 / (3.0 * n);
}

}