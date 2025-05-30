#include "analytical.hpp"
#include "../core/utils.hpp"
#include <complex>
#include <algorithm>

namespace math::pricing {

BarrierEngine::BarrierEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real BarrierEngine::knockOut(const core::MarketData& market, core::Real barrier,
                                  core::Real rebate, bool isUp) const {
    const core::Real S = market.spot;
    const core::Real K = market.strike;
    const core::Real B = barrier;
    const core::Real T = market.timeToExpiry;
    const core::Real r = market.riskFreeRate;
    const core::Real v = market.volatility;
    const core::Real q = market.dividendYield;
    
    const core::Real mu = (r - q - 0.5 * v * v) / (v * v);
    const core::Real lambda = std::sqrt(mu * mu + 2.0 * r / (v * v));
    
    const core::Real discount = std::exp(-r * T);
    
    const core::Real sqrtT = std::sqrt(T);
    const core::Real x1 = (std::log(S / B) + (r - q + 0.5 * v * v) * T) / (v * sqrtT);
    const core::Real y1 = (std::log(B * B / (S * K)) + (r - q + 0.5 * v * v) * T) / (v * sqrtT);
    
    core::Real price = 0.0;
    
    if (isUp && S >= B) {
        return rebate * discount;
    }
    if (!isUp && S <= B) {
        return rebate * discount;
    }
    
    const core::Real standardPrice = engine_->blackScholes(market).price;
    
    if (isUp) {
        if (K >= B) {
            price = rebate * discount;
        } else {
            const core::Real eta = market.isCall ? 1.0 : -1.0;
            const core::Real barrierTerm = std::pow(B / S, 2.0 * mu);
            
            const core::Real A = eta * S * std::exp(-q * T) * engine_->normalCDF(eta * x1);
            const core::Real B_term = eta * K * discount * engine_->normalCDF(eta * x1 - eta * v * sqrtT);
            const core::Real C = eta * S * std::exp(-q * T) * std::pow(B / S, 2.0 * (mu + 1.0)) * 
                                engine_->normalCDF(eta * y1);
            const core::Real D = eta * K * discount * barrierTerm * engine_->normalCDF(eta * y1 - eta * v * sqrtT);
            
            price = standardPrice - (A - B_term - C + D) + rebate * discount;
        }
    } else {
        if (K <= B) {
            price = rebate * discount;
        } else {
            const core::Real eta = market.isCall ? 1.0 : -1.0;
            const core::Real barrierTerm = std::pow(B / S, 2.0 * mu);
            
            const core::Real A = eta * S * std::exp(-q * T) * engine_->normalCDF(eta * x1);
            const core::Real B_term = eta * K * discount * engine_->normalCDF(eta * x1 - eta * v * sqrtT);
            const core::Real C = eta * S * std::exp(-q * T) * std::pow(B / S, 2.0 * (mu + 1.0)) * 
                                engine_->normalCDF(eta * y1);
            const core::Real D = eta * K * discount * barrierTerm * engine_->normalCDF(eta * y1 - eta * v * sqrtT);
            
            price = standardPrice - (A - B_term - C + D) + rebate * discount;
        }
    }
    
    return std::max(price, rebate * discount);
}

core::Real BarrierEngine::knockIn(const core::MarketData& market, core::Real barrier,
                                 core::Real rebate, bool isUp) const {
    const core::Real knockOutPrice = knockOut(market, barrier, 0.0, isUp);
    const core::Real standardPrice = engine_->blackScholes(market).price;
    const core::Real discount = std::exp(-market.riskFreeRate * market.timeToExpiry);
    
    return standardPrice - knockOutPrice + rebate * discount;
}

core::Greeks BarrierEngine::knockOutGreeks(const core::MarketData& market, core::Real barrier,
                                          core::Real rebate, bool isUp) const {
    core::Greeks greeks;
    
    const core::Real basePrice = knockOut(market, barrier, rebate, isUp);
    greeks.price = basePrice;
    
    const core::Real spotBump = market.spot * 0.01;
    const core::Real volBump = 0.01;
    const core::Real timeBump = 1.0 / 365.0;
    const core::Real rateBump = 0.01;
    
    core::MarketData upMarket = market;
    upMarket.spot += spotBump;
    const core::Real upPrice = knockOut(upMarket, barrier, rebate, isUp);
    
    core::MarketData downMarket = market;
    downMarket.spot -= spotBump;
    const core::Real downPrice = knockOut(downMarket, barrier, rebate, isUp);
    
    greeks.delta = (upPrice - downPrice) / (2.0 * spotBump);
    greeks.gamma = (upPrice - 2.0 * basePrice + downPrice) / (spotBump * spotBump);
    
    if (market.timeToExpiry > timeBump) {
        core::MarketData timeMarket = market;
        timeMarket.timeToExpiry -= timeBump;
        greeks.theta = knockOut(timeMarket, barrier, rebate, isUp) - basePrice;
    }
    
    core::MarketData volMarket = market;
    volMarket.volatility += volBump;
    greeks.vega = knockOut(volMarket, barrier, rebate, isUp) - basePrice;
    
    core::MarketData rateMarket = market;
    rateMarket.riskFreeRate += rateBump;
    greeks.rho = knockOut(rateMarket, barrier, rebate, isUp) - basePrice;
    
    greeks.impliedVol = market.volatility;
    
    return greeks;
}

AsianEngine::AsianEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real AsianEngine::geometricAverage(const core::MarketData& market,
                                        const core::Vector& fixingTimes) const {
    if (fixingTimes.empty()) return 0.0;
    
    const core::Real T = market.timeToExpiry;
    const core::Real n = static_cast<core::Real>(fixingTimes.size());
    
    core::Real sumTimes = 0.0;
    core::Real sumProducts = 0.0;
    
    for (core::Size i = 0; i < fixingTimes.size(); ++i) {
        sumTimes += fixingTimes[i];
        for (core::Size j = i; j < fixingTimes.size(); ++j) {
            const core::Real weight = (i == j) ? 1.0 : 2.0;
            sumProducts += weight * std::min(fixingTimes[i], fixingTimes[j]);
        }
    }
    
    const core::Real avgTime = sumTimes / n;
    const core::Real adjVol = market.volatility * std::sqrt(sumProducts / (n * n * T));
    const core::Real adjRate = 0.5 * (market.riskFreeRate + market.dividendYield + 
                                     adjVol * adjVol - market.volatility * market.volatility);
    
    core::MarketData adjMarket = market;
    adjMarket.volatility = adjVol;
    adjMarket.riskFreeRate = adjRate;
    adjMarket.timeToExpiry = avgTime;
    
    return engine_->blackScholes(adjMarket).price;
}

core::Real AsianEngine::arithmeticAverage(const core::MarketData& market,
                                         const core::Vector& fixingTimes,
                                         core::Integer numSimulations) const {
    if (fixingTimes.empty()) return 0.0;
    
    const core::Real dt = market.timeToExpiry / fixingTimes.size();
    const core::Real drift = (market.riskFreeRate - market.dividendYield - 
                             0.5 * market.volatility * market.volatility) * dt;
    const core::Real diffusion = market.volatility * std::sqrt(dt);
    
    auto& rng = engine_->getRNG();
    std::normal_distribution<core::Real> normal(0.0, 1.0);
    
    core::Real totalPayoff = 0.0;
    
    for (core::Integer sim = 0; sim < numSimulations; ++sim) {
        core::Real spot = market.spot;
        core::Real sumPrices = 0.0;
        
        for (core::Size i = 0; i < fixingTimes.size(); ++i) {
            const core::Real z = normal(rng);
            spot *= std::exp(drift + diffusion * z);
            sumPrices += spot;
        }
        
        const core::Real avgPrice = sumPrices / fixingTimes.size();
        core::Real payoff = 0.0;
        
        if (market.isCall) {
            payoff = std::max(0.0, avgPrice - market.strike);
        } else {
            payoff = std::max(0.0, market.strike - avgPrice);
        }
        
        totalPayoff += payoff;
    }
    
    const core::Real avgPayoff = totalPayoff / numSimulations;
    return avgPayoff * std::exp(-market.riskFreeRate * market.timeToExpiry);
}

core::Greeks AsianEngine::geometricAverageGreeks(const core::MarketData& market,
                                                const core::Vector& fixingTimes) const {
    core::Greeks greeks;
    
    const core::Real basePrice = geometricAverage(market, fixingTimes);
    greeks.price = basePrice;
    
    const core::Real spotBump = market.spot * 0.01;
    const core::Real volBump = 0.01;
    const core::Real timeBump = 1.0 / 365.0;
    const core::Real rateBump = 0.01;
    
    core::MarketData upMarket = market;
    upMarket.spot += spotBump;
    const core::Real upPrice = geometricAverage(upMarket, fixingTimes);
    
    core::MarketData downMarket = market;
    downMarket.spot -= spotBump;
    const core::Real downPrice = geometricAverage(downMarket, fixingTimes);
    
    greeks.delta = (upPrice - downPrice) / (2.0 * spotBump);
    greeks.gamma = (upPrice - 2.0 * basePrice + downPrice) / (spotBump * spotBump);
    
    if (market.timeToExpiry > timeBump) {
        core::MarketData timeMarket = market;
        timeMarket.timeToExpiry -= timeBump;
        greeks.theta = geometricAverage(timeMarket, fixingTimes) - basePrice;
    }
    
    core::MarketData volMarket = market;
    volMarket.volatility += volBump;
    greeks.vega = geometricAverage(volMarket, fixingTimes) - basePrice;
    
    core::MarketData rateMarket = market;
    rateMarket.riskFreeRate += rateBump;
    greeks.rho = geometricAverage(rateMarket, fixingTimes) - basePrice;
    
    greeks.impliedVol = market.volatility;
    
    return greeks;
}

LookbackEngine::LookbackEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real LookbackEngine::fixedStrike(const core::MarketData& market, core::Real extremeValue,
                                      bool isFloatingMax) const {
    const core::Real S = market.spot;
    const core::Real K = market.strike;
    const core::Real T = market.timeToExpiry;
    const core::Real r = market.riskFreeRate;
    const core::Real v = market.volatility;
    const core::Real q = market.dividendYield;
    
    const core::Real sqrtT = std::sqrt(T);
    const core::Real M = extremeValue;
    
    const core::Real a1 = (std::log(S / K) + (r - q + 0.5 * v * v) * T) / (v * sqrtT);
    const core::Real a2 = a1 - v * sqrtT;
    const core::Real a3 = (std::log(S / M) + (r - q + 0.5 * v * v) * T) / (v * sqrtT);
    
    const core::Real lambda = (-r + q + 0.5 * v * v) / (v * v);
    const core::Real discount = std::exp(-r * T);
    const core::Real forwardFactor = std::exp(-q * T);
    
    core::Real price = 0.0;
    
    if (market.isCall && isFloatingMax) {
        if (M > K) {
            price = S * forwardFactor * engine_->normalCDF(a1) - 
                   K * discount * engine_->normalCDF(a2) +
                   S * forwardFactor * (v * v) / (2.0 * (r - q)) * 
                   (std::pow(S / M, -2.0 * lambda) * engine_->normalCDF(-a3) - 
                    std::exp((r - q) * T) * engine_->normalCDF(-a1));
        } else {
            price = (M - K) * discount;
        }
    } else if (!market.isCall && !isFloatingMax) {
        if (M < K) {
            price = K * discount * engine_->normalCDF(-a2) - 
                   S * forwardFactor * engine_->normalCDF(-a1) +
                   S * forwardFactor * (v * v) / (2.0 * (r - q)) * 
                   (std::pow(S / M, -2.0 * lambda) * engine_->normalCDF(a3) - 
                    std::exp((r - q) * T) * engine_->normalCDF(a1));
        } else {
            price = (K - M) * discount;
        }
    }
    
    return std::max(price, 0.0);
}

core::Real LookbackEngine::floatingStrike(const core::MarketData& market, core::Real extremeValue,
                                         bool isCall) const {
    const core::Real S = market.spot;
    const core::Real T = market.timeToExpiry;
    const core::Real r = market.riskFreeRate;
    const core::Real v = market.volatility;
    const core::Real q = market.dividendYield;
    const core::Real M = extremeValue;
    
    const core::Real sqrtT = std::sqrt(T);
    const core::Real lambda = (-r + q + 0.5 * v * v) / (v * v);
    const core::Real discount = std::exp(-r * T);
    const core::Real forwardFactor = std::exp(-q * T);
    
    const core::Real b1 = (std::log(S / M) + (r - q + 0.5 * v * v) * T) / (v * sqrtT);
    const core::Real b2 = b1 - v * sqrtT;
    const core::Real b3 = (std::log(S / M) + (r - q - 0.5 * v * v) * T) / (v * sqrtT);
    
    core::Real price = 0.0;
    
    if (isCall) {
        price = S * forwardFactor * engine_->normalCDF(b1) - 
               M * discount * engine_->normalCDF(b2) +
               S * forwardFactor * (v * v) / (2.0 * (r - q)) * 
               (-std::pow(S / M, -2.0 * lambda) * engine_->normalCDF(-b3) + 
                std::exp((r - q) * T) * engine_->normalCDF(b1));
    } else {
        price = M * discount * engine_->normalCDF(-b2) - 
               S * forwardFactor * engine_->normalCDF(-b1) +
               S * forwardFactor * (v * v) / (2.0 * (r - q)) * 
               (std::pow(S / M, -2.0 * lambda) * engine_->normalCDF(b3) - 
                std::exp((r - q) * T) * engine_->normalCDF(-b1));
    }
    
    return std::max(price, 0.0);
}

CompoundEngine::CompoundEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real CompoundEngine::callOnCall(const core::MarketData& outerOption,
                                     const core::MarketData& innerOption) const {
    const core::Real S = outerOption.spot;
    const core::Real K1 = outerOption.strike;
    const core::Real K2 = innerOption.strike;
    const core::Real T1 = outerOption.timeToExpiry;
    const core::Real T2 = innerOption.timeToExpiry;
    const core::Real r = outerOption.riskFreeRate;
    const core::Real v = outerOption.volatility;
    
    if (T1 >= T2) return 0.0;
    
    const core::Real sqrtT1 = std::sqrt(T1);
    const core::Real sqrtT2 = std::sqrt(T2);
    const core::Real rho = std::sqrt(T1 / T2);
    
    const core::Real critS = criticalPrice(innerOption);
    
    const core::Real y1 = (std::log(S / critS) + (r + 0.5 * v * v) * T1) / (v * sqrtT1);
    const core::Real y2 = (std::log(S / K2) + (r + 0.5 * v * v) * T2) / (v * sqrtT2);
    
    const core::Real discount1 = std::exp(-r * T1);
    const core::Real discount2 = std::exp(-r * T2);
    
    const core::Real M1 = bivariateCDF(y1, y2, rho);
    const core::Real M2 = bivariateCDF(y1 - v * sqrtT1, y2 - v * sqrtT2, rho);
    
    return S * M1 - K2 * discount2 * M2 - K1 * discount1 * engine_->normalCDF(y1);
}

core::Real CompoundEngine::bivariateCDF(core::Real x, core::Real y, core::Real rho) const {
    const core::Integer n = 20;
    const core::Real h = 0.5;
    
    core::Real sum = 0.0;
    for (core::Integer i = 0; i < n; ++i) {
        const core::Real xi = -3.0 + 6.0 * i / (n - 1);
        const core::Real weight = (i == 0 || i == n - 1) ? 1.0 : ((i % 2 == 0) ? 2.0 : 4.0);
        
        const core::Real integrand = std::exp(-0.5 * (xi * xi - 2.0 * rho * xi * y + y * y) / (1.0 - rho * rho));
        sum += weight * integrand;
    }
    
    return engine_->normalCDF(x) * engine_->normalCDF(y) + 
           rho / (2.0 * core::PI * std::sqrt(1.0 - rho * rho)) * sum * h / 3.0;
}

core::Real CompoundEngine::criticalPrice(const core::MarketData& innerOption) const {
    auto objectiveFunc = [&](core::Real S) -> core::Real {
        core::MarketData tempMarket = innerOption;
        tempMarket.spot = S;
        return engine_->blackScholes(tempMarket).price - innerOption.strike;
    };
    
    return core::NumericalMethods::brent(objectiveFunc, innerOption.spot * 0.5, innerOption.spot * 2.0);
}

ExchangeEngine::ExchangeEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real ExchangeEngine::exchangeOption(core::Real S1, core::Real S2, core::Real T,
                                         core::Real r, core::Real v1, core::Real v2,
                                         core::Real rho, core::Real q1, core::Real q2) const {
    const core::Real sigma = std::sqrt(v1 * v1 + v2 * v2 - 2.0 * rho * v1 * v2);
    const core::Real sqrtT = std::sqrt(T);
    
    const core::Real d1 = (std::log(S1 / S2) + (r - q1 + q2 + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const core::Real d2 = d1 - sigma * sqrtT;
    
    const core::Real discount1 = std::exp(-q1 * T);
    const core::Real discount2 = std::exp(-q2 * T);
    
    return S1 * discount1 * engine_->normalCDF(d1) - S2 * discount2 * engine_->normalCDF(d2);
}

core::Real ExchangeEngine::basketCall(const core::Vector& spots, const core::Vector& weights,
                                     core::Real strike, core::Real T, core::Real r,
                                     const core::Matrix& correlations, const core::Vector& vols,
                                     const core::Vector& dividends) const {
    if (spots.size() != weights.size() || spots.size() != vols.size()) {
        return 0.0;
    }
    
    const core::Real basketVol = basketVolatility(weights, correlations, vols);
    
    core::Real basketSpot = 0.0;
    for (core::Size i = 0; i < spots.size(); ++i) {
        const core::Real q = i < dividends.size() ? dividends[i] : 0.0;
        basketSpot += weights[i] * spots[i] * std::exp(-q * T);
    }
    
    const core::Real sqrtT = std::sqrt(T);
    const core::Real d1 = (std::log(basketSpot / strike) + (r + 0.5 * basketVol * basketVol) * T) / 
                          (basketVol * sqrtT);
    const core::Real d2 = d1 - basketVol * sqrtT;
    
    const core::Real discount = std::exp(-r * T);
    
    return basketSpot * engine_->normalCDF(d1) - strike * discount * engine_->normalCDF(d2);
}

core::Real ExchangeEngine::basketVolatility(const core::Vector& weights, const core::Matrix& correlations,
                                           const core::Vector& vols) const {
    core::Real variance = 0.0;
    
    for (core::Size i = 0; i < weights.size(); ++i) {
        for (core::Size j = 0; j < weights.size(); ++j) {
            if (i < correlations.size() && j < correlations[i].size()) {
                variance += weights[i] * weights[j] * vols[i] * vols[j] * correlations[i][j];
            }
        }
    }
    
    return std::sqrt(variance);
}

}