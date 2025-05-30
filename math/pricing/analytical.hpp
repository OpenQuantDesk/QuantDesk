#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include <memory>

namespace math::pricing {

class BlackScholesEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit BlackScholesEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Greeks price(const core::MarketData& market) const;
    core::Real priceOnly(const core::MarketData& market) const;
    
    void priceBatch(const core::MarketDataVector& markets, core::GreeksVector& results) const;
    void priceOnlyBatch(const core::MarketDataVector& markets, core::Vector& prices) const;
    
    core::Real impliedVolatility(core::Real marketPrice, const core::MarketData& market) const;
    void impliedVolatilityBatch(const core::Vector& marketPrices, 
                               const core::MarketDataVector& markets,
                               core::Vector& impliedVols) const;
    
    core::Greeks digitalOption(const core::MarketData& market, core::Real cashPayoff = 1.0) const;
    core::Greeks assetOrNothingOption(const core::MarketData& market) const;
    
private:
    core::Real calculateD1(core::Real S, core::Real K, core::Real T, 
                          core::Real r, core::Real v, core::Real q = 0.0) const;
    core::Real calculateD2(core::Real d1, core::Real v, core::Real T) const;
};

class BarrierEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit BarrierEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real knockOut(const core::MarketData& market, core::Real barrier, 
                       core::Real rebate = 0.0, bool isUp = true) const;
    
    core::Real knockIn(const core::MarketData& market, core::Real barrier, 
                      core::Real rebate = 0.0, bool isUp = true) const;
    
    core::Greeks knockOutGreeks(const core::MarketData& market, core::Real barrier,
                               core::Real rebate = 0.0, bool isUp = true) const;
    
private:
    core::Real barrierFactor(core::Real S, core::Real B, core::Real T, 
                            core::Real r, core::Real v, core::Real q = 0.0) const;
};

class AsianEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit AsianEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real geometricAverage(const core::MarketData& market, 
                               const core::Vector& fixingTimes) const;
    
    core::Real arithmeticAverage(const core::MarketData& market,
                                const core::Vector& fixingTimes,
                                core::Integer numSimulations = 100000) const;
    
    core::Greeks geometricAverageGreeks(const core::MarketData& market,
                                       const core::Vector& fixingTimes) const;
    
private:
    core::Real adjustedVolatility(core::Real vol, const core::Vector& fixingTimes) const;
    core::Real adjustedDrift(core::Real r, core::Real vol, 
                            const core::Vector& fixingTimes) const;
};

class LookbackEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit LookbackEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real fixedStrike(const core::MarketData& market, core::Real extremeValue,
                          bool isFloatingMax = true) const;
    
    core::Real floatingStrike(const core::MarketData& market, core::Real extremeValue,
                             bool isCall = true) const;
    
    core::Greeks fixedStrikeGreeks(const core::MarketData& market, core::Real extremeValue,
                                  bool isFloatingMax = true) const;
    
private:
    core::Real d1Lookback(core::Real S, core::Real E, core::Real T, 
                         core::Real r, core::Real v, core::Real q = 0.0) const;
    core::Real d2Lookback(core::Real S, core::Real E, core::Real T, 
                         core::Real r, core::Real v, core::Real q = 0.0) const;
};

class CompoundEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit CompoundEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real callOnCall(const core::MarketData& outerOption, 
                         const core::MarketData& innerOption) const;
    
    core::Real putOnPut(const core::MarketData& outerOption,
                       const core::MarketData& innerOption) const;
    
    core::Real callOnPut(const core::MarketData& outerOption,
                        const core::MarketData& innerOption) const;
    
    core::Real putOnCall(const core::MarketData& outerOption,
                        const core::MarketData& innerOption) const;
    
private:
    core::Real bivariateCDF(core::Real x, core::Real y, core::Real rho) const;
    core::Real criticalPrice(const core::MarketData& innerOption) const;
};

class ExchangeEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit ExchangeEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real exchangeOption(core::Real S1, core::Real S2, core::Real T,
                             core::Real r, core::Real v1, core::Real v2, 
                             core::Real rho, core::Real q1 = 0.0, core::Real q2 = 0.0) const;
    
    core::Real basketCall(const core::Vector& spots, const core::Vector& weights,
                         core::Real strike, core::Real T, core::Real r,
                         const core::Matrix& correlations, const core::Vector& vols,
                         const core::Vector& dividends = {}) const;
    
    core::Real basketPut(const core::Vector& spots, const core::Vector& weights,
                        core::Real strike, core::Real T, core::Real r,
                        const core::Matrix& correlations, const core::Vector& vols,
                        const core::Vector& dividends = {}) const;
    
    core::Real rainbowOption(const core::Vector& spots, core::Real T, core::Real r,
                            const core::Matrix& correlations, const core::Vector& vols,
                            bool isCall = true, bool isMax = true) const;
    
private:
    core::Real basketVolatility(const core::Vector& weights, const core::Matrix& correlations,
                               const core::Vector& vols) const;
    core::Real effectiveStrike(const core::Vector& spots, const core::Vector& weights,
                              core::Real strike) const;
};

class SpreadEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit SpreadEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real spreadOption(core::Real S1, core::Real S2, core::Real K,
                           core::Real T, core::Real r, core::Real v1, core::Real v2,
                           core::Real rho, bool isCall = true) const;
    
    core::Real ratioSpread(core::Real S1, core::Real S2, core::Real ratio,
                          core::Real T, core::Real r, core::Real v1, core::Real v2,
                          core::Real rho) const;
    
    core::Real crackSpread(core::Real crude, core::Real gasoline, core::Real heating,
                          const core::Vector& ratios, core::Real T, core::Real r,
                          const core::Vector& vols, const core::Matrix& correlations) const;
    
private:
    core::Real spreadVolatility(core::Real v1, core::Real v2, core::Real rho) const;
    core::Real approximateSpreadPrice(core::Real F1, core::Real F2, core::Real K,
                                     core::Real T, core::Real spreadVol, bool isCall) const;
};

class PowerEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit PowerEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real powerOption(const core::MarketData& market, core::Real power) const;
    core::Real logOption(const core::MarketData& market) const;
    core::Real expOption(const core::MarketData& market) const;
    
    core::Greeks powerOptionGreeks(const core::MarketData& market, core::Real power) const;
    
private:
    core::Real adjustedForward(core::Real S, core::Real r, core::Real q, 
                              core::Real T, core::Real power) const;
    core::Real adjustedVolatility(core::Real vol, core::Real power) const;
};

class QuantoEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit QuantoEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real quantoOption(const core::MarketData& market, core::Real fxVol,
                           core::Real correlation, core::Real foreignRate) const;
    
    core::Greeks quantoGreeks(const core::MarketData& market, core::Real fxVol,
                             core::Real correlation, core::Real foreignRate) const;
    
private:
    core::Real adjustedDrift(core::Real domesticRate, core::Real foreignRate,
                            core::Real assetVol, core::Real fxVol, core::Real correlation) const;
};

inline BlackScholesEngine::BlackScholesEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

inline core::Greeks BlackScholesEngine::price(const core::MarketData& market) const {
    return engine_->blackScholes(market);
}

inline core::Real BlackScholesEngine::priceOnly(const core::MarketData& market) const {
    return engine_->blackScholes(market).price;
}

inline void BlackScholesEngine::priceBatch(const core::MarketDataVector& markets, 
                                          core::GreeksVector& results) const {
    engine_->blackScholesBatch(markets, results);
}

inline core::Real BlackScholesEngine::impliedVolatility(core::Real marketPrice, 
                                                        const core::MarketData& market) const {
    return engine_->impliedVolatility(marketPrice, market);
}

inline core::Real BlackScholesEngine::calculateD1(core::Real S, core::Real K, core::Real T,
                                                  core::Real r, core::Real v, core::Real q) const {
    return (std::log(S / K) + (r - q + 0.5 * v * v) * T) / (v * std::sqrt(T));
}

inline core::Real BlackScholesEngine::calculateD2(core::Real d1, core::Real v, core::Real T) const {
    return d1 - v * std::sqrt(T);
}

inline core::Greeks BlackScholesEngine::digitalOption(const core::MarketData& market, 
                                                      core::Real cashPayoff) const {
    const core::Real d2 = calculateD2(calculateD1(market.spot, market.strike, market.timeToExpiry,
                                                  market.riskFreeRate, market.volatility), 
                                     market.volatility, market.timeToExpiry);
    
    const core::Real discountFactor = std::exp(-market.riskFreeRate * market.timeToExpiry);
    const core::Real nd2 = engine_->normalCDF(market.isCall ? d2 : -d2);
    
    core::Greeks greeks;
    greeks.price = cashPayoff * discountFactor * nd2;
    
    const core::Real pdf = engine_->normalPDF(d2);
    const core::Real sqrtT = std::sqrt(market.timeToExpiry);
    
    greeks.delta = (market.isCall ? 1.0 : -1.0) * cashPayoff * discountFactor * pdf / 
                   (market.spot * market.volatility * sqrtT);
    greeks.gamma = -greeks.delta * (d2 / (market.volatility * sqrtT) + 1.0 / (market.volatility * sqrtT)) / market.spot;
    greeks.theta = -market.riskFreeRate * greeks.price - 
                   (market.isCall ? 1.0 : -1.0) * cashPayoff * discountFactor * pdf / 
                   (2.0 * market.volatility * sqrtT);
    greeks.vega = -(market.isCall ? 1.0 : -1.0) * cashPayoff * discountFactor * pdf * d2 / market.volatility;
    greeks.rho = (market.isCall ? 1.0 : -1.0) * cashPayoff * market.timeToExpiry * discountFactor * pdf / 
                 (market.volatility * sqrtT);
    
    return greeks;
}

}