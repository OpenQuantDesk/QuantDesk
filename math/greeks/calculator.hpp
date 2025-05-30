#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include "../core/threading.hpp"
#include <memory>
#include <functional>

namespace math::greeks {

class GreeksCalculator {
private:
    std::shared_ptr<core::MathEngine> engine_;
    mutable core::GreeksPool greeksPool_;
    
public:
    explicit GreeksCalculator(std::shared_ptr<core::MathEngine> engine);
    
    core::Greeks analytical(const core::MarketData& market) const;
    core::ExtendedGreeks analyticalExtended(const core::MarketData& market, 
                                           core::Real quantity = 1.0) const;
    
    core::Greeks numerical(const core::MarketData& market, 
                          core::Real bumpSize = 0.01) const;
    core::ExtendedGreeks numericalExtended(const core::MarketData& market,
                                          core::Real bumpSize = 0.01,
                                          core::Real quantity = 1.0) const;
    
    void analyticalBatch(const core::MarketDataVector& markets,
                        core::GreeksVector& results) const;
    void extendedBatch(const core::MarketDataVector& markets,
                      const core::Vector& quantities,
                      core::ExtendedGreeksVector& results) const;
    
    core::Greeks aggregate(const core::GreeksVector& positions,
                          const core::Vector& quantities) const;
    core::ExtendedGreeks aggregateExtended(const core::ExtendedGreeksVector& positions) const;
    
    core::Real calculatePinRisk(core::Real spot, core::Real strike, 
                               core::Real timeToExpiry, core::Real quantity = 1.0) const;
    
    core::Real calculateDollarDelta(const core::Greeks& greeks, core::Real spot,
                                   core::Real quantity, core::Real multiplier = 100.0) const;
    
    core::Real calculateDollarGamma(const core::Greeks& greeks, core::Real spot,
                                   core::Real quantity, core::Real multiplier = 100.0) const;
    
    void calculateRiskLadder(const core::MarketDataVector& markets,
                            const core::Vector& quantities,
                            const core::Vector& strikeLadder,
                            core::Matrix& deltaLadder,
                            core::Matrix& gammaLadder) const;
    
    struct StressTestParams {
        core::Vector spotShifts;
        core::Vector volShifts;
        core::Vector rateShifts;
        core::Vector timeDecays;
    };
    
    struct StressTestResults {
        core::Matrix spotStress;
        core::Matrix volStress;
        core::Matrix rateStress;
        core::Matrix timeStress;
        core::Vector worstCase;
        core::Vector bestCase;
    };
    
    StressTestResults stressTest(const core::MarketDataVector& markets,
                                const core::Vector& quantities,
                                const StressTestParams& params) const;
    
private:
    core::ExtendedGreeks calculateSecondOrder(const core::MarketData& market,
                                             core::Real quantity,
                                             core::Real bumpSize) const;
    
    core::Real calculateVolga(const core::MarketData& market, core::Real bumpSize) const;
    core::Real calculateVanna(const core::MarketData& market, core::Real bumpSize) const;
    core::Real calculateCharm(const core::MarketData& market, core::Real bumpSize) const;
    core::Real calculateSpeed(const core::MarketData& market, core::Real bumpSize) const;
    core::Real calculateZomma(const core::MarketData& market, core::Real bumpSize) const;
    core::Real calculateColor(const core::MarketData& market, core::Real bumpSize) const;
};

class RiskEngine {
private:
    std::shared_ptr<GreeksCalculator> calculator_;
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit RiskEngine(std::shared_ptr<GreeksCalculator> calculator);
    
    struct PortfolioRisk {
        core::Real totalDelta;
        core::Real totalGamma;
        core::Real totalVega;
        core::Real totalTheta;
        core::Real totalRho;
        core::Real var95;
        core::Real var99;
        core::Real expectedShortfall95;
        core::Real expectedShortfall99;
        core::Real maxDrawdown;
        core::Vector concentrationRisk;
        core::Matrix correlationMatrix;
    };
    
    PortfolioRisk analyzePortfolio(const core::MarketDataVector& markets,
                                  const core::Vector& quantities,
                                  const core::Vector& historicalReturns) const;
    
    core::Real calculateVaR(const core::GreeksVector& greeks,
                           const core::Vector& quantities,
                           const core::Matrix& covarianceMatrix,
                           core::Real confidence = 0.95) const;
    
    core::Real calculateExpectedShortfall(const core::GreeksVector& greeks,
                                         const core::Vector& quantities,
                                         const core::Matrix& covarianceMatrix,
                                         core::Real confidence = 0.95) const;
    
    core::Vector calculateMarginalVaR(const core::GreeksVector& greeks,
                                     const core::Vector& quantities,
                                     const core::Matrix& covarianceMatrix) const;
    
    core::Vector calculateComponentVaR(const core::GreeksVector& greeks,
                                      const core::Vector& quantities,
                                      const core::Matrix& covarianceMatrix) const;
    
    core::Matrix calculateCorrelationMatrix(const core::GreeksVector& greeks) const;
    
    struct HedgeRecommendation {
        core::Vector deltaHedge;
        core::Vector gammaHedge;
        core::Vector vegaHedge;
        core::Real hedgingCost;
        core::Real residualRisk;
        core::Vector optimalQuantities;
    };
    
    HedgeRecommendation recommendHedge(const core::GreeksVector& portfolioGreeks,
                                      const core::GreeksVector& hedgeInstruments,
                                      const core::Vector& hedgeCosts) const;
    
private:
    core::Matrix buildCovarianceMatrix(const core::Vector& returns, core::Size lag = 252) const;
    core::Real portfolioVariance(const core::Vector& weights, const core::Matrix& covariance) const;
    core::Vector optimizeHedge(const core::GreeksVector& portfolio,
                              const core::GreeksVector& instruments,
                              const core::Vector& costs) const;
};

class SensitivityAnalyzer {
private:
    std::shared_ptr<GreeksCalculator> calculator_;
    
public:
    explicit SensitivityAnalyzer(std::shared_ptr<GreeksCalculator> calculator);
    
    struct SensitivityProfile {
        core::Vector spotSensitivity;
        core::Vector volSensitivity;
        core::Vector rateSensitivity;
        core::Vector timeSensitivity;
        core::Matrix crossSensitivities;
        core::Vector keyRiskFactors;
        core::Real totalSensitivity;
    };
    
    SensitivityProfile analyzeSensitivities(const core::MarketDataVector& markets,
                                           const core::Vector& quantities,
                                           const core::Vector& bumpSizes) const;
    
    core::Matrix calculateHessian(const core::MarketData& market,
                                 const std::vector<std::function<core::Real(const core::MarketData&)>>& functions,
                                 core::Real bumpSize = 0.01) const;
    
    core::Vector calculateJacobian(const core::MarketData& market,
                                  const std::function<core::Vector(const core::MarketData&)>& function,
                                  core::Real bumpSize = 0.01) const;
    
    struct ScenarioAnalysis {
        core::Vector scenarios;
        core::Vector probabilities;
        core::Vector expectedValues;
        core::Vector worstCase;
        core::Vector bestCase;
        core::Real expectedReturn;
        core::Real variance;
        core::Real skewness;
        core::Real kurtosis;
    };
    
    ScenarioAnalysis runScenarioAnalysis(const core::MarketDataVector& markets,
                                        const core::Vector& quantities,
                                        const core::Matrix& scenarioShocks,
                                        const core::Vector& scenarioProbabilities) const;
    
private:
    core::Real calculateNumericalDerivative(const std::function<core::Real(core::Real)>& func,
                                           core::Real x, core::Real h, core::Integer order = 1) const;
    
    core::Matrix calculateCrossDerivatives(const core::MarketData& market,
                                          const std::function<core::Real(const core::MarketData&)>& func,
                                          core::Real bumpSize) const;
};

inline GreeksCalculator::GreeksCalculator(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

inline core::Greeks GreeksCalculator::analytical(const core::MarketData& market) const {
    return engine_->blackScholes(market);
}

inline core::ExtendedGreeks GreeksCalculator::analyticalExtended(const core::MarketData& market,
                                                                 core::Real quantity) const {
    const core::Greeks basic = analytical(market);
    core::ExtendedGreeks extended;
    
    extended.price = basic.price * quantity;
    extended.delta = basic.delta * quantity;
    extended.gamma = basic.gamma * quantity;
    extended.theta = basic.theta * quantity;
    extended.vega = basic.vega * quantity;
    extended.rho = basic.rho * quantity;
    extended.impliedVol = basic.impliedVol;
    
    const core::Real bumpSize = 0.01;
    extended.volga = calculateVolga(market, bumpSize) * quantity;
    extended.vanna = calculateVanna(market, bumpSize) * quantity;
    extended.charm = calculateCharm(market, bumpSize) * quantity;
    extended.speed = calculateSpeed(market, bumpSize) * quantity;
    extended.zomma = calculateZomma(market, bumpSize) * quantity;
    extended.color = calculateColor(market, bumpSize) * quantity;
    
    extended.dollarDelta = calculateDollarDelta(basic, market.spot, quantity);
    extended.dollarGamma = calculateDollarGamma(basic, market.spot, quantity);
    extended.pinRisk = calculatePinRisk(market.spot, market.strike, market.timeToExpiry, quantity);
    
    return extended;
}

inline core::Greeks GreeksCalculator::numerical(const core::MarketData& market,
                                                core::Real bumpSize) const {
    core::Greeks greeks;
    
    const core::Real basePrice = engine_->blackScholes(market).price;
    greeks.price = basePrice;
    
    core::MarketData upMarket = market;
    upMarket.spot += bumpSize;
    const core::Real upPrice = engine_->blackScholes(upMarket).price;
    
    core::MarketData downMarket = market;
    downMarket.spot -= bumpSize;
    const core::Real downPrice = engine_->blackScholes(downMarket).price;
    
    greeks.delta = (upPrice - downPrice) / (2.0 * bumpSize);
    greeks.gamma = (upPrice - 2.0 * basePrice + downPrice) / (bumpSize * bumpSize);
    
    const core::Real timeBump = 1.0 / 365.0;
    if (market.timeToExpiry > timeBump) {
        core::MarketData timeMarket = market;
        timeMarket.timeToExpiry -= timeBump;
        greeks.theta = engine_->blackScholes(timeMarket).price - basePrice;
    }
    
    const core::Real volBump = 0.01;
    core::MarketData volMarket = market;
    volMarket.volatility += volBump;
    greeks.vega = engine_->blackScholes(volMarket).price - basePrice;
    
    const core::Real rateBump = 0.01;
    core::MarketData rateMarket = market;
    rateMarket.riskFreeRate += rateBump;
    greeks.rho = engine_->blackScholes(rateMarket).price - basePrice;
    
    greeks.impliedVol = market.volatility;
    
    return greeks;
}

inline void GreeksCalculator::analyticalBatch(const core::MarketDataVector& markets,
                                             core::GreeksVector& results) const {
    engine_->blackScholesBatch(markets, results);
}

inline core::Greeks GreeksCalculator::aggregate(const core::GreeksVector& positions,
                                               const core::Vector& quantities) const {
    return engine_->aggregateGreeks(positions, quantities);
}

inline core::ExtendedGreeks GreeksCalculator::aggregateExtended(const core::ExtendedGreeksVector& positions) const {
    core::ExtendedGreeks total;
    
    for (const auto& pos : positions) {
        total.price += pos.price;
        total.delta += pos.delta;
        total.gamma += pos.gamma;
        total.theta += pos.theta;
        total.vega += pos.vega;
        total.rho += pos.rho;
        total.volga += pos.volga;
        total.vanna += pos.vanna;
        total.charm += pos.charm;
        total.speed += pos.speed;
        total.zomma += pos.zomma;
        total.color += pos.color;
        total.dollarDelta += pos.dollarDelta;
        total.dollarGamma += pos.dollarGamma;
        total.pinRisk += pos.pinRisk;
    }
    
    return total;
}

inline core::Real GreeksCalculator::calculatePinRisk(core::Real spot, core::Real strike,
                                                     core::Real timeToExpiry, core::Real quantity) const {
    if (timeToExpiry > 7.0 / 365.0) return 0.0;
    
    const core::Real distance = std::abs(spot - strike) / spot;
    const core::Real timeWeight = 1.0 - timeToExpiry / (7.0 / 365.0);
    const core::Real pinDistance = 0.02;
    
    if (distance < pinDistance) {
        return (pinDistance - distance) / pinDistance * timeWeight * std::abs(quantity);
    }
    
    return 0.0;
}

inline core::Real GreeksCalculator::calculateDollarDelta(const core::Greeks& greeks, core::Real spot,
                                                        core::Real quantity, core::Real multiplier) const {
    return greeks.delta * spot * quantity * multiplier * 0.01;
}

inline core::Real GreeksCalculator::calculateDollarGamma(const core::Greeks& greeks, core::Real spot,
                                                        core::Real quantity, core::Real multiplier) const {
    return greeks.gamma * spot * spot * quantity * multiplier * 0.01 * 0.01;
}

inline core::Real GreeksCalculator::calculateVolga(const core::MarketData& market, core::Real bumpSize) const {
    core::MarketData volUp = market;
    volUp.volatility += bumpSize;
    const core::Real vegaUp = engine_->blackScholes(volUp).vega;
    
    core::MarketData volDown = market;
    volDown.volatility -= bumpSize;
    const core::Real vegaDown = engine_->blackScholes(volDown).vega;
    
    return (vegaUp - vegaDown) / (2.0 * bumpSize);
}

inline core::Real GreeksCalculator::calculateVanna(const core::MarketData& market, core::Real bumpSize) const {
    const core::Real spotBump = market.spot * 0.01;
    const core::Real volBump = 0.01;
    
    core::MarketData upUp = market;
    upUp.spot += spotBump;
    upUp.volatility += volBump;
    const core::Real deltaUpUp = engine_->blackScholes(upUp).delta;
    
    core::MarketData upDown = market;
    upDown.spot += spotBump;
    upDown.volatility -= volBump;
    const core::Real deltaUpDown = engine_->blackScholes(upDown).delta;
    
    core::MarketData downUp = market;
    downUp.spot -= spotBump;
    downUp.volatility += volBump;
    const core::Real deltaDownUp = engine_->blackScholes(downUp).delta;
    
    core::MarketData downDown = market;
    downDown.spot -= spotBump;
    downDown.volatility -= volBump;
    const core::Real deltaDownDown = engine_->blackScholes(downDown).delta;
    
    return (deltaUpUp - deltaUpDown - deltaDownUp + deltaDownDown) / (4.0 * spotBump * volBump);
}

inline RiskEngine::RiskEngine(std::shared_ptr<GreeksCalculator> calculator)
    : calculator_(std::move(calculator)) {}

inline SensitivityAnalyzer::SensitivityAnalyzer(std::shared_ptr<GreeksCalculator> calculator)
    : calculator_(std::move(calculator)) {}

}