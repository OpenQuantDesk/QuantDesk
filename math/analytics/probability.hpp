#pragma once

#include "../core/types.hpp"
#include "../metrics/probability.hpp"
#include <memory>

namespace math {

class MathEngine;

class ProbabilityEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit ProbabilityEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::ProbabilityMetrics analyzeOption(const core::MarketData& market, core::Real premium) const;
    
    core::ProbabilityMetrics analyzeStrategy(const std::vector<core::MarketData>& legs,
                                           const core::Vector& quantities,
                                           const core::Vector& premiums) const;
    
    core::ProbabilityMetrics calculateHistoricalMetrics(const core::Vector& prices,
                                                       const core::MarketData& market) const;
    
    core::Real calculateTouchProbability(core::Real spot, core::Real barrier,
                                        core::Real vol, core::Real time,
                                        bool isUp = true) const;
    
    core::Real calculateBreakevenProbability(const core::Vector& strikes,
                                            const core::Vector& quantities,
                                            const core::Vector& premiums,
                                            const core::MarketData& market) const;
    
    core::Real calculateMaxPain(const core::Vector& strikes,
                               const core::Vector& callOI,
                               const core::Vector& putOI) const;
    
    core::Real calculateExpectedMove(const core::MarketData& market) const;
    
    struct MonteCarloParams {
        core::Integer numSimulations = 100000;
        core::Integer timeSteps = 252;
        bool useAntitheticVariates = true;
        bool useControlVariates = false;
        uint32_t seed = 12345;
    };
    
    core::ProbabilityMetrics monteCarloAnalysis(const core::MarketData& market,
                                              const std::function<core::Real(core::Real)>& payoffFunc,
                                              const MonteCarloParams& params = {}) const;
    
    struct ScenarioAnalysis {
        core::Vector spotScenarios;
        core::Vector probabilities;
        core::Vector payoffs;
        core::Real expectedPayoff;
        core::Real variance;
        core::Real probabilityProfit;
        core::Real var95;
        core::Real var99;
        core::Real expectedShortfall95;
        core::Real expectedShortfall99;
    };
    
    ScenarioAnalysis runScenarioAnalysis(const core::Vector& spotScenarios,
                                        const core::Vector& probabilities,
                                        const std::function<core::Real(core::Real)>& payoffFunc) const;
    
private:
    core::Real calculateKellyCriterion(core::Real winRate, core::Real avgWin, core::Real avgLoss) const;
    
    core::Real calculateProbabilityITM(const core::MarketData& market) const;
    
    core::Vector generatePriceScenarios(const core::MarketData& market,
                                       const MonteCarloParams& params) const;
    
    core::Real bivariateNormalCDF(core::Real x, core::Real y, core::Real rho) const;
};

inline ProbabilityEngine::ProbabilityEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

inline core::Real ProbabilityEngine::calculateExpectedMove(const core::MarketData& market) const {
    return market.spot * market.volatility * std::sqrt(market.timeToExpiry);
}

inline core::Real ProbabilityEngine::calculateProbabilityITM(const core::MarketData& market) const {
    if (market.timeToExpiry <= 0.0 || market.volatility <= 0.0) return 0.0;
    
    const core::Real sqrtT = std::sqrt(market.timeToExpiry);
    const core::Real d2 = (std::log(market.spot / market.strike) + 
                          (market.riskFreeRate - 0.5 * market.volatility * market.volatility) * market.timeToExpiry) / 
                         (market.volatility * sqrtT);
    
    const core::Real nd2 = engine_->normalCDF(d2);
    return market.isCall ? nd2 : (1.0 - nd2);
}

inline core::Real ProbabilityEngine::calculateTouchProbability(core::Real spot, core::Real barrier,
                                                              core::Real vol, core::Real time, bool isUp) const {
    if (time <= 0.0 || vol <= 0.0 || spot <= 0.0) return 0.0;
    
    if ((isUp && spot >= barrier) || (!isUp && spot <= barrier)) return 1.0;
    
    const core::Real mu = -0.5 * vol * vol;
    const core::Real logRatio = std::log(barrier / spot);
    const core::Real sqrtT = std::sqrt(time);
    
    const core::Real lambda1 = (logRatio - mu * time) / (vol * sqrtT);
    const core::Real lambda2 = (logRatio + mu * time) / (vol * sqrtT);
    
    const core::Real barrierFactor = std::pow(barrier / spot, 2.0 * mu / (vol * vol));
    
    core::Real touchProb = engine_->normalCDF(lambda1) + barrierFactor * engine_->normalCDF(lambda2);
    
    if (!isUp) {
        touchProb = 1.0 - touchProb;
    }
    
    return std::clamp(touchProb, 0.0, 1.0);
}

inline core::Real ProbabilityEngine::calculateKellyCriterion(core::Real winRate, core::Real avgWin, core::Real avgLoss) const {
    if (avgLoss <= 0.0 || winRate <= 0.0 || winRate >= 1.0) return 0.0;
    
    const core::Real lossRate = 1.0 - winRate;
    const core::Real payoffRatio = avgWin / avgLoss;
    const core::Real kelly = (winRate * payoffRatio - lossRate) / payoffRatio;
    
    return std::clamp(kelly, 0.0, 1.0);
}

}