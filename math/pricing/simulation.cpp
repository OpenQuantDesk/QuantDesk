#include "simulation.hpp"
#include "../core/utils.hpp"
#include <algorithm>
#include <numeric>

namespace math::pricing {

MonteCarloEngine::MonteCarloEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)),
      pathGenerator_(std::make_shared<simulation::PathGenerator>(engine)) {}

MonteCarloEngine::SimulationResult MonteCarloEngine::priceEuropean(
    const core::MarketData& market, const SimulationParams& params) const {
    
    SimulationResult result;
    
    simulation::PathGenerationParams pathParams;
    pathParams.spot = market.spot;
    pathParams.volatility = market.volatility;
    pathParams.riskFreeRate = market.riskFreeRate;
    pathParams.dividendYield = market.dividendYield;
    pathParams.timeToMaturity = market.timeToExpiry;
    pathParams.timeSteps = 1;
    pathParams.numPaths = params.numSimulations;
    pathParams.useAntitheticVariates = params.useAntitheticVariates;
    pathParams.seed = params.seed;
    
    const core::Matrix paths = pathGenerator_->generatePaths(pathParams);
    
    result.payoffs.resize(params.numSimulations);
    
    for (core::Integer i = 0; i < params.numSimulations; ++i) {
        const core::Real finalSpot = paths[i].back();
        
        if (market.isCall) {
            result.payoffs[i] = std::max(0.0, finalSpot - market.strike);
        } else {
            result.payoffs[i] = std::max(0.0, market.strike - finalSpot);
        }
    }
    
    const core::Real discountFactor = std::exp(-market.riskFreeRate * market.timeToExpiry);
    
    result.price = core::Statistics::mean(result.payoffs.begin(), result.payoffs.end()) * discountFactor;
    
    const core::Real variance = core::Statistics::variance(result.payoffs.begin(), result.payoffs.end(), 
                                                          result.price / discountFactor);
    result.standardError = std::sqrt(variance / params.numSimulations) * discountFactor;
    
    const core::Real confidenceInterval = 1.96 * result.standardError;
    result.confidence95Lower = result.price - confidenceInterval;
    result.confidence95Upper = result.price + confidenceInterval;
    
    const core::Size profitableOutcomes = std::count_if(result.payoffs.begin(), result.payoffs.end(),
                                                        [](core::Real payoff) { return payoff > 0.0; });
    result.probabilityProfit = static_cast<core::Real>(profitableOutcomes) / params.numSimulations;
    
    return result;
}

MonteCarloEngine::SimulationResult MonteCarloEngine::priceAsian(
    const core::MarketData& market, const SimulationParams& params, bool isArithmetic) const {
    
    SimulationResult result;
    
    simulation::PathGenerationParams pathParams;
    pathParams.spot = market.spot;
    pathParams.volatility = market.volatility;
    pathParams.riskFreeRate = market.riskFreeRate;
    pathParams.dividendYield = market.dividendYield;
    pathParams.timeToMaturity = market.timeToExpiry;
    pathParams.timeSteps = params.timeSteps;
    pathParams.numPaths = params.numSimulations;
    pathParams.useAntitheticVariates = params.useAntitheticVariates;
    pathParams.seed = params.seed;
    
    const core::Matrix paths = pathGenerator_->generatePaths(pathParams);
    
    result.payoffs.resize(params.numSimulations);
    
    for (core::Integer i = 0; i < params.numSimulations; ++i) {
        core::Real average;
        
        if (isArithmetic) {
            average = core::Statistics::mean(paths[i].begin(), paths[i].end());
        } else {
            core::Real logSum = 0.0;
            for (core::Real price : paths[i]) {
                logSum += std::log(std::max(price, core::EPSILON));
            }
            average = std::exp(logSum / paths[i].size());
        }
        
        if (market.isCall) {
            result.payoffs[i] = std::max(0.0, average - market.strike);
        } else {
            result.payoffs[i] = std::max(0.0, market.strike - average);
        }
    }
    
    const core::Real discountFactor = std::exp(-market.riskFreeRate * market.timeToExpiry);
    
    result.price = core::Statistics::mean(result.payoffs.begin(), result.payoffs.end()) * discountFactor;
    
    const core::Real variance = core::Statistics::variance(result.payoffs.begin(), result.payoffs.end(),
                                                          result.price / discountFactor);
    result.standardError = std::sqrt(variance / params.numSimulations) * discountFactor;
    
    const core::Real confidenceInterval = 1.96 * result.standardError;
    result.confidence95Lower = result.price - confidenceInterval;
    result.confidence95Upper = result.price + confidenceInterval;
    
    const core::Size profitableOutcomes = std::count_if(result.payoffs.begin(), result.payoffs.end(),
                                                        [](core::Real payoff) { return payoff > 0.0; });
    result.probabilityProfit = static_cast<core::Real>(profitableOutcomes) / params.numSimulations;
    
    return result;
}

MonteCarloEngine::SimulationResult MonteCarloEngine::priceBarrier(
    const core::MarketData& market, const SimulationParams& params,
    core::Real barrier, bool isKnockOut, bool isUp, core::Real rebate) const {
    
    SimulationResult result;
    
    simulation::PathGenerationParams pathParams;
    pathParams.spot = market.spot;
    pathParams.volatility = market.volatility;
    pathParams.riskFreeRate = market.riskFreeRate;
    pathParams.dividendYield = market.dividendYield;
    pathParams.timeToMaturity = market.timeToExpiry;
    pathParams.timeSteps = params.timeSteps;
    pathParams.numPaths = params.numSimulations;
    pathParams.useAntitheticVariates = params.useAntitheticVariates;
    pathParams.seed = params.seed;
    
    const core::Matrix paths = pathGenerator_->generatePaths(pathParams);
    
    result.payoffs.resize(params.numSimulations);
    
    for (core::Integer i = 0; i < params.numSimulations; ++i) {
        bool barrierHit = false;
        
        if (isUp) {
            for (core::Real price : paths[i]) {
                if (price >= barrier) {
                    barrierHit = true;
                    break;
                }
            }
        } else {
            for (core::Real price : paths[i]) {
                if (price <= barrier) {
                    barrierHit = true;
                    break;
                }
            }
        }
        
        const core::Real finalSpot = paths[i].back();
        core::Real intrinsicValue = 0.0;
        
        if (market.isCall) {
            intrinsicValue = std::max(0.0, finalSpot - market.strike);
        } else {
            intrinsicValue = std::max(0.0, market.strike - finalSpot);
        }
        
        if (isKnockOut) {
            result.payoffs[i] = barrierHit ? rebate : intrinsicValue;
        } else {
            result.payoffs[i] = barrierHit ? intrinsicValue : rebate;
        }
    }
    
    const core::Real discountFactor = std::exp(-market.riskFreeRate * market.timeToExpiry);
    
    result.price = core::Statistics::mean(result.payoffs.begin(), result.payoffs.end()) * discountFactor;
    
    const core::Real variance = core::Statistics::variance(result.payoffs.begin(), result.payoffs.end(),
                                                          result.price / discountFactor);
    result.standardError = std::sqrt(variance / params.numSimulations) * discountFactor;
    
    const core::Real confidenceInterval = 1.96 * result.standardError;
    result.confidence95Lower = result.price - confidenceInterval;
    result.confidence95Upper = result.price + confidenceInterval;
    
    const core::Size profitableOutcomes = std::count_if(result.payoffs.begin(), result.payoffs.end(),
                                                        [](core::Real payoff) { return payoff > 0.0; });
    result.probabilityProfit = static_cast<core::Real>(profitableOutcomes) / params.numSimulations;
    
    return result;
}

MonteCarloEngine::SimulationResult MonteCarloEngine::priceLookback(
    const core::MarketData& market, const SimulationParams& params, bool isFloating) const {
    
    SimulationResult result;
    
    simulation::PathGenerationParams pathParams;
    pathParams.spot = market.spot;
    pathParams.volatility = market.volatility;
    pathParams.riskFreeRate = market.riskFreeRate;
    pathParams.dividendYield = market.dividendYield;
    pathParams.timeToMaturity = market.timeToExpiry;
    pathParams.timeSteps = params.timeSteps;
    pathParams.numPaths = params.numSimulations;
    pathParams.useAntitheticVariates = params.useAntitheticVariates;
    pathParams.seed = params.seed;
    
    const core::Matrix paths = pathGenerator_->generatePaths(pathParams);
    
    result.payoffs.resize(params.numSimulations);
    
    for (core::Integer i = 0; i < params.numSimulations; ++i) {
        const core::Real maxPrice = *std::max_element(paths[i].begin(), paths[i].end());
        const core::Real minPrice = *std::min_element(paths[i].begin(), paths[i].end());
        const core::Real finalPrice = paths[i].back();
        
        if (isFloating) {
            if (market.isCall) {
                result.payoffs[i] = finalPrice - minPrice;
            } else {
                result.payoffs[i] = maxPrice - finalPrice;
            }
        } else {
            if (market.isCall) {
                result.payoffs[i] = std::max(0.0, maxPrice - market.strike);
            } else {
                result.payoffs[i] = std::max(0.0, market.strike - minPrice);
            }
        }
    }
    
    const core::Real discountFactor = std::exp(-market.riskFreeRate * market.timeToExpiry);
    
    result.price = core::Statistics::mean(result.payoffs.begin(), result.payoffs.end()) * discountFactor;
    
    const core::Real variance = core::Statistics::variance(result.payoffs.begin(), result.payoffs.end(),
                                                          result.price / discountFactor);
    result.standardError = std::sqrt(variance / params.numSimulations) * discountFactor;
    
    const core::Real confidenceInterval = 1.96 * result.standardError;
    result.confidence95Lower = result.price - confidenceInterval;
    result.confidence95Upper = result.price + confidenceInterval;
    
    const core::Size profitableOutcomes = std::count_if(result.payoffs.begin(), result.payoffs.end(),
                                                        [](core::Real payoff) { return payoff > 0.0; });
    result.probabilityProfit = static_cast<core::Real>(profitableOutcomes) / params.numSimulations;
    
    return result;
}

MonteCarloEngine::VarianceReductionResult MonteCarloEngine::compareVarianceReduction(
    const core::MarketData& market, const SimulationParams& params) const {
    
    VarianceReductionResult result;
    
    SimulationParams standardParams = params;
    standardParams.useAntitheticVariates = false;
    standardParams.useControlVariates = false;
    
    SimulationParams antitheticParams = params;
    antitheticParams.useAntitheticVariates = true;
    antitheticParams.useControlVariates = false;
    
    SimulationParams controlParams = params;
    controlParams.useAntitheticVariates = false;
    controlParams.useControlVariates = true;
    
    const auto standardResult = priceEuropean(market, standardParams);
    const auto antitheticResult = priceEuropean(market, antitheticParams);
    const auto controlResult = priceEuropean(market, controlParams);
    
    result.standardMC = standardResult.standardError;
    result.antitheticMC = antitheticResult.standardError;
    result.controlVariateMC = controlResult.standardError;
    result.importanceSamplingMC = standardResult.standardError;
    
    if (result.standardMC > 0.0) {
        result.varianceReduction = 1.0 - (std::min({result.antitheticMC, result.controlVariateMC}) / result.standardMC);
    }
    
    return result;
}

QuasiMonteCarloEngine::QuasiMonteCarloEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

QuasiMonteCarloEngine::QMCResult QuasiMonteCarloEngine::priceSobol(
    const core::MarketData& market, const QMCParams& params) const {
    
    QMCResult result;
    
    const core::Size numDimensions = params.timeSteps;
    core::Matrix sobolSequence = generateSobolSequence(params.numSimulations, numDimensions);
    
    result.payoffs.resize(params.numSimulations);
    
    const core::Real dt = market.timeToExpiry / params.timeSteps;
    const core::Real drift = (market.riskFreeRate - market.dividendYield - 
                             0.5 * market.volatility * market.volatility) * dt;
    const core::Real diffusion = market.volatility * std::sqrt(dt);
    
    for (core::Size i = 0; i < params.numSimulations; ++i) {
        core::Real spot = market.spot;
        
        for (core::Size step = 0; step < params.timeSteps; ++step) {
            const core::Real u = sobolSequence[i][step];
            const core::Real z = engine_->normalInverse(u);
            spot *= std::exp(drift + diffusion * z);
        }
        
        if (market.isCall) {
            result.payoffs[i] = std::max(0.0, spot - market.strike);
        } else {
            result.payoffs[i] = std::max(0.0, market.strike - spot);
        }
    }
    
    const core::Real discountFactor = std::exp(-market.riskFreeRate * market.timeToExpiry);
    
    result.price = core::Statistics::mean(result.payoffs.begin(), result.payoffs.end()) * discountFactor;
    
    const core::Real variance = core::Statistics::variance(result.payoffs.begin(), result.payoffs.end(),
                                                          result.price / discountFactor);
    result.standardError = std::sqrt(variance / params.numSimulations) * discountFactor;
    
    result.convergenceRate = calculateConvergenceRate(result.payoffs);
    
    return result;
}

QuasiMonteCarloEngine::QMCResult QuasiMonteCarloEngine::priceHalton(
    const core::MarketData& market, const QMCParams& params) const {
    
    QMCResult result;
    
    const core::Size numDimensions = params.timeSteps;
    core::Matrix haltonSequence = generateHaltonSequence(params.numSimulations, numDimensions);
    
    result.payoffs.resize(params.numSimulations);
    
    const core::Real dt = market.timeToExpiry / params.timeSteps;
    const core::Real drift = (market.riskFreeRate - market.dividendYield - 
                             0.5 * market.volatility * market.volatility) * dt;
    const core::Real diffusion = market.volatility * std::sqrt(dt);
    
    for (core::Size i = 0; i < params.numSimulations; ++i) {
        core::Real spot = market.spot;
        
        for (core::Size step = 0; step < params.timeSteps; ++step) {
            const core::Real u = haltonSequence[i][step];
            const core::Real z = engine_->normalInverse(u);
            spot *= std::exp(drift + diffusion * z);
        }
        
        if (market.isCall) {
            result.payoffs[i] = std::max(0.0, spot - market.strike);
        } else {
            result.payoffs[i] = std::max(0.0, market.strike - spot);
        }
    }
    
    const core::Real discountFactor = std::exp(-market.riskFreeRate * market.timeToExpiry);
    
    result.price = core::Statistics::mean(result.payoffs.begin(), result.payoffs.end()) * discountFactor;
    
    const core::Real variance = core::Statistics::variance(result.payoffs.begin(), result.payoffs.end(),
                                                          result.price / discountFactor);
    result.standardError = std::sqrt(variance / params.numSimulations) * discountFactor;
    
    result.convergenceRate = calculateConvergenceRate(result.payoffs);
    
    return result;
}

core::Matrix QuasiMonteCarloEngine::generateSobolSequence(core::Size numPoints, core::Size dimensions) const {
    core::Matrix sequence(numPoints, core::Vector(dimensions));
    
    for (core::Size i = 0; i < numPoints; ++i) {
        for (core::Size d = 0; d < dimensions; ++d) {
            sequence[i][d] = sobolPoint(i, d);
        }
    }
    
    return sequence;
}

core::Matrix QuasiMonteCarloEngine::generateHaltonSequence(core::Size numPoints, core::Size dimensions) const {
    core::Matrix sequence(numPoints, core::Vector(dimensions));
    
    const core::Vector primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    
    for (core::Size i = 0; i < numPoints; ++i) {
        for (core::Size d = 0; d < dimensions && d < primes.size(); ++d) {
            sequence[i][d] = haltonPoint(i + 1, static_cast<core::Integer>(primes[d]));
        }
    }
    
    return sequence;
}

core::Real QuasiMonteCarloEngine::sobolPoint(core::Size index, core::Size dimension) const {
    if (dimension == 0) {
        core::Real result = 0.0;
        core::Real fraction = 0.5;
        core::Size n = index;
        
        while (n > 0) {
            if (n & 1) {
                result += fraction;
            }
            fraction *= 0.5;
            n >>= 1;
        }
        
        return result;
    }
    
    return static_cast<core::Real>(index) / static_cast<core::Real>(1ULL << 32);
}

core::Real QuasiMonteCarloEngine::haltonPoint(core::Size index, core::Integer base) const {
    core::Real result = 0.0;
    core::Real fraction = 1.0 / base;
    core::Size n = index;
    
    while (n > 0) {
        result += (n % base) * fraction;
        n /= base;
        fraction /= base;
    }
    
    return result;
}

core::Real QuasiMonteCarloEngine::calculateConvergenceRate(const core::Vector& payoffs) const {
    if (payoffs.size() < 100) return 0.0;
    
    const core::Size n1 = payoffs.size() / 4;
    const core::Size n2 = payoffs.size() / 2;
    const core::Size n3 = payoffs.size();
    
    const core::Real mean1 = core::Statistics::mean(payoffs.begin(), payoffs.begin() + n1);
    const core::Real mean2 = core::Statistics::mean(payoffs.begin(), payoffs.begin() + n2);
    const core::Real mean3 = core::Statistics::mean(payoffs.begin(), payoffs.begin() + n3);
    
    const core::Real error1 = std::abs(mean1 - mean3);
    const core::Real error2 = std::abs(mean2 - mean3);
    
    if (error1 > 0.0 && error2 > 0.0) {
        return std::log(error1 / error2) / std::log(static_cast<core::Real>(n2) / n1);
    }
    
    return 0.5;
}

}