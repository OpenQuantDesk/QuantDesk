#include "variance_reduction.hpp"
#include "../core/utils.hpp"
#include <algorithm>
#include <numeric>
#include <future>
#include <thread>

namespace math::simulation {

VarianceReducer::VarianceReducer(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::SimulationResult VarianceReducer::antitheticVariates(
    const core::SimulationParams& params,
    const std::function<core::Real(core::Real)>& payoffFunc) const {
    
    core::SimulationResult result;
    
    if (params.numSimulations % 2 != 0) {
        core::SimulationParams adjustedParams = params;
        adjustedParams.numSimulations += 1;
        return antitheticVariates(adjustedParams, payoffFunc);
    }
    
    const core::Size numPairs = params.numSimulations / 2;
    const core::Size numThreads = engine_->getThreadPool().getNumThreads();
    const core::Size pairsPerThread = (numPairs + numThreads - 1) / numThreads;
    
    std::vector<std::future<std::pair<core::Vector, core::Vector>>> futures;
    futures.reserve(numThreads);
    
    const core::Real dt = params.timeToExpiry / params.timeSteps;
    const core::Real drift = (params.riskFreeRate - 0.5 * params.volatility * params.volatility) * dt;
    const core::Real diffusion = params.volatility * std::sqrt(dt);
    
    for (core::Size t = 0; t < numThreads; ++t) {
        const core::Size startPair = t * pairsPerThread;
        const core::Size endPair = std::min(startPair + pairsPerThread, numPairs);
        
        if (startPair >= numPairs) break;
        
        futures.emplace_back(std::async(std::launch::async, [=, &payoffFunc]() {
            auto& rng = engine_->getRNG();
            rng.seed(params.seed + t);
            std::normal_distribution<core::Real> normal(0.0, 1.0);
            
            core::Vector regularPayoffs;
            core::Vector antitheticPayoffs;
            regularPayoffs.reserve(endPair - startPair);
            antitheticPayoffs.reserve(endPair - startPair);
            
            for (core::Size pair = startPair; pair < endPair; ++pair) {
                core::Real spot1 = params.spot;
                core::Real spot2 = params.spot;
                
                for (core::Integer step = 0; step < params.timeSteps; ++step) {
                    const core::Real z = normal(rng);
                    spot1 *= std::exp(drift + diffusion * z);
                    spot2 *= std::exp(drift + diffusion * (-z));
                }
                
                regularPayoffs.push_back(payoffFunc(spot1));
                antitheticPayoffs.push_back(payoffFunc(spot2));
            }
            
            return std::make_pair(regularPayoffs, antitheticPayoffs);
        }));
    }
    
    core::Vector allPayoffs;
    allPayoffs.reserve(params.numSimulations);
    
    for (auto& future : futures) {
        auto [regular, antithetic] = future.get();
        allPayoffs.insert(allPayoffs.end(), regular.begin(), regular.end());
        allPayoffs.insert(allPayoffs.end(), antithetic.begin(), antithetic.end());
    }
    
    result.payoffs = allPayoffs;
    result.expectedValue = core::Statistics::mean(allPayoffs.begin(), allPayoffs.end());
    
    const core::Real variance = core::Statistics::variance(allPayoffs.begin(), allPayoffs.end(), 
                                                          result.expectedValue);
    result.standardError = std::sqrt(variance / allPayoffs.size());
    
    const core::Integer profitable = std::count_if(allPayoffs.begin(), allPayoffs.end(),
                                                  [](core::Real p) { return p > 0; });
    result.probabilityProfit = static_cast<core::Real>(profitable) / allPayoffs.size();
    
    std::sort(allPayoffs.begin(), allPayoffs.end());
    result.valueAtRisk95 = core::Statistics::percentile(allPayoffs.begin(), allPayoffs.end(), 0.05);
    result.valueAtRisk99 = core::Statistics::percentile(allPayoffs.begin(), allPayoffs.end(), 0.01);
    
    return result;
}

core::SimulationResult VarianceReducer::controlVariates(
    const core::SimulationParams& params,
    const std::function<core::Real(core::Real)>& payoffFunc,
    const std::function<core::Real(core::Real)>& controlFunc) const {
    
    core::SimulationResult result;
    
    const core::Size numThreads = engine_->getThreadPool().getNumThreads();
    const core::Size simsPerThread = (params.numSimulations + numThreads - 1) / numThreads;
    
    std::vector<std::future<std::tuple<core::Vector, core::Vector, core::Vector>>> futures;
    futures.reserve(numThreads);
    
    const core::Real dt = params.timeToExpiry / params.timeSteps;
    const core::Real drift = (params.riskFreeRate - 0.5 * params.volatility * params.volatility) * dt;
    const core::Real diffusion = params.volatility * std::sqrt(dt);
    
    for (core::Size t = 0; t < numThreads; ++t) {
        const core::Size startSim = t * simsPerThread;
        const core::Size endSim = std::min(startSim + simsPerThread, 
                                          static_cast<core::Size>(params.numSimulations));
        
        if (startSim >= static_cast<core::Size>(params.numSimulations)) break;
        
        futures.emplace_back(std::async(std::launch::async, [=, &payoffFunc, &controlFunc]() {
            auto& rng = engine_->getRNG();
            rng.seed(params.seed + t);
            std::normal_distribution<core::Real> normal(0.0, 1.0);
            
            core::Vector payoffs;
            core::Vector controlValues;
            core::Vector finalSpots;
            
            payoffs.reserve(endSim - startSim);
            controlValues.reserve(endSim - startSim);
            finalSpots.reserve(endSim - startSim);
            
            for (core::Size sim = startSim; sim < endSim; ++sim) {
                core::Real spot = params.spot;
                
                for (core::Integer step = 0; step < params.timeSteps; ++step) {
                    const core::Real z = normal(rng);
                    spot *= std::exp(drift + diffusion * z);
                }
                
                finalSpots.push_back(spot);
                payoffs.push_back(payoffFunc(spot));
                controlValues.push_back(controlFunc(spot));
            }
            
            return std::make_tuple(payoffs, controlValues, finalSpots);
        }));
    }
    
    core::Vector allPayoffs;
    core::Vector allControlValues;
    core::Vector allFinalSpots;
    
    allPayoffs.reserve(params.numSimulations);
    allControlValues.reserve(params.numSimulations);
    allFinalSpots.reserve(params.numSimulations);
    
    for (auto& future : futures) {
        auto [payoffs, controls, spots] = future.get();
        allPayoffs.insert(allPayoffs.end(), payoffs.begin(), payoffs.end());
        allControlValues.insert(allControlValues.end(), controls.begin(), controls.end());
        allFinalSpots.insert(allFinalSpots.end(), spots.begin(), spots.end());
    }
    
    const core::Real payoffMean = core::Statistics::mean(allPayoffs.begin(), allPayoffs.end());
    const core::Real controlMean = core::Statistics::mean(allControlValues.begin(), allControlValues.end());
    
    core::Real covariance = 0.0;
    core::Real controlVariance = 0.0;
    
    for (core::Size i = 0; i < allPayoffs.size(); ++i) {
        const core::Real payoffDev = allPayoffs[i] - payoffMean;
        const core::Real controlDev = allControlValues[i] - controlMean;
        
        covariance += payoffDev * controlDev;
        controlVariance += controlDev * controlDev;
    }
    
    covariance /= (allPayoffs.size() - 1);
    controlVariance /= (allPayoffs.size() - 1);
    
    const core::Real expectedControlValue = params.spot * std::exp(params.riskFreeRate * params.timeToExpiry);
    const core::Real optimalBeta = controlVariance > 0.0 ? covariance / controlVariance : 0.0;
    
    core::Vector adjustedPayoffs(allPayoffs.size());
    for (core::Size i = 0; i < allPayoffs.size(); ++i) {
        adjustedPayoffs[i] = allPayoffs[i] - optimalBeta * (allControlValues[i] - expectedControlValue);
    }
    
    result.payoffs = adjustedPayoffs;
    result.expectedValue = core::Statistics::mean(adjustedPayoffs.begin(), adjustedPayoffs.end());
    
    const core::Real variance = core::Statistics::variance(adjustedPayoffs.begin(), adjustedPayoffs.end(),
                                                          result.expectedValue);
    result.standardError = std::sqrt(variance / adjustedPayoffs.size());
    
    const core::Integer profitable = std::count_if(adjustedPayoffs.begin(), adjustedPayoffs.end(),
                                                  [](core::Real p) { return p > 0; });
    result.probabilityProfit = static_cast<core::Real>(profitable) / adjustedPayoffs.size();
    
    std::sort(adjustedPayoffs.begin(), adjustedPayoffs.end());
    result.valueAtRisk95 = core::Statistics::percentile(adjustedPayoffs.begin(), adjustedPayoffs.end(), 0.05);
    result.valueAtRisk99 = core::Statistics::percentile(adjustedPayoffs.begin(), adjustedPayoffs.end(), 0.01);
    
    return result;
}

core::SimulationResult VarianceReducer::importanceSampling(
    const core::SimulationParams& params,
    const std::function<core::Real(core::Real)>& payoffFunc,
    core::Real optimalDrift) const {
    
    core::SimulationResult result;
    
    const core::Size numThreads = engine_->getThreadPool().getNumThreads();
    const core::Size simsPerThread = (params.numSimulations + numThreads - 1) / numThreads;
    
    std::vector<std::future<core::Vector>> futures;
    futures.reserve(numThreads);
    
    const core::Real dt = params.timeToExpiry / params.timeSteps;
    const core::Real originalDrift = (params.riskFreeRate - 0.5 * params.volatility * params.volatility) * dt;
    const core::Real importanceDrift = originalDrift + optimalDrift * dt;
    const core::Real diffusion = params.volatility * std::sqrt(dt);
    
    for (core::Size t = 0; t < numThreads; ++t) {
        const core::Size startSim = t * simsPerThread;
        const core::Size endSim = std::min(startSim + simsPerThread,
                                          static_cast<core::Size>(params.numSimulations));
        
        if (startSim >= static_cast<core::Size>(params.numSimulations)) break;
        
        futures.emplace_back(std::async(std::launch::async, [=, &payoffFunc]() {
            auto& rng = engine_->getRNG();
            rng.seed(params.seed + t);
            std::normal_distribution<core::Real> normal(0.0, 1.0);
            
            core::Vector weightedPayoffs;
            weightedPayoffs.reserve(endSim - startSim);
            
            for (core::Size sim = startSim; sim < endSim; ++sim) {
                core::Real spot = params.spot;
                core::Real logLikelihoodRatio = 0.0;
                
                for (core::Integer step = 0; step < params.timeSteps; ++step) {
                    const core::Real z = normal(rng);
                    spot *= std::exp(importanceDrift + diffusion * z);
                    
                    const core::Real originalZ = z - optimalDrift * std::sqrt(dt) / diffusion;
                    logLikelihoodRatio += -0.5 * (z * z - originalZ * originalZ);
                }
                
                const core::Real likelihoodRatio = std::exp(logLikelihoodRatio);
                const core::Real payoff = payoffFunc(spot);
                weightedPayoffs.push_back(payoff * likelihoodRatio);
            }
            
            return weightedPayoffs;
        }));
    }
    
    core::Vector allPayoffs;
    allPayoffs.reserve(params.numSimulations);
    
    for (auto& future : futures) {
        auto payoffs = future.get();
        allPayoffs.insert(allPayoffs.end(), payoffs.begin(), payoffs.end());
    }
    
    result.payoffs = allPayoffs;
    result.expectedValue = core::Statistics::mean(allPayoffs.begin(), allPayoffs.end());
    
    const core::Real variance = core::Statistics::variance(allPayoffs.begin(), allPayoffs.end(),
                                                          result.expectedValue);
    result.standardError = std::sqrt(variance / allPayoffs.size());
    
    const core::Integer profitable = std::count_if(allPayoffs.begin(), allPayoffs.end(),
                                                  [](core::Real p) { return p > 0; });
    result.probabilityProfit = static_cast<core::Real>(profitable) / allPayoffs.size();
    
    std::sort(allPayoffs.begin(), allPayoffs.end());
    result.valueAtRisk95 = core::Statistics::percentile(allPayoffs.begin(), allPayoffs.end(), 0.05);
    result.valueAtRisk99 = core::Statistics::percentile(allPayoffs.begin(), allPayoffs.end(), 0.01);
    
    return result;
}

VarianceReducer::VarianceMetrics VarianceReducer::compareVarianceReduction(
    const core::SimulationParams& params,
    const std::function<core::Real(core::Real)>& payoffFunc) const {
    
    VarianceMetrics metrics;
    
    auto standardPayoff = [&](core::Real finalSpot) -> core::Real {
        return payoffFunc(finalSpot);
    };
    
    auto controlPayoff = [&](core::Real finalSpot) -> core::Real {
        return finalSpot;
    };
    
    core::SimulationParams testParams = params;
    testParams.numSimulations = std::min(params.numSimulations, core::Integer(10000));
    
    const auto standardResult = standardMonteCarlo(testParams, standardPayoff);
    const auto antitheticResult = antitheticVariates(testParams, standardPayoff);
    const auto controlResult = controlVariates(testParams, standardPayoff, controlPayoff);
    
    metrics.standardMC = standardResult.standardError * standardResult.standardError;
    metrics.antitheticMC = antitheticResult.standardError * antitheticResult.standardError;
    metrics.controlVariateMC = controlResult.standardError * controlResult.standardError;
    
    if (metrics.standardMC > 0.0) {
        metrics.varianceReduction = 1.0 - metrics.antitheticMC / metrics.standardMC;
        metrics.efficiency = metrics.standardMC / metrics.controlVariateMC;
    }
    
    return metrics;
}

core::SimulationResult VarianceReducer::standardMonteCarlo(
    const core::SimulationParams& params,
    const std::function<core::Real(core::Real)>& payoffFunc) const {
    
    core::SimulationResult result;
    
    const core::Size numThreads = engine_->getThreadPool().getNumThreads();
    const core::Size simsPerThread = (params.numSimulations + numThreads - 1) / numThreads;
    
    std::vector<std::future<core::Vector>> futures;
    futures.reserve(numThreads);
    
    const core::Real dt = params.timeToExpiry / params.timeSteps;
    const core::Real drift = (params.riskFreeRate - 0.5 * params.volatility * params.volatility) * dt;
    const core::Real diffusion = params.volatility * std::sqrt(dt);
    
    for (core::Size t = 0; t < numThreads; ++t) {
        const core::Size startSim = t * simsPerThread;
        const core::Size endSim = std::min(startSim + simsPerThread,
                                          static_cast<core::Size>(params.numSimulations));
        
        if (startSim >= static_cast<core::Size>(params.numSimulations)) break;
        
        futures.emplace_back(std::async(std::launch::async, [=, &payoffFunc]() {
            auto& rng = engine_->getRNG();
            rng.seed(params.seed + t);
            std::normal_distribution<core::Real> normal(0.0, 1.0);
            
            core::Vector payoffs;
            payoffs.reserve(endSim - startSim);
            
            for (core::Size sim = startSim; sim < endSim; ++sim) {
                core::Real spot = params.spot;
                
                for (core::Integer step = 0; step < params.timeSteps; ++step) {
                    const core::Real z = normal(rng);
                    spot *= std::exp(drift + diffusion * z);
                }
                
                payoffs.push_back(payoffFunc(spot));
            }
            
            return payoffs;
        }));
    }
    
    core::Vector allPayoffs;
    allPayoffs.reserve(params.numSimulations);
    
    for (auto& future : futures) {
        auto payoffs = future.get();
        allPayoffs.insert(allPayoffs.end(), payoffs.begin(), payoffs.end());
    }
    
    result.payoffs = allPayoffs;
    result.expectedValue = core::Statistics::mean(allPayoffs.begin(), allPayoffs.end());
    
    const core::Real variance = core::Statistics::variance(allPayoffs.begin(), allPayoffs.end(),
                                                          result.expectedValue);
    result.standardError = std::sqrt(variance / allPayoffs.size());
    
    const core::Integer profitable = std::count_if(allPayoffs.begin(), allPayoffs.end(),
                                                  [](core::Real p) { return p > 0; });
    result.probabilityProfit = static_cast<core::Real>(profitable) / allPayoffs.size();
    
    std::sort(allPayoffs.begin(), allPayoffs.end());
    result.valueAtRisk95 = core::Statistics::percentile(allPayoffs.begin(), allPayoffs.end(), 0.05);
    result.valueAtRisk99 = core::Statistics::percentile(allPayoffs.begin(), allPayoffs.end(), 0.01);
    
    return result;
}

}