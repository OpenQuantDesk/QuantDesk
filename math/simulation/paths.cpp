#include "paths.hpp"
#include "../core/utils.hpp"
#include <algorithm>
#include <numeric>

namespace math::simulation {

core::Vector PathGenerator::generateSinglePath(const PathGenerationParams& params) const {
    core::Vector path(params.timeSteps + 1);
    path[0] = params.spot;
    
    const core::Real dt = params.timeToMaturity / params.timeSteps;
    const core::Real drift = (params.riskFreeRate - params.dividendYield - 0.5 * params.volatility * params.volatility) * dt;
    const core::Real diffusion = params.volatility * std::sqrt(dt);
    
    auto& rng = rngStorage_.get();
    std::normal_distribution<core::Real> normal(0.0, 1.0);
    
    for (core::Integer step = 1; step <= params.timeSteps; ++step) {
        const core::Real z = normal(rng);
        path[step] = path[step - 1] * std::exp(drift + diffusion * z);
    }
    
    return path;
}

core::Matrix PathGenerator::generatePaths(const PathGenerationParams& params) const {
    core::Matrix paths(params.numPaths, core::Vector(params.timeSteps + 1));
    
    const core::Size numThreads = engine_->getThreadPool().getNumThreads();
    const core::Size pathsPerThread = (params.numPaths + numThreads - 1) / numThreads;
    
    engine_->getThreadPool().parallelFor(
        core::Vector(numThreads).begin(),
        core::Vector(numThreads).end(),
        [this, &paths, &params, pathsPerThread, numThreads](core::Real threadIndex) {
            const core::Size threadIdx = static_cast<core::Size>(threadIndex);
            const core::Size startPath = threadIdx * pathsPerThread;
            const core::Size endPath = std::min(startPath + pathsPerThread, static_cast<core::Size>(params.numPaths));
            
            auto& rng = rngStorage_.get();
            rng.seed(params.seed + threadIdx);
            
            const core::Real dt = params.timeToMaturity / params.timeSteps;
            const core::Real drift = (params.riskFreeRate - params.dividendYield - 0.5 * params.volatility * params.volatility) * dt;
            const core::Real diffusion = params.volatility * std::sqrt(dt);
            
            std::normal_distribution<core::Real> normal(0.0, 1.0);
            
            for (core::Size pathIdx = startPath; pathIdx < endPath; ++pathIdx) {
                paths[pathIdx][0] = params.spot;
                
                for (core::Integer step = 1; step <= params.timeSteps; ++step) {
                    const core::Real z = normal(rng);
                    paths[pathIdx][step] = paths[pathIdx][step - 1] * std::exp(drift + diffusion * z);
                }
            }
        });
    
    return paths;
}

core::Matrix PathGenerator::generateAntitheticPaths(const PathGenerationParams& params) const {
    if (params.numPaths % 2 != 0) {
        PathGenerationParams adjustedParams = params;
        adjustedParams.numPaths += 1;
        return generateAntitheticPaths(adjustedParams);
    }
    
    core::Matrix paths(params.numPaths, core::Vector(params.timeSteps + 1));
    
    const core::Real dt = params.timeToMaturity / params.timeSteps;
    const core::Real drift = (params.riskFreeRate - params.dividendYield - 0.5 * params.volatility * params.volatility) * dt;
    const core::Real diffusion = params.volatility * std::sqrt(dt);
    
    auto& rng = rngStorage_.get();
    rng.seed(params.seed);
    std::normal_distribution<core::Real> normal(0.0, 1.0);
    
    for (core::Integer pathIdx = 0; pathIdx < params.numPaths; pathIdx += 2) {
        paths[pathIdx][0] = params.spot;
        paths[pathIdx + 1][0] = params.spot;
        
        for (core::Integer step = 1; step <= params.timeSteps; ++step) {
            const core::Real z = normal(rng);
            
            paths[pathIdx][step] = paths[pathIdx][step - 1] * std::exp(drift + diffusion * z);
            paths[pathIdx + 1][step] = paths[pathIdx + 1][step - 1] * std::exp(drift + diffusion * (-z));
        }
    }
    
    return paths;
}

std::vector<core::Matrix> PathGenerator::generateMultiAssetPaths(const MultiAssetParams& params) const {
    const core::Size numAssets = params.spots.size();
    std::vector<core::Matrix> allPaths(numAssets);
    
    for (core::Size i = 0; i < numAssets; ++i) {
        allPaths[i].resize(params.numPaths, core::Vector(params.timeSteps + 1));
    }
    
    const core::Matrix cholesky = choleskyDecomposition(params.correlationMatrix);
    const core::Real dt = params.timeToMaturity / params.timeSteps;
    
    auto& rng = rngStorage_.get();
    rng.seed(params.seed);
    std::normal_distribution<core::Real> normal(0.0, 1.0);
    
    for (core::Integer pathIdx = 0; pathIdx < params.numPaths; ++pathIdx) {
        for (core::Size assetIdx = 0; assetIdx < numAssets; ++assetIdx) {
            allPaths[assetIdx][pathIdx][0] = params.spots[assetIdx];
        }
        
        for (core::Integer step = 1; step <= params.timeSteps; ++step) {
            core::Vector independentNormals(numAssets);
            for (core::Size i = 0; i < numAssets; ++i) {
                independentNormals[i] = normal(rng);
            }
            
            core::Vector correlatedNormals(numAssets, 0.0);
            for (core::Size i = 0; i < numAssets; ++i) {
                for (core::Size j = 0; j <= i; ++j) {
                    correlatedNormals[i] += cholesky[i][j] * independentNormals[j];
                }
            }
            
            for (core::Size assetIdx = 0; assetIdx < numAssets; ++assetIdx) {
                const core::Real vol = assetIdx < params.volatilities.size() ? params.volatilities[assetIdx] : 0.2;
                const core::Real rate = assetIdx < params.riskFreeRates.size() ? params.riskFreeRates[assetIdx] : 0.05;
                const core::Real div = assetIdx < params.dividendYields.size() ? params.dividendYields[assetIdx] : 0.0;
                
                const core::Real drift = (rate - div - 0.5 * vol * vol) * dt;
                const core::Real diffusion = vol * std::sqrt(dt);
                
                allPaths[assetIdx][pathIdx][step] = allPaths[assetIdx][pathIdx][step - 1] * 
                                                   std::exp(drift + diffusion * correlatedNormals[assetIdx]);
            }
        }
    }
    
    return allPaths;
}

core::Vector PathGenerator::generateJumpDiffusionPath(const PathGenerationParams& params,
                                                     core::Real jumpIntensity,
                                                     core::Real jumpMean,
                                                     core::Real jumpStddev) const {
    core::Vector path(params.timeSteps + 1);
    path[0] = params.spot;
    
    const core::Real dt = params.timeToMaturity / params.timeSteps;
    const core::Real drift = (params.riskFreeRate - params.dividendYield - 0.5 * params.volatility * params.volatility) * dt;
    const core::Real diffusion = params.volatility * std::sqrt(dt);
    
    auto& rng = rngStorage_.get();
    std::normal_distribution<core::Real> normal(0.0, 1.0);
    std::normal_distribution<core::Real> jumpNormal(jumpMean, jumpStddev);
    std::exponential_distribution<core::Real> exponential(jumpIntensity);
    
    core::Real nextJumpTime = exponential(rng);
    core::Real currentTime = 0.0;
    
    for (core::Integer step = 1; step <= params.timeSteps; ++step) {
        const core::Real z = normal(rng);
        const core::Real stepTime = step * dt;
        
        core::Real jumpComponent = 0.0;
        while (nextJumpTime <= stepTime) {
            jumpComponent += jumpNormal(rng);
            nextJumpTime += exponential(rng);
        }
        
        path[step] = path[step - 1] * std::exp(drift + diffusion * z + jumpComponent);
    }
    
    return path;
}

core::Vector PathGenerator::generateHestonPath(const PathGenerationParams& params,
                                              core::Real kappa, core::Real theta,
                                              core::Real xi, core::Real rho, core::Real v0) const {
    core::Vector spotPath(params.timeSteps + 1);
    core::Vector varPath(params.timeSteps + 1);
    
    spotPath[0] = params.spot;
    varPath[0] = v0;
    
    const core::Real dt = params.timeToMaturity / params.timeSteps;
    const core::Real sqrtDt = std::sqrt(dt);
    const core::Real rhoComp = std::sqrt(1.0 - rho * rho);
    
    auto& rng = rngStorage_.get();
    std::normal_distribution<core::Real> normal(0.0, 1.0);
    
    for (core::Integer step = 1; step <= params.timeSteps; ++step) {
        const core::Real z1 = normal(rng);
        const core::Real z2 = normal(rng);
        const core::Real zCorr = rho * z1 + rhoComp * z2;
        
        const core::Real currentVar = std::max(varPath[step - 1], 0.0);
        const core::Real sqrtVar = std::sqrt(currentVar);
        
        varPath[step] = currentVar + kappa * (theta - currentVar) * dt + xi * sqrtVar * sqrtDt * z1;
        varPath[step] = std::max(varPath[step], 0.0);
        
        const core::Real drift = (params.riskFreeRate - params.dividendYield - 0.5 * currentVar) * dt;
        spotPath[step] = spotPath[step - 1] * std::exp(drift + sqrtVar * sqrtDt * zCorr);
    }
    
    return spotPath;
}

core::Vector PathGenerator::generateVarianceGammaPath(const PathGenerationParams& params,
                                                     core::Real sigma, core::Real nu, core::Real theta) const {
    core::Vector path(params.timeSteps + 1);
    path[0] = params.spot;
    
    const core::Real dt = params.timeToMaturity / params.timeSteps;
    
    auto& rng = rngStorage_.get();
    std::normal_distribution<core::Real> normal(0.0, 1.0);
    
    for (core::Integer step = 1; step <= params.timeSteps; ++step) {
        const core::Real gammaVariate = generateGammaVariate(dt / nu, nu, rng);
        const core::Real normalVariate = normal(rng);
        
        const core::Real vgIncrement = theta * gammaVariate + sigma * std::sqrt(gammaVariate) * normalVariate;
        const core::Real drift = (params.riskFreeRate - params.dividendYield) * dt;
        
        path[step] = path[step - 1] * std::exp(drift + vgIncrement);
    }
    
    return path;
}

PathGenerator::PathStatistics PathGenerator::analyzePathStatistics(const core::Matrix& paths) const {
    PathStatistics stats;
    
    if (paths.empty() || paths[0].empty()) return stats;
    
    const core::Size numPaths = paths.size();
    const core::Size numSteps = paths[0].size();
    
    core::Vector finalValues(numPaths);
    for (core::Size i = 0; i < numPaths; ++i) {
        finalValues[i] = paths[i].back();
    }
    
    stats.finalMean = core::Statistics::mean(finalValues.begin(), finalValues.end());
    stats.finalStddev = core::Statistics::standardDeviation(finalValues.begin(), finalValues.end(), stats.finalMean);
    
    stats.maxValue = *std::max_element(finalValues.begin(), finalValues.end());
    stats.minValue = *std::min_element(finalValues.begin(), finalValues.end());
    
    stats.percentiles.resize(9);
    const core::Vector percentileLevels = {0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.999};
    
    core::Vector sortedValues = finalValues;
    std::sort(sortedValues.begin(), sortedValues.end());
    
    for (core::Size i = 0; i < percentileLevels.size(); ++i) {
        stats.percentiles[i] = core::Statistics::percentile(sortedValues.begin(), sortedValues.end(), percentileLevels[i]);
    }
    
    core::Real totalRealizedVar = 0.0;
    for (core::Size pathIdx = 0; pathIdx < numPaths; ++pathIdx) {
        core::Real pathRealizedVar = 0.0;
        for (core::Size step = 1; step < numSteps; ++step) {
            const core::Real logReturn = std::log(paths[pathIdx][step] / paths[pathIdx][step - 1]);
            pathRealizedVar += logReturn * logReturn;
        }
        totalRealizedVar += pathRealizedVar;
    }
    
    stats.realizedVolatility = std::sqrt(totalRealizedVar / (numPaths * (numSteps - 1))) * std::sqrt(252.0);
    
    if (numSteps > 1) {
        const core::Real initialValue = paths[0][0];
        const core::Real timeToMaturity = 1.0;
        stats.driftRate = std::log(stats.finalMean / initialValue) / timeToMaturity;
    }
    
    return stats;
}

core::Matrix PathGenerator::choleskyDecomposition(const core::Matrix& correlation) const {
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

core::Real PathGenerator::generateGammaVariate(core::Real shape, core::Real scale, std::mt19937& rng) const {
    std::gamma_distribution<core::Real> gamma(shape, scale);
    return gamma(rng);
}

core::Real PathGenerator::generatePoissonVariate(core::Real lambda, std::mt19937& rng) const {
    std::poisson_distribution<int> poisson(lambda);
    return static_cast<core::Real>(poisson(rng));
}

MonteCarloEngine::MCResult MonteCarloEngine::priceEuropeanOption(const core::MarketData& market,
                                                                const PathGenerationParams& params) const {
    MCResult result;
    
    const core::Matrix paths = pathGenerator_->generatePaths(params);
    result.allPayoffs.resize(params.numPaths);
    
    for (core::Integer i = 0; i < params.numPaths; ++i) {
        const core::Real finalSpot = paths[i].back();
        result.allPayoffs[i] = calculateEuropeanPayoff(finalSpot, market);
    }
    
    result.price = core::Statistics::mean(result.allPayoffs.begin(), result.allPayoffs.end());
    result.price *= std::exp(-market.riskFreeRate * market.timeToExpiry);
    
    const core::Real variance = core::Statistics::variance(result.allPayoffs.begin(), result.allPayoffs.end(), result.price);
    result.standardError = std::sqrt(variance / params.numPaths);
    
    const core::Real confidenceInterval = 1.96 * result.standardError;
    result.confidence95Lower = result.price - confidenceInterval;
    result.confidence95Upper = result.price + confidenceInterval;
    
    return result;
}

MonteCarloEngine::MCResult MonteCarloEngine::priceAsianOption(const core::MarketData& market,
                                                             const PathGenerationParams& params,
                                                             bool isArithmetic) const {
    MCResult result;
    
    const core::Matrix paths = pathGenerator_->generatePaths(params);
    result.allPayoffs.resize(params.numPaths);
    
    for (core::Integer i = 0; i < params.numPaths; ++i) {
        result.allPayoffs[i] = calculateAsianPayoff(paths[i], market, isArithmetic);
    }
    
    result.price = core::Statistics::mean(result.allPayoffs.begin(), result.allPayoffs.end());
    result.price *= std::exp(-market.riskFreeRate * market.timeToExpiry);
    
    const core::Real variance = core::Statistics::variance(result.allPayoffs.begin(), result.allPayoffs.end(), result.price);
    result.standardError = std::sqrt(variance / params.numPaths);
    
    const core::Real confidenceInterval = 1.96 * result.standardError;
    result.confidence95Lower = result.price - confidenceInterval;
    result.confidence95Upper = result.price + confidenceInterval;
    
    return result;
}

MonteCarloEngine::MCResult MonteCarloEngine::priceBarrierOption(const core::MarketData& market,
                                                               const PathGenerationParams& params,
                                                               core::Real barrier, bool isKnockOut,
                                                               bool isUp, core::Real rebate) const {
    MCResult result;
    
    const core::Matrix paths = pathGenerator_->generatePaths(params);
    result.allPayoffs.resize(params.numPaths);
    
    for (core::Integer i = 0; i < params.numPaths; ++i) {
        result.allPayoffs[i] = calculateBarrierPayoff(paths[i], market, barrier, isKnockOut, isUp, rebate);
    }
    
    result.price = core::Statistics::mean(result.allPayoffs.begin(), result.allPayoffs.end());
    result.price *= std::exp(-market.riskFreeRate * market.timeToExpiry);
    
    const core::Real variance = core::Statistics::variance(result.allPayoffs.begin(), result.allPayoffs.end(), result.price);
    result.standardError = std::sqrt(variance / params.numPaths);
    
    const core::Real confidenceInterval = 1.96 * result.standardError;
    result.confidence95Lower = result.price - confidenceInterval;
    result.confidence95Upper = result.price + confidenceInterval;
    
    return result;
}

core::Real MonteCarloEngine::calculateBarrierPayoff(const core::Vector& path, const core::MarketData& market,
                                                   core::Real barrier, bool isKnockOut, bool isUp, core::Real rebate) const {
    if (path.empty()) return 0.0;
    
    bool barrierHit = false;
    
    if (isUp) {
        for (core::Real price : path) {
            if (price >= barrier) {
                barrierHit = true;
                break;
            }
        }
    } else {
        for (core::Real price : path) {
            if (price <= barrier) {
                barrierHit = true;
                break;
            }
        }
    }
    
    const core::Real finalSpot = path.back();
    const core::Real intrinsicValue = calculateEuropeanPayoff(finalSpot, market);
    
    if (isKnockOut) {
        return barrierHit ? rebate : intrinsicValue;
    } else {
        return barrierHit ? intrinsicValue : rebate;
    }
}

bool MonteCarloEngine::checkConvergence(const core::Vector& runningMeans, core::Real tolerance) const {
    if (runningMeans.size() < 100) return false;
    
    const core::Size windowSize = std::min(static_cast<core::Size>(50), runningMeans.size() / 2);
    const core::Real recentMean = core::Statistics::mean(runningMeans.end() - windowSize, runningMeans.end());
    const core::Real earlierMean = core::Statistics::mean(runningMeans.end() - 2 * windowSize, runningMeans.end() - windowSize);
    
    return std::abs(recentMean - earlierMean) < tolerance;
}

}