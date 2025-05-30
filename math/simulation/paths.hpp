#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include "../core/threading.hpp"
#include <random>
#include <memory>

namespace math::simulation {

struct PathGenerationParams {
    core::Real spot = 100.0;
    core::Real volatility = 0.2;
    core::Real riskFreeRate = 0.05;
    core::Real dividendYield = 0.0;
    core::Real timeToMaturity = 1.0;
    core::Integer timeSteps = 252;
    core::Integer numPaths = 10000;
    bool useAntitheticVariates = true;
    bool useControlVariates = false;
    uint32_t seed = 12345;
};

struct MultiAssetParams {
    core::Vector spots;
    core::Vector volatilities;
    core::Vector riskFreeRates;
    core::Vector dividendYields;
    core::Matrix correlationMatrix;
    core::Real timeToMaturity = 1.0;
    core::Integer timeSteps = 252;
    core::Integer numPaths = 10000;
    bool useAntitheticVariates = true;
    uint32_t seed = 12345;
};

class PathGenerator {
private:
    std::shared_ptr<core::MathEngine> engine_;
    mutable core::ThreadLocalStorage<std::mt19937> rngStorage_;
    
public:
    explicit PathGenerator(std::shared_ptr<core::MathEngine> engine);
    
    core::Vector generateSinglePath(const PathGenerationParams& params) const;
    
    core::Matrix generatePaths(const PathGenerationParams& params) const;
    
    core::Matrix generateAntitheticPaths(const PathGenerationParams& params) const;
    
    std::vector<core::Matrix> generateMultiAssetPaths(const MultiAssetParams& params) const;
    
    core::Vector generateJumpDiffusionPath(const PathGenerationParams& params,
                                          core::Real jumpIntensity,
                                          core::Real jumpMean,
                                          core::Real jumpStddev) const;
    
    core::Vector generateHestonPath(const PathGenerationParams& params,
                                   core::Real kappa,
                                   core::Real theta,
                                   core::Real xi,
                                   core::Real rho,
                                   core::Real v0) const;
    
    core::Vector generateVarianceGammaPath(const PathGenerationParams& params,
                                          core::Real sigma,
                                          core::Real nu,
                                          core::Real theta) const;
    
    struct PathStatistics {
        core::Real finalMean = 0.0;
        core::Real finalStddev = 0.0;
        core::Real driftRate = 0.0;
        core::Real realizedVolatility = 0.0;
        core::Real maxValue = 0.0;
        core::Real minValue = 0.0;
        core::Vector percentiles;
    };
    
    PathStatistics analyzePathStatistics(const core::Matrix& paths) const;
    
private:
    void generateCorrelatedNormals(core::Matrix& normals, const core::Matrix& cholesky) const;
    
    core::Matrix choleskyDecomposition(const core::Matrix& correlation) const;
    
    core::Real generateGammaVariate(core::Real shape, core::Real scale, std::mt19937& rng) const;
    
    core::Real generatePoissonVariate(core::Real lambda, std::mt19937& rng) const;
};

class MonteCarloEngine {
private:
    std::shared_ptr<PathGenerator> pathGenerator_;
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit MonteCarloEngine(std::shared_ptr<core::MathEngine> engine);
    
    struct MCResult {
        core::Real price = 0.0;
        core::Real standardError = 0.0;
        core::Real confidence95Lower = 0.0;
        core::Real confidence95Upper = 0.0;
        core::Vector allPayoffs;
        core::Integer convergenceIterations = 0;
    };
    
    MCResult priceEuropeanOption(const core::MarketData& market,
                                const PathGenerationParams& params) const;
    
    MCResult priceAmericanOption(const core::MarketData& market,
                                const PathGenerationParams& params) const;
    
    MCResult priceAsianOption(const core::MarketData& market,
                             const PathGenerationParams& params,
                             bool isArithmetic = true) const;
    
    MCResult priceLookbackOption(const core::MarketData& market,
                                const PathGenerationParams& params,
                                bool isFloating = true) const;
    
    MCResult priceBarrierOption(const core::MarketData& market,
                               const PathGenerationParams& params,
                               core::Real barrier,
                               bool isKnockOut = true,
                               bool isUp = true,
                               core::Real rebate = 0.0) const;
    
    MCResult priceMultiAssetOption(const MultiAssetParams& params,
                                  const std::function<core::Real(const core::Vector&)>& payoffFunc) const;
    
    struct VarianceReductionResult {
        core::Real standardMC = 0.0;
        core::Real antitheticMC = 0.0;
        core::Real controlVariateMC = 0.0;
        core::Real importanceSamplingMC = 0.0;
        core::Real varianceReduction = 0.0;
    };
    
    VarianceReductionResult compareVarianceReduction(const core::MarketData& market,
                                                    const PathGenerationParams& params) const;
    
private:
    core::Real calculateEuropeanPayoff(core::Real finalSpot, const core::MarketData& market) const;
    
    core::Real calculateAsianPayoff(const core::Vector& path, const core::MarketData& market,
                                   bool isArithmetic) const;
    
    core::Real calculateLookbackPayoff(const core::Vector& path, const core::MarketData& market,
                                      bool isFloating) const;
    
    core::Real calculateBarrierPayoff(const core::Vector& path, const core::MarketData& market,
                                     core::Real barrier, bool isKnockOut, bool isUp, core::Real rebate) const;
    
    core::Real calculateAmericanPayoff(const core::Vector& path, const core::MarketData& market) const;
    
    bool checkConvergence(const core::Vector& runningMeans, core::Real tolerance = 1e-4) const;
};

inline PathGenerator::PathGenerator(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)),
      rngStorage_([]() { return std::make_unique<std::mt19937>(std::random_device{}()); }) {}

inline MonteCarloEngine::MonteCarloEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)),
      pathGenerator_(std::make_shared<PathGenerator>(engine)) {}

inline core::Real MonteCarloEngine::calculateEuropeanPayoff(core::Real finalSpot, const core::MarketData& market) const {
    if (market.isCall) {
        return std::max(0.0, finalSpot - market.strike);
    } else {
        return std::max(0.0, market.strike - finalSpot);
    }
}

inline core::Real MonteCarloEngine::calculateAsianPayoff(const core::Vector& path, const core::MarketData& market,
                                                        bool isArithmetic) const {
    if (path.empty()) return 0.0;
    
    core::Real average;
    if (isArithmetic) {
        average = core::Statistics::mean(path.begin(), path.end());
    } else {
        core::Real logSum = 0.0;
        for (core::Real price : path) {
            logSum += std::log(std::max(price, core::EPSILON));
        }
        average = std::exp(logSum / path.size());
    }
    
    if (market.isCall) {
        return std::max(0.0, average - market.strike);
    } else {
        return std::max(0.0, market.strike - average);
    }
}

inline core::Real MonteCarloEngine::calculateLookbackPayoff(const core::Vector& path, const core::MarketData& market,
                                                           bool isFloating) const {
    if (path.empty()) return 0.0;
    
    const core::Real maxPrice = *std::max_element(path.begin(), path.end());
    const core::Real minPrice = *std::min_element(path.begin(), path.end());
    const core::Real finalPrice = path.back();
    
    if (isFloating) {
        if (market.isCall) {
            return finalPrice - minPrice;
        } else {
            return maxPrice - finalPrice;
        }
    } else {
        if (market.isCall) {
            return std::max(0.0, maxPrice - market.strike);
        } else {
            return std::max(0.0, market.strike - minPrice);
        }
    }
}

}