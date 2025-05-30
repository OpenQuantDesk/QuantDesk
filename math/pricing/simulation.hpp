#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include "../simulation/paths.hpp"
#include <memory>
#include <functional>

namespace math::pricing {

class MonteCarloEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    std::shared_ptr<simulation::PathGenerator> pathGenerator_;
    
public:
    explicit MonteCarloEngine(std::shared_ptr<core::MathEngine> engine);
    
    struct SimulationParams {
        core::Integer numSimulations = 100000;
        core::Integer timeSteps = 252;
        bool useAntitheticVariates = true;
        bool useControlVariates = false;
        core::Real controlVariateCorrelation = 0.9;
        uint32_t seed = 12345;
        bool enableProgressReporting = false;
    };
    
    struct SimulationResult {
        core::Real price = 0.0;
        core::Real standardError = 0.0;
        core::Real confidence95Lower = 0.0;
        core::Real confidence95Upper = 0.0;
        core::Real probabilityProfit = 0.0;
        core::Vector payoffs;
        core::Integer convergenceIterations = 0;
        core::Real convergenceError = 0.0;
    };
    
    SimulationResult priceEuropean(const core::MarketData& market,
                                  const SimulationParams& params = SimulationParams{}) const;
    
    SimulationResult priceAmerican(const core::MarketData& market,
                                  const SimulationParams& params = SimulationParams{}) const;
    
    SimulationResult priceAsian(const core::MarketData& market,
                               const SimulationParams& params = SimulationParams{},
                               bool isArithmetic = true) const;
    
    SimulationResult priceLookback(const core::MarketData& market,
                                  const SimulationParams& params = SimulationParams{},
                                  bool isFloating = true) const;
    
    SimulationResult priceBarrier(const core::MarketData& market,
                                 const SimulationParams& params,
                                 core::Real barrier, bool isKnockOut = true,
                                 bool isUp = true, core::Real rebate = 0.0) const;
    
    SimulationResult priceCustomPayoff(const core::MarketData& market,
                                      const std::function<core::Real(const core::Vector&)>& payoffFunc,
                                      const SimulationParams& params = SimulationParams{}) const;
    
    struct VarianceReductionResult {
        core::Real standardMC = 0.0;
        core::Real antitheticMC = 0.0;
        core::Real controlVariateMC = 0.0;
        core::Real importanceSamplingMC = 0.0;
        core::Real varianceReduction = 0.0;
        core::Real optimalAlpha = 0.0;
    };
    
    VarianceReductionResult compareVarianceReduction(const core::MarketData& market,
                                                    const SimulationParams& params) const;
    
    struct AdaptiveParams {
        core::Real targetAccuracy = 1e-4;
        core::Integer minSimulations = 10000;
        core::Integer maxSimulations = 1000000;
        core::Integer batchSize = 10000;
        core::Real confidenceLevel = 0.95;
    };
    
    SimulationResult adaptivePricing(const core::MarketData& market,
                                    const std::function<core::Real(const core::Vector&)>& payoffFunc,
                                    const AdaptiveParams& params = AdaptiveParams{}) const;
    
    struct SensitivityResult {
        core::Real delta = 0.0;
        core::Real gamma = 0.0;
        core::Real vega = 0.0;
        core::Real theta = 0.0;
        core::Real rho = 0.0;
        core::Vector deltaStandardErrors;
        core::Vector gammaStandardErrors;
    };
    
    SensitivityResult calculateSensitivities(const core::MarketData& market,
                                            const SimulationParams& params = SimulationParams{}) const;
    
private:
    core::Real calculateControlVariateAdjustment(const core::Vector& payoffs,
                                                 const core::Vector& controlPayoffs,
                                                 core::Real controlPrice) const;
    
    bool checkConvergence(const core::Vector& runningMeans, core::Real tolerance) const;
    
    core::Real calculateImportanceSamplingWeight(core::Real spot, core::Real drift,
                                                 core::Real newDrift) const;
};

class QuasiMonteCarloEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit QuasiMonteCarloEngine(std::shared_ptr<core::MathEngine> engine);
    
    struct QMCParams {
        core::Size numSimulations = 65536;
        core::Size timeSteps = 252;
        bool useScrambling = true;
        bool useBrownianBridge = true;
        core::Integer seed = 12345;
    };
    
    struct QMCResult {
        core::Real price = 0.0;
        core::Real standardError = 0.0;
        core::Vector payoffs;
        core::Real convergenceRate = 0.0;
        core::Real discrepancy = 0.0;
    };
    
    QMCResult priceSobol(const core::MarketData& market,
                        const QMCParams& params = QMCParams{}) const;
    
    QMCResult priceHalton(const core::MarketData& market,
                         const QMCParams& params = QMCParams{}) const;
    
    QMCResult priceFaure(const core::MarketData& market,
                        const QMCParams& params = QMCParams{}) const;
    
    QMCResult priceNiederreiter(const core::MarketData& market,
                               const QMCParams& params = QMCParams{}) const;
    
    struct QMCComparison {
        QMCResult sobolResult;
        QMCResult haltonResult;
        QMCResult faureResult;
        std::string bestMethod;
        core::Real improvementRatio = 0.0;
    };
    
    QMCComparison compareQMCMethods(const core::MarketData& market,
                                   const QMCParams& params = QMCParams{}) const;
    
    QMCResult priceMultiAsset(const std::vector<core::MarketData>& markets,
                             const std::function<core::Real(const core::Vector&)>& payoffFunc,
                             const QMCParams& params = QMCParams{}) const;
    
private:
    core::Matrix generateSobolSequence(core::Size numPoints, core::Size dimensions) const;
    core::Matrix generateHaltonSequence(core::Size numPoints, core::Size dimensions) const;
    core::Matrix generateFaureSequence(core::Size numPoints, core::Size dimensions) const;
    
    core::Real sobolPoint(core::Size index, core::Size dimension) const;
    core::Real haltonPoint(core::Size index, core::Integer base) const;
    core::Real faurePoint(core::Size index, core::Size dimension, core::Integer base) const;
    
    core::Matrix applyScrambling(const core::Matrix& sequence) const;
    core::Matrix applyBrownianBridge(const core::Matrix& sequence) const;
    
    core::Real calculateConvergenceRate(const core::Vector& payoffs) const;
    core::Real calculateDiscrepancy(const core::Matrix& sequence) const;
};

class AccelerationEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit AccelerationEngine(std::shared_ptr<core::MathEngine> engine);
    
    struct AccelerationParams {
        bool useMultilevelMC = true;
        bool useMultiindexMC = false;
        core::Integer numLevels = 5;
        core::Real levelRatio = 2.0;
        core::Real toleranceRatio = 0.1;
        core::Integer minSamplesPerLevel = 1000;
    };
    
    struct MLMCResult {
        core::Real price = 0.0;
        core::Real variance = 0.0;
        core::Vector levelResults;
        core::Vector levelVariances;
        core::Vector optimalSamples;
        core::Real totalCost = 0.0;
        core::Real speedup = 0.0;
    };
    
    MLMCResult multilevelMonteCarlo(const core::MarketData& market,
                                   const std::function<core::Real(const core::Vector&, core::Size)>& payoffFunc,
                                   const AccelerationParams& params = AccelerationParams{}) const;
    
    struct ControlVariateParams {
        core::Vector controlPrices;
        std::vector<std::function<core::Real(const core::Vector&)>> controlPayoffs;
        bool useOptimalWeights = true;
        core::Real correlationThreshold = 0.5;
    };
    
    MonteCarloEngine::SimulationResult enhancedControlVariates(
        const core::MarketData& market,
        const std::function<core::Real(const core::Vector&)>& payoffFunc,
        const ControlVariateParams& cvParams,
        const MonteCarloEngine::SimulationParams& simParams = MonteCarloEngine::SimulationParams{}) const;
    
    struct ImportanceSamplingParams {
        core::Real shiftParameter = 0.1;
        bool useOptimalShift = true;
        bool useExponentialTilting = true;
        core::Real adaptiveRate = 0.01;
    };
    
    MonteCarloEngine::SimulationResult importanceSampling(
        const core::MarketData& market,
        const std::function<core::Real(const core::Vector&)>& payoffFunc,
        const ImportanceSamplingParams& isParams,
        const MonteCarloEngine::SimulationParams& simParams = MonteCarloEngine::SimulationParams{}) const;
    
    struct StratificationParams {
        core::Size numStrata = 10;
        bool useOptimalAllocation = true;
        bool useProportionalAllocation = false;
        core::Real stratificationVariable = 0.0;
    };
    
    MonteCarloEngine::SimulationResult stratifiedSampling(
        const core::MarketData& market,
        const std::function<core::Real(const core::Vector&)>& payoffFunc,
        const StratificationParams& stratParams,
        const MonteCarloEngine::SimulationParams& simParams = MonteCarloEngine::SimulationParams{}) const;
    
private:
    core::Vector calculateOptimalAllocation(const core::Vector& stratumVariances,
                                           const core::Vector& stratumCosts) const;
    
    core::Real calculateOptimalShift(const core::MarketData& market,
                                    const std::function<core::Real(const core::Vector&)>& payoffFunc) const;
    
    core::Vector calculateOptimalWeights(const core::Matrix& controlPayoffs,
                                        const core::Vector& mainPayoffs) const;
    
    core::Real estimateLevelVariance(const core::MarketData& market,
                                    const std::function<core::Real(const core::Vector&, core::Size)>& payoffFunc,
                                    core::Size level, core::Integer numSamples) const;
};

class RegressionEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit RegressionEngine(std::shared_ptr<core::MathEngine> engine);
    
    struct RegressionParams {
        core::Integer polynomialOrder = 3;
        bool useLaguerre = true;
        bool useHermite = false;
        bool useChebyshev = false;
        core::Real regularizationParameter = 0.0;
        bool useRidgeRegression = false;
    };
    
    struct LSMResult {
        core::Real price = 0.0;
        core::Real standardError = 0.0;
        core::Vector exerciseBoundary;
        core::Vector exerciseTimes;
        core::Real earlyExercisePremium = 0.0;
        core::Matrix regressionCoefficients;
    };
    
    LSMResult leastSquaresMonteCarlo(const core::MarketData& market,
                                    const core::Vector& exerciseTimes,
                                    const RegressionParams& regParams = RegressionParams{},
                                    const MonteCarloEngine::SimulationParams& simParams = MonteCarloEngine::SimulationParams{}) const;
    
    LSMResult americanOption(const core::MarketData& market,
                            const RegressionParams& regParams = RegressionParams{},
                            const MonteCarloEngine::SimulationParams& simParams = MonteCarloEngine::SimulationParams{}) const;
    
    LSMResult bermudanOption(const core::MarketData& market,
                            const core::Vector& exerciseTimes,
                            const RegressionParams& regParams = RegressionParams{},
                            const MonteCarloEngine::SimulationParams& simParams = MonteCarloEngine::SimulationParams{}) const;
    
    struct PathwiseDerivativeResult {
        core::Real delta = 0.0;
        core::Real gamma = 0.0;
        core::Real vega = 0.0;
        core::Real deltaStandardError = 0.0;
        core::Real gammaStandardError = 0.0;
        core::Real vegaStandardError = 0.0;
    };
    
    PathwiseDerivativeResult pathwiseDerivatives(const core::MarketData& market,
                                                const std::function<core::Real(const core::Vector&)>& payoffFunc,
                                                const MonteCarloEngine::SimulationParams& simParams = MonteCarloEngine::SimulationParams{}) const;
    
    PathwiseDerivativeResult likelihoodRatioMethod(const core::MarketData& market,
                                                  const std::function<core::Real(const core::Vector&)>& payoffFunc,
                                                  const MonteCarloEngine::SimulationParams& simParams = MonteCarloEngine::SimulationParams{}) const;
    
private:
    core::Vector laguerreBasis(core::Real x, core::Integer order) const;
    core::Vector hermiteBasis(core::Real x, core::Integer order) const;
    core::Vector chebyshevBasis(core::Real x, core::Integer order) const;
    
    core::Matrix buildRegressionMatrix(const core::Vector& spots, const RegressionParams& params) const;
    core::Vector solveRegression(const core::Matrix& X, const core::Vector& y, 
                                const RegressionParams& params) const;
    
    core::Real calculateContinuationValue(core::Real spot, const core::Vector& coefficients,
                                         const RegressionParams& params) const;
};

inline MonteCarloEngine::MonteCarloEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)),
      pathGenerator_(std::make_shared<simulation::PathGenerator>(engine)) {}

inline QuasiMonteCarloEngine::QuasiMonteCarloEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

inline AccelerationEngine::AccelerationEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

inline RegressionEngine::RegressionEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

}