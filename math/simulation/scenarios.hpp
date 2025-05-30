#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include "../core/threading.hpp"
#include <memory>
#include <random>
#include <vector>

namespace math::simulation {

struct ScenarioParams {
    core::Real baseSpot = 100.0;
    core::Real baseVol = 0.2;
    core::Real baseRate = 0.05;
    core::Real timeHorizon = 1.0;
    core::Integer numScenarios = 10000;
    bool useAntitheticVariates = true;
    uint32_t seed = 12345;
};

struct MacroFactors {
    core::Real interestRate = 0.05;
    core::Real inflation = 0.02;
    core::Real gdpGrowth = 0.03;
    core::Real unemployment = 0.04;
    core::Real volatilityIndex = 0.15;
    core::Real creditSpread = 0.01;
    core::Real fxRate = 1.0;
};

struct MarketScenario {
    core::Real spotPrice = 100.0;
    core::Real volatility = 0.2;
    core::Real riskFreeRate = 0.05;
    core::Real dividendYield = 0.0;
    MacroFactors macroFactors;
    core::Real scenarioProbability = 1.0;
    std::string description;
};

class ScenarioGenerator {
private:
    std::shared_ptr<core::MathEngine> engine_;
    mutable core::ThreadLocalStorage<std::mt19937> rngStorage_;
    
public:
    explicit ScenarioGenerator(std::shared_ptr<core::MathEngine> engine);
    
    std::vector<MarketScenario> generateHistoricalScenarios(
        const ScenarioParams& params) const;
    
    std::vector<MarketScenario> generateMonteCarloScenarios(
        const ScenarioParams& params) const;
    
    std::vector<MarketScenario> generateStressScenarios(
        const ScenarioParams& params) const;
    
    std::vector<MarketScenario> generateTailRiskScenarios(
        const ScenarioParams& params, core::Real confidence = 0.01) const;
    
    struct VolatilityRegime {
        core::Real meanVol = 0.2;
        core::Real volOfVol = 0.05;
        core::Real persistence = 0.9;
        core::Real meanReversion = 0.1;
        core::Real probability = 1.0;
        std::string description;
    };
    
    std::vector<MarketScenario> generateRegimeSwitchingScenarios(
        const std::vector<VolatilityRegime>& regimes,
        const ScenarioParams& params) const;
    
    struct CopulaParams {
        core::Matrix correlationMatrix;
        std::string copulaType = "gaussian";
        core::Vector degreesFreedom;
    };
    
    std::vector<std::vector<MarketScenario>> generateMultiAssetScenarios(
        const std::vector<ScenarioParams>& assetParams,
        const CopulaParams& copulaParams) const;
    
    struct ClimateScenario {
        core::Real temperatureChange = 0.0;
        core::Real carbonPrice = 0.0;
        core::Real transitionRisk = 0.0;
        core::Real physicalRisk = 0.0;
        core::Real regulatoryRisk = 0.0;
    };
    
    std::vector<MarketScenario> generateClimateScenarios(
        const ScenarioParams& params,
        const std::vector<ClimateScenario>& climateInputs) const;
    
    struct ScenarioStatistics {
        core::Real meanSpot = 0.0;
        core::Real meanVol = 0.0;
        core::Real meanRate = 0.0;
        core::Real spotVolatility = 0.0;
        core::Real volVolatility = 0.0;
        core::Real rateVolatility = 0.0;
        core::Real skewness = 0.0;
        core::Real kurtosis = 0.0;
        core::Matrix correlationMatrix;
    };
    
    ScenarioStatistics analyzeScenarios(
        const std::vector<MarketScenario>& scenarios) const;
    
    void exportScenarios(const std::vector<MarketScenario>& scenarios,
                        const std::string& filename,
                        const std::string& format = "csv") const;
    
    std::vector<MarketScenario> importScenarios(
        const std::string& filename,
        const std::string& format = "csv") const;
    
    std::vector<MarketScenario> filterScenarios(
        const std::vector<MarketScenario>& scenarios,
        const std::function<bool(const MarketScenario&)>& filter) const;
    
    void calibrateToHistoricalData(const core::Matrix& historicalReturns,
                                  ScenarioParams& params) const;
    
private:
    MarketScenario generateBaseScenario(const ScenarioParams& params) const;
    
    core::Real generateCorrelatedNormal(core::Real correlation,
                                       core::Real independent,
                                       std::mt19937& rng) const;
    
    core::Matrix generateCopulaSample(const CopulaParams& params,
                                     core::Size numSamples) const;
    
    core::Real calculateRegimeProbability(const VolatilityRegime& regime,
                                         core::Real currentVol,
                                         core::Real timeStep) const;
    
    void applyMacroFactorShocks(MarketScenario& scenario,
                               const MacroFactors& shocks) const;
    
    core::Real calculateClimateImpact(const ClimateScenario& climate,
                                     const std::string& sector) const;
};

class ScenarioAggregator {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit ScenarioAggregator(std::shared_ptr<core::MathEngine> engine);
    
    struct AggregationResult {
        core::Vector expectedValues;
        core::Vector standardErrors;
        core::Vector percentiles5;
        core::Vector percentiles95;
        core::Vector worstCase;
        core::Vector bestCase;
        core::Matrix scenarioContributions;
    };
    
    AggregationResult aggregateScenarios(
        const std::vector<MarketScenario>& scenarios,
        const std::function<core::Real(const MarketScenario&)>& valueFunction) const;
    
    struct PortfolioResult {
        core::Real expectedPnl = 0.0;
        core::Real var95 = 0.0;
        core::Real var99 = 0.0;
        core::Real expectedShortfall95 = 0.0;
        core::Real expectedShortfall99 = 0.0;
        core::Vector scenarioPnl;
        core::Vector scenarioProbabilities;
    };
    
    PortfolioResult evaluatePortfolio(
        const std::vector<MarketScenario>& scenarios,
        const std::function<core::Real(const MarketScenario&)>& portfolioValuation) const;
    
    core::Matrix calculateScenarioWeights(
        const std::vector<MarketScenario>& scenarios,
        const std::string& weightingScheme = "equal") const;
    
    std::vector<MarketScenario> reduceScenarios(
        const std::vector<MarketScenario>& scenarios,
        core::Size targetSize,
        const std::string& method = "kmeans") const;
    
    core::Real calculateScenarioDistance(const MarketScenario& s1,
                                        const MarketScenario& s2) const;
    
private:
    std::vector<core::Size> performKMeansClustering(
        const std::vector<MarketScenario>& scenarios,
        core::Size numClusters) const;
    
    MarketScenario calculateClusterCentroid(
        const std::vector<MarketScenario>& scenarios,
        const std::vector<core::Size>& clusterIndices) const;
};

inline ScenarioGenerator::ScenarioGenerator(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)),
      rngStorage_([]() { return std::make_unique<std::mt19937>(std::random_device{}()); }) {}

inline ScenarioAggregator::ScenarioAggregator(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

}