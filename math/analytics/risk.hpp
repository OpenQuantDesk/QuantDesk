#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include <memory>

namespace math::analytics {

struct StressTestScenarios {
    core::Vector shocks;
    core::Vector probabilities;
    std::vector<std::string> descriptions;
};

struct StressTestResults {
    core::Vector scenarioResults;
    core::Vector scenarioProbabilities;
    core::Real worstCase = 0.0;
    core::Real bestCase = 0.0;
    core::Real expectedOutcome = 0.0;
    core::Real var95 = 0.0;
    core::Real var99 = 0.0;
};

struct CorrelationAnalysis {
    core::Matrix correlationMatrix;
    core::Real averageCorrelation = 0.0;
    core::Real maxCorrelation = 0.0;
    core::Real minCorrelation = 0.0;
    core::Vector eigenvalues;
    core::Matrix eigenvectors;
};

struct LiquidityRisk {
    core::Real averageVolume = 0.0;
    core::Real averageSpread = 0.0;
    core::Real volumeVolatility = 0.0;
    core::Real spreadVolatility = 0.0;
    core::Real illiquidityRatio = 0.0;
    core::Real liquidityScore = 0.0;
    core::Real timeToLiquidate = 0.0;
    core::Real liquidityCost = 0.0;
};

class RiskAnalyzer {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit RiskAnalyzer(std::shared_ptr<core::MathEngine> engine);
    
    core::RiskMetrics analyze(const core::Vector& returns) const;
    
    core::RiskMetrics analyzeWithBenchmark(const core::Vector& returns,
                                          const core::Vector& benchmarkReturns) const;
    
    core::Real calculateVaR(const core::Vector& returns, core::Real confidence = 0.95) const;
    
    core::Real calculateExpectedShortfall(const core::Vector& returns, 
                                         core::Real confidence = 0.95) const;
    
    core::Real calculateMaxDrawdown(const core::Vector& returns) const;
    
    core::Real calculateParametricVaR(const core::Vector& returns, 
                                     core::Real confidence = 0.95) const;
    
    core::Real calculateHistoricalVaR(const core::Vector& returns,
                                     core::Real confidence = 0.95) const;
    
    core::Real calculateModifiedVaR(const core::Vector& returns,
                                   core::Real confidence = 0.95) const;
    
    StressTestResults stressTest(const core::Vector& returns,
                                const StressTestScenarios& scenarios) const;
    
    CorrelationAnalysis analyzeCorrelations(const core::Matrix& returns) const;
    
    LiquidityRisk assessLiquidityRisk(const core::Vector& volumes,
                                     const core::Vector& spreads,
                                     const core::Vector& prices) const;
    
    struct ComponentVaR {
        core::Vector marginalVaR;
        core::Vector componentVaR;
        core::Vector contributions;
        core::Real portfolioVaR = 0.0;
    };
    
    ComponentVaR calculateComponentVaR(const core::Matrix& returns,
                                      const core::Vector& weights,
                                      core::Real confidence = 0.95) const;
    
    struct BacktestResults {
        core::Vector violations;
        core::Real violationRate = 0.0;
        core::Real expectedViolationRate = 0.0;
        core::Real kupiecLR = 0.0;
        core::Real christoffersenLR = 0.0;
        bool passesKupiec = false;
        bool passesChristoffersen = false;
    };
    
    BacktestResults backtestVaR(const core::Vector& returns,
                               const core::Vector& varForecasts,
                               core::Real confidence = 0.95) const;
    
    struct TailRisk {
        core::Real tailRatio = 0.0;
        core::Real tailExpectation = 0.0;
        core::Real tailVariance = 0.0;
        core::Real extremeValueIndex = 0.0;
    };
    
    TailRisk analyzeTailRisk(const core::Vector& returns) const;
    
private:
    core::Real calculatePortfolioVariance(const core::Matrix& returns,
                                         const core::Vector& weights) const;
    
    core::Vector calculateMarginalContributions(const core::Matrix& returns,
                                               const core::Vector& weights) const;
};

inline RiskAnalyzer::RiskAnalyzer(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

}