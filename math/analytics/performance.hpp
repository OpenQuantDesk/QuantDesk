#pragma once

#include "../core/engine.hpp"
#include "../core/types.hpp"
#include <memory>

namespace math::analytics {

    struct PerformanceMetrics {
        core::Real totalReturn = 0.0;
        core::Real annualizedReturn = 0.0;
        core::Real annualizedVolatility = 0.0;
        core::Real sharpeRatio = 0.0;
        core::Real sortinoRatio = 0.0;
        core::Real calmarRatio = 0.0;
        core::Real maxDrawdown = 0.0;
        core::Real skewness = 0.0;
        core::Real kurtosis = 0.0;
        core::Real winRate = 0.0;
        core::Real averageWin = 0.0;
        core::Real averageLoss = 0.0;
        core::Real profitFactor = 0.0;
        core::Real var95 = 0.0;
        core::Real var99 = 0.0;
        core::Real expectedShortfall95 = 0.0;
        core::Real expectedShortfall99 = 0.0;
        core::Real alpha = 0.0;
        core::Real beta = 0.0;
        core::Real trackingError = 0.0;
        core::Real informationRatio = 0.0;
    };

    struct RollingMetrics {
        core::Vector rollingReturns;
        core::Vector rollingVolatility;
        core::Vector rollingSharpe;
        core::Vector rollingDrawdown;
    };

    struct RiskAttribution {
        core::Vector factorContributions;
        core::Vector factorRisks;
        core::Real specificRisk = 0.0;
        core::Real specificContribution = 0.0;
        core::Real totalRisk = 0.0;
    };

    class PerformanceAnalyzer {
    private:
        std::shared_ptr<core::MathEngine> engine_;
    public:
        explicit PerformanceAnalyzer(std::shared_ptr<core::MathEngine> engine);

        PerformanceMetrics analyze(const core::Vector& returns) const;

        PerformanceMetrics
        compareToBenchmark(const core::Vector& returns,
                           const core::Vector& benchmarkReturns) const;

        RollingMetrics calculateRollingMetrics(const core::Vector& returns,
                                               core::Size windowSize
                                               = 252) const;

        RiskAttribution
        analyzeRiskAttribution(const core::Vector& portfolioReturns,
                               const core::Matrix& factorReturns,
                               const core::Vector& factorLoadings) const;

        core::Real calculateMaxDrawdown(const core::Vector& returns) const;

        core::Real
        calculateUpsideCapture(const core::Vector& returns,
                               const core::Vector& benchmarkReturns) const;

        core::Real
        calculateDownsideCapture(const core::Vector& returns,
                                 const core::Vector& benchmarkReturns) const;

        core::Real calculateTreynorRatio(const core::Vector& returns,
                                         const core::Vector& benchmarkReturns,
                                         core::Real riskFreeRate = 0.0) const;

        core::Real calculateJensenAlpha(const core::Vector& returns,
                                        const core::Vector& benchmarkReturns,
                                        core::Real riskFreeRate = 0.0) const;

        struct StyleAnalysis {
            core::Vector exposures;
            core::Real rSquared = 0.0;
            core::Real trackingError = 0.0;
            core::Vector residuals;
        };

        StyleAnalysis
        performStyleAnalysis(const core::Vector& returns,
                             const core::Matrix& styleFactors) const;

        struct DrawdownAnalysis {
            core::Vector drawdowns;
            core::Vector durations;
            core::Real averageDrawdown = 0.0;
            core::Real averageDuration = 0.0;
            core::Real maxDuration = 0.0;
            core::Real recoveryTime = 0.0;
        };

        DrawdownAnalysis analyzeDrawdowns(const core::Vector& returns) const;
    private:
        core::Vector
        calculateCumulativeReturns(const core::Vector& returns) const;
        core::Vector
        calculateDrawdownSeries(const core::Vector& cumulativeReturns) const;
        core::Real calculateBeta(const core::Vector& returns,
                                 const core::Vector& benchmark) const;
    };

    inline PerformanceAnalyzer::PerformanceAnalyzer(
        std::shared_ptr<core::MathEngine> engine)
        : engine_(std::move(engine))
    {}

inline core::Real PerformanceAnalyzer::calculateUpsideCapture(