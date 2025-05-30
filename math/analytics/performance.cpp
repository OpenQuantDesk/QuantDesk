#include "../analytics/performance.hpp"
#include "../core/utils.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace math::analytics {

PerformanceAnalyzer::PerformanceAnalyzer(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

PerformanceMetrics PerformanceAnalyzer::analyze(const core::Vector& returns) const {
    PerformanceMetrics metrics;
    
    if (returns.empty()) return metrics;
    
    const core::Real mean = core::Statistics::mean(returns.begin(), returns.end());
    const core::Real variance = core::Statistics::variance(returns.begin(), returns.end(), mean);
    const core::Real stddev = std::sqrt(variance);
    
    metrics.totalReturn = std::accumulate(returns.begin(), returns.end(), 0.0);
    metrics.annualizedReturn = mean * 252.0;
    metrics.annualizedVolatility = stddev * std::sqrt(252.0);
    
    if (metrics.annualizedVolatility > 0) {
        metrics.sharpeRatio = metrics.annualizedReturn / metrics.annualizedVolatility;
    }
    
    core::Vector downsideReturns;
    std::copy_if(returns.begin(), returns.end(), std::back_inserter(downsideReturns),
                [](core::Real r) { return r < 0; });
    
    if (!downsideReturns.empty()) {
        const core::Real downsideVariance = core::Statistics::variance(downsideReturns.begin(), 
                                                                      downsideReturns.end(), 0.0);
        const core::Real downsideStddev = std::sqrt(downsideVariance) * std::sqrt(252.0);
        
        if (downsideStddev > 0) {
            metrics.sortinoRatio = metrics.annualizedReturn / downsideStddev;
        }
    }
    
    metrics.maxDrawdown = calculateMaxDrawdown(returns);
    
    if (std::abs(metrics.maxDrawdown) > 1e-10) {
        metrics.calmarRatio = metrics.annualizedReturn / std::abs(metrics.maxDrawdown);
    }
    
    metrics.skewness = core::Statistics::skewness(returns.begin(), returns.end(), mean, stddev);
    metrics.kurtosis = core::Statistics::kurtosis(returns.begin(), returns.end(), mean, stddev);
    
    metrics.winRate = static_cast<core::Real>(std::count_if(returns.begin(), returns.end(),
                                                           [](core::Real r) { return r > 0; })) / returns.size();
    
    if (metrics.winRate > 0) {
        core::Real sumWins = 0.0;
        core::Integer winCount = 0;
        for (core::Real r : returns) {
            if (r > 0) {
                sumWins += r;
                ++winCount;
            }
        }
        metrics.averageWin = winCount > 0 ? sumWins / winCount : 0.0;
    }
    
    if (metrics.winRate < 1.0) {
        core::Real sumLosses = 0.0;
        core::Integer lossCount = 0;
        for (core::Real r : returns) {
            if (r < 0) {
                sumLosses += std::abs(r);
                ++lossCount;
            }
        }
        metrics.averageLoss = lossCount > 0 ? sumLosses / lossCount : 0.0;
    }
    
    if (metrics.averageLoss > 0) {
        metrics.profitFactor = metrics.averageWin / metrics.averageLoss;
    }
    
    core::Vector sortedReturns = returns;
    std::sort(sortedReturns.begin(), sortedReturns.end());
    
    metrics.var95 = core::Statistics::percentile(sortedReturns.begin(), sortedReturns.end(), 0.05);
    metrics.var99 = core::Statistics::percentile(sortedReturns.begin(), sortedReturns.end(), 0.01);
    
    const core::Size tail95Count = static_cast<core::Size>(sortedReturns.size() * 0.05);
    if (tail95Count > 0) {
        metrics.expectedShortfall95 = core::Statistics::mean(sortedReturns.begin(), 
                                                            sortedReturns.begin() + tail95Count);
    }
    
    const core::Size tail99Count = static_cast<core::Size>(sortedReturns.size() * 0.01);
    if (tail99Count > 0) {
        metrics.expectedShortfall99 = core::Statistics::mean(sortedReturns.begin(),
                                                            sortedReturns.begin() + tail99Count);
    }
    
    return metrics;
}

PerformanceMetrics PerformanceAnalyzer::compareToBenchmark(const core::Vector& returns,
                                                          const core::Vector& benchmarkReturns) const {
    PerformanceMetrics metrics = analyze(returns);
    
    if (returns.size() != benchmarkReturns.size() || returns.empty()) {
        return metrics;
    }
    
    core::Vector excessReturns(returns.size());
    std::transform(returns.begin(), returns.end(), benchmarkReturns.begin(),
                  excessReturns.begin(), std::minus<core::Real>());
    
    const core::Real excessMean = core::Statistics::mean(excessReturns.begin(), excessReturns.end());
    const core::Real excessStddev = core::Statistics::standardDeviation(excessReturns.begin(),
                                                                       excessReturns.end(), excessMean);
    
    metrics.alpha = excessMean * 252.0;
    metrics.trackingError = excessStddev * std::sqrt(252.0);
    
    if (metrics.trackingError > 0) {
        metrics.informationRatio = metrics.alpha / metrics.trackingError;
    }
    
    const core::Real benchmarkMean = core::Statistics::mean(benchmarkReturns.begin(), benchmarkReturns.end());
    const core::Real benchmarkVariance = core::Statistics::variance(benchmarkReturns.begin(),
                                                                   benchmarkReturns.end(), benchmarkMean);
    
    if (benchmarkVariance > 0) {
        core::Real covariance = 0.0;
        for (core::Size i = 0; i < returns.size(); ++i) {
            covariance += (returns[i] - core::Statistics::mean(returns.begin(), returns.end())) *
                         (benchmarkReturns[i] - benchmarkMean);
        }
        covariance /= (returns.size() - 1);
        
        metrics.beta = covariance / benchmarkVariance;
    }
    
    return metrics;
}

core::Real PerformanceAnalyzer::calculateMaxDrawdown(const core::Vector& returns) const {
    if (returns.empty()) return 0.0;
    
    core::Vector cumulativeReturns(returns.size());
    std::partial_sum(returns.begin(), returns.end(), cumulativeReturns.begin());
    
    core::Real peak = cumulativeReturns[0];
    core::Real maxDrawdown = 0.0;
    
    for (core::Real value : cumulativeReturns) {
        if (value > peak) {
            peak = value;
        }
        
        const core::Real drawdown = peak - value;
        if (drawdown > maxDrawdown) {
            maxDrawdown = drawdown;
        }
    }
    
    return -maxDrawdown;
}

RollingMetrics PerformanceAnalyzer::calculateRollingMetrics(const core::Vector& returns,
                                                          core::Size windowSize) const {
    RollingMetrics metrics;
    
    if (returns.size() < windowSize || windowSize == 0) {
        return metrics;
    }
    
    const core::Size numWindows = returns.size() - windowSize + 1;
    metrics.rollingReturns.reserve(numWindows);
    metrics.rollingVolatility.reserve(numWindows);
    metrics.rollingSharpe.reserve(numWindows);
    metrics.rollingDrawdown.reserve(numWindows);
    
    for (core::Size i = 0; i < numWindows; ++i) {
        core::Vector window(returns.begin() + i, returns.begin() + i + windowSize);
        
        const core::Real windowReturn = std::accumulate(window.begin(), window.end(), 0.0);
        const core::Real windowMean = windowReturn / windowSize;
        const core::Real windowVariance = core::Statistics::variance(window.begin(), window.end(), windowMean);
        const core::Real windowStddev = std::sqrt(windowVariance);
        
        metrics.rollingReturns.push_back(windowReturn);
        metrics.rollingVolatility.push_back(windowStddev * std::sqrt(252.0));
        
        if (windowStddev > 0) {
            metrics.rollingSharpe.push_back((windowMean * 252.0) / (windowStddev * std::sqrt(252.0)));
        } else {
            metrics.rollingSharpe.push_back(0.0);
        }
        
        metrics.rollingDrawdown.push_back(calculateMaxDrawdown(window));
    }
    
    return metrics;
}

RiskAttribution PerformanceAnalyzer::analyzeRiskAttribution(const core::Vector& portfolioReturns,
                                                          const core::Matrix& factorReturns,
                                                          const core::Vector& factorLoadings) const {
    RiskAttribution attribution;
    
    if (portfolioReturns.empty() || factorReturns.empty() || factorLoadings.empty()) {
        return attribution;
    }
    
    const core::Size numObs = portfolioReturns.size();
    const core::Size numFactors = factorLoadings.size();
    
    if (factorReturns.size() != numFactors || factorReturns[0].size() != numObs) {
        return attribution;
    }
    
    attribution.factorContributions.resize(numFactors);
    attribution.factorRisks.resize(numFactors);
    
    const core::Real portfolioVariance = core::Statistics::variance(portfolioReturns.begin(),
                                                                   portfolioReturns.end(),
                                                                   core::Statistics::mean(portfolioReturns.begin(),
                                                                                         portfolioReturns.end()));
    
    core::Real totalFactorRisk = 0.0;
    
    for (core::Size i = 0; i < numFactors; ++i) {
        const core::Real factorMean = core::Statistics::mean(factorReturns[i].begin(), factorReturns[i].end());
        const core::Real factorVariance = core::Statistics::variance(factorReturns[i].begin(),
                                                                    factorReturns[i].end(), factorMean);
        
        const core::Real factorRisk = factorLoadings[i] * factorLoadings[i] * factorVariance;
        attribution.factorRisks[i] = factorRisk;
        totalFactorRisk += factorRisk;
        
        attribution.factorContributions[i] = factorRisk / portfolioVariance;
    }
    
    attribution.specificRisk = portfolioVariance - totalFactorRisk;
    attribution.specificContribution = attribution.specificRisk / portfolioVariance;
    attribution.totalRisk = std::sqrt(portfolioVariance);
    
    return attribution;
}

}