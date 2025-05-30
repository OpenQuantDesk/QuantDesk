#include "../analytics/risk.hpp"
#include "../core/utils.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace math::analytics {

RiskAnalyzer::RiskAnalyzer(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::RiskMetrics RiskAnalyzer::analyze(const core::Vector& returns) const {
    core::RiskMetrics metrics;
    
    if (returns.empty()) return metrics;
    
    core::Vector sortedReturns = returns;
    std::sort(sortedReturns.begin(), sortedReturns.end());
    
    metrics.var95 = core::Statistics::percentile(sortedReturns.begin(), sortedReturns.end(), 0.05);
    metrics.var99 = core::Statistics::percentile(sortedReturns.begin(), sortedReturns.end(), 0.01);
    
    const core::Size tail95Count = static_cast<core::Size>(sortedReturns.size() * 0.05);
    if (tail95Count > 0) {
        metrics.cvar95 = core::Statistics::mean(sortedReturns.begin(), 
                                              sortedReturns.begin() + tail95Count);
    }
    
    const core::Size tail99Count = static_cast<core::Size>(sortedReturns.size() * 0.01);
    if (tail99Count > 0) {
        metrics.cvar99 = core::Statistics::mean(sortedReturns.begin(),
                                              sortedReturns.begin() + tail99Count);
    }
    
    metrics.maxDrawdown = calculateMaxDrawdown(returns);
    
    const core::Real mean = core::Statistics::mean(returns.begin(), returns.end());
    const core::Real stddev = core::Statistics::standardDeviation(returns.begin(), returns.end(), mean);
    
    metrics.sharpeRatio = stddev > 0 ? (mean * 252.0) / (stddev * std::sqrt(252.0)) : 0.0;
    
    core::Vector downsideReturns;
    std::copy_if(returns.begin(), returns.end(), std::back_inserter(downsideReturns),
                [](core::Real r) { return r < 0; });
    
    if (!downsideReturns.empty()) {
        const core::Real downsideStddev = core::Statistics::standardDeviation(downsideReturns.begin(),
                                                                             downsideReturns.end(), 0.0);
        metrics.sortinoRatio = downsideStddev > 0 ? (mean * 252.0) / (downsideStddev * std::sqrt(252.0)) : 0.0;
    }
    
    if (std::abs(metrics.maxDrawdown) > 1e-10) {
        metrics.calmarRatio = (mean * 252.0) / std::abs(metrics.maxDrawdown);
    }
    
    return metrics;
}

core::RiskMetrics RiskAnalyzer::analyzeWithBenchmark(const core::Vector& returns,
                                                    const core::Vector& benchmarkReturns) const {
    core::RiskMetrics metrics = analyze(returns);
    
    if (returns.size() != benchmarkReturns.size() || returns.empty()) {
        return metrics;
    }
    
    const core::Real portfolioMean = core::Statistics::mean(returns.begin(), returns.end());
    const core::Real benchmarkMean = core::Statistics::mean(benchmarkReturns.begin(), benchmarkReturns.end());
    const core::Real benchmarkVariance = core::Statistics::variance(benchmarkReturns.begin(),
                                                                   benchmarkReturns.end(), benchmarkMean);
    
    if (benchmarkVariance > 0) {
        core::Real covariance = 0.0;
        for (core::Size i = 0; i < returns.size(); ++i) {
            covariance += (returns[i] - portfolioMean) * (benchmarkReturns[i] - benchmarkMean);
        }
        covariance /= (returns.size() - 1);
        
        metrics.beta = covariance / benchmarkVariance;
        metrics.alpha = (portfolioMean - benchmarkMean) * 252.0;
    }
    
    core::Vector excessReturns(returns.size());
    std::transform(returns.begin(), returns.end(), benchmarkReturns.begin(),
                  excessReturns.begin(), std::minus<core::Real>());
    
    const core::Real trackingErrorVariance = core::Statistics::variance(excessReturns.begin(),
                                                                       excessReturns.end(),
                                                                       core::Statistics::mean(excessReturns.begin(),
                                                                                             excessReturns.end()));
    metrics.trackingError = std::sqrt(trackingErrorVariance) * std::sqrt(252.0);
    
    if (metrics.trackingError > 0) {
        metrics.informationRatio = metrics.alpha / metrics.trackingError;
    }
    
    return metrics;
}

core::Real RiskAnalyzer::calculateVaR(const core::Vector& returns, core::Real confidence) const {
    if (returns.empty() || confidence <= 0.0 || confidence >= 1.0) return 0.0;
    
    core::Vector sortedReturns = returns;
    std::sort(sortedReturns.begin(), sortedReturns.end());
    
    return core::Statistics::percentile(sortedReturns.begin(), sortedReturns.end(), 1.0 - confidence);
}

core::Real RiskAnalyzer::calculateExpectedShortfall(const core::Vector& returns, core::Real confidence) const {
    if (returns.empty() || confidence <= 0.0 || confidence >= 1.0) return 0.0;
    
    core::Vector sortedReturns = returns;
    std::sort(sortedReturns.begin(), sortedReturns.end());
    
    const core::Size tailCount = static_cast<core::Size>(sortedReturns.size() * (1.0 - confidence));
    if (tailCount == 0) return sortedReturns[0];
    
    return core::Statistics::mean(sortedReturns.begin(), sortedReturns.begin() + tailCount);
}

core::Real RiskAnalyzer::calculateMaxDrawdown(const core::Vector& returns) const {
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

core::Real RiskAnalyzer::calculateParametricVaR(const core::Vector& returns, core::Real confidence) const {
    if (returns.empty()) return 0.0;
    
    const core::Real mean = core::Statistics::mean(returns.begin(), returns.end());
    const core::Real stddev = core::Statistics::standardDeviation(returns.begin(), returns.end(), mean);
    
    const core::Real zScore = engine_->normalInverse(1.0 - confidence);
    return mean + zScore * stddev;
}

core::Real RiskAnalyzer::calculateHistoricalVaR(const core::Vector& returns, core::Real confidence) const {
    return calculateVaR(returns, confidence);
}

core::Real RiskAnalyzer::calculateModifiedVaR(const core::Vector& returns, core::Real confidence) const {
    if (returns.empty()) return 0.0;
    
    const core::Real mean = core::Statistics::mean(returns.begin(), returns.end());
    const core::Real stddev = core::Statistics::standardDeviation(returns.begin(), returns.end(), mean);
    const core::Real skew = core::Statistics::skewness(returns.begin(), returns.end(), mean, stddev);
    const core::Real kurt = core::Statistics::kurtosis(returns.begin(), returns.end(), mean, stddev);
    
    const core::Real zScore = engine_->normalInverse(1.0 - confidence);
    
    const core::Real modifiedZ = zScore + 
                                (std::pow(zScore, 2) - 1.0) * skew / 6.0 +
                                (std::pow(zScore, 3) - 3.0 * zScore) * kurt / 24.0 -
                                (2.0 * std::pow(zScore, 3) - 5.0 * zScore) * std::pow(skew, 2) / 36.0;
    
    return mean + modifiedZ * stddev;
}

StressTestResults RiskAnalyzer::stressTest(const core::Vector& returns,
                                         const StressTestScenarios& scenarios) const {
    StressTestResults results;
    
    if (returns.empty() || scenarios.shocks.empty()) return results;
    
    const core::Real baseValue = std::accumulate(returns.begin(), returns.end(), 0.0);
    
    results.scenarioResults.reserve(scenarios.shocks.size());
    results.scenarioProbabilities = scenarios.probabilities;
    
    for (core::Size i = 0; i < scenarios.shocks.size(); ++i) {
        const core::Real shock = scenarios.shocks[i];
        const core::Real stressedValue = baseValue * (1.0 + shock);
        const core::Real pnl = stressedValue - baseValue;
        
        results.scenarioResults.push_back(pnl);
    }
    
    if (!results.scenarioResults.empty()) {
        results.worstCase = *std::min_element(results.scenarioResults.begin(), results.scenarioResults.end());
        results.bestCase = *std::max_element(results.scenarioResults.begin(), results.scenarioResults.end());
        
        if (!scenarios.probabilities.empty() && scenarios.probabilities.size() == scenarios.shocks.size()) {
            results.expectedOutcome = std::inner_product(results.scenarioResults.begin(),
                                                        results.scenarioResults.end(),
                                                        scenarios.probabilities.begin(), 0.0);
        } else {
            results.expectedOutcome = core::Statistics::mean(results.scenarioResults.begin(),
                                                            results.scenarioResults.end());
        }
        
        core::Vector sortedResults = results.scenarioResults;
        std::sort(sortedResults.begin(), sortedResults.end());
        
        results.var95 = core::Statistics::percentile(sortedResults.begin(), sortedResults.end(), 0.05);
        results.var99 = core::Statistics::percentile(sortedResults.begin(), sortedResults.end(), 0.01);
    }
    
    return results;
}

CorrelationAnalysis RiskAnalyzer::analyzeCorrelations(const core::Matrix& returns) const {
    CorrelationAnalysis analysis;
    
    if (returns.empty() || returns[0].empty()) return analysis;
    
    const core::Size numAssets = returns.size();
    const core::Size numObs = returns[0].size();
    
    analysis.correlationMatrix.resize(numAssets, core::Vector(numAssets, 0.0));
    analysis.averageCorrelation = 0.0;
    analysis.maxCorrelation = -1.0;
    analysis.minCorrelation = 1.0;
    
    core::Vector means(numAssets);
    for (core::Size i = 0; i < numAssets; ++i) {
        means[i] = core::Statistics::mean(returns[i].begin(), returns[i].end());
    }
    
    core::Size correlationCount = 0;
    core::Real correlationSum = 0.0;
    
    for (core::Size i = 0; i < numAssets; ++i) {
        for (core::Size j = 0; j < numAssets; ++j) {
            if (i == j) {
                analysis.correlationMatrix[i][j] = 1.0;
                continue;
            }
            
            core::Real numerator = 0.0;
            core::Real denomI = 0.0;
            core::Real denomJ = 0.0;
            
            for (core::Size k = 0; k < numObs; ++k) {
                const core::Real devI = returns[i][k] - means[i];
                const core::Real devJ = returns[j][k] - means[j];
                
                numerator += devI * devJ;
                denomI += devI * devI;
                denomJ += devJ * devJ;
            }
            
            const core::Real correlation = (denomI > 0 && denomJ > 0) ? 
                                         numerator / std::sqrt(denomI * denomJ) : 0.0;
            
            analysis.correlationMatrix[i][j] = correlation;
            
            if (i < j) {
                correlationSum += correlation;
                ++correlationCount;
                
                if (correlation > analysis.maxCorrelation) {
                    analysis.maxCorrelation = correlation;
                }
                if (correlation < analysis.minCorrelation) {
                    analysis.minCorrelation = correlation;
                }
            }
        }
    }
    
    if (correlationCount > 0) {
        analysis.averageCorrelation = correlationSum / correlationCount;
    }
    
    return analysis;
}

LiquidityRisk RiskAnalyzer::assessLiquidityRisk(const core::Vector& volumes,
                                               const core::Vector& spreads,
                                               const core::Vector& prices) const {
    LiquidityRisk risk;
    
    if (volumes.empty() || spreads.empty() || prices.empty()) return risk;
    
    const core::Size n = std::min({volumes.size(), spreads.size(), prices.size()});
    
    risk.averageVolume = core::Statistics::mean(volumes.begin(), volumes.begin() + n);
    risk.averageSpread = core::Statistics::mean(spreads.begin(), spreads.begin() + n);
    
    const core::Real volumeStddev = core::Statistics::standardDeviation(volumes.begin(), volumes.begin() + n,
                                                                       risk.averageVolume);
    const core::Real spreadStddev = core::Statistics::standardDeviation(spreads.begin(), spreads.begin() + n,
                                                                       risk.averageSpread);
    
    risk.volumeVolatility = volumeStddev / risk.averageVolume;
    risk.spreadVolatility = spreadStddev / risk.averageSpread;
    
    core::Vector returns(n - 1);
    for (core::Size i = 1; i < n; ++i) {
        returns[i - 1] = (prices[i] - prices[i - 1]) / prices[i - 1];
    }
    
    core::Vector absReturns(returns.size());
    std::transform(returns.begin(), returns.end(), absReturns.begin(),
                  [](core::Real r) { return std::abs(r); });
    
    if (!absReturns.empty() && risk.averageVolume > 0) {
        const core::Real avgAbsReturn = core::Statistics::mean(absReturns.begin(), absReturns.end());
        risk.illiquidityRatio = avgAbsReturn / risk.averageVolume;
    }
    
    const core::Real volumeThreshold = risk.averageVolume * 0.5;
    const auto lowVolumeCount = std::count_if(volumes.begin(), volumes.begin() + n,
                                            [volumeThreshold](core::Real v) { return v < volumeThreshold; });
    risk.liquidityScore = 1.0 - static_cast<core::Real>(lowVolumeCount) / n;
    
    return risk;
}

}