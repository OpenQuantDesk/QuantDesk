#include "sensitivities.hpp"
#include "../core/utils.hpp"
#include <algorithm>
#include <numeric>

namespace math::greeks {

SensitivityEngine::SensitivityEngine(std::shared_ptr<GreeksCalculator> calculator)
    : calculator_(std::move(calculator)) {}

SensitivityEngine::SensitivityProfile SensitivityEngine::calculateSensitivities(
    const core::MarketDataVector& markets, const core::Vector& quantities,
    const SensitivityParams& params) const {
    
    SensitivityProfile profile;
    
    if (markets.empty() || quantities.empty()) {
        return profile;
    }
    
    const core::Size numPositions = markets.size();
    
    profile.spotSensitivities.resize(numPositions);
    profile.volSensitivities.resize(numPositions);
    profile.rateSensitivities.resize(numPositions);
    profile.timeSensitivities.resize(numPositions);
    
    for (core::Size i = 0; i < numPositions; ++i) {
        const core::Real quantity = i < quantities.size() ? quantities[i] : 1.0;
        
        profile.spotSensitivities[i] = calculateSpotSensitivity(markets[i], quantity, params);
        profile.volSensitivities[i] = calculateVolSensitivity(markets[i], quantity, params);
        profile.rateSensitivities[i] = calculateRateSensitivity(markets[i], quantity, params);
        profile.timeSensitivities[i] = calculateTimeSensitivity(markets[i], quantity, params);
    }
    
    profile.crossSensitivities = calculateCrossSensitivities(markets, quantities, params);
    
    profile.totalSpotSensitivity = std::accumulate(profile.spotSensitivities.begin(),
                                                  profile.spotSensitivities.end(), 0.0);
    profile.totalVolSensitivity = std::accumulate(profile.volSensitivities.begin(),
                                                 profile.volSensitivities.end(), 0.0);
    profile.totalRateSensitivity = std::accumulate(profile.rateSensitivities.begin(),
                                                  profile.rateSensitivities.end(), 0.0);
    profile.totalTimeSensitivity = std::accumulate(profile.timeSensitivities.begin(),
                                                  profile.timeSensitivities.end(), 0.0);
    
    return profile;
}

SensitivityEngine::StressTestResult SensitivityEngine::runStressTest(
    const core::MarketDataVector& markets, const core::Vector& quantities,
    const StressTestParams& params) const {
    
    StressTestResult result;
    
    if (markets.empty() || quantities.empty()) {
        return result;
    }
    
    const core::Real baseValue = calculatePortfolioValue(markets, quantities);
    
    result.spotStresses.reserve(params.spotShocks.size());
    result.volStresses.reserve(params.volShocks.size());
    result.rateStresses.reserve(params.rateShocks.size());
    result.timeStresses.reserve(params.timeShocks.size());
    
    for (core::Real spotShock : params.spotShocks) {
        core::MarketDataVector stressedMarkets = markets;
        for (auto& market : stressedMarkets) {
            market.spot *= (1.0 + spotShock);
        }
        
        const core::Real stressedValue = calculatePortfolioValue(stressedMarkets, quantities);
        result.spotStresses.push_back(stressedValue - baseValue);
    }
    
    for (core::Real volShock : params.volShocks) {
        core::MarketDataVector stressedMarkets = markets;
        for (auto& market : stressedMarkets) {
            market.volatility += volShock;
            market.volatility = std::max(0.01, market.volatility);
        }
        
        const core::Real stressedValue = calculatePortfolioValue(stressedMarkets, quantities);
        result.volStresses.push_back(stressedValue - baseValue);
    }
    
    for (core::Real rateShock : params.rateShocks) {
        core::MarketDataVector stressedMarkets = markets;
        for (auto& market : stressedMarkets) {
            market.riskFreeRate += rateShock;
            market.riskFreeRate = std::max(-0.10, market.riskFreeRate);
        }
        
        const core::Real stressedValue = calculatePortfolioValue(stressedMarkets, quantities);
        result.rateStresses.push_back(stressedValue - baseValue);
    }
    
    for (core::Real timeShock : params.timeShocks) {
        core::MarketDataVector stressedMarkets = markets;
        for (auto& market : stressedMarkets) {
            market.timeToExpiry = std::max(0.001, market.timeToExpiry - timeShock);
        }
        
        const core::Real stressedValue = calculatePortfolioValue(stressedMarkets, quantities);
        result.timeStresses.push_back(stressedValue - baseValue);
    }
    
    core::Vector allStresses;
    allStresses.insert(allStresses.end(), result.spotStresses.begin(), result.spotStresses.end());
    allStresses.insert(allStresses.end(), result.volStresses.begin(), result.volStresses.end());
    allStresses.insert(allStresses.end(), result.rateStresses.begin(), result.rateStresses.end());
    allStresses.insert(allStresses.end(), result.timeStresses.begin(), result.timeStresses.end());
    
    if (!allStresses.empty()) {
        std::sort(allStresses.begin(), allStresses.end());
        
        result.worstCase = allStresses.front();
        result.bestCase = allStresses.back();
        result.var95 = core::Statistics::percentile(allStresses.begin(), allStresses.end(), 0.05);
        result.var99 = core::Statistics::percentile(allStresses.begin(), allStresses.end(), 0.01);
        
        const core::Size tail95Count = static_cast<core::Size>(allStresses.size() * 0.05);
        if (tail95Count > 0) {
            result.expectedShortfall95 = core::Statistics::mean(allStresses.begin(),
                                                               allStresses.begin() + tail95Count);
        }
        
        const core::Size tail99Count = static_cast<core::Size>(allStresses.size() * 0.01);
        if (tail99Count > 0) {
            result.expectedShortfall99 = core::Statistics::mean(allStresses.begin(),
                                                               allStresses.begin() + tail99Count);
        }
    }
    
    return result;
}

SensitivityEngine::ScenarioAnalysis SensitivityEngine::runScenarioAnalysis(
    const core::MarketDataVector& markets, const core::Vector& quantities,
    const core::Matrix& scenarios, const core::Vector& probabilities) const {
    
    ScenarioAnalysis analysis;
    
    if (markets.empty() || scenarios.empty() || scenarios[0].size() < 4) {
        return analysis;
    }
    
    const core::Real baseValue = calculatePortfolioValue(markets, quantities);
    const core::Size numScenarios = scenarios.size();
    
    analysis.scenarioValues.reserve(numScenarios);
    analysis.pnlValues.reserve(numScenarios);
    
    for (core::Size s = 0; s < numScenarios; ++s) {
        if (scenarios[s].size() < 4) continue;
        
        core::MarketDataVector scenarioMarkets = markets;
        
        const core::Real spotShock = scenarios[s][0];
        const core::Real volShock = scenarios[s][1];
        const core::Real rateShock = scenarios[s][2];
        const core::Real timeShock = scenarios[s][3];
        
        for (auto& market : scenarioMarkets) {
            market.spot *= (1.0 + spotShock);
            market.volatility += volShock;
            market.volatility = std::max(0.01, market.volatility);
            market.riskFreeRate += rateShock;
            market.riskFreeRate = std::max(-0.10, market.riskFreeRate);
            market.timeToExpiry = std::max(0.001, market.timeToExpiry - timeShock);
        }
        
        const core::Real scenarioValue = calculatePortfolioValue(scenarioMarkets, quantities);
        const core::Real pnl = scenarioValue - baseValue;
        
        analysis.scenarioValues.push_back(scenarioValue);
        analysis.pnlValues.push_back(pnl);
    }
    
    if (!analysis.pnlValues.empty()) {
        if (probabilities.size() == analysis.pnlValues.size()) {
            analysis.expectedPnl = std::inner_product(analysis.pnlValues.begin(),
                                                     analysis.pnlValues.end(),
                                                     probabilities.begin(), 0.0);
        } else {
            analysis.expectedPnl = core::Statistics::mean(analysis.pnlValues.begin(),
                                                         analysis.pnlValues.end());
        }
        
        const core::Real meanPnl = analysis.expectedPnl;
        core::Real variance = 0.0;
        
        if (probabilities.size() == analysis.pnlValues.size()) {
            for (core::Size i = 0; i < analysis.pnlValues.size(); ++i) {
                const core::Real diff = analysis.pnlValues[i] - meanPnl;
                variance += probabilities[i] * diff * diff;
            }
        } else {
            variance = core::Statistics::variance(analysis.pnlValues.begin(),
                                                 analysis.pnlValues.end(), meanPnl);
        }
        
        analysis.pnlVolatility = std::sqrt(variance);
        
        core::Vector sortedPnl = analysis.pnlValues;
        std::sort(sortedPnl.begin(), sortedPnl.end());
        
        analysis.worstCase = sortedPnl.front();
        analysis.bestCase = sortedPnl.back();
        
        const core::Size numNegative = std::count_if(analysis.pnlValues.begin(),
                                                    analysis.pnlValues.end(),
                                                    [](core::Real pnl) { return pnl < 0.0; });
        analysis.probabilityOfLoss = static_cast<core::Real>(numNegative) / analysis.pnlValues.size();
    }
    
    return analysis;
}

SensitivityEngine::CorrelationAnalysis SensitivityEngine::analyzeCorrelations(
    const core::MarketDataVector& markets, const core::Matrix& historicalReturns) const {
    
    CorrelationAnalysis analysis;
    
    if (markets.empty() || historicalReturns.empty() || historicalReturns[0].empty()) {
        return analysis;
    }
    
    const core::Size numAssets = markets.size();
    const core::Size numObs = historicalReturns[0].size();
    
    if (historicalReturns.size() != numAssets) {
        return analysis;
    }
    
    analysis.correlationMatrix.resize(numAssets, core::Vector(numAssets, 0.0));
    
    core::Vector means(numAssets);
    for (core::Size i = 0; i < numAssets; ++i) {
        means[i] = core::Statistics::mean(historicalReturns[i].begin(),
                                         historicalReturns[i].end());
    }
    
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
                const core::Real devI = historicalReturns[i][k] - means[i];
                const core::Real devJ = historicalReturns[j][k] - means[j];
                
                numerator += devI * devJ;
                denomI += devI * devI;
                denomJ += devJ * devJ;
            }
            
            if (denomI > 0.0 && denomJ > 0.0) {
                analysis.correlationMatrix[i][j] = numerator / std::sqrt(denomI * denomJ);
            }
        }
    }
    
    core::Real sumCorr = 0.0;
    core::Size corrCount = 0;
    core::Real maxCorr = -1.0;
    core::Real minCorr = 1.0;
    
    for (core::Size i = 0; i < numAssets; ++i) {
        for (core::Size j = i + 1; j < numAssets; ++j) {
            const core::Real corr = analysis.correlationMatrix[i][j];
            sumCorr += corr;
            ++corrCount;
            maxCorr = std::max(maxCorr, corr);
            minCorr = std::min(minCorr, corr);
        }
    }
    
    analysis.averageCorrelation = corrCount > 0 ? sumCorr / corrCount : 0.0;
    analysis.maxCorrelation = maxCorr;
    analysis.minCorrelation = minCorr;
    
    return analysis;
}

core::Real SensitivityEngine::calculateSpotSensitivity(const core::MarketData& market,
                                                      core::Real quantity,
                                                      const SensitivityParams& params) const {
    const core::Real spotShift = market.spot * params.spotBumpSize;
    
    core::MarketData upMarket = market;
    upMarket.spot += spotShift;
    const core::Real upPrice = calculator_->analytical(upMarket).price;
    
    core::MarketData downMarket = market;
    downMarket.spot -= spotShift;
    const core::Real downPrice = calculator_->analytical(downMarket).price;
    
    return quantity * (upPrice - downPrice) / (2.0 * spotShift);
}

core::Real SensitivityEngine::calculateVolSensitivity(const core::MarketData& market,
                                                     core::Real quantity,
                                                     const SensitivityParams& params) const {
    const core::Real basePrice = calculator_->analytical(market).price;
    
    core::MarketData volMarket = market;
    volMarket.volatility += params.volBumpSize;
    const core::Real volPrice = calculator_->analytical(volMarket).price;
    
    return quantity * (volPrice - basePrice) / params.volBumpSize;
}

core::Real SensitivityEngine::calculateRateSensitivity(const core::MarketData& market,
                                                      core::Real quantity,
                                                      const SensitivityParams& params) const {
    const core::Real basePrice = calculator_->analytical(market).price;
    
    core::MarketData rateMarket = market;
    rateMarket.riskFreeRate += params.rateBumpSize;
    const core::Real ratePrice = calculator_->analytical(rateMarket).price;
    
    return quantity * (ratePrice - basePrice) / params.rateBumpSize;
}

core::Real SensitivityEngine::calculateTimeSensitivity(const core::MarketData& market,
                                                      core::Real quantity,
                                                      const SensitivityParams& params) const {
    if (market.timeToExpiry <= params.timeBumpSize) {
        return 0.0;
    }
    
    const core::Real basePrice = calculator_->analytical(market).price;
    
    core::MarketData timeMarket = market;
    timeMarket.timeToExpiry -= params.timeBumpSize;
    const core::Real timePrice = calculator_->analytical(timeMarket).price;
    
    return quantity * (timePrice - basePrice) / params.timeBumpSize;
}

core::Matrix SensitivityEngine::calculateCrossSensitivities(const core::MarketDataVector& markets,
                                                           const core::Vector& quantities,
                                                           const SensitivityParams& params) const {
    const core::Size numMarkets = markets.size();
    core::Matrix crossSens(numMarkets, core::Vector(numMarkets, 0.0));
    
    if (numMarkets < 2) {
        return crossSens;
    }
    
    const core::Real baseValue = calculatePortfolioValue(markets, quantities);
    
    for (core::Size i = 0; i < numMarkets; ++i) {
        for (core::Size j = i; j < numMarkets; ++j) {
            if (i == j) {
                crossSens[i][j] = calculateSecondOrderSensitivity(markets[i], quantities[i], params);
            } else {
                crossSens[i][j] = calculateCrossDerivative(markets, quantities, i, j, params, baseValue);
                crossSens[j][i] = crossSens[i][j];
            }
        }
    }
    
    return crossSens;
}

core::Real SensitivityEngine::calculateSecondOrderSensitivity(const core::MarketData& market,
                                                             core::Real quantity,
                                                             const SensitivityParams& params) const {
    const core::Real spotShift = market.spot * params.spotBumpSize;
    
    const core::Real basePrice = calculator_->analytical(market).price;
    
    core::MarketData upMarket = market;
    upMarket.spot += spotShift;
    const core::Real upPrice = calculator_->analytical(upMarket).price;
    
    core::MarketData downMarket = market;
    downMarket.spot -= spotShift;
    const core::Real downPrice = calculator_->analytical(downMarket).price;
    
    return quantity * (upPrice - 2.0 * basePrice + downPrice) / (spotShift * spotShift);
}

core::Real SensitivityEngine::calculateCrossDerivative(const core::MarketDataVector& markets,
                                                      const core::Vector& quantities,
                                                      core::Size i, core::Size j,
                                                      const SensitivityParams& params,
                                                      core::Real baseValue) const {
    const core::Real spotShift1 = markets[i].spot * params.spotBumpSize;
    const core::Real spotShift2 = markets[j].spot * params.spotBumpSize;
    
    core::MarketDataVector upUpMarkets = markets;
    upUpMarkets[i].spot += spotShift1;
    upUpMarkets[j].spot += spotShift2;
    const core::Real upUpValue = calculatePortfolioValue(upUpMarkets, quantities);
    
    core::MarketDataVector upDownMarkets = markets;
    upDownMarkets[i].spot += spotShift1;
    upDownMarkets[j].spot -= spotShift2;
    const core::Real upDownValue = calculatePortfolioValue(upDownMarkets, quantities);
    
    core::MarketDataVector downUpMarkets = markets;
    downUpMarkets[i].spot -= spotShift1;
    downUpMarkets[j].spot += spotShift2;
    const core::Real downUpValue = calculatePortfolioValue(downUpMarkets, quantities);
    
    core::MarketDataVector downDownMarkets = markets;
    downDownMarkets[i].spot -= spotShift1;
    downDownMarkets[j].spot -= spotShift2;
    const core::Real downDownValue = calculatePortfolioValue(downDownMarkets, quantities);
    
    return (upUpValue - upDownValue - downUpValue + downDownValue) / 
           (4.0 * spotShift1 * spotShift2);
}

core::Real SensitivityEngine::calculatePortfolioValue(const core::MarketDataVector& markets,
                                                     const core::Vector& quantities) const {
    core::GreeksVector greeks;
    calculator_->analyticalBatch(markets, greeks);
    
    core::Real totalValue = 0.0;
    for (core::Size i = 0; i < greeks.size() && i < quantities.size(); ++i) {
        totalValue += greeks[i].price * quantities[i];
    }
    
    return totalValue;
}