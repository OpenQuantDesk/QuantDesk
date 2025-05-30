#include "calculator.hpp"
#include "../core/utils.hpp"
#include <algorithm>

namespace math::greeks {

core::ExtendedGreeks GreeksCalculator::numericalExtended(const core::MarketData& market,
                                                        core::Real bumpSize, core::Real quantity) const {
    core::ExtendedGreeks extended;
    
    const core::Greeks base = numerical(market, bumpSize);
    extended.price = base.price * quantity;
    extended.delta = base.delta * quantity;
    extended.gamma = base.gamma * quantity;
    extended.theta = base.theta * quantity;
    extended.vega = base.vega * quantity;
    extended.rho = base.rho * quantity;
    extended.impliedVol = base.impliedVol;
    
    extended.volga = calculateVolga(market, bumpSize) * quantity;
    extended.vanna = calculateVanna(market, bumpSize) * quantity;
    extended.charm = calculateCharm(market, bumpSize) * quantity;
    extended.speed = calculateSpeed(market, bumpSize) * quantity;
    extended.zomma = calculateZomma(market, bumpSize) * quantity;
    extended.color = calculateColor(market, bumpSize) * quantity;
    
    extended.dollarDelta = calculateDollarDelta(base, market.spot, quantity);
    extended.dollarGamma = calculateDollarGamma(base, market.spot, quantity);
    extended.pinRisk = calculatePinRisk(market.spot, market.strike, market.timeToExpiry, quantity);
    
    return extended;
}

void GreeksCalculator::extendedBatch(const core::MarketDataVector& markets,
                                    const core::Vector& quantities,
                                    core::ExtendedGreeksVector& results) const {
    results.resize(markets.size());
    
    if (markets.size() < engine_->getThreadPool().getNumThreads() * 4) {
        for (core::Size i = 0; i < markets.size(); ++i) {
            const core::Real qty = i < quantities.size() ? quantities[i] : 1.0;
            results[i] = analyticalExtended(markets[i], qty);
        }
        return;
    }
    
    engine_->getThreadPool().parallelFor(markets.begin(), markets.end(),
        [this, &markets, &quantities, &results](const core::MarketData& market) {
            core::Size index = &market - &markets[0];
            const core::Real qty = index < quantities.size() ? quantities[index] : 1.0;
            results[index] = analyticalExtended(market, qty);
        });
}

void GreeksCalculator::calculateRiskLadder(const core::MarketDataVector& markets,
                                          const core::Vector& quantities,
                                          const core::Vector& strikeLadder,
                                          core::Matrix& deltaLadder,
                                          core::Matrix& gammaLadder) const {
    const core::Size numPositions = markets.size();
    const core::Size numStrikes = strikeLadder.size();
    
    deltaLadder.resize(numPositions, core::Vector(numStrikes, 0.0));
    gammaLadder.resize(numPositions, core::Vector(numStrikes, 0.0));
    
    for (core::Size i = 0; i < numPositions; ++i) {
        const core::Real quantity = i < quantities.size() ? quantities[i] : 1.0;
        
        for (core::Size j = 0; j < numStrikes; ++j) {
            core::MarketData shiftedMarket = markets[i];
            shiftedMarket.spot = strikeLadder[j];
            
            const core::Greeks greeks = analytical(shiftedMarket);
            deltaLadder[i][j] = greeks.delta * quantity;
            gammaLadder[i][j] = greeks.gamma * quantity;
        }
    }
}

GreeksCalculator::StressTestResults GreeksCalculator::stressTest(const core::MarketDataVector& markets,
                                                               const core::Vector& quantities,
                                                               const StressTestParams& params) const {
    StressTestResults results;
    
    const core::Size numPositions = markets.size();
    const core::Size numSpotShifts = params.spotShifts.size();
    const core::Size numVolShifts = params.volShifts.size();
    const core::Size numRateShifts = params.rateShifts.size();
    const core::Size numTimeDecays = params.timeDecays.size();
    
    results.spotStress.resize(numPositions, core::Vector(numSpotShifts, 0.0));
    results.volStress.resize(numPositions, core::Vector(numVolShifts, 0.0));
    results.rateStress.resize(numPositions, core::Vector(numRateShifts, 0.0));
    results.timeStress.resize(numPositions, core::Vector(numTimeDecays, 0.0));
    
    for (core::Size i = 0; i < numPositions; ++i) {
        const core::Real quantity = i < quantities.size() ? quantities[i] : 1.0;
        const core::Real basePrice = analytical(markets[i]).price * quantity;
        
        for (core::Size j = 0; j < numSpotShifts; ++j) {
            core::MarketData stressedMarket = markets[i];
            stressedMarket.spot *= (1.0 + params.spotShifts[j]);
            const core::Real stressedPrice = analytical(stressedMarket).price * quantity;
            results.spotStress[i][j] = stressedPrice - basePrice;
        }
        
        for (core::Size j = 0; j < numVolShifts; ++j) {
            core::MarketData stressedMarket = markets[i];
            stressedMarket.volatility += params.volShifts[j];
            const core::Real stressedPrice = analytical(stressedMarket).price * quantity;
            results.volStress[i][j] = stressedPrice - basePrice;
        }
        
        for (core::Size j = 0; j < numRateShifts; ++j) {
            core::MarketData stressedMarket = markets[i];
            stressedMarket.riskFreeRate += params.rateShifts[j];
            const core::Real stressedPrice = analytical(stressedMarket).price * quantity;
            results.rateStress[i][j] = stressedPrice - basePrice;
        }
        
        for (core::Size j = 0; j < numTimeDecays; ++j) {
            core::MarketData stressedMarket = markets[i];
            stressedMarket.timeToExpiry = std::max(0.001, stressedMarket.timeToExpiry - params.timeDecays[j]);
            const core::Real stressedPrice = analytical(stressedMarket).price * quantity;
            results.timeStress[i][j] = stressedPrice - basePrice;
        }
    }
    
    const core::Size totalScenarios = numSpotShifts + numVolShifts + numRateShifts + numTimeDecays;
    results.worstCase.resize(numPositions);
    results.bestCase.resize(numPositions);
    
    for (core::Size i = 0; i < numPositions; ++i) {
        core::Vector allStresses;
        allStresses.reserve(totalScenarios);
        
        allStresses.insert(allStresses.end(), results.spotStress[i].begin(), results.spotStress[i].end());
        allStresses.insert(allStresses.end(), results.volStress[i].begin(), results.volStress[i].end());
        allStresses.insert(allStresses.end(), results.rateStress[i].begin(), results.rateStress[i].end());
        allStresses.insert(allStresses.end(), results.timeStress[i].begin(), results.timeStress[i].end());
        
        if (!allStresses.empty()) {
            results.worstCase[i] = *std::min_element(allStresses.begin(), allStresses.end());
            results.bestCase[i] = *std::max_element(allStresses.begin(), allStresses.end());
        }
    }
    
    return results;
}

core::Real GreeksCalculator::calculateCharm(const core::MarketData& market, core::Real bumpSize) const {
    const core::Real timeBump = 1.0 / 365.0;
    
    if (market.timeToExpiry <= timeBump) return 0.0;
    
    const core::Real baseDelta = engine_->blackScholes(market).delta;
    
    core::MarketData timeMarket = market;
    timeMarket.timeToExpiry -= timeBump;
    const core::Real timeDelta = engine_->blackScholes(timeMarket).delta;
    
    return (baseDelta - timeDelta) / timeBump;
}

core::Real GreeksCalculator::calculateSpeed(const core::MarketData& market, core::Real bumpSize) const {
    const core::Real spotBump = market.spot * bumpSize;
    
    core::MarketData upMarket = market;
    upMarket.spot += spotBump;
    const core::Real upGamma = engine_->blackScholes(upMarket).gamma;
    
    core::MarketData downMarket = market;
    downMarket.spot -= spotBump;
    const core::Real downGamma = engine_->blackScholes(downMarket).gamma;
    
    const core::Real baseGamma = engine_->blackScholes(market).gamma;
    
    return (upGamma - 2.0 * baseGamma + downGamma) / (spotBump * spotBump);
}

core::Real GreeksCalculator::calculateZomma(const core::MarketData& market, core::Real bumpSize) const {
    const core::Real volBump = 0.01;
    
    core::MarketData upMarket = market;
    upMarket.volatility += volBump;
    const core::Real upGamma = engine_->blackScholes(upMarket).gamma;
    
    core::MarketData downMarket = market;
    downMarket.volatility -= volBump;
    const core::Real downGamma = engine_->blackScholes(downMarket).gamma;
    
    return (upGamma - downGamma) / (2.0 * volBump);
}

core::Real GreeksCalculator::calculateColor(const core::MarketData& market, core::Real bumpSize) const {
    const core::Real timeBump = 1.0 / 365.0;
    
    if (market.timeToExpiry <= timeBump) return 0.0;
    
    const core::Real baseGamma = engine_->blackScholes(market).gamma;
    
    core::MarketData timeMarket = market;
    timeMarket.timeToExpiry -= timeBump;
    const core::Real timeGamma = engine_->blackScholes(timeMarket).gamma;
    
    return (baseGamma - timeGamma) / timeBump;
}

RiskEngine::PortfolioRisk RiskEngine::analyzePortfolio(const core::MarketDataVector& markets,
                                                      const core::Vector& quantities,
                                                      const core::Vector& historicalReturns) const {
    PortfolioRisk risk;
    
    core::GreeksVector greeks;
    calculator_->analyticalBatch(markets, greeks);
    
    for (core::Size i = 0; i < greeks.size() && i < quantities.size(); ++i) {
        const core::Real qty = quantities[i];
        risk.totalDelta += greeks[i].delta * qty;
        risk.totalGamma += greeks[i].gamma * qty;
        risk.totalVega += greeks[i].vega * qty;
        risk.totalTheta += greeks[i].theta * qty;
        risk.totalRho += greeks[i].rho * qty;
    }
    
    if (!historicalReturns.empty()) {
        core::Matrix covarianceMatrix = buildCovarianceMatrix(historicalReturns);
        
        risk.var95 = calculateVaR(greeks, quantities, covarianceMatrix, 0.95);
        risk.var99 = calculateVaR(greeks, quantities, covarianceMatrix, 0.99);
        risk.expectedShortfall95 = calculateExpectedShortfall(greeks, quantities, covarianceMatrix, 0.95);
        risk.expectedShortfall99 = calculateExpectedShortfall(greeks, quantities, covarianceMatrix, 0.99);
        
        risk.correlationMatrix = covarianceMatrix;
    }
    
    core::Vector cumulativeReturns(historicalReturns.size());
    if (!historicalReturns.empty()) {
        std::partial_sum(historicalReturns.begin(), historicalReturns.end(), cumulativeReturns.begin());
        
        core::Real peak = cumulativeReturns[0];
        core::Real maxDrawdown = 0.0;
        
        for (core::Real value : cumulativeReturns) {
            if (value > peak) peak = value;
            const core::Real drawdown = peak - value;
            if (drawdown > maxDrawdown) maxDrawdown = drawdown;
        }
        
        risk.maxDrawdown = -maxDrawdown;
    }
    
    risk.concentrationRisk.resize(quantities.size());
    const core::Real totalNotional = std::accumulate(quantities.begin(), quantities.end(), 0.0,
                                                    [](core::Real sum, core::Real q) { return sum + std::abs(q); });
    
    if (totalNotional > 0) {
        for (core::Size i = 0; i < quantities.size(); ++i) {
            risk.concentrationRisk[i] = std::abs(quantities[i]) / totalNotional;
        }
    }
    
    return risk;
}

core::Real RiskEngine::calculateVaR(const core::GreeksVector& greeks,
                                   const core::Vector& quantities,
                                   const core::Matrix& covarianceMatrix,
                                   core::Real confidence) const {
    if (greeks.empty() || quantities.empty() || covarianceMatrix.empty()) return 0.0;
    
    core::Vector weights(quantities.size());
    for (core::Size i = 0; i < weights.size() && i < quantities.size(); ++i) {
        weights[i] = quantities[i];
    }
    
    const core::Real portfolioVariance = portfolioVariance(weights, covarianceMatrix);
    const core::Real portfolioStddev = std::sqrt(portfolioVariance);
    
    const core::Real zScore = engine_->normalInverse(1.0 - confidence);
    return -zScore * portfolioStddev;
}

core::Real RiskEngine::portfolioVariance(const core::Vector& weights, const core::Matrix& covariance) const {
    if (weights.empty() || covariance.empty() || weights.size() != covariance.size()) return 0.0;
    
    core::Real variance = 0.0;
    const core::Size n = weights.size();
    
    for (core::Size i = 0; i < n; ++i) {
        for (core::Size j = 0; j < n && j < covariance[i].size(); ++j) {
            variance += weights[i] * weights[j] * covariance[i][j];
        }
    }
    
    return variance;
}

core::Matrix RiskEngine::buildCovarianceMatrix(const core::Vector& returns, core::Size lag) const {
    if (returns.size() < lag + 1) return {};
    
    const core::Size n = returns.size() - lag;
    core::Matrix covariance(n, core::Vector(n, 0.0));
    
    for (core::Size i = 0; i < n; ++i) {
        for (core::Size j = 0; j < n; ++j) {
            if (i == j) {
                const core::Real mean = core::Statistics::mean(returns.begin() + i, returns.begin() + i + lag);
                covariance[i][j] = core::Statistics::variance(returns.begin() + i, returns.begin() + i + lag, mean);
            } else {
                core::Real covar = 0.0;
                const core::Real meanI = core::Statistics::mean(returns.begin() + i, returns.begin() + i + lag);
                const core::Real meanJ = core::Statistics::mean(returns.begin() + j, returns.begin() + j + lag);
                
                for (core::Size k = 0; k < lag; ++k) {
                    covar += (returns[i + k] - meanI) * (returns[j + k] - meanJ);
                }
                covariance[i][j] = covar / (lag - 1);
            }
        }
    }
    
    return covariance;
}

}