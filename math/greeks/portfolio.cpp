#include "portfolio.hpp"
#include "../core/utils.hpp"
#include <algorithm>
#include <numeric>

namespace math::greeks {

PortfolioManager::PortfolioManager(std::shared_ptr<GreeksCalculator> calculator)
    : calculator_(std::move(calculator)) {}

void PortfolioManager::addPosition(const Position& position) {
    core::UniqueLockGuard<core::ReadWriteLock> lock(positionsLock_);
    positions_.push_back(position);
    invalidateCache();
}

void PortfolioManager::removePosition(core::Size index) {
    core::UniqueLockGuard<core::ReadWriteLock> lock(positionsLock_);
    if (index < positions_.size()) {
        positions_.erase(positions_.begin() + index);
        invalidateCache();
    }
}

void PortfolioManager::updatePosition(core::Size index, const Position& position) {
    core::UniqueLockGuard<core::ReadWriteLock> lock(positionsLock_);
    if (index < positions_.size()) {
        positions_[index] = position;
        invalidateCache();
    }
}

PortfolioManager::PortfolioGreeks PortfolioManager::calculatePortfolioGreeks() const {
    {
        core::SharedLockGuard<core::ReadWriteLock> lock(cacheLock_);
        if (cacheValid_) {
            return cachedPortfolioGreeks_;
        }
    }
    
    PortfolioGreeks portfolio;
    
    {
        core::SharedLockGuard<core::ReadWriteLock> lock(positionsLock_);
        
        if (positions_.empty()) {
            core::UniqueLockGuard<core::ReadWriteLock> cacheLock(cacheLock_);
            cachedPortfolioGreeks_ = portfolio;
            cacheValid_ = true;
            return portfolio;
        }
        
        core::MarketDataVector markets;
        core::Vector quantities;
        
        markets.reserve(positions_.size());
        quantities.reserve(positions_.size());
        
        for (const auto& position : positions_) {
            markets.push_back(position.market);
            quantities.push_back(position.quantity);
        }
        
        core::GreeksVector positionGreeks;
        calculator_->analyticalBatch(markets, positionGreeks);
        
        portfolio.totalGreeks = calculator_->aggregate(positionGreeks, quantities);
        
        portfolio.positionGreeks.reserve(positions_.size());
        for (core::Size i = 0; i < positions_.size(); ++i) {
            PositionGreeks posGreeks;
            posGreeks.greeks = positionGreeks[i];
            posGreeks.scaledGreeks.price = positionGreeks[i].price * quantities[i];
            posGreeks.scaledGreeks.delta = positionGreeks[i].delta * quantities[i];
            posGreeks.scaledGreeks.gamma = positionGreeks[i].gamma * quantities[i];
            posGreeks.scaledGreeks.theta = positionGreeks[i].theta * quantities[i];
            posGreeks.scaledGreeks.vega = positionGreeks[i].vega * quantities[i];
            posGreeks.scaledGreeks.rho = positionGreeks[i].rho * quantities[i];
            posGreeks.scaledGreeks.impliedVol = positionGreeks[i].impliedVol;
            
            posGreeks.dollarDelta = calculator_->calculateDollarDelta(positionGreeks[i], 
                                                                     positions_[i].market.spot,
                                                                     quantities[i]);
            posGreeks.dollarGamma = calculator_->calculateDollarGamma(positionGreeks[i],
                                                                     positions_[i].market.spot,
                                                                     quantities[i]);
            posGreeks.pinRisk = calculator_->calculatePinRisk(positions_[i].market.spot,
                                                             positions_[i].market.strike,
                                                             positions_[i].market.timeToExpiry,
                                                             quantities[i]);
            
            portfolio.positionGreeks.push_back(posGreeks);
        }
        
        portfolio.totalDollarDelta = 0.0;
        portfolio.totalDollarGamma = 0.0;
        portfolio.totalPinRisk = 0.0;
        
        for (const auto& posGreeks : portfolio.positionGreeks) {
            portfolio.totalDollarDelta += posGreeks.dollarDelta;
            portfolio.totalDollarGamma += posGreeks.dollarGamma;
            portfolio.totalPinRisk += posGreeks.pinRisk;
        }
    }
    
    core::UniqueLockGuard<core::ReadWriteLock> cacheLock(cacheLock_);
    cachedPortfolioGreeks_ = portfolio;
    cacheValid_ = true;
    
    return portfolio;
}

PortfolioManager::RiskBreakdown PortfolioManager::analyzeRisk() const {
    RiskBreakdown breakdown;
    
    const auto portfolio = calculatePortfolioGreeks();
    
    breakdown.deltaRisk = std::abs(portfolio.totalGreeks.delta);
    breakdown.gammaRisk = std::abs(portfolio.totalGreeks.gamma);
    breakdown.vegaRisk = std::abs(portfolio.totalGreeks.vega);
    breakdown.thetaRisk = std::abs(portfolio.totalGreeks.theta);
    breakdown.rhoRisk = std::abs(portfolio.totalGreeks.rho);
    
    breakdown.concentrationRisk.resize(positions_.size());
    const core::Real totalNotional = std::accumulate(positions_.begin(), positions_.end(), 0.0,
        [](core::Real sum, const Position& pos) {
            return sum + std::abs(pos.quantity * pos.market.spot);
        });
    
    if (totalNotional > 0.0) {
        for (core::Size i = 0; i < positions_.size(); ++i) {
            const core::Real notional = std::abs(positions_[i].quantity * positions_[i].market.spot);
            breakdown.concentrationRisk[i] = notional / totalNotional;
        }
    }
    
    breakdown.sectorRisk = calculateSectorRisk();
    breakdown.maturityRisk = calculateMaturityRisk();
    breakdown.strikeRisk = calculateStrikeRisk();
    
    return breakdown;
}

core::Matrix PortfolioManager::calculateRiskLadder(const core::Vector& spotLadder) const {
    core::SharedLockGuard<core::ReadWriteLock> lock(positionsLock_);
    
    if (positions_.empty() || spotLadder.empty()) {
        return {};
    }
    
    const core::Size numPositions = positions_.size();
    const core::Size numSpots = spotLadder.size();
    
    core::Matrix deltaLadder(numPositions, core::Vector(numSpots, 0.0));
    core::Matrix gammaLadder(numPositions, core::Vector(numSpots, 0.0));
    
    core::MarketDataVector markets;
    core::Vector quantities;
    
    for (const auto& position : positions_) {
        markets.push_back(position.market);
        quantities.push_back(position.quantity);
    }
    
    calculator_->calculateRiskLadder(markets, quantities, spotLadder, deltaLadder, gammaLadder);
    
    return deltaLadder;
}

PortfolioManager::ScenarioResults PortfolioManager::runScenarioAnalysis(const ScenarioParams& params) const {
    ScenarioResults results;
    
    core::SharedLockGuard<core::ReadWriteLock> lock(positionsLock_);
    
    if (positions_.empty()) {
        return results;
    }
    
    core::MarketDataVector baseMarkets;
    core::Vector quantities;
    
    for (const auto& position : positions_) {
        baseMarkets.push_back(position.market);
        quantities.push_back(position.quantity);
    }
    
    const core::Real basePortfolioValue = calculatePortfolioValue(baseMarkets, quantities);
    
    const core::Size numSpotScenarios = params.spotShifts.size();
    const core::Size numVolScenarios = params.volShifts.size();
    const core::Size numTimeScenarios = params.timeDecays.size();
    const core::Size totalScenarios = numSpotScenarios * numVolScenarios * numTimeScenarios;
    
    results.scenarioValues.reserve(totalScenarios);
    results.pnlValues.reserve(totalScenarios);
    
    for (core::Real spotShift : params.spotShifts) {
        for (core::Real volShift : params.volShifts) {
            for (core::Real timeDecay : params.timeDecays) {
                core::MarketDataVector scenarioMarkets = baseMarkets;
                
                for (auto& market : scenarioMarkets) {
                    market.spot *= (1.0 + spotShift);
                    market.volatility += volShift;
                    market.timeToExpiry = std::max(0.001, market.timeToExpiry - timeDecay);
                }
                
                const core::Real scenarioValue = calculatePortfolioValue(scenarioMarkets, quantities);
                const core::Real pnl = scenarioValue - basePortfolioValue;
                
                results.scenarioValues.push_back(scenarioValue);
                results.pnlValues.push_back(pnl);
            }
        }
    }
    
    if (!results.pnlValues.empty()) {
        auto sortedPnl = results.pnlValues;
        std::sort(sortedPnl.begin(), sortedPnl.end());
        
        results.worstCase = sortedPnl.front();
        results.bestCase = sortedPnl.back();
        results.expectedPnl = core::Statistics::mean(sortedPnl.begin(), sortedPnl.end());
        
        const core::Size var95Index = static_cast<core::Size>(sortedPnl.size() * 0.05);
        const core::Size var99Index = static_cast<core::Size>(sortedPnl.size() * 0.01);
        
        results.var95 = sortedPnl[var95Index];
        results.var99 = sortedPnl[var99Index];
        
        if (var95Index > 0) {
            results.expectedShortfall95 = core::Statistics::mean(sortedPnl.begin(), 
                                                                sortedPnl.begin() + var95Index);
        }
        
        if (var99Index > 0) {
            results.expectedShortfall99 = core::Statistics::mean(sortedPnl.begin(),
                                                                sortedPnl.begin() + var99Index);
        }
    }
    
    return results;
}

PortfolioManager::HedgeRecommendation PortfolioManager::generateHedgeRecommendation(
    const HedgeParams& params) const {
    
    HedgeRecommendation recommendation;
    
    const auto portfolio = calculatePortfolioGreeks();
    
    recommendation.currentExposure.delta = portfolio.totalGreeks.delta;
    recommendation.currentExposure.gamma = portfolio.totalGreeks.gamma;
    recommendation.currentExposure.vega = portfolio.totalGreeks.vega;
    recommendation.currentExposure.theta = portfolio.totalGreeks.theta;
    
    if (std::abs(portfolio.totalGreeks.delta) > params.deltaThreshold) {
        HedgeInstrument deltaHedge;
        deltaHedge.instrumentType = "Stock";
        deltaHedge.quantity = -portfolio.totalGreeks.delta;
        deltaHedge.cost = std::abs(deltaHedge.quantity) * 0.01;
        deltaHedge.purpose = "Delta neutrality";
        recommendation.hedgeInstruments.push_back(deltaHedge);
    }
    
    if (std::abs(portfolio.totalGreeks.gamma) > params.gammaThreshold) {
        HedgeInstrument gammaHedge;
        gammaHedge.instrumentType = "ATM Options";
        gammaHedge.quantity = -portfolio.totalGreeks.gamma / 0.1;
        gammaHedge.cost = std::abs(gammaHedge.quantity) * 0.05;
        gammaHedge.purpose = "Gamma neutrality";
        recommendation.hedgeInstruments.push_back(gammaHedge);
    }
    
    if (std::abs(portfolio.totalGreeks.vega) > params.vegaThreshold) {
        HedgeInstrument vegaHedge;
        vegaHedge.instrumentType = "Long-dated Options";
        vegaHedge.quantity = -portfolio.totalGreeks.vega / 0.2;
        vegaHedge.cost = std::abs(vegaHedge.quantity) * 0.03;
        vegaHedge.purpose = "Vega neutrality";
        recommendation.hedgeInstruments.push_back(vegaHedge);
    }
    
    recommendation.totalHedgeCost = std::accumulate(recommendation.hedgeInstruments.begin(),
                                                   recommendation.hedgeInstruments.end(), 0.0,
                                                   [](core::Real sum, const HedgeInstrument& inst) {
                                                       return sum + inst.cost;
                                                   });
    
    recommendation.projectedExposure.delta = portfolio.totalGreeks.delta;
    recommendation.projectedExposure.gamma = portfolio.totalGreeks.gamma;
    recommendation.projectedExposure.vega = portfolio.totalGreeks.vega;
    recommendation.projectedExposure.theta = portfolio.totalGreeks.theta;
    
    for (const auto& hedge : recommendation.hedgeInstruments) {
        if (hedge.instrumentType == "Stock") {
            recommendation.projectedExposure.delta += hedge.quantity;
        } else if (hedge.instrumentType == "ATM Options") {
            recommendation.projectedExposure.gamma += hedge.quantity * 0.1;
        } else if (hedge.instrumentType == "Long-dated Options") {
            recommendation.projectedExposure.vega += hedge.quantity * 0.2;
        }
    }
    
    recommendation.hedgeEffectiveness = 1.0 - (
        std::abs(recommendation.projectedExposure.delta) / std::max(std::abs(portfolio.totalGreeks.delta), 1.0) +
        std::abs(recommendation.projectedExposure.gamma) / std::max(std::abs(portfolio.totalGreeks.gamma), 1.0) +
        std::abs(recommendation.projectedExposure.vega) / std::max(std::abs(portfolio.totalGreeks.vega), 1.0)
    ) / 3.0;
    
    return recommendation;
}

void PortfolioManager::setRealTimeUpdates(bool enabled) {
    realTimeUpdates_ = enabled;
    if (enabled) {
        invalidateCache();
    }
}

std::vector<PortfolioManager::Position> PortfolioManager::getPositions() const {
    core::SharedLockGuard<core::ReadWriteLock> lock(positionsLock_);
    return positions_;
}

core::Size PortfolioManager::getPositionCount() const {
    core::SharedLockGuard<core::ReadWriteLock> lock(positionsLock_);
    return positions_.size();
}

void PortfolioManager::clearPositions() {
    core::UniqueLockGuard<core::ReadWriteLock> lock(positionsLock_);
    positions_.clear();
    invalidateCache();
}

void PortfolioManager::invalidateCache() const {
    core::UniqueLockGuard<core::ReadWriteLock> lock(cacheLock_);
    cacheValid_ = false;
}

core::Vector PortfolioManager::calculateSectorRisk() const {
    std::map<std::string, core::Real> sectorExposure;
    
    for (const auto& position : positions_) {
        const core::Real notional = std::abs(position.quantity * position.market.spot);
        sectorExposure[position.sector] += notional;
    }
    
    core::Vector sectorRisks;
    sectorRisks.reserve(sectorExposure.size());
    
    const core::Real totalNotional = std::accumulate(sectorExposure.begin(), sectorExposure.end(), 0.0,
        [](core::Real sum, const auto& pair) { return sum + pair.second; });
    
    if (totalNotional > 0.0) {
        for (const auto& [sector, exposure] : sectorExposure) {
            sectorRisks.push_back(exposure / totalNotional);
        }
    }
    
    return sectorRisks;
}

core::Vector PortfolioManager::calculateMaturityRisk() const {
    std::map<core::Real, core::Real> maturityBuckets;
    
    for (const auto& position : positions_) {
        const core::Real maturityBucket = std::floor(position.market.timeToExpiry * 12.0) / 12.0;
        const core::Real notional = std::abs(position.quantity * position.market.spot);
        maturityBuckets[maturityBucket] += notional;
    }
    
    core::Vector maturityRisks;
    maturityRisks.reserve(maturityBuckets.size());
    
    const core::Real totalNotional = std::accumulate(maturityBuckets.begin(), maturityBuckets.end(), 0.0,
        [](core::Real sum, const auto& pair) { return sum + pair.second; });
    
    if (totalNotional > 0.0) {
        for (const auto& [maturity, exposure] : maturityBuckets) {
            maturityRisks.push_back(exposure / totalNotional);
        }
    }
    
    return maturityRisks;
}

core::Vector PortfolioManager::calculateStrikeRisk() const {
    std::map<core::Real, core::Real> strikeBuckets;
    
    for (const auto& position : positions_) {
        const core::Real moneyness = position.market.strike / position.market.spot;
        const core::Real strikeBucket = std::floor(moneyness * 20.0) / 20.0;
        const core::Real notional = std::abs(position.quantity * position.market.spot);
        strikeBuckets[strikeBucket] += notional;
    }
    
    core::Vector strikeRisks;
    strikeRisks.reserve(strikeBuckets.size());
    
    const core::Real totalNotional = std::accumulate(strikeBuckets.begin(), strikeBuckets.end(), 0.0,
        [](core::Real sum, const auto& pair) { return sum + pair.second; });
    
    if (totalNotional > 0.0) {
        for (const auto& [strike, exposure] : strikeBuckets) {
            strikeRisks.push_back(exposure / totalNotional);
        }
    }
    
    return strikeRisks;
}

core::Real PortfolioManager::calculatePortfolioValue(const core::MarketDataVector& markets,
                                                    const core::Vector& quantities) const {
    core::GreeksVector greeks;
    calculator_->analyticalBatch(markets, greeks);
    
    core::Real totalValue = 0.0;
    for (core::Size i = 0; i < greeks.size() && i < quantities.size(); ++i) {
        totalValue += greeks[i].price * quantities[i];
    }
    
    return totalValue;
}

}