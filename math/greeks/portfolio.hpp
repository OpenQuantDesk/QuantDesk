#pragma once

#include "../core/types.hpp"
#include "../core/threading.hpp"
#include "calculator.hpp"
#include <string>
#include <map>

namespace math::greeks {

class PortfolioManager {
public:
    struct Position {
        core::MarketData market;
        core::Real quantity = 0.0;
        std::string instrumentId;
        std::string sector = "Default";
        core::Real costBasis = 0.0;
        core::Real marketValue = 0.0;
        std::string strategy = "Outright";
        core::Integer entryDate = 0;
        bool isHedge = false;
    };
    
    struct PositionGreeks {
        core::Greeks greeks;
        core::Greeks scaledGreeks;
        core::Real dollarDelta = 0.0;
        core::Real dollarGamma = 0.0;
        core::Real pinRisk = 0.0;
        core::Real percentOfPortfolio = 0.0;
        core::Real pnl = 0.0;
    };
    
    struct PortfolioGreeks {
        core::Greeks totalGreeks;
        core::Real totalDollarDelta = 0.0;
        core::Real totalDollarGamma = 0.0;
        core::Real totalPinRisk = 0.0;
        core::Real totalValue = 0.0;
        core::Real totalPnl = 0.0;
        std::vector<PositionGreeks> positionGreeks;
    };
    
    struct RiskBreakdown {
        core::Real deltaRisk = 0.0;
        core::Real gammaRisk = 0.0;
        core::Real vegaRisk = 0.0;
        core::Real thetaRisk = 0.0;
        core::Real rhoRisk = 0.0;
        core::Vector concentrationRisk;
        core::Vector sectorRisk;
        core::Vector maturityRisk;
        core::Vector strikeRisk;
        core::Real liquidityRisk = 0.0;
        core::Real correlationRisk = 0.0;
    };
    
    struct ScenarioParams {
        core::Vector spotShifts = {-0.20, -0.10, -0.05, 0.0, 0.05, 0.10, 0.20};
        core::Vector volShifts = {-0.05, -0.02, -0.01, 0.0, 0.01, 0.02, 0.05};
        core::Vector timeDecays = {0.0, 1.0/365.0, 7.0/365.0, 30.0/365.0};
        bool includeRateShocks = false;
        core::Vector rateShifts = {-0.01, 0.0, 0.01};
    };
    
    struct ScenarioResults {
        core::Vector scenarioValues;
        core::Vector pnlValues;
        core::Real worstCase = 0.0;
        core::Real bestCase = 0.0;
        core::Real expectedPnl = 0.0;
        core::Real var95 = 0.0;
        core::Real var99 = 0.0;
        core::Real expectedShortfall95 = 0.0;
        core::Real expectedShortfall99 = 0.0;
        core::Real probabilityOfLoss = 0.0;
    };
    
    struct HedgeParams {
        core::Real deltaThreshold = 100.0;
        core::Real gammaThreshold = 50.0;
        core::Real vegaThreshold = 1000.0;
        core::Real thetaThreshold = 500.0;
        bool hedgeDelta = true;
        bool hedgeGamma = true;
        bool hedgeVega = true;
        bool allowPartialHedges = true;
        core::Real maxHedgeCost = 10000.0;
    };
    
    struct HedgeInstrument {
        std::string instrumentType;
        core::Real quantity = 0.0;
        core::Real cost = 0.0;
        std::string purpose;
        core::Real effectiveness = 0.0;
    };
    
    struct HedgeRecommendation {
        PositionGreeks currentExposure;
        PositionGreeks projectedExposure;
        std::vector<HedgeInstrument> hedgeInstruments;
        core::Real totalHedgeCost = 0.0;
        core::Real hedgeEffectiveness = 0.0;
        core::Real residualRisk = 0.0;
        std::string rationale;
    };
    
private:
    std::shared_ptr<GreeksCalculator> calculator_;
    std::vector<Position> positions_;
    mutable core::ReadWriteLock positionsLock_;
    mutable core::ReadWriteLock cacheLock_;
    mutable bool cacheValid_ = false;
    mutable PortfolioGreeks cachedPortfolioGreeks_;
    bool realTimeUpdates_ = false;
    
public:
    explicit PortfolioManager(std::shared_ptr<GreeksCalculator> calculator);
    
    void addPosition(const Position& position);
    void removePosition(core::Size index);
    void updatePosition(core::Size index, const Position& position);
    void clearPositions();
    
    PortfolioGreeks calculatePortfolioGreeks() const;
    RiskBreakdown analyzeRisk() const;
    
    core::Matrix calculateRiskLadder(const core::Vector& spotLadder) const;
    core::Matrix calculateVolatilityLadder(const core::Vector& volShifts) const;
    core::Matrix calculateTimeLadder(const core::Vector& timeDecays) const;
    
    ScenarioResults runScenarioAnalysis(const ScenarioParams& params = ScenarioParams{}) const;
    HedgeRecommendation generateHedgeRecommendation(const HedgeParams& params = HedgeParams{}) const;
    
    struct PerformanceMetrics {
        core::Real totalReturn = 0.0;
        core::Real sharpeRatio = 0.0;
        core::Real maxDrawdown = 0.0;
        core::Real winRate = 0.0;
        core::Real profitFactor = 0.0;
        core::Real volatility = 0.0;
        core::Vector dailyPnl;
        core::Vector cumulativePnl;
    };
    
    PerformanceMetrics calculatePerformance(const core::Vector& historicalPnl) const;
    
    struct ExposureAnalysis {
        std::map<std::string, core::Real> sectorExposure;
        std::map<core::Real, core::Real> maturityExposure;
        std::map<core::Real, core::Real> strikeExposure;
        std::map<std::string, core::Real> strategyExposure;
        core::Real leverageRatio = 0.0;
        core::Real netExposure = 0.0;
        core::Real grossExposure = 0.0;
    };
    
    ExposureAnalysis analyzeExposure() const;
    
    struct MarginRequirement {
        core::Real initialMargin = 0.0;
        core::Real maintenanceMargin = 0.0;
        core::Real excessLiquidity = 0.0;
        core::Real marginUtilization = 0.0;
        bool marginCall = false;
        core::Vector positionMargins;
    };
    
    MarginRequirement calculateMarginRequirement() const;
    
    void setRealTimeUpdates(bool enabled);
    bool getRealTimeUpdates() const { return realTimeUpdates_; }
    
    std::vector<Position> getPositions() const;
    Position getPosition(core::Size index) const;
    core::Size getPositionCount() const;
    
    struct PortfolioSummary {
        core::Size totalPositions = 0;
        core::Real totalValue = 0.0;
        core::Real totalPnl = 0.0;
        core::Real totalDelta = 0.0;
        core::Real totalGamma = 0.0;
        core::Real totalVega = 0.0;
        core::Real totalTheta = 0.0;
        core::Real avgTimeToExpiry = 0.0;
        core::Real avgImpliedVol = 0.0;
        std::string largestPosition;
        std::string riskiestPosition;
    };
    
    PortfolioSummary getSummary() const;
    
    void exportToCSV(const std::string& filename) const;
    void importFromCSV(const std::string& filename);
    
    struct AlertCondition {
        std::string type;
        core::Real threshold = 0.0;
        bool enabled = true;
        std::string description;
    };
    
    struct Alert {
        std::string message;
        std::string severity;
        core::Integer timestamp = 0;
        std::string positionId;
    };
    
    void addAlertCondition(const AlertCondition& condition);
    std::vector<Alert> checkAlerts() const;
    
private:
    void invalidateCache() const;
    
    core::Vector calculateSectorRisk() const;
    core::Vector calculateMaturityRisk() const;
    core::Vector calculateStrikeRisk() const;
    
    core::Real calculatePortfolioValue(const core::MarketDataVector& markets,
                                      const core::Vector& quantities) const;
    
    std::vector<AlertCondition> alertConditions_;
    mutable std::vector<Alert> activeAlerts_;
    
    struct PositionLimits {
        core::Real maxPositionSize = 1000000.0;
        core::Real maxDelta = 10000.0;
        core::Real maxGamma = 5000.0;
        core::Real maxVega = 50000.0;
        core::Real maxConcentration = 0.25;
    };
    
    PositionLimits limits_;
    
public:
    void setPositionLimits(const PositionLimits& limits) { limits_ = limits; }
    const PositionLimits& getPositionLimits() const { return limits_; }
    
    bool validatePosition(const Position& position) const;
    std::vector<std::string> getPositionViolations(const Position& position) const;
};

class PortfolioOptimizer {
private:
    std::shared_ptr<PortfolioManager> portfolioManager_;
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit PortfolioOptimizer(std::shared_ptr<PortfolioManager> portfolioManager);
    
    struct OptimizationObjective {
        core::Real returnWeight = 1.0;
        core::Real riskWeight = -1.0;
        core::Real deltaWeight = 0.0;
        core::Real gammaWeight = 0.0;
        core::Real vegaWeight = 0.0;
        core::Real thetaWeight = 0.0;
        core::Real liquidityWeight = 0.0;
    };
    
    struct OptimizationConstraints {
        core::Real maxDelta = 1000.0;
        core::Real maxGamma = 500.0;
        core::Real maxVega = 10000.0;
        core::Real maxConcentration = 0.3;
        core::Vector sectorLimits;
        core::Vector maturityLimits;
        bool allowShorts = true;
        core::Real maxLeverage = 2.0;
    };
    
    struct OptimizationResult {
        bool success = false;
        core::Vector optimalWeights;
        core::Real expectedReturn = 0.0;
        core::Real expectedRisk = 0.0;
        core::Real sharpeRatio = 0.0;
        std::string errorMessage;
        core::Integer iterations = 0;
    };
    
    OptimizationResult optimizePortfolio(const std::vector<PortfolioManager::Position>& candidates,
                                        const OptimizationObjective& objective,
                                        const OptimizationConstraints& constraints) const;
    
    struct EfficientFrontier {
        core::Vector returns;
        core::Vector risks;
        core::Vector sharpeRatios;
        core::Matrix weights;
        core::Real maxSharpeReturn = 0.0;
        core::Real maxSharpeRisk = 0.0;
    };
    
    EfficientFrontier generateEfficientFrontier(const std::vector<PortfolioManager::Position>& candidates,
                                               core::Integer numPoints = 20) const;
    
    OptimizationResult deltaHedge(core::Real targetDelta = 0.0) const;
    OptimizationResult gammaHedge(core::Real targetGamma = 0.0) const;
    OptimizationResult vegaHedge(core::Real targetVega = 0.0) const;
    OptimizationResult neutralizeRisk() const;
    
private:
    core::Real calculateObjectiveFunction(const core::Vector& weights,
                                         const std::vector<PortfolioManager::Position>& positions,
                                         const OptimizationObjective& objective) const;
    
    bool satisfiesConstraints(const core::Vector& weights,
                             const std::vector<PortfolioManager::Position>& positions,
                             const OptimizationConstraints& constraints) const;
    
    core::Vector calculatePortfolioGreeks(const core::Vector& weights,
                                         const std::vector<PortfolioManager::Position>& positions) const;
};

inline PortfolioManager::PortfolioManager(std::shared_ptr<GreeksCalculator> calculator)
    : calculator_(std::move(calculator)) {}

inline PortfolioOptimizer::PortfolioOptimizer(std::shared_ptr<PortfolioManager> portfolioManager)
    : portfolioManager_(std::move(portfolioManager)) {}

inline PortfolioManager::Position PortfolioManager::getPosition(core::Size index) const {
    core::SharedLockGuard<core::ReadWriteLock> lock(positionsLock_);
    if (index < positions_.size()) {
        return positions_[index];
    }
    return Position{};
}

}