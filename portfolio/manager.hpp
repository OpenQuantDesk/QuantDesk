#pragma once

#include "common/types.hpp"
#include "math/engine.hpp"
#include "math/metrics/probability.hpp"
#include <map>
#include <vector>
#include <memory>
#include <shared_mutex>
#include <atomic>
#include <functional>
#include <future>
#include <thread>
#include <chrono>

namespace portfolio {

struct EnhancedPosition {
    std::string symbol;
    std::string optionSymbol;
    std::string underlying;
    double strike = 0.0;
    std::string expiration;
    bool isCall = true;
    double quantity = 0.0;
    double avgPrice = 0.0;
    double currentPrice = 0.0;
    double underlyingPrice = 0.0;
    
    common::Greeks greeks;
    math::ExtendedGreeks extendedGreeks;
    math::ProbabilityMetrics probMetrics;
    
    double theoreticalValue = 0.0;
    double impliedVol = 0.0;
    double unrealizedPnL = 0.0;
    double dayChange = 0.0;
    double deltaAdjustedExposure = 0.0;
    double gammaRisk = 0.0;
    double vegaRisk = 0.0;
    double thetaDecay = 0.0;
    
    std::chrono::system_clock::time_point entryTime;
    std::chrono::system_clock::time_point lastUpdate;
    
    void updateAnalytics(double underlyingPx, double vol, double riskFreeRate,
                        const math::MathEngine& engine);
    double getDaysToExpiry() const;
    double getNotionalValue() const;
    double getImpliedProbability() const;
    double getMoneyness() const;
    bool isExpiring(int days = 7) const;
    
private:
    std::chrono::system_clock::time_point parseExpirationDate(const std::string& expiration) const;
};

struct PortfolioMetrics {
    double totalValue = 0.0;
    double totalPnL = 0.0;
    double dayPnL = 0.0;
    double cash = 0.0;
    
    double netDelta = 0.0;
    double netGamma = 0.0;
    double netTheta = 0.0;
    double netVega = 0.0;
    double netRho = 0.0;
    
    double netVolga = 0.0;
    double netVanna = 0.0;
    double dollarDelta = 0.0;
    double dollarGamma = 0.0;
    
    double var95 = 0.0;
    double var99 = 0.0;
    double expectedShortfall = 0.0;
    double maxDrawdown = 0.0;
    double leverageRatio = 0.0;
    
    double portfolioProbProfit = 0.0;
    double sharpeRatio = 0.0;
    double kellyOptimal = 0.0;
    double correlationRisk = 0.0;
    
    std::map<std::string, double> underlyingExposures;
    std::map<std::string, double> sectorExposures;
    std::map<std::string, double> strategyAllocations;
    
    std::chrono::system_clock::time_point timestamp;
};

struct RiskLimits {
    double maxDelta = 100.0;
    double maxGamma = 10.0;
    double maxTheta = -100.0;
    double maxVega = 500.0;
    double maxPositionSize = 10000.0;
    double maxPortfolioRisk = 50000.0;
    double maxLeverage = 2.0;
    double maxConcentration = 0.2;
    double maxVaR = -5000.0;
    double maxDrawdown = -0.15;
};

struct RiskAssessment {
    bool withinLimits = true;
    std::vector<std::string> violations;
    double riskScore = 0.0;
    std::map<std::string, double> recommendations;
    std::vector<common::OrderRequest> hedgeOrders;
    double capitalEfficiency = 0.0;
    double portfolioHeat = 0.0;
};

class PortfolioManager {
private:
    std::map<std::string, EnhancedPosition> positions_;
    std::vector<PortfolioMetrics> metricsHistory_;
    std::map<std::string, std::vector<double>> pnlHistory_;
    
    double cash_ = 0.0;
    double initialCapital_ = 100000.0;
    
    std::shared_ptr<math::MathEngine> mathEngine_;
    
    mutable std::shared_mutex positionsLock_;
    mutable std::shared_mutex metricsLock_;
    
    std::atomic<double> totalValue_{0.0};
    std::atomic<double> totalPnL_{0.0};
    std::atomic<double> riskScore_{0.0};
    
    std::vector<std::function<void(const EnhancedPosition&)>> positionCallbacks_;
    std::vector<std::function<void(const PortfolioMetrics&)>> metricsCallbacks_;
    std::vector<std::function<void(const RiskAssessment&)>> riskCallbacks_;
    
    std::thread analyticsThread_;
    std::atomic<bool> running_{false};
    
public:
    explicit PortfolioManager(std::shared_ptr<math::MathEngine> engine);
    ~PortfolioManager();
    
    PortfolioManager(const PortfolioManager&) = delete;
    PortfolioManager& operator=(const PortfolioManager&) = delete;
    
    void start();
    void stop();
    
    void updatePosition(const common::Position& position);
    void updatePosition(const std::string& symbol, double quantity, double price);
    void removePosition(const std::string& symbol);
    
    void updatePositionPrice(const std::string& symbol, double price, double underlyingPrice);
    void updatePositionGreeks(const std::string& symbol, const common::Greeks& greeks);
    void updateAllPositions(const std::map<std::string, common::Quote>& quotes,
                           const std::map<std::string, common::OptionChain>& chains);
    
    std::future<PortfolioMetrics> calculateMetrics() const;
    std::future<std::vector<EnhancedPosition>> getPositions() const;
    std::future<EnhancedPosition> getPosition(const std::string& symbol) const;
    
    std::future<std::vector<EnhancedPosition>> getPositionsByUnderlying(const std::string& underlying) const;
    std::future<std::vector<EnhancedPosition>> getExpiringPositions(int days = 7) const;
    std::future<std::vector<EnhancedPosition>> getHighRiskPositions(double threshold = 0.8) const;
    
    double getTotalValue() const { return totalValue_.load(); }
    double getTotalPnL() const { return totalPnL_.load(); }
    double getRiskScore() const { return riskScore_.load(); }
    
    double getCash() const;
    void setCash(double cash);
    void addCash(double amount);
    
    common::Greeks getPortfolioGreeks() const;
    math::ExtendedGreeks getExtendedGreeks() const;
    
    std::vector<std::string> getUnderlyings() const;
    std::map<std::string, double> getUnderlyingExposures() const;
    std::map<std::string, double> getSectorExposures() const;
    
    std::future<std::vector<double>> calculatePnLAttribution(const std::string& period = "1D") const;
    std::future<double> calculateBeta(const std::string& benchmark = "SPY") const;
    std::future<std::map<std::string, double>> calculateCorrelations() const;
    
    void addPositionCallback(std::function<void(const EnhancedPosition&)> callback);
    void addMetricsCallback(std::function<void(const PortfolioMetrics&)> callback);
    void addRiskCallback(std::function<void(const RiskAssessment&)> callback);
    
    std::vector<PortfolioMetrics> getMetricsHistory(int days = 30) const;
    std::map<std::string, std::vector<double>> getPnLHistory(int days = 30) const;
    
private:
    void analyticsLoop();
    void calculatePortfolioRisk(PortfolioMetrics& metrics) const;
    void calculateAdvancedMetrics(PortfolioMetrics& metrics) const;
    void updateMetricsHistory(const PortfolioMetrics& metrics);
    void updatePnLHistory();
    
    void notifyPositionUpdated(const EnhancedPosition& position);
    void notifyMetricsUpdated(const PortfolioMetrics& metrics);
    void notifyRiskUpdated(const RiskAssessment& assessment);
    
    void updateTotalValues();
    double calculateConcentrationRisk() const;
    double calculateLeverageRatio() const;
    double calculatePositionRisk(const EnhancedPosition& position) const;
};

class RiskManager {
private:
    std::shared_ptr<PortfolioManager> portfolio_;
    std::shared_ptr<math::MathEngine> mathEngine_;
    RiskLimits limits_;
    
    std::vector<RiskAssessment> assessmentHistory_;
    mutable std::shared_mutex assessmentLock_;
    
    std::thread monitoringThread_;
    std::atomic<bool> monitoring_{false};
    
    std::vector<std::function<void(const RiskAssessment&)>> alertCallbacks_;
    
public:
    RiskManager(std::shared_ptr<PortfolioManager> portfolio,
               std::shared_ptr<math::MathEngine> engine);
    ~RiskManager();
    
    RiskManager(const RiskManager&) = delete;
    RiskManager& operator=(const RiskManager&) = delete;
    
    void startMonitoring();
    void stopMonitoring();
    
    std::future<RiskAssessment> assessRisk() const;
    std::future<bool> validateOrder(const common::OrderRequest& order) const;
    std::future<std::vector<common::OrderRequest>> generateHedgeOrders() const;
    
    void setRiskLimits(const RiskLimits& limits);
    RiskLimits getRiskLimits() const { return limits_; }
    
    std::future<double> calculatePortfolioVaR(double confidence = 0.95, int horizon = 1) const;
    std::future<double> calculateExpectedShortfall(double confidence = 0.95) const;
    std::future<double> calculateMaxDrawdown() const;
    
    std::future<std::map<std::string, double>> calculateComponentVaR() const;
    std::future<std::map<std::string, double>> calculateMarginalVaR() const;
    
    void addAlertCallback(std::function<void(const RiskAssessment&)> callback);
    std::vector<RiskAssessment> getAssessmentHistory(int days = 7) const;
    
private:
    void monitoringLoop();
    double calculateOrderImpact(const common::OrderRequest& order) const;
    double calculateRiskScore(const PortfolioMetrics& metrics) const;
    std::vector<common::OrderRequest> generateDeltaHedge(double targetDelta = 0.0) const;
    std::vector<common::OrderRequest> generateGammaHedge() const;
    std::vector<common::OrderRequest> generateVegaHedge() const;
    
    void notifyRiskAlert(const RiskAssessment& assessment);
    bool isViolation(const std::string& metric, double value, double limit) const;
};

}