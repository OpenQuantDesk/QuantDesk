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
    common::ExtendedGreeks extendedGreeks;
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
    
    double getCash() const;
    void setCash(double cash);
    void addCash(double amount);
    
    common::Greeks getPortfolioGreeks() const;
    std::vector<std::string> getUnderlyings() const;
    std::map<std::string, double> getUnderlyingExposures() const;
    
    void addPositionCallback(std::function<void(const EnhancedPosition&)> callback);
    void addMetricsCallback(std::function<void(const PortfolioMetrics&)> callback);
    
private:
    void analyticsLoop();
    void calculatePortfolioRisk(PortfolioMetrics& metrics) const;
    void calculateAdvancedMetrics(PortfolioMetrics& metrics) const;
    void updateMetricsHistory(const PortfolioMetrics& metrics);
    void updatePnLHistory();
    void notifyPositionUpdated(const EnhancedPosition& position);
    void updateTotalValues();
    double calculateConcentrationRisk() const;
    double calculateLeverageRatio() const;
    double calculatePositionRisk(const EnhancedPosition& position) const;
};

}