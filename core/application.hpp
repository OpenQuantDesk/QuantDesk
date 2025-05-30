#pragma once

#include "common/types.hpp"
#include "broker/interface.hpp"
#include "math/engine.hpp"
#include "data/economic.hpp"
#include "portfolio/manager.hpp"
#include "strategy/engine.hpp"
#include <memory>
#include <atomic>
#include <thread>
#include <vector>
#include <map>
#include <functional>
#include <shared_mutex>
#include <chrono>

namespace core {

struct ApplicationConfig {
    std::vector<common::BrokerConfig> brokers;
    std::map<std::string, std::string> economicDataKeys;
    std::vector<std::string> watchedSymbols = {"SPY", "QQQ", "IWM"};
    int updateIntervalSeconds = 10;
    bool enableBackgroundAnalysis = true;
    bool enableRiskManagement = true;
    std::string configFile = "config.json";
    
    struct PerformanceConfig {
        bool enableHardwareOptimization = true;
        size_t maxThreads = 0;
        bool enableCaching = true;
        size_t cacheSize = 1000;
    } performance;
    
    struct LoggingConfig {
        bool enableDebugLogging = false;
        std::string logLevel = "INFO";
        std::string logFile = "application.log";
    } logging;
};

class Application {
private:
    ApplicationConfig config_;
    
    std::shared_ptr<math::MathEngine> mathEngine_;
    std::shared_ptr<data::EconomicDataManager> economicData_;
    std::shared_ptr<portfolio::PortfolioManager> portfolio_;
    std::shared_ptr<portfolio::RiskManager> riskManager_;
    std::shared_ptr<strategy::StrategyEngine> strategyEngine_;
    std::unique_ptr<broker::BrokerManager> brokerManager_;
    std::unique_ptr<broker::PluginLoader> pluginLoader_;
    
    std::atomic<bool> running_{false};
    std::vector<std::thread> threads_;
    
    std::map<std::string, common::Quote> currentQuotes_;
    std::map<std::string, common::OptionChain> currentOptionChains_;
    std::vector<strategy::StrategyOpportunity> currentOpportunities_;
    
    mutable std::shared_mutex dataLock_;
    
    std::vector<std::function<void(const std::map<std::string, common::Quote>&)>> quoteCallbacks_;
    std::vector<std::function<void(const std::vector<strategy::StrategyOpportunity>&)>> opportunityCallbacks_;
    std::vector<std::function<void(const std::string&)>> statusCallbacks_;
    
    std::atomic<std::chrono::system_clock::time_point> lastDataUpdate_;
    std::atomic<std::chrono::system_clock::time_point> lastAnalysisUpdate_;
    
public:
    Application();
    ~Application();
    
    Application(const Application&) = delete;
    Application& operator=(const Application&) = delete;
    
    void initialize();
    void initialize(const ApplicationConfig& config);
    void run();
    void shutdown();
    
    bool isRunning() const { return running_; }
    
    std::shared_ptr<math::MathEngine> getMathEngine() const { return mathEngine_; }
    std::shared_ptr<data::EconomicDataManager> getEconomicData() const { return economicData_; }
    std::shared_ptr<portfolio::PortfolioManager> getPortfolio() const { return portfolio_; }
    std::shared_ptr<portfolio::RiskManager> getRiskManager() const { return riskManager_; }
    std::shared_ptr<strategy::StrategyEngine> getStrategyEngine() const { return strategyEngine_; }
    broker::BrokerManager* getBrokerManager() const { return brokerManager_.get(); }
    
    std::map<std::string, common::Quote> getCurrentQuotes() const;
    std::map<std::string, common::OptionChain> getCurrentOptionChains() const;
    std::vector<strategy::StrategyOpportunity> getCurrentOpportunities() const;
    
    void addQuoteCallback(std::function<void(const std::map<std::string, common::Quote>&)> callback);
    void addOpportunityCallback(std::function<void(const std::vector<strategy::StrategyOpportunity>&)> callback);
    void addStatusCallback(std::function<void(const std::string&)> callback);
    
    std::future<broker::Result<common::OrderResponse>> placeOrder(const common::OrderRequest& request);
    std::future<broker::Result<common::OrderResponse>> cancelOrder(const std::string& accountId, 
                                                                  const std::string& orderId);
    
    void reloadWatchlist();
    void addSymbol(const std::string& symbol);
    void removeSymbol(const std::string& symbol);
    std::vector<std::string> getWatchedSymbols() const;
    
    struct SystemMetrics {
        double cpuUsage = 0.0;
        double memoryUsage = 0.0;
        size_t threadCount = 0;
        std::chrono::system_clock::time_point lastUpdate;
        std::map<std::string, double> componentTiming;
    };
    
    SystemMetrics getSystemMetrics() const;
    
private:
    void initializeComponents();
    void loadConfiguration();
    void saveConfiguration() const;
    void loadPlugins();
    void initializeBrokers();
    
    void startBackgroundThreads();
    void stopBackgroundThreads();
    
    void dataCollectionLoop();
    void analysisLoop();
    void riskMonitoringLoop();
    void systemMonitoringLoop();
    
    void updateQuotes();
    void updateOptionChains();
    void updatePortfolioAnalytics();
    void runStrategyAnalysis();
    void performRiskCheck();
    
    void notifyQuoteUpdate(const std::map<std::string, common::Quote>& quotes);
    void notifyOpportunityUpdate(const std::vector<strategy::StrategyOpportunity>& opportunities);
    void notifyStatusUpdate(const std::string& status);
    
    ApplicationConfig loadDefaultConfig() const;
    std::vector<std::string> loadWatchlistFromFile(const std::string& filename) const;
    void saveWatchlistToFile(const std::string& filename) const;
    
    void optimizePerformance();
    void setupLogging();
};

class ConfigManager {
public:
    static ApplicationConfig loadFromFile(const std::string& filename);
    static void saveToFile(const ApplicationConfig& config, const std::string& filename);
    static ApplicationConfig loadFromEnvironment();
    
private:
    static common::BrokerConfig parseBrokerConfig(const std::map<std::string, std::string>& params);
    static std::map<std::string, std::string> parseEconomicDataConfig(const std::map<std::string, std::string>& params);
};

class PerformanceMonitor {
private:
    std::atomic<double> cpuUsage_{0.0};
    std::atomic<double> memoryUsage_{0.0};
    std::map<std::string, std::chrono::high_resolution_clock::time_point> timingPoints_;
    mutable std::shared_mutex timingLock_;
    
public:
    void startTiming(const std::string& component);
    void endTiming(const std::string& component);
    double getComponentTiming(const std::string& component) const;
    
    void updateSystemMetrics();
    double getCpuUsage() const { return cpuUsage_.load(); }
    double getMemoryUsage() const { return memoryUsage_.load(); }
};

}