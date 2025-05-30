#pragma once

#include <atomic>
#include <chrono>
#include <future>
#include <map>
#include <memory>
#include <shared_mutex>
#include <string>
#include <vector>
#include <thread>

namespace data {

struct EconomicIndicator {
    std::string symbol;
    std::string name;
    double value = 0.0;
    std::string unit;
    std::chrono::system_clock::time_point date;
    std::string source;
    double previousValue = 0.0;
    double change = 0.0;
    double changePercent = 0.0;
    double historicalAverage = 0.0;
    double volatility = 0.0;
};

struct MarketSentiment {
    double vixLevel = 0.0;
    double putCallRatio = 0.0;
    double fearGreedIndex = 0.0;
    std::string marketRegime;
    double correlationSPY_VIX = 0.0;
    double skewIndex = 0.0;
    double termStructure = 0.0;
    std::chrono::system_clock::time_point timestamp;
};

struct TreasuryYieldCurve {
    std::map<int, double> yields;
    double slope = 0.0;
    double curvature = 0.0;
    bool inverted = false;
    double parallelShift = 0.0;
    double twistRisk = 0.0;
    std::chrono::system_clock::time_point timestamp;
};

struct StrategySignals {
    bool highVolatilityEnvironment = false;
    bool yieldCurveInversion = false;
    bool recessionary = false;
    double optimalIVRank = 0.5;
    std::string recommendedStrategy;
    std::vector<std::string> avoidStrategies;
    double confidenceLevel = 0.0;
    std::map<std::string, double> strategyWeights;
};

class IDataProvider {
public:
    virtual ~IDataProvider() = default;
    virtual std::future<EconomicIndicator> getIndicator(const std::string& symbol) = 0;
    virtual std::string getName() const = 0;
    virtual bool isAvailable() const = 0;
    virtual std::chrono::minutes getRateLimit() const { return std::chrono::minutes(1); }
};

class FREDProvider : public IDataProvider {
private:
    std::string apiKey_;
    std::string baseUrl_;
    std::atomic<std::chrono::system_clock::time_point> lastRequest_;
    
public:
    explicit FREDProvider(const std::string& apiKey);
    
    std::future<EconomicIndicator> getIndicator(const std::string& symbol) override;
    std::future<TreasuryYieldCurve> getYieldCurve();
    std::future<std::vector<EconomicIndicator>> getBulkIndicators(const std::vector<std::string>& symbols);
    
    std::string getName() const override { return "FRED"; }
    bool isAvailable() const override;
    
private:
    std::string makeRequest(const std::string& endpoint, 
                          const std::map<std::string, std::string>& params);
    void enforceRateLimit();
};

class YahooFinanceProvider : public IDataProvider {
private:
    std::string baseUrl_;
    std::atomic<std::chrono::system_clock::time_point> lastRequest_;
    
public:
    YahooFinanceProvider();
    
    std::future<EconomicIndicator> getIndicator(const std::string& symbol) override;
    std::future<double> getVIX();
    std::future<std::map<std::string, double>> getMajorIndices();
    std::future<MarketSentiment> getMarketSentiment();
    
    std::string getName() const override { return "Yahoo Finance"; }
    bool isAvailable() const override;
    
private:
    std::string makeRequest(const std::string& symbol);
    void enforceRateLimit();
};

class AlphaVantageProvider : public IDataProvider {
private:
    std::string apiKey_;
    std::string baseUrl_;
    std::atomic<std::chrono::system_clock::time_point> lastRequest_;
    
public:
    explicit AlphaVantageProvider(const std::string& apiKey);
    
    std::future<EconomicIndicator> getIndicator(const std::string& symbol) override;
    std::future<std::vector<EconomicIndicator>> getEconomicSuite();
    
    std::string getName() const override { return "Alpha Vantage"; }
    bool isAvailable() const override;
    std::chrono::minutes getRateLimit() const override { return std::chrono::minutes(1); }
    
private:
    std::string makeRequest(const std::string& endpoint,
                          const std::map<std::string, std::string>& params);
    void enforceRateLimit();
};

class EconomicDataManager {
private:
    std::vector<std::unique_ptr<IDataProvider>> providers_;
    std::map<std::string, EconomicIndicator> cachedIndicators_;
    MarketSentiment cachedSentiment_;
    TreasuryYieldCurve cachedYieldCurve_;
    
    std::chrono::system_clock::time_point lastUpdate_;
    std::chrono::minutes updateInterval_{30};
    
    mutable std::shared_mutex dataLock_;
    std::atomic<bool> updating_{false};
    std::thread updateThread_;
    std::atomic<bool> running_{false};
    
    std::vector<std::function<void(const EconomicIndicator&)>> indicatorCallbacks_;
    std::vector<std::function<void(const MarketSentiment&)>> sentimentCallbacks_;
    std::vector<std::function<void(const TreasuryYieldCurve&)>> yieldCurveCallbacks_;
    
public:
    EconomicDataManager();
    ~EconomicDataManager();
    
    EconomicDataManager(const EconomicDataManager&) = delete;
    EconomicDataManager& operator=(const EconomicDataManager&) = delete;
    
    void start();
    void stop();
    
    void addProvider(std::unique_ptr<IDataProvider> provider);
    void removeProvider(const std::string& name);
    
    std::future<void> updateAll();
    std::future<void> updateIndicators();
    std::future<void> updateSentiment();
    std::future<void> updateYieldCurve();
    
    EconomicIndicator getIndicator(const std::string& symbol) const;
    MarketSentiment getSentiment() const;
    TreasuryYieldCurve getYieldCurve() const;
    
    std::string analyzeMarketRegime() const;
    double calculateMarketStress() const;
    bool isRecessionary() const;
    StrategySignals getStrategySignals() const;
    
    std::vector<EconomicIndicator> getAllIndicators() const;
    std::map<std::string, double> getCorrelationMatrix() const;
    
    void addIndicatorCallback(std::function<void(const EconomicIndicator&)> callback);
    void addSentimentCallback(std::function<void(const MarketSentiment&)> callback);
    void addYieldCurveCallback(std::function<void(const TreasuryYieldCurve&)> callback);
    
    void setUpdateInterval(std::chrono::minutes interval);
    std::chrono::minutes getUpdateInterval() const { return updateInterval_; }
    
    struct DataQuality {
        double completeness = 0.0;
        double freshness = 0.0;
        double reliability = 0.0;
        std::chrono::system_clock::time_point lastAssessment;
    };
    
    DataQuality assessDataQuality() const;
    
private:
    void updateLoop();
    bool needsUpdate() const;
    void processIndicatorUpdate(const std::string& symbol, const EconomicIndicator& indicator);
    void calculateDerivedMetrics();
    void notifyCallbacks();
    
    double calculateIndicatorVolatility(const std::string& symbol) const;
    double calculateCorrelation(const std::string& symbol1, const std::string& symbol2) const;
};

class EconomicDataFactory {
public:
    static std::unique_ptr<EconomicDataManager> create(const std::map<std::string, std::string>& config);
    
private:
    static std::unique_ptr<IDataProvider> createFREDProvider(const std::string& apiKey);
    static std::unique_ptr<IDataProvider> createYahooProvider();
    static std::unique_ptr<IDataProvider> createAlphaVantageProvider(const std::string& apiKey);
};

}