#pragma once

#include <atomic>
#include <chrono>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

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
    };

    struct MarketSentiment {
        double vixLevel = 0.0;
        double putCallRatio = 0.0;
        double fearGreedIndex = 0.0;
        std::string marketRegime;
        double correlationSPY_VIX = 0.0;
    };

    struct TreasuryYieldCurve {
        std::map<int, double> yields;
        double slope = 0.0;
        double curvature = 0.0;
        bool inverted = false;
    };

    struct StrategySignals {
        bool highVolatilityEnvironment = false;
        bool yieldCurveInversion = false;
        bool recessionary = false;
        double optimalIVRank = 0.5;
        std::string recommendedStrategy;
        std::vector<std::string> avoidStrategies;
    };

    class IDataProvider {
    public:
        virtual ~IDataProvider() = default;
        virtual std::future<EconomicIndicator>
        getIndicator(const std::string& symbol) = 0;
        virtual std::string getName() const = 0;
    };

    class FREDProvider : public IDataProvider {
    private:
        std::string apiKey_;
        std::string baseUrl_;
    public:
        explicit FREDProvider(const std::string& apiKey);
        std::future<EconomicIndicator>
        getIndicator(const std::string& symbol) override;
        std::future<TreasuryYieldCurve> getYieldCurve();
        std::string getName() const override { return "FRED"; }
    private:
        std::string
        makeRequest(const std::string& endpoint,
                    const std::map<std::string, std::string>& params);
    };

    class YahooFinanceProvider : public IDataProvider {
    private:
        std::string baseUrl_;
    public:
        YahooFinanceProvider();
        std::future<EconomicIndicator>
        getIndicator(const std::string& symbol) override;
        std::future<double> getVIX();
        std::future<std::map<std::string, double>> getMajorIndices();
        std::string getName() const override { return "Yahoo Finance"; }
    private:
        std::string makeRequest(const std::string& symbol);
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
    public:
        EconomicDataManager();
        ~EconomicDataManager() = default;

        EconomicDataManager(const EconomicDataManager&) = delete;
        EconomicDataManager& operator=(const EconomicDataManager&) = delete;

        void addProvider(std::unique_ptr<IDataProvider> provider);

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
    private:
        bool needsUpdate() const;
        void processIndicatorUpdate(const std::string& symbol,
                                    const EconomicIndicator& indicator);
    };

    class EconomicDataFactory {
    public:
        static std::unique_ptr<EconomicDataManager>
        create(const std::map<std::string, std::string>& config);
    private:
        static std::unique_ptr<IDataProvider>
        createFREDProvider(const std::string& apiKey);
        static std::unique_ptr<IDataProvider> createYahooProvider();
    };

} // namespace data