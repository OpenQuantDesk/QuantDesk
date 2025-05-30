#pragma once

#include "common/types.hpp"
#include "math/engine.hpp"
#include "data/economic.hpp"
#include <vector>
#include <memory>
#include <future>
#include <functional>
#include <string>
#include <map>
#include <atomic>
#include <shared_mutex>

namespace strategy {

struct StrategyOpportunity {
    std::string name;
    std::string underlying;
    std::vector<std::string> legs;
    double expectedProfit = 0.0;
    double maxRisk = 0.0;
    double probabilityProfit = 0.0;
    double sharpeRatio = 0.0;
    std::string reasoning;
    double confidence = 0.0;
    bool suitableForCurrentVol = false;
    bool suitableForCurrentRegime = false;
    double optimalIVRank = 0.0;
    std::vector<common::OrderRequest> orderRequests;
};

}

struct MarketConditions {
    double vixLevel = 0.0;
    double ivRank = 0.0;
    double hvRank = 0.0;
    std::string regime;
    bool highVolEnvironment = false;
    double correlations = 0.0;
};

class IStrategy {
public:
    virtual ~IStrategy() = default;
    virtual std::string getName() const = 0;
    virtual StrategyOpportunity analyze(
        const std::string& underlying,
        double spot,
        const common::OptionChain& chain,
        const MarketConditions& conditions) const = 0;
    virtual double getMinConfidence() const { return 0.6; }
};

class IronCondorStrategy : public IStrategy {
public:
    std::string getName() const override { return "Iron Condor"; }
    StrategyOpportunity analyze(
        const std::string& underlying,
        double spot,
        const common::OptionChain& chain,
        const MarketConditions& conditions) const override;
};

class ShortStrangleStrategy : public IStrategy {
public:
    std::string getName() const override { return "Short Strangle"; }
    StrategyOpportunity analyze(
        const std::string& underlying,
        double spot,
        const common::OptionChain& chain,
        const MarketConditions& conditions) const override;
};

class CalendarSpreadStrategy : public IStrategy {
public:
    std::string getName() const override { return "Calendar Spread"; }
    StrategyOpportunity analyze(
        const std::string& underlying,
        double spot,
        const common::OptionChain& chain,
        const MarketConditions& conditions) const override;
};

class LongStraddleStrategy : public IStrategy {
public:
    std::string getName() const override { return "Long Straddle"; }
    StrategyOpportunity analyze(
        const std::string& underlying,
        double spot,
        const common::OptionChain& chain,
        const MarketConditions& conditions) const override;
};

class StrategyEngine {
private:
    std::vector<std::unique_ptr<IStrategy>> strategies_;
    std::shared_ptr<math::MathEngine> mathEngine_;
    std::shared_ptr<data::EconomicDataManager> economicData_;
    
    mutable std::shared_mutex strategiesLock_;
    std::atomic<bool> analysisRunning_{false};
    
public:
    StrategyEngine(std::shared_ptr<math::MathEngine> mathEngine,
                  std::shared_ptr<data::EconomicDataManager> economicData);
    ~StrategyEngine() = default;
    
    StrategyEngine(const StrategyEngine&) = delete;
    StrategyEngine& operator=(const StrategyEngine&) = delete;
    
    void addStrategy(std::unique_ptr<IStrategy> strategy);
    void removeStrategy(const std::string& name);
    
    std::future<std::vector<StrategyOpportunity>> scanOpportunities(
        const std::map<std::string, common::OptionChain>& optionChains,
        const std::map<std::string, common::Quote>& underlyingQuotes) const;
    
    std::future<StrategyOpportunity> analyzeSpecificStrategy(
        const std::string& strategyName,
        const std::string& underlying,
        const common::OptionChain& chain,
        const common::Quote& underlyingQuote) const;
    
    MarketConditions assessMarketConditions(
        const std::map<std::string, common::Quote>& quotes) const;
    
    std::vector<std::string> getAvailableStrategies() const;
    
private:
    void initializeDefaultStrategies();
    double calculateIVRank(const std::vector<double>& impliedVols) const;
    double calculateHVRank(const std::vector<double>& prices) const;
    std::vector<double> extractPrices(const std::map<std::string, common::Quote>& quotes) const;
    double calculateCorrelations(const std::map<std::string, common::Quote>& quotes) const;
};

class StrategyOptimizer {
private:
    std::shared_ptr<math::MathEngine> mathEngine_;
    
    struct OptimizationTarget {
        double maxRisk = 1000.0;
        double minProbability = 0.60;
        double minExpectedReturn = 0.10;
        double maxTimeDecay = -50.0;
        bool preferHighProbability = true;
    };
    
public:
    explicit StrategyOptimizer(std::shared_ptr<math::MathEngine> mathEngine);
    
    struct OptimizationResult {
        std::vector<common::OrderRequest> optimalOrders;
        double expectedProfit = 0.0;
        double maxRisk = 0.0;
        double probabilityProfit = 0.0;
        double score = 0.0;
    };
    
    std::future<OptimizationResult> optimizeStrategy(
        const StrategyOpportunity& opportunity,
        const OptimizationTarget& target = {}) const;
    
    std::future<std::vector<OptimizationResult>> optimizePortfolio(
        const std::vector<StrategyOpportunity>& opportunities,
        double availableCapital,
        const OptimizationTarget& target = {}) const;
    
private:
    double scoreStrategy(const StrategyOpportunity& opportunity,
                        const OptimizationTarget& target) const;
};

class BacktestEngine {
private:
    std::shared_ptr<math::MathEngine> mathEngine_;
    
public:
    explicit BacktestEngine(std::shared_ptr<math::MathEngine> mathEngine);
    
    struct BacktestParams {
        std::chrono::system_clock::time_point startDate;
        std::chrono::system_clock::time_point endDate;
        double initialCapital = 100000.0;
        std::vector<std::string> strategies;
        std::vector<std::string> underlyings;
    };
    
    struct BacktestResult {
        std::string strategyName;
        double totalReturn = 0.0;
        double sharpeRatio = 0.0;
        double maxDrawdown = 0.0;
        double winRate = 0.0;
        double avgWin = 0.0;
        double avgLoss = 0.0;
        double profitFactor = 0.0;
        std::vector<double> equityCurve;
        std::vector<std::string> trades;
    };
    
    std::future<std::vector<BacktestResult>> runBacktest(
        const BacktestParams& params) const;
    
private:
    BacktestResult simulateStrategy(
        const std::string& strategy,
        const std::vector<std::string>& underlyings,
        const BacktestParams& params) const;
};