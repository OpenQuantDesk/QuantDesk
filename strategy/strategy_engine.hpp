#pragma once

#include "strategy/strategy_opportunity.hpp"
#include "strategy/market_conditions.hpp"
#include "strategy/strategy_interface.hpp"
#include "common/quote.hpp"
#include "common/option_chain.hpp"
#include "math/math_engine.hpp"
#include "data/economic_data_manager.hpp"
#include <vector>
#include <memory>
#include <future>
#include <string>
#include <map>
#include <atomic>
#include <shared_mutex>

namespace strategy {

class StrategyEngine {
private:
    std::vector<std::unique_ptr<IStrategy>> strategies_;
    std::shared_ptr<math::MathEngine> mathEngine_;
    std::shared_ptr<data::EconomicDataManager> economicData_;
    
    mutable std::shared_mutex strategiesLock_;
    std::atomic<bool> analysisRunning_;

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
        const std::string& strategyName, const std::string& underlying,
        const common::OptionChain& chain, const common::Quote& underlyingQuote) const;
    
    MarketConditions assessMarketConditions(const std::map<std::string, common::Quote>& quotes) const;
    std::vector<std::string> getAvailableStrategies() const;
    
    bool isAnalysisRunning() const { return analysisRunning_.load(); }
    
private:
    void initializeDefaultStrategies();
    double calculateIVRank(const std::vector<double>& impliedVols) const;
    double calculateHVRank(const std::vector<double>& prices) const;
    std::vector<double> extractPrices(const std::map<std::string, common::Quote>& quotes) const;
    double calculateCorrelations(const std::map<std::string, common::Quote>& quotes) const;
    void setAnalysisRunning(bool running) const { analysisRunning_.store(running); }
};

}