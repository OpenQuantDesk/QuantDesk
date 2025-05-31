/*
 * Filename: strategy_engine.hpp
 * Developer: Benjamin Cance
 * Date: 5/31/2025
 * 
 * Copyright 2025 Open Quant Desk, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include "strategy/strategy_opportunity.hpp"
#include "strategy/market_conditions.hpp"
#include "strategy/strategy_interface.hpp"
#include "common/quote.hpp"
#include "common/option_chain.hpp"
#include "OptionsQuantLib/core/engine.hpp"
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
    std::shared_ptr<math::core::MathEngine> mathEngine_;
    std::shared_ptr<data::EconomicDataManager> economicData_;
    
    mutable std::shared_mutex strategiesLock_;
    std::atomic<bool> analysisRunning_;

public:
    StrategyEngine(std::shared_ptr<math::core::MathEngine> mathEngine,
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