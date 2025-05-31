/*
 * Filename: strategy_engine.cpp
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

// To do#include "strategy/engine.hpp"
#include <algorithm>
#include <future>

namespace strategy {

StrategyEngine::StrategyEngine(std::shared_ptr<math::core::MathEngine> mathEngine,
                              std::shared_ptr<data::EconomicDataManager> economicData)
    : mathEngine_(mathEngine), economicData_(economicData) {
    initializeDefaultStrategies();
}

void StrategyEngine::addStrategy(std::unique_ptr<IStrategy> strategy) {
    strategies_.push_back(std::move(strategy));
}

void StrategyEngine::removeStrategy(const std::string& name) {
    strategies_.erase(
        std::remove_if(strategies_.begin(), strategies_.end(),
                      [&name](const auto& strategy) { return strategy->getName() == name; }),
        strategies_.end());
}

std::future<std::vector<StrategyOpportunity>> StrategyEngine::scanOpportunities(
    const std::map<std::string, common::OptionChain>& optionChains,
    const std::map<std::string, common::Quote>& underlyingQuotes) const {
    
    return std::async(std::launch::async, [this, &optionChains, &underlyingQuotes]() {
        std::vector<StrategyOpportunity> opportunities;
        MarketConditions conditions = assessMarketConditions(underlyingQuotes);
        
        std::mutex opportunitiesMutex;
        std::vector<std::future<void>> futures;
        
        for (const auto& [underlying, chain] : optionChains) {
            auto quoteIt = underlyingQuotes.find(underlying);
            if (quoteIt == underlyingQuotes.end() || !quoteIt->second.last.has_value()) {
                continue;
            }
            
            double spot = *quoteIt->second.last;
            
            futures.emplace_back(std::async(std::launch::async, [this, &opportunities, &opportunitiesMutex,
                                                               underlying, spot, &chain, &conditions]() {
                std::vector<StrategyOpportunity> localOpportunities;
                
                for (const auto& strategy : strategies_) {
                    try {
                        auto opportunity = strategy->analyze(underlying, spot, chain, conditions);
                        if (opportunity.confidence >= strategy->getMinConfidence()) {
                            localOpportunities.push_back(opportunity);
                        }
                    } catch (const std::exception& e) {
                        
                    }
                }
                
                std::lock_guard<std::mutex> lock(opportunitiesMutex);
                opportunities.insert(opportunities.end(), localOpportunities.begin(), localOpportunities.end());
            }));
        }
        
        for (auto& future : futures) {
            future.wait();
        }
        
        std::sort(opportunities.begin(), opportunities.end(),
                 [](const auto& a, const auto& b) {
                     return (a.expectedProfit / a.maxRisk * a.confidence) > 
                            (b.expectedProfit / b.maxRisk * b.confidence);
                 });
        
        return opportunities;
    });
}

std::future<StrategyOpportunity> StrategyEngine::analyzeSpecificStrategy(
    const std::string& strategyName,
    const std::string& underlying,
    const common::OptionChain& chain,
    const common::Quote& underlyingQuote) const {
    
    return std::async(std::launch::async, [this, strategyName, underlying, &chain, &underlyingQuote]() {
        if (!underlyingQuote.last.has_value()) {
            return StrategyOpportunity{};
        }
        
        double spot = *underlyingQuote.last;
        MarketConditions conditions = assessMarketConditions({{underlying, underlyingQuote}});
        
        for (const auto& strategy : strategies_) {
            if (strategy->getName() == strategyName) {
                return strategy->analyze(underlying, spot, chain, conditions);
            }
        }
        
        return StrategyOpportunity{};
    });
}

MarketConditions StrategyEngine::assessMarketConditions(
    const std::map<std::string, common::Quote>& quotes) const {
    
    MarketConditions conditions;
    
    auto vixIt = quotes.find("^VIX");
    if (vixIt != quotes.end() && vixIt->second.last.has_value()) {
        conditions.vixLevel = *vixIt->second.last;
        conditions.highVolEnvironment = conditions.vixLevel > 25.0;
    } else {
        conditions.vixLevel = 20.0;
        conditions.highVolEnvironment = false;
    }
    
    std::vector<double> impliedVols;
    for (const auto& [symbol, quote] : quotes) {
        if (quote.greeks.has_value()) {
            impliedVols.push_back(quote.greeks->impliedVol);
        }
    }
    
    if (!impliedVols.empty()) {
        conditions.ivRank = calculateIVRank(impliedVols);
        conditions.hvRank = calculateHVRank(extractPrices(quotes));
    }
    
    if (conditions.vixLevel > 30) {
        conditions.regime = "Fear";
    } else if (conditions.vixLevel < 15) {
        conditions.regime = "Greed";
    } else {
        conditions.regime = "Neutral";
    }
    
    conditions.correlations = calculateCorrelations(quotes);
    
    return conditions;
}

std::vector<std::string> StrategyEngine::getAvailableStrategies() const {
    std::vector<std::string> names;
    names.reserve(strategies_.size());
    
    for (const auto& strategy : strategies_) {
        names.push_back(strategy->getName());
    }
    
    return names;
}

void StrategyEngine::initializeDefaultStrategies() {
    strategies_.push_back(std::make_unique<IronCondorStrategy>());
    strategies_.push_back(std::make_unique<ShortStrangleStrategy>());
    strategies_.push_back(std::make_unique<CalendarSpreadStrategy>());
    strategies_.push_back(std::make_unique<LongStraddleStrategy>());
}

double StrategyEngine::calculateIVRank(const std::vector<double>& impliedVols) const {
    if (impliedVols.empty()) return 0.5;
    
    auto minMax = std::minmax_element(impliedVols.begin(), impliedVols.end());
    double current = impliedVols.back();
    double range = *minMax.second - *minMax.first;
    
    if (range <= 0) return 0.5;
    
    return (current - *minMax.first) / range;
}

double StrategyEngine::calculateHVRank(const std::vector<double>& prices) const {
    if (prices.size() < 20) return 0.5;
    
    return mathEngine_->realizedVolatility(prices) / 0.3;
}

std::vector<double> StrategyEngine::extractPrices(const std::map<std::string, common::Quote>& quotes) const {
    std::vector<double> prices;
    
    for (const auto& [symbol, quote] : quotes) {
        if (quote.last.has_value()) {
            prices.push_back(*quote.last);
        }
    }
    
    return prices;
}

double StrategyEngine::calculateCorrelations(const std::map<std::string, common::Quote>& quotes) const {
    return 0.7;
}

StrategyOpportunity IronCondorStrategy::analyze(
    const std::string& underlying,
    double spot,
    const common::OptionChain& chain,
    const MarketConditions& conditions) const {
    
    StrategyOpportunity opportunity;
    opportunity.name = getName();
    opportunity.underlying = underlying;
    
    if (!conditions.highVolEnvironment || conditions.ivRank < 0.6) {
        opportunity.confidence = 0.3;
        opportunity.reasoning = "Low IV environment not suitable for premium selling";
        return opportunity;
    }
    
    std::vector<common::Quote> calls = chain.calls;
    std::vector<common::Quote> puts = chain.puts;
    
    auto findClosestStrike = [spot](const std::vector<common::Quote>& options, double target) -> common::Quote {
        auto it = std::min_element(options.begin(), options.end(),
            [target](const auto& a, const auto& b) {
                return std::abs(a.symbol.find(std::to_string(static_cast<int>(target))) - target) <
                       std::abs(b.symbol.find(std::to_string(static_cast<int>(target))) - target);
            });
        return (it != options.end()) ? *it : common::Quote{};
    };
    
    double putStrike1 = spot * 0.95;
    double putStrike2 = spot * 0.90;
    double callStrike1 = spot * 1.05;
    double callStrike2 = spot * 1.10;
    
    auto shortPut = findClosestStrike(puts, putStrike1);
    auto longPut = findClosestStrike(puts, putStrike2);
    auto shortCall = findClosestStrike(calls, callStrike1);
    auto longCall = findClosestStrike(calls, callStrike2);
    
    double totalCredit = 0.0;
    if (shortPut.bid.has_value()) totalCredit += *shortPut.bid;
    if (shortCall.bid.has_value()) totalCredit += *shortCall.bid;
    if (longPut.ask.has_value()) totalCredit -= *longPut.ask;
    if (longCall.ask.has_value()) totalCredit -= *longCall.ask;
    
    opportunity.expectedProfit = totalCredit * 100;
    opportunity.maxRisk = (500 - totalCredit) * 100;
    opportunity.probabilityProfit = 0.68;
    opportunity.sharpeRatio = opportunity.expectedProfit / opportunity.maxRisk;
    opportunity.confidence = conditions.ivRank > 0.75 ? 0.85 : 0.70;
    opportunity.suitableForCurrentVol = true;
    opportunity.suitableForCurrentRegime = conditions.regime != "Fear";
    opportunity.optimalIVRank = 0.75;
    
    opportunity.reasoning = "High IV rank (" + std::to_string(static_cast<int>(conditions.ivRank * 100)) + 
                           "%) favors premium selling strategies";
    
    opportunity.legs = {
        shortPut.symbol + " SELL",
        longPut.symbol + " BUY", 
        shortCall.symbol + " SELL",
        longCall.symbol + " BUY"
    };
    
    return opportunity;
}

StrategyOpportunity ShortStrangleStrategy::analyze(
    const std::string& underlying,
    double spot,
    const common::OptionChain& chain,
    const MarketConditions& conditions) const {
    
    StrategyOpportunity opportunity;
    opportunity.name = getName();
    opportunity.underlying = underlying;
    
    if (!conditions.highVolEnvironment) {
        opportunity.confidence = 0.4;
        opportunity.reasoning = "Normal volatility environment";
        return opportunity;
    }
    
    auto findOTMOption = [](const std::vector<common::Quote>& options, double spot, bool isCall) -> common::Quote {
        double targetDelta = isCall ? 0.3 : -0.3;
        
        auto it = std::min_element(options.begin(), options.end(),
            [targetDelta](const auto& a, const auto& b) {
                double deltaA = a.greeks.has_value() ? a.greeks->delta : (a.symbol.find("C") != std::string::npos ? 0.5 : -0.5);
                double deltaB = b.greeks.has_value() ? b.greeks->delta : (b.symbol.find("C") != std::string::npos ? 0.5 : -0.5);
                return std::abs(deltaA - targetDelta) < std::abs(deltaB - targetDelta);
            });
        
        return (it != options.end()) ? *it : common::Quote{};
    };
    
    auto shortPut = findOTMOption(chain.puts, spot, false);
    auto shortCall = findOTMOption(chain.calls, spot, true);
    
    double totalCredit = 0.0;
    if (shortPut.bid.has_value()) totalCredit += *shortPut.bid;
    if (shortCall.bid.has_value()) totalCredit += *shortCall.bid;
    
    opportunity.expectedProfit = totalCredit * 100;
    opportunity.maxRisk = (spot * 0.2 - totalCredit) * 100;
    opportunity.probabilityProfit = 0.60;
    opportunity.sharpeRatio = opportunity.expectedProfit / opportunity.maxRisk;
    opportunity.confidence = conditions.ivRank > 0.8 ? 0.80 : 0.65;
    opportunity.suitableForCurrentVol = true;
    opportunity.suitableForCurrentRegime = conditions.regime == "Neutral";
    opportunity.optimalIVRank = 0.80;
    
    opportunity.reasoning = "High volatility with neutral market outlook";
    
    opportunity.legs = {
        shortPut.symbol + " SELL",
        shortCall.symbol + " SELL"
    };
    
    return opportunity;
}

StrategyOpportunity CalendarSpreadStrategy::analyze(
    const std::string& underlying,
    double spot,
    const common::OptionChain& chain,
    const MarketConditions& conditions) const {
    
    StrategyOpportunity opportunity;
    opportunity.name = getName();
    opportunity.underlying = underlying;
    opportunity.confidence = 0.65;
    opportunity.expectedProfit = 150.0;
    opportunity.maxRisk = 300.0;
    opportunity.probabilityProfit = 0.55;
    opportunity.reasoning = "Time decay advantage in neutral environment";
    opportunity.suitableForCurrentVol = conditions.ivRank > 0.4 && conditions.ivRank < 0.8;
    opportunity.suitableForCurrentRegime = conditions.regime == "Neutral";
    opportunity.optimalIVRank = 0.60;
    
    return opportunity;
}

StrategyOpportunity LongStraddleStrategy::analyze(
    const std::string& underlying,
    double spot,
    const common::OptionChain& chain,
    const MarketConditions& conditions) const {
    
    StrategyOpportunity opportunity;
    opportunity.name = getName();
    opportunity.underlying = underlying;
    
    if (conditions.ivRank > 0.5) {
        opportunity.confidence = 0.3;
        opportunity.reasoning = "High IV makes long premium expensive";
        return opportunity;
    }
    
    opportunity.confidence = 0.75;
    opportunity.expectedProfit = 300.0;
    opportunity.maxRisk = 250.0;
    opportunity.probabilityProfit = 0.45;
    opportunity.reasoning = "Low IV with volatility expansion expected";
    opportunity.suitableForCurrentVol = conditions.ivRank < 0.3;
    opportunity.suitableForCurrentRegime = conditions.regime == "Fear" || conditions.regime == "Greed";
    opportunity.optimalIVRank = 0.25;
    
    return opportunity;
}