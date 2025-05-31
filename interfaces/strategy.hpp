/*
 * Filename: strategy.hpp
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

#include "strategy/opportunity.hpp"
#include "market/conditions.hpp"
#include "common/option_chain.hpp"
#include <string>

namespace strategy {

class IStrategy {
public:
    virtual ~IStrategy() = default;
    
    virtual std::string getName() const = 0;
    virtual StrategyOpportunity analyze(const std::string& underlying, 
                                      double spot,
                                      const common::OptionChain& chain,
                                      const MarketConditions& conditions) const = 0;
    virtual double getMinConfidence() const { return 0.6; }
    virtual bool isApplicable(const MarketConditions& conditions) const { return true; }
    virtual std::string getDescription() const { return ""; }
    virtual std::vector<std::string> getRequiredMarketData() const { return {}; }
};

}