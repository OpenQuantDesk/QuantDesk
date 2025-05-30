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