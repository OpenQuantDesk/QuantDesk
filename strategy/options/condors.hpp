#pragma once

#include "strategy/strategy_interface.hpp"

namespace strategy {

class IronCondorStrategy : public IStrategy {
private:
    double minIVRank_;
    double maxDaysToExpiration_;
    double targetProbabilityOTM_;

public:
    IronCondorStrategy();
    ~IronCondorStrategy() = default;
    
    std::string getName() const override { return "Iron Condor"; }
    
    StrategyOpportunity analyze(const std::string& underlying, 
                              double spot,
                              const common::OptionChain& chain,
                              const MarketConditions& conditions) const override;
    
    double getMinConfidence() const override { return 0.65; }
    bool isApplicable(const MarketConditions& conditions) const override;
    std::string getDescription() const override;
    std::vector<std::string> getRequiredMarketData() const override;
    
    void setMinIVRank(double minIVRank) { minIVRank_ = minIVRank; }
    void setMaxDaysToExpiration(double maxDaysToExpiration) { maxDaysToExpiration_ = maxDaysToExpiration; }
    void setTargetProbabilityOTM(double targetProbabilityOTM) { targetProbabilityOTM_ = targetProbabilityOTM; }
    
private:
    common::Quote findOptionByStrike(const std::vector<common::Quote>& options, double targetStrike) const;
    double calculateCreditReceived(const common::Quote& shortPut, const common::Quote& longPut,
                                 const common::Quote& shortCall, const common::Quote& longCall) const;
    double calculateMaxRisk(double creditReceived, double strikeWidth) const;
    double estimateProbabilityProfit(double spot, double shortPutStrike, double shortCallStrike,
                                   double timeToExpiry, double volatility) const;
};

}