#pragma once

#include "common/types.hpp"
#include <memory>
#include <vector>

namespace math {

class MathEngine;

class FirstOrderGreeks {
private:
    std::shared_ptr<MathEngine> mathEngine_;
    
public:
    explicit FirstOrderGreeks(std::shared_ptr<MathEngine> engine);
    ~FirstOrderGreeks() = default;
    
    FirstOrderGreeks(const FirstOrderGreeks&) = delete;
    FirstOrderGreeks& operator=(const FirstOrderGreeks&) = delete;
    
    common::Greeks calculate(double spot, double strike, double timeToExpiry,
                           double riskFreeRate, double volatility, bool isCall) const;
    
    common::Greeks calculateWithBumps(double spot, double strike, double timeToExpiry,
                                    double riskFreeRate, double volatility, bool isCall,
                                    double bumpSize = 0.01) const;
    
    void calculateBatch(const std::vector<double>& spots,
                       const std::vector<double>& strikes,
                       const std::vector<double>& timeToExpiries,
                       const std::vector<double>& riskFreeRates,
                       const std::vector<double>& volatilities,
                       const std::vector<bool>& isCall,
                       std::vector<common::Greeks>& results) const;
    
    common::Greeks aggregatePortfolio(const std::vector<common::Greeks>& positions,
                                    const std::vector<double>& quantities) const;
    
    double calculatePinRisk(double spot, double strike, double timeToExpiry) const;
    
    double calculateDollarDelta(const common::Greeks& greeks, double spot, 
                              double quantity, double multiplier = 100.0) const;
    
    double calculateDollarGamma(const common::Greeks& greeks, double spot,
                              double quantity, double multiplier = 100.0) const;
};

}