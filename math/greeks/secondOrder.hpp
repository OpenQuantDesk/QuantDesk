#pragma once

#include "common/types.hpp"
#include "math/engine.hpp"
#include <memory>
#include <vector>
#include <thread>

namespace math {

class SecondOrderGreeks {
private:
    std::shared_ptr<MathEngine> mathEngine_;
    
public:
    explicit SecondOrderGreeks(std::shared_ptr<MathEngine> engine);
    ~SecondOrderGreeks() = default;
    
    SecondOrderGreeks(const SecondOrderGreeks&) = delete;
    SecondOrderGreeks& operator=(const SecondOrderGreeks&) = delete;
    
    ExtendedGreeks calculate(double spot, double strike, double timeToExpiry,
                           double riskFreeRate, double volatility, bool isCall,
                           double quantity = 1.0) const;
    
    ExtendedGreeks aggregatePortfolio(const std::vector<ExtendedGreeks>& positions) const;
    
    void calculateBatch(const std::vector<double>& spots,
                       const std::vector<double>& strikes,
                       const std::vector<double>& timeToExpiries,
                       const std::vector<double>& riskFreeRates,
                       const std::vector<double>& volatilities,
                       const std::vector<bool>& isCall,
                       const std::vector<double>& quantities,
                       std::vector<ExtendedGreeks>& results) const;
    
private:
    double calculatePinRisk(double spot, double strike, double timeToExpiry) const;
};

}