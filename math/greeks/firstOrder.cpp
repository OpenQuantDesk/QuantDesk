#include "math/greeks/firstOrder.hpp"
#include "math/engine.hpp"
#include <cmath>

namespace math {

FirstOrderGreeks::FirstOrderGreeks(std::shared_ptr<MathEngine> engine) 
    : mathEngine_(engine) {}

Greeks FirstOrderGreeks::calculate(double spot, double strike, double timeToExpiry,
                                 double riskFreeRate, double volatility, bool isCall) const {
    return mathEngine_->blackScholes(spot, strike, timeToExpiry, riskFreeRate, volatility, isCall);
}

Greeks FirstOrderGreeks::calculateWithBumps(double spot, double strike, double timeToExpiry,
                                           double riskFreeRate, double volatility, bool isCall,
                                           double bumpSize) const {
    Greeks greeks;
    
    auto baseOption = mathEngine_->blackScholes(spot, strike, timeToExpiry, riskFreeRate, volatility, isCall);
    greeks.price = baseOption.price;
    
    auto spotUp = mathEngine_->blackScholes(spot + bumpSize, strike, timeToExpiry, riskFreeRate, volatility, isCall);
    auto spotDown = mathEngine_->blackScholes(spot - bumpSize, strike, timeToExpiry, riskFreeRate, volatility, isCall);
    greeks.delta = (spotUp.price - spotDown.price) / (2.0 * bumpSize);
    
    greeks.gamma = (spotUp.price - 2.0 * baseOption.price + spotDown.price) / (bumpSize * bumpSize);
    
    auto timeDown = mathEngine_->blackScholes(spot, strike, timeToExpiry - 1.0/365.0, riskFreeRate, volatility, isCall);
    greeks.theta = (timeDown.price - baseOption.price);
    
    auto volUp = mathEngine_->blackScholes(spot, strike, timeToExpiry, riskFreeRate, volatility + 0.01, isCall);
    greeks.vega = volUp.price - baseOption.price;
    
    auto rateUp = mathEngine_->blackScholes(spot, strike, timeToExpiry, riskFreeRate + 0.01, isCall);
    greeks.rho = rateUp.price - baseOption.price;
    
    greeks.impliedVol = volatility;
    
    return greeks;
}

void FirstOrderGreeks::calculateBatch(const std::vector<double>& spots,
                                    const std::vector<double>& strikes,
                                    const std::vector<double>& timeToExpiries,
                                    const std::vector<double>& riskFreeRates,
                                    const std::vector<double>& volatilities,
                                    const std::vector<bool>& isCall,
                                    std::vector<Greeks>& results) const {
    mathEngine_->blackScholesBatch(spots, strikes, timeToExpiries, riskFreeRates, volatilities, isCall, results);
}

Greeks FirstOrderGreeks::aggregatePortfolio(const std::vector<Greeks>& positions,
                                           const std::vector<double>& quantities) const {
    return mathEngine_->aggregateGreeks(positions, quantities);
}

double FirstOrderGreeks::calculatePinRisk(double spot, double strike, double timeToExpiry) const {
    if (timeToExpiry > 7.0/365.0) return 0.0;
    
    double distance = std::abs(spot - strike) / spot;
    double timeWeight = 1.0 - timeToExpiry / (7.0/365.0);
    
    return std::max(0.0, (0.02 - distance) * timeWeight);
}

double FirstOrderGreeks::calculateDollarDelta(const Greeks& greeks, double spot, 
                                             double quantity, double multiplier) const {
    return greeks.delta * spot * quantity * multiplier * 0.01;
}

double FirstOrderGreeks::calculateDollarGamma(const Greeks& greeks, double spot,
                                             double quantity, double multiplier) const {
    return greeks.gamma * spot * spot * quantity * multiplier * 0.01 * 0.01;
}

}