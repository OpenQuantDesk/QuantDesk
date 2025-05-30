#include "portfolio/manager.hpp"
#include "math/greeks/secondOrder.hpp"
#include "math/metrics/probability.hpp"

namespace portfolio {

void Position::updateAnalytics(double underlyingPx, double vol, double riskFreeRate,
                                     const math::MathEngine& engine) {
    underlyingPrice = underlyingPx;
    impliedVol = vol;
    lastUpdate = std::chrono::system_clock::now();
    
    double timeToExpiry = getDaysToExpiry() / 365.0;
    
    if (timeToExpiry > 0 && vol > 0 && underlyingPrice > 0 && strike > 0) {
        greeks = engine.blackScholes(underlyingPrice, strike, timeToExpiry, 
                                   riskFreeRate, vol, isCall);
        theoreticalValue = greeks.price;
        
        math::SecondOrderGreeks advGreeksCalc(std::make_shared<math::MathEngine>(engine));
        extendedGreeks = advGreeksCalc.calculate(underlyingPrice, strike, timeToExpiry,
                                               riskFreeRate, vol, isCall, quantity);
        
        math::ProbabilityAnalyzer probAnalyzer;
        probMetrics = probAnalyzer.analyzeOption(underlyingPrice, strike, vol,
                                                timeToExpiry, isCall, currentPrice);
    }
    
    unrealizedPnL = (currentPrice - avgPrice) * quantity * 100;
    deltaAdjustedExposure = greeks.delta * quantity * underlyingPrice;
    gammaRisk = std::abs(greeks.gamma * quantity * underlyingPrice * underlyingPrice * 0.01);
    vegaRisk = std::abs(greeks.vega * quantity);
    thetaDecay = greeks.theta * quantity;
}

double Position::getDaysToExpiry() const {
    auto now = std::chrono::system_clock::now();
    auto expiry = parseExpirationDate(expiration);
    auto duration = std::chrono::duration_cast<std::chrono::hours>(expiry - now);
    return duration.count() / 24.0;
}

double Position::getNotionalValue() const {
    return currentPrice * std::abs(quantity) * 100;
}

double Position::getImpliedProbability() const {
    return probMetrics.probabilityITM;
}

double Position::getMoneyness() const {
    if (underlyingPrice <= 0 || strike <= 0) return 1.0;
    
    if (isCall) {
        return underlyingPrice / strike;
    } else {
        return strike / underlyingPrice;
    }
}

bool Position::isExpiring(int days) const {
    return getDaysToExpiry() <= days;
}

std::chrono::system_clock::time_point Position::parseExpirationDate(const std::string& expiration) const {
    std::tm tm = {};
    std::istringstream ss(expiration);
    ss >> std::get_time(&tm, "%Y-%m-%d");
    
    if (ss.fail()) {
        return std::chrono::system_clock::now() + std::chrono::hours(24 * 30);
    }
    
    return std::chrono::system_clock::from_time_t(std::mktime(&tm));
}

}