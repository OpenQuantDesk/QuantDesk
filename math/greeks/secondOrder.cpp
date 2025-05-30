#include "math/greeks/secondOrder.hpp"
#include "math/engine.hpp"
#include <cmath>

namespace math {

SecondOrderGreeks::SecondOrderGreeks(std::shared_ptr<MathEngine> engine) 
    : mathEngine_(engine) {}

ExtendedGreeks SecondOrderGreeks::calculate(double spot, double strike, double timeToExpiry,
                                          double riskFreeRate, double volatility, bool isCall,
                                          double quantity) const {
    ExtendedGreeks greeks;
    
    const double spotBump = spot * 0.01;
    const double volBump = 0.01;
    const double timeBump = 1.0 / 365.0;
    
    auto base = mathEngine_->blackScholes(spot, strike, timeToExpiry, riskFreeRate, volatility, isCall);
    greeks.price = base.price * quantity;
    greeks.delta = base.delta * quantity;
    greeks.gamma = base.gamma * quantity;
    greeks.theta = base.theta * quantity;
    greeks.vega = base.vega * quantity;
    greeks.rho = base.rho * quantity;
    
    auto volUp = mathEngine_->blackScholes(spot, strike, timeToExpiry, riskFreeRate, volatility + volBump, isCall);
    auto volDown = mathEngine_->blackScholes(spot, strike, timeToExpiry, riskFreeRate, volatility - volBump, isCall);
    greeks.volga = (volUp.vega - volDown.vega) / (2.0 * volBump) * quantity;
    
    auto spotUpVolUp = mathEngine_->blackScholes(spot + spotBump, strike, timeToExpiry, riskFreeRate, volatility + volBump, isCall);
    auto spotDownVolDown = mathEngine_->blackScholes(spot - spotBump, strike, timeToExpiry, riskFreeRate, volatility - volBump, isCall);
    auto spotUpVolDown = mathEngine_->blackScholes(spot + spotBump, strike, timeToExpiry, riskFreeRate, volatility - volBump, isCall);
    auto spotDownVolUp = mathEngine_->blackScholes(spot - spotBump, strike, timeToExpiry, riskFreeRate, volatility + volBump, isCall);
    
    greeks.vanna = (spotUpVolUp.delta - spotUpVolDown.delta - spotDownVolUp.delta + spotDownVolDown.delta) / 
                   (4.0 * spotBump * volBump) * quantity;
    
    auto timeDown = mathEngine_->blackScholes(spot, strike, timeToExpiry - timeBump, riskFreeRate, volatility, isCall);
    auto spotUpTimeDown = mathEngine_->blackScholes(spot + spotBump, strike, timeToExpiry - timeBump, riskFreeRate, volatility, isCall);
    auto spotDownTimeDown = mathEngine_->blackScholes(spot - spotBump, strike, timeToExpiry - timeBump, riskFreeRate, volatility, isCall);
    
    double deltaTimeDown = (spotUpTimeDown.price - spotDownTimeDown.price) / (2.0 * spotBump);
    greeks.charm = (base.delta - deltaTimeDown) / timeBump * quantity;
    
    auto spotUp = mathEngine_->blackScholes(spot + spotBump, strike, timeToExpiry, riskFreeRate, volatility, isCall);
    auto spotDown = mathEngine_->blackScholes(spot - spotBump, strike, timeToExpiry, riskFreeRate, volatility, isCall);
    auto spotUpUp = mathEngine_->blackScholes(spot + 2.0 * spotBump, strike, timeToExpiry, riskFreeRate, volatility, isCall);
    
    greeks.speed = (spotUpUp.gamma - 2.0 * spotUp.gamma + base.gamma) / (spotBump * spotBump) * quantity;
    
    auto volUpGamma = (spotUp.gamma + volUp.gamma - base.gamma) / volBump;
    auto volDownGamma = (spotUp.gamma + volDown.gamma - base.gamma) / (-volBump);
    greeks.zomma = (volUpGamma - volDownGamma) / (2.0 * volBump) * quantity;
    
    auto timeDownGamma = (spotUpTimeDown.gamma + spotDownTimeDown.gamma - 2.0 * timeDown.gamma) / (spotBump * spotBump);
    greeks.color = (base.gamma - timeDownGamma) / timeBump * quantity;
    
    greeks.dollarDelta = greeks.delta * spot * 0.01;
    greeks.dollarGamma = greeks.gamma * spot * spot * 0.01 * 0.01;
    
    greeks.pinRisk = calculatePinRisk(spot, strike, timeToExpiry) * std::abs(quantity);
    
    return greeks;
}

ExtendedGreeks SecondOrderGreeks::aggregatePortfolio(const std::vector<ExtendedGreeks>& positions) const {
    ExtendedGreeks total;
    
    for (const auto& pos : positions) {
        total.price += pos.price;
        total.delta += pos.delta;
        total.gamma += pos.gamma;
        total.theta += pos.theta;
        total.vega += pos.vega;
        total.rho += pos.rho;
        total.volga += pos.volga;
        total.vanna += pos.vanna;
        total.charm += pos.charm;
        total.speed += pos.speed;
        total.zomma += pos.zomma;
        total.color += pos.color;
        total.dollarDelta += pos.dollarDelta;
        total.dollarGamma += pos.dollarGamma;
        total.pinRisk += pos.pinRisk;
    }
    
    return total;
}

double SecondOrderGreeks::calculatePinRisk(double spot, double strike, double timeToExpiry) const {
    if (timeToExpiry > 7.0/365.0) return 0.0;
    
    double distance = std::abs(spot - strike) / spot;
    double timeWeight = 1.0 - timeToExpiry / (7.0/365.0);
    double pinDistance = 0.02;
    
    if (distance < pinDistance) {
        return (pinDistance - distance) / pinDistance * timeWeight;
    }
    
    return 0.0;
}

void SecondOrderGreeks::calculateBatch(const std::vector<double>& spots,
                                     const std::vector<double>& strikes,
                                     const std::vector<double>& timeToExpiries,
                                     const std::vector<double>& riskFreeRates,
                                     const std::vector<double>& volatilities,
                                     const std::vector<bool>& isCall,
                                     const std::vector<double>& quantities,
                                     std::vector<ExtendedGreeks>& results) const {
    const size_t n = spots.size();
    results.resize(n);
    
    const size_t numThreads = std::thread::hardware_concurrency();
    const size_t chunkSize = (n + numThreads - 1) / numThreads;
    
    std::vector<std::thread> threads;
    threads.reserve(numThreads);
    
    for (size_t i = 0; i < numThreads; ++i) {
        const size_t start = i * chunkSize;
        const size_t end = std::min(start + chunkSize, n);
        
        if (start >= n) break;
        
        threads.emplace_back([this, &spots, &strikes, &timeToExpiries, &riskFreeRates,
                             &volatilities, &isCall, &quantities, &results, start, end]() {
            for (size_t j = start; j < end; ++j) {
                results[j] = calculate(spots[j], strikes[j], timeToExpiries[j],
                                     riskFreeRates[j], volatilities[j], isCall[j], quantities[j]);
            }
        });
    }
    
    for (auto& thread : threads) {
        thread.join();
    }
}