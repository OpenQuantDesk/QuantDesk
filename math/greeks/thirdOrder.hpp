/* These are separate because they were not in the original design */

#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include <memory>

namespace math::greeks {

struct ThirdOrderGreeks {
    core::Real speed = 0.0;        // d³V/dS³ - third derivative with respect to spot
    core::Real zomma = 0.0;        // d²V/dS/dσ - gamma sensitivity to volatility 
    core::Real color = 0.0;        // d²V/dS²/dt - gamma decay with time
    core::Real volga = 0.0;        // d²V/dσ² - vega convexity
    core::Real vanna = 0.0;        // d²V/dS/dσ - delta sensitivity to volatility
    core::Real charm = 0.0;        // d²V/dS/dt - delta decay with time
    core::Real veta = 0.0;         // d²V/dσ/dt - vega decay with time
    core::Real vera = 0.0;         // d²V/dσ/dr - vega sensitivity to interest rate
    core::Real lambda = 0.0;       // d²V/dS/dr - delta sensitivity to interest rate
    
    // True third-order terms
    core::Real ultima = 0.0;       // d³V/dσ³ - third derivative with respect to volatility
    core::Real totto = 0.0;        // d³V/dS²/dσ - speed sensitivity to volatility
    core::Real dvannaDvol = 0.0;   // d³V/dS/dσ² - vanna sensitivity to volatility  
    core::Real dcharmDspot = 0.0;  // d³V/dS²/dt - charm sensitivity to spot
    core::Real dzommaDvol = 0.0;   // d³V/dS/dσ² - zomma sensitivity to volatility
    core::Real dcolorDspot = 0.0;  // d³V/dS³/dt - color sensitivity to spot
    core::Real dvolgaDvol = 0.0;   // d³V/dσ³ - volga sensitivity to volatility (same as ultima)
    core::Real dvetaDvol = 0.0;    // d³V/dσ²/dt - veta sensitivity to volatility
    core::Real dveraDrate = 0.0;   // d³V/dσ/dr² - vera sensitivity to interest rate
    core::Real dlambdaDspot = 0.0; // d³V/dS²/dr - lambda sensitivity to spot
    
    // Cross third-order terms
    core::Real dspeedDvol = 0.0;   // d³V/dS³/dσ - speed sensitivity to volatility
    core::Real dspeedDtime = 0.0;  // d³V/dS³/dt - speed sensitivity to time
    core::Real dzommaDtime = 0.0;  // d³V/dS/dσ/dt - zomma sensitivity to time
    core::Real dcolorDvol = 0.0;   // d³V/dS²/dt/dσ - color sensitivity to volatility
    core::Real dvannaDtime = 0.0;  // d³V/dS/dσ/dt - vanna sensitivity to time
    core::Real dcharmDvol = 0.0;   // d³V/dS/dt/dσ - charm sensitivity to volatility
};

class ThirdOrderCalculator {
private:
    std::shared_ptr<core::MathEngine> engine_;
    core::Real defaultBumpSize_;
    
public:
    explicit ThirdOrderCalculator(std::shared_ptr<core::MathEngine> engine, 
                                 core::Real bumpSize = 0.01);
    
    ThirdOrderGreeks calculate(const core::MarketData& market, 
                              core::Real quantity = 1.0) const;
    
    ThirdOrderGreeks calculateAnalytical(const core::MarketData& market,
                                        core::Real quantity = 1.0) const;
    
    void calculateBatch(const core::MarketDataVector& markets,
                       const core::Vector& quantities,
                       std::vector<ThirdOrderGreeks>& results) const;
    
    // Individual third-order Greek calculations
    core::Real calculateSpeed(const core::MarketData& market) const;
    core::Real calculateZomma(const core::MarketData& market) const;
    core::Real calculateColor(const core::MarketData& market) const;
    core::Real calculateUltima(const core::MarketData& market) const;
    core::Real calculateTotto(const core::MarketData& market) const;
    core::Real calculateDvannaDvol(const core::MarketData& market) const;
    core::Real calculateDcharmDspot(const core::MarketData& market) const;
    core::Real calculateDzommaDvol(const core::MarketData& market) const;
    core::Real calculateDcolorDspot(const core::MarketData& market) const;
    core::Real calculateDvetaDvol(const core::MarketData& market) const;
    core::Real calculateDveraDrate(const core::MarketData& market) const;
    core::Real calculateDlambdaDspot(const core::MarketData& market) const;
    
    // Cross derivatives
    core::Real calculateDspeedDvol(const core::MarketData& market) const;
    core::Real calculateDspeedDtime(const core::MarketData& market) const;
    core::Real calculateDzommaDtime(const core::MarketData& market) const;
    core::Real calculateDcolorDvol(const core::MarketData& market) const;
    core::Real calculateDvannaDtime(const core::MarketData& market) const;
    core::Real calculateDcharmDvol(const core::MarketData& market) const;
    
    struct ThirdOrderRiskMetrics {
        core::Real convexityRisk = 0.0;        // Overall convexity exposure
        core::Real volOfVolRisk = 0.0;         // Volatility-of-volatility risk
        core::Real timeDecayAcceleration = 0.0; // Rate of time decay acceleration
        core::Real crossRisk = 0.0;            // Cross-derivative risk
        core::Real tailRisk = 0.0;             // Tail risk from third-order effects
        core::Vector riskContributions;        // Individual risk contributions
    };
    
    ThirdOrderRiskMetrics analyzePortfolioRisk(const std::vector<ThirdOrderGreeks>& positions,
                                              const core::Vector& quantities) const;
    
    struct ConvexityProfile {
        core::Vector spotLadder;
        core::Vector volLadder;
        core::Vector timeLadder;
        core::Matrix speedProfile;    // Speed across spot/vol grid
        core::Matrix zommaProfile;    // Zomma across spot/vol grid  
        core::Matrix colorProfile;    // Color across spot/time grid
        core::Matrix ultimaProfile;   // Ultima across vol/time grid
    };
    
    ConvexityProfile buildConvexityProfile(const core::MarketData& market,
                                          const core::Vector& spotRange,
                                          const core::Vector& volRange,
                                          const core::Vector& timeRange) const;
    
    struct HedgingRecommendation {
        core::Vector thirdOrderExposures;
        core::Vector hedgeInstruments;
        core::Vector hedgeQuantities;
        core::Real residualRisk;
        core::Real hedgingCost;
        core::Matrix hedgeEffectiveness;
    };
    
    HedgingRecommendation recommendThirdOrderHedge(const std::vector<ThirdOrderGreeks>& portfolio,
                                                  const std::vector<ThirdOrderGreeks>& hedgeInstruments,
                                                  const core::Vector& hedgeCosts) const;
    
private:
    // Analytical formulas for Black-Scholes third-order Greeks
    core::Real analyticalSpeed(const core::MarketData& market) const;
    core::Real analyticalZomma(const core::MarketData& market) const;
    core::Real analyticalColor(const core::MarketData& market) const;
    core::Real analyticalUltima(const core::MarketData& market) const;
    
    // Finite difference implementations
    core::Real numericalThirdDerivative(const std::function<core::Real(core::Real)>& func,
                                       core::Real x, core::Real h) const;
    
    core::Real numericalMixedDerivative(const std::function<core::Real(core::Real, core::Real)>& func,
                                       core::Real x, core::Real y, 
                                       core::Real hx, core::Real hy,
                                       bool isThirdOrder = true) const;
    
    // Helper functions for Black-Scholes parameters
    core::Real calculateD1(const core::MarketData& market) const;
    core::Real calculateD2(const core::MarketData& market) const;
    core::Real normalPDF(core::Real x) const;
    core::Real normalCDF(core::Real x) const;
    
    // Risk aggregation utilities
    core::Real calculatePortfolioConvexity(const std::vector<ThirdOrderGreeks>& positions,
                                          const core::Vector& quantities) const;
    
    core::Vector calculateRiskContributions(const std::vector<ThirdOrderGreeks>& positions,
                                           const core::Vector& quantities) const;
};

// Inline implementations for performance-critical functions
inline ThirdOrderCalculator::ThirdOrderCalculator(std::shared_ptr<core::MathEngine> engine, 
                                                  core::Real bumpSize)
    : engine_(std::move(engine)), defaultBumpSize_(bumpSize) {}

inline core::Real ThirdOrderCalculator::calculateD1(const core::MarketData& market) const {
    const core::Real sqrtT = std::sqrt(market.timeToExpiry);
    return (std::log(market.spot / market.strike) + 
            (market.riskFreeRate + 0.5 * market.volatility * market.volatility) * market.timeToExpiry) /
           (market.volatility * sqrtT);
}

inline core::Real ThirdOrderCalculator::calculateD2(const core::MarketData& market) const {
    return calculateD1(market) - market.volatility * std::sqrt(market.timeToExpiry);
}

inline core::Real ThirdOrderCalculator::normalPDF(core::Real x) const {
    return engine_->normalPDF(x);
}

inline core::Real ThirdOrderCalculator::normalCDF(core::Real x) const {
    return engine_->normalCDF(x);
}

}