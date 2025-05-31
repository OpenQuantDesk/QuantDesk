#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include "../curves/yield_curve.hpp"
#include <memory>

namespace math::pricing {

class EnhancedBlackScholesEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    std::shared_ptr<curves::YieldCurve> domesticCurve_;
    std::shared_ptr<curves::YieldCurve> foreignCurve_;
    std::shared_ptr<curves::YieldCurve> dividendCurve_;
    
public:
    explicit EnhancedBlackScholesEngine(std::shared_ptr<core::MathEngine> engine);
    
    void setDomesticCurve(std::shared_ptr<curves::YieldCurve> curve);
    void setForeignCurve(std::shared_ptr<curves::YieldCurve> curve);
    void setDividendCurve(std::shared_ptr<curves::YieldCurve> curve);
    
    struct CurveMarketData {
        core::Real spot = 100.0;
        core::Real strike = 100.0;
        core::Real timeToExpiry = 1.0;
        core::Real volatility = 0.2;
        bool isCall = true;
        bool useYieldCurve = true;
        bool useDividendCurve = false;
        bool useForeignCurve = false;
    };
    
    core::Greeks price(const CurveMarketData& market) const;
    core::Real priceOnly(const CurveMarketData& market) const;
    
    void priceBatch(const std::vector<CurveMarketData>& markets,
                   core::GreeksVector& results) const;
    
    core::Real impliedVolatility(core::Real marketPrice, const CurveMarketData& market) const;
    
    struct TermStructureGreeks {
        core::Greeks spotGreeks;
        core::Vector rateGreeks;
        core::Vector divGreeks;
        core::Matrix bucketGreeks;
        core::Real parallelDelta01 = 0.0;
        core::Real curveDuration = 0.0;
        core::Real curveConvexity = 0.0;
    };
    
    TermStructureGreeks calculateTermStructureGreeks(const CurveMarketData& market) const;
    
    core::Real americanPrice(const CurveMarketData& market, core::Integer timeSteps = 100) const;
    
private:
    core::Real getDiscountFactor(core::Real timeToExpiry) const;
    core::Real getDividendDiscountFactor(core::Real timeToExpiry) const;
    core::Real getForeignDiscountFactor(core::Real timeToExpiry) const;
    
    core::Real calculateForwardPrice(const CurveMarketData& market) const;
    core::Real calculateAdjustedRiskFreeRate(const CurveMarketData& market) const;
    
    core::Greeks calculateBucketSensitivities(const CurveMarketData& market,
                                             const core::Vector& bucketTenors,
                                             core::Real bumpSize = 0.0001) const;
};

}