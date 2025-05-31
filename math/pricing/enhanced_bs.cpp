#include "enhanced_bs.hpp"
#include "../core/utils.hpp"
#include <algorithm>

namespace math::pricing {

EnhancedBlackScholesEngine::EnhancedBlackScholesEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

void EnhancedBlackScholesEngine::setDomesticCurve(std::shared_ptr<curves::YieldCurve> curve) {
    domesticCurve_ = std::move(curve);
}

void EnhancedBlackScholesEngine::setForeignCurve(std::shared_ptr<curves::YieldCurve> curve) {
    foreignCurve_ = std::move(curve);
}

void EnhancedBlackScholesEngine::setDividendCurve(std::shared_ptr<curves::YieldCurve> curve) {
    dividendCurve_ = std::move(curve);
}

core::Greeks EnhancedBlackScholesEngine::price(const CurveMarketData& market) const {
    if (!market.useYieldCurve || !domesticCurve_) {
        core::MarketData standardMarket;
        standardMarket.spot = market.spot;
        standardMarket.strike = market.strike;
        standardMarket.timeToExpiry = market.timeToExpiry;
        standardMarket.riskFreeRate = 0.05;
        standardMarket.volatility = market.volatility;
        standardMarket.isCall = market.isCall;
        
        return engine_->blackScholes(standardMarket);
    }
    
    const core::Real discountFactor = getDiscountFactor(market.timeToExpiry);
    const core::Real dividendDiscount = getDividendDiscountFactor(market.timeToExpiry);
    const core::Real forwardPrice = calculateForwardPrice(market);
    
    const core::Real adjustedRate = calculateAdjustedRiskFreeRate(market);
    
    core::MarketData adjustedMarket;
    adjustedMarket.spot = forwardPrice * dividendDiscount;
    adjustedMarket.strike = market.strike;
    adjustedMarket.timeToExpiry = market.timeToExpiry;
    adjustedMarket.riskFreeRate = adjustedRate;
    adjustedMarket.volatility = market.volatility;
    adjustedMarket.isCall = market.isCall;
    
    core::Greeks greeks = engine_->blackScholes(adjustedMarket);
    
    greeks.price *= discountFactor;
    greeks.delta *= dividendDiscount;
    greeks.gamma *= dividendDiscount * dividendDiscount;
    greeks.vega *= discountFactor;
    greeks.rho *= market.timeToExpiry * discountFactor;
    
    return greeks;
}

core::Real EnhancedBlackScholesEngine::priceOnly(const CurveMarketData& market) const {
    return price(market).price;
}

void EnhancedBlackScholesEngine::priceBatch(const std::vector<CurveMarketData>& markets,
                                           core::GreeksVector& results) const {
    results.resize(markets.size());
    
    if (markets.size() < engine_->getThreadPool().getNumThreads() * 4) {
        for (core::Size i = 0; i < markets.size(); ++i) {
            results[i] = price(markets[i]);
        }
        return;
    }
    
    engine_->getThreadPool().parallelFor(markets.begin(), markets.end(),
        [this, &markets, &results](const CurveMarketData& market) {
            core::Size index = &market - &markets[0];
            results[index] = price(market);
        });
}

EnhancedBlackScholesEngine::TermStructureGreeks 
EnhancedBlackScholesEngine::calculateTermStructureGreeks(const CurveMarketData& market) const {
    TermStructureGreeks tsGreeks;
    
    tsGreeks.spotGreeks = price(market);
    
    if (!domesticCurve_) return tsGreeks;
    
    const core::Vector bucketTenors = {0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0};
    tsGreeks.bucketGreeks = calculateBucketSensitivities(market, bucketTenors).price;
    
    const core::Real basePrice = price(market).price;
    
    auto shiftedCurve = domesticCurve_->createShiftedCurve(0.0001);
    auto originalCurve = domesticCurve_;
    
    domesticCurve_ = std::make_shared<curves::YieldCurve>(shiftedCurve);
    const core::Real shiftedPrice = price(market).price;
    domesticCurve_ = originalCurve;
    
    tsGreeks.parallelDelta01 = (shiftedPrice - basePrice) / 0.0001;
    
    return tsGreeks;
}

core::Real EnhancedBlackScholesEngine::getDiscountFactor(core::Real timeToExpiry) const {
    if (!domesticCurve_) return std::exp(-0.05 * timeToExpiry);
    return domesticCurve_->getDiscountFactor(timeToExpiry);
}

core::Real EnhancedBlackScholesEngine::getDividendDiscountFactor(core::Real timeToExpiry) const {
    if (!dividendCurve_) return 1.0;
    return dividendCurve_->getDiscountFactor(timeToExpiry);
}

core::Real EnhancedBlackScholesEngine::getForeignDiscountFactor(core::Real timeToExpiry) const {
    if (!foreignCurve_) return 1.0;
    return foreignCurve_->getDiscountFactor(timeToExpiry);
}

core::Real EnhancedBlackScholesEngine::calculateForwardPrice(const CurveMarketData& market) const {
    const core::Real domesticDiscount = getDiscountFactor(market.timeToExpiry);
    const core::Real dividendDiscount = getDividendDiscountFactor(market.timeToExpiry);
    const core::Real foreignDiscount = getForeignDiscountFactor(market.timeToExpiry);
    
    if (market.useForeignCurve) {
        return market.spot * foreignDiscount / domesticDiscount;
    }
    
    return market.spot / (dividendDiscount * domesticDiscount);
}

core::Real EnhancedBlackScholesEngine::calculateAdjustedRiskFreeRate(const CurveMarketData& market) const {
    if (!domesticCurve_) return 0.05;
    
    const core::Real discountFactor = getDiscountFactor(market.timeToExpiry);
    return -std::log(discountFactor) / market.timeToExpiry;
}

core::Greeks EnhancedBlackScholesEngine::calculateBucketSensitivities(const CurveMarketData& market,
                                                                     const core::Vector& bucketTenors,
                                                                     core::Real bumpSize) const {
    core::Greeks bucketGreeks;
    
    if (!domesticCurve_) return bucketGreeks;
    
    const core::Real basePrice = price(market).price;
    
    for (core::Real tenor : bucketTenors) {
        auto tempCurve = *domesticCurve_;
        tempCurve.addRate(tenor, tempCurve.getRate(tenor) + bumpSize);
        
        auto originalCurve = domesticCurve_;
        domesticCurve_ = std::make_shared<curves::YieldCurve>(tempCurve);
        
        const core::Real bumpedPrice = price(market).price;
        const core::Real bucketSensitivity = (bumpedPrice - basePrice) / bumpSize;
        
        domesticCurve_ = originalCurve;
        
        bucketGreeks.price += bucketSensitivity;
    }
    
    return bucketGreeks;
}

}