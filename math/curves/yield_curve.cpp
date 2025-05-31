#include "yield_curve.hpp"
#include "../core/utils.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace math::curves {

YieldCurve::YieldCurve(const std::string& currency) 
    : currency_(currency), valuationDate_(0) {}

void YieldCurve::addRate(core::Real maturity, core::Real rate) {
    rates_[maturity] = rate;
    clearCache();
}

void YieldCurve::addRates(const std::vector<RateTenor>& tenors) {
    for (const auto& tenor : tenors) {
        rates_[tenor.maturity] = tenor.rate;
    }
    clearCache();
}

core::Real YieldCurve::getRate(core::Real maturity, InterpolationType type) const {
    if (rates_.empty()) return 0.05;
    
    auto it = rates_.find(maturity);
    if (it != rates_.end()) {
        return it->second;
    }
    
    switch (type) {
        case InterpolationType::Linear:
            return linearInterpolation(maturity);
        case InterpolationType::LogLinear:
            return logLinearInterpolation(maturity);
        case InterpolationType::CubicSpline:
            return cubicSplineInterpolation(maturity);
        case InterpolationType::Nelson:
            return nelsonSiegelRate(maturity);
        case InterpolationType::Svensson:
            return svenssonRate(maturity);
        default:
            return linearInterpolation(maturity);
    }
}

core::Real YieldCurve::getDiscountFactor(core::Real maturity) const {
    {
        core::SharedLockGuard<core::ReadWriteLock> lock(cacheLock_);
        auto it = discountCache_.find(maturity);
        if (it != discountCache_.end()) {
            return it->second;
        }
    }
    
    const core::Real rate = getRate(maturity);
    const core::Real discount = std::exp(-rate * maturity);
    
    {
        core::UniqueLockGuard<core::ReadWriteLock> lock(cacheLock_);
        discountCache_[maturity] = discount;
    }
    
    return discount;
}

core::Real YieldCurve::getForwardRate(core::Real startTime, core::Real endTime) const {
    if (startTime >= endTime) return getRate(startTime);
    
    const core::Real df1 = getDiscountFactor(startTime);
    const core::Real df2 = getDiscountFactor(endTime);
    
    return std::log(df1 / df2) / (endTime - startTime);
}

core::Real YieldCurve::getInstantaneousForward(core::Real maturity) const {
    const core::Real h = 0.001;
    const core::Real rate1 = getRate(maturity);
    const core::Real rate2 = getRate(maturity + h);
    
    return rate1 + maturity * (rate2 - rate1) / h;
}

core::Real YieldCurve::linearInterpolation(core::Real maturity) const {
    if (rates_.size() < 2) return rates_.empty() ? 0.05 : rates_.begin()->second;
    
    auto upper = rates_.upper_bound(maturity);
    if (upper == rates_.begin()) return upper->second;
    if (upper == rates_.end()) return rates_.rbegin()->second;
    
    auto lower = std::prev(upper);
    
    const core::Real t1 = lower->first;
    const core::Real t2 = upper->first;
    const core::Real r1 = lower->second;
    const core::Real r2 = upper->second;
    
    return r1 + (r2 - r1) * (maturity - t1) / (t2 - t1);
}

core::Real YieldCurve::logLinearInterpolation(core::Real maturity) const {
    if (rates_.size() < 2) return rates_.empty() ? 0.05 : rates_.begin()->second;
    
    auto upper = rates_.upper_bound(maturity);
    if (upper == rates_.begin()) return upper->second;
    if (upper == rates_.end()) return rates_.rbegin()->second;
    
    auto lower = std::prev(upper);
    
    const core::Real t1 = lower->first;
    const core::Real t2 = upper->first;
    const core::Real df1 = std::exp(-lower->second * t1);
    const core::Real df2 = std::exp(-upper->second * t2);
    
    const core::Real logDf = std::log(df1) + (std::log(df2) - std::log(df1)) * (maturity - t1) / (t2 - t1);
    
    return -logDf / maturity;
}

core::Real YieldCurve::cubicSplineInterpolation(core::Real maturity) const {
    if (rates_.size() < 4) return linearInterpolation(maturity);
    
    core::Vector tenors;
    core::Vector rates;
    
    for (const auto& [tenor, rate] : rates_) {
        tenors.push_back(tenor);
        rates.push_back(rate);
    }
    
    const auto interpolated = core::Interpolation::cubicSpline(tenors, rates, {maturity});
    return interpolated.empty() ? 0.05 : interpolated[0];
}

core::Real YieldCurve::nelsonSiegelRate(core::Real maturity) const {
    const core::Real tau = maturity;
    const core::Real lambda = nsParams_.lambda;
    
    const core::Real term1 = nsParams_.beta0;
    const core::Real term2 = nsParams_.beta1 * (1.0 - std::exp(-lambda * tau)) / (lambda * tau);
    const core::Real term3 = nsParams_.beta2 * ((1.0 - std::exp(-lambda * tau)) / (lambda * tau) - std::exp(-lambda * tau));
    
    return std::max(0.0001, term1 + term2 + term3);
}

core::Real YieldCurve::svenssonRate(core::Real maturity) const {
    const core::Real nsRate = nelsonSiegelRate(maturity);
    const core::Real tau = maturity;
    const core::Real lambda2 = svParams_.lambda2;
    
    const core::Real additionalTerm = svParams_.beta3 * ((1.0 - std::exp(-lambda2 * tau)) / (lambda2 * tau) - std::exp(-lambda2 * tau));
    
    return std::max(0.0001, nsRate + additionalTerm);
}

void YieldCurve::clearCache() const {
    core::UniqueLockGuard<core::ReadWriteLock> lock(cacheLock_);
    discountCache_.clear();
}

void YieldCurve::shift(core::Real parallelShift) {
    for (auto& [maturity, rate] : rates_) {
        rate += parallelShift;
    }
    clearCache();
}

YieldCurve YieldCurve::createShiftedCurve(core::Real parallelShift) const {
    YieldCurve shifted(currency_);
    for (const auto& [maturity, rate] : rates_) {
        shifted.addRate(maturity, rate + parallelShift);
    }
    return shifted;
}

YieldCurveBuilder::YieldCurveBuilder(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

YieldCurve YieldCurveBuilder::bootstrapCurve(const std::vector<InstrumentQuote>& instruments) const {
    YieldCurve curve;
    
    std::vector<InstrumentQuote> sortedInstruments = instruments;
    std::sort(sortedInstruments.begin(), sortedInstruments.end(),
              [](const InstrumentQuote& a, const InstrumentQuote& b) {
                  return a.maturity < b.maturity;
              });
    
    for (const auto& instrument : sortedInstruments) {
        if (instrument.type == "deposit" || instrument.type == "libor") {
            const core::Real rate = instrument.quote;
            curve.addRate(instrument.maturity, rate);
        } else if (instrument.type == "swap") {
            const core::Real solvedRate = solveForRate(instrument, curve);
            curve.addRate(instrument.maturity, solvedRate);
        }
    }
    
    return curve;
}

core::Real YieldCurveBuilder::solveForRate(const InstrumentQuote& instrument, 
                                          const YieldCurve& partialCurve) const {
    auto objectiveFunc = [&](core::Real rate) -> core::Real {
        YieldCurve tempCurve = partialCurve;
        tempCurve.addRate(instrument.maturity, rate);
        
        const core::Real modelRate = calculateSwapRate(tempCurve, instrument.maturity);
        return modelRate - instrument.quote;
    };
    
    return core::NumericalMethods::brent(objectiveFunc, 0.001, 0.20);
}

core::Real YieldCurveBuilder::calculateSwapRate(const YieldCurve& curve, core::Real maturity, 
                                               core::Real frequency) const {
    const core::Integer numPayments = static_cast<core::Integer>(maturity * frequency);
    const core::Real paymentInterval = 1.0 / frequency;
    
    core::Real pvFixedLeg = 0.0;
    core::Real pvFloatingLeg = 1.0 - curve.getDiscountFactor(maturity);
    
    for (core::Integer i = 1; i <= numPayments; ++i) {
        const core::Real paymentTime = i * paymentInterval;
        pvFixedLeg += paymentInterval * curve.getDiscountFactor(paymentTime);
    }
    
    return pvFloatingLeg / pvFixedLeg;
}

}