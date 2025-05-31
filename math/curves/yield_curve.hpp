#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include "../core/utils.hpp"
#include <map>
#include <memory>

namespace math::curves {

struct RateTenor {
    core::Real maturity;
    core::Real rate;
    std::string tenor;
};

class YieldCurve {
private:
    std::map<core::Real, core::Real> rates_;
    std::string currency_;
    core::Integer valuationDate_;
    mutable std::map<core::Real, core::Real> discountCache_;
    mutable core::ReadWriteLock cacheLock_;
    
public:
    enum class InterpolationType {
        Linear,
        LogLinear,
        CubicSpline,
        Nelson,
        Svensson
    };
    
    explicit YieldCurve(const std::string& currency = "USD");
    
    void addRate(core::Real maturity, core::Real rate);
    void addRates(const std::vector<RateTenor>& tenors);
    
    core::Real getRate(core::Real maturity, InterpolationType type = InterpolationType::CubicSpline) const;
    core::Real getDiscountFactor(core::Real maturity) const;
    core::Real getForwardRate(core::Real startTime, core::Real endTime) const;
    
    core::Real getInstantaneousForward(core::Real maturity) const;
    core::Real getZeroRate(core::Real maturity) const;
    
    void calibrateNelsonSiegel(const std::vector<RateTenor>& marketRates);
    void calibrateSvensson(const std::vector<RateTenor>& marketRates);
    
    struct NelsonSiegelParams {
        core::Real beta0 = 0.05;
        core::Real beta1 = -0.02;
        core::Real beta2 = 0.01;
        core::Real lambda = 1.0;
    };
    
    struct SvenssonParams : NelsonSiegelParams {
        core::Real beta3 = 0.0;
        core::Real lambda2 = 2.0;
    };
    
    void setNelsonSiegelParams(const NelsonSiegelParams& params);
    void setSvenssonParams(const SvenssonParams& params);
    
    core::Vector getTenors() const;
    core::Vector getRates() const;
    
    void shift(core::Real parallelShift);
    void twist(core::Real shortEndShift, core::Real longEndShift);
    void butterfly(core::Real shortShift, core::Real midShift, core::Real longShift);
    
    YieldCurve createShiftedCurve(core::Real parallelShift) const;
    
    core::Real duration(core::Real cashFlow, core::Real maturity) const;
    core::Real convexity(core::Real cashFlow, core::Real maturity) const;
    
private:
    core::Real linearInterpolation(core::Real maturity) const;
    core::Real logLinearInterpolation(core::Real maturity) const;
    core::Real cubicSplineInterpolation(core::Real maturity) const;
    core::Real nelsonSiegelRate(core::Real maturity) const;
    core::Real svenssonRate(core::Real maturity) const;
    
    void clearCache() const;
    
    NelsonSiegelParams nsParams_;
    SvenssonParams svParams_;
    InterpolationType defaultInterpolation_ = InterpolationType::CubicSpline;
};

class YieldCurveBuilder {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit YieldCurveBuilder(std::shared_ptr<core::MathEngine> engine);
    
    struct InstrumentQuote {
        std::string type;
        core::Real maturity;
        core::Real quote;
        core::Real spread = 0.0;
    };
    
    YieldCurve buildFromDepositsSwaps(const std::vector<InstrumentQuote>& quotes,
                                     const std::string& currency = "USD") const;
    
    YieldCurve buildFromBonds(const std::vector<InstrumentQuote>& bondQuotes,
                             const std::string& currency = "USD") const;
    
    YieldCurve bootstrapCurve(const std::vector<InstrumentQuote>& instruments) const;
    
private:
    core::Real calculateSwapRate(const YieldCurve& curve, core::Real maturity, 
                                core::Real frequency = 2.0) const;
    core::Real solveForRate(const InstrumentQuote& instrument, 
                           const YieldCurve& partialCurve) const;
};

}