#pragma once

#include "core/engine.hpp"
#include "core/types.hpp"
#include "optimization/hardware.hpp"
#include "greeks/calculator.hpp"
#include "pricing/analytical.hpp"
#include "pricing/numerical.hpp"
#include "volatility/surface.hpp"
#include "simulation/paths.hpp"
#include "analytics/probability.hpp"
#include "analytics/performance.hpp"
#include "analytics/risk.hpp"
#include <memory>

namespace math {

class MathEngine {
private:
    std::shared_ptr<core::MathEngine> coreEngine_;
    std::shared_ptr<HardwareCapabilities> hardware_;
    std::shared_ptr<greeks::GreeksCalculator> greeksCalculator_;
    std::shared_ptr<pricing::BlackScholesEngine> bsEngine_;
    std::shared_ptr<volatility::VolatilitySurface> volSurface_;
    std::shared_ptr<simulation::PathGenerator> pathGen_;
    std::shared_ptr<analytics::PerformanceAnalyzer> perfAnalyzer_;
    std::shared_ptr<analytics::RiskAnalyzer> riskAnalyzer_;
    std::shared_ptr<ProbabilityEngine> probEngine_;

public:
    MathEngine();
    ~MathEngine() = default;
    
    MathEngine(const MathEngine&) = delete;
    MathEngine& operator=(const MathEngine&) = delete;
    
    std::shared_ptr<core::MathEngine> getCoreEngine() const { return coreEngine_; }
    std::shared_ptr<greeks::GreeksCalculator> getGreeksCalculator() const { return greeksCalculator_; }
    std::shared_ptr<pricing::BlackScholesEngine> getPricingEngine() const { return bsEngine_; }
    std::shared_ptr<volatility::VolatilitySurface> getVolatilitySurface() const { return volSurface_; }
    std::shared_ptr<simulation::PathGenerator> getPathGenerator() const { return pathGen_; }
    std::shared_ptr<analytics::PerformanceAnalyzer> getPerformanceAnalyzer() const { return perfAnalyzer_; }
    std::shared_ptr<analytics::RiskAnalyzer> getRiskAnalyzer() const { return riskAnalyzer_; }
    std::shared_ptr<ProbabilityEngine> getProbabilityEngine() const { return probEngine_; }
    
    const HardwareCapabilities& getHardwareInfo() const { return *hardware_; }
    
    core::Greeks blackScholes(const core::MarketData& market) const;
    void blackScholesBatch(const core::MarketDataVector& markets, 
                          core::GreeksVector& results) const;
    
    core::Real impliedVolatility(core::Real marketPrice, 
                                const core::MarketData& market) const;
    
    core::Real normalCDF(core::Real x) const;
    core::Real normalPDF(core::Real x) const;
    core::Real normalInverse(core::Real p) const;
    
    void normalCDF_vectorized(const core::Real* input, core::Real* output, 
                             core::Size count) const;
    
    std::mt19937& getRNG() const;
    
    template<typename T>
    T* allocate(core::Size count) {
        return coreEngine_->allocate<T>(count);
    }
    
    template<typename T>
    void deallocate(T* ptr, core::Size count) {
        coreEngine_->deallocate<T>(ptr, count);
    }
};

inline MathEngine::MathEngine() 
    : coreEngine_(std::make_shared<core::MathEngine>()),
      hardware_(std::make_shared<HardwareCapabilities>(HardwareCapabilities::detect())),
      greeksCalculator_(std::make_shared<greeks::GreeksCalculator>(coreEngine_)),
      bsEngine_(std::make_shared<pricing::BlackScholesEngine>(coreEngine_)),
      volSurface_(std::make_shared<volatility::VolatilitySurface>(coreEngine_)),
      pathGen_(std::make_shared<simulation::PathGenerator>(coreEngine_)),
      perfAnalyzer_(std::make_shared<analytics::PerformanceAnalyzer>(coreEngine_)),
      riskAnalyzer_(std::make_shared<analytics::RiskAnalyzer>(coreEngine_)),
      probEngine_(std::make_shared<ProbabilityEngine>(coreEngine_)) {}

inline core::Greeks MathEngine::blackScholes(const core::MarketData& market) const {
    return coreEngine_->blackScholes(market);
}

inline void MathEngine::blackScholesBatch(const core::MarketDataVector& markets, 
                                         core::GreeksVector& results) const {
    coreEngine_->blackScholesBatch(markets, results);
}

inline core::Real MathEngine::impliedVolatility(core::Real marketPrice, 
                                               const core::MarketData& market) const {
    return coreEngine_->impliedVolatility(marketPrice, market);
}

inline core::Real MathEngine::normalCDF(core::Real x) const {
    return coreEngine_->normalCDF(x);
}

inline core::Real MathEngine::normalPDF(core::Real x) const {
    return coreEngine_->normalPDF(x);
}

inline core::Real MathEngine::normalInverse(core::Real p) const {
    return coreEngine_->normalInverse(p);
}

inline void MathEngine::normalCDF_vectorized(const core::Real* input, core::Real* output, 
                                            core::Size count) const {
    coreEngine_->normalCDF_vectorized(input, output, count);
}

inline std::mt19937& MathEngine::getRNG() const {
    return coreEngine_->getRNG();
}

}