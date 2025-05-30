#pragma once

#include "common/types.hpp"
#include "math/optimization/hardware.hpp"
#include "math/metrics/volatility.hpp"
#include <vector>
#include <memory>
#include <thread>
#include <cmath>

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace math {

struct ExtendedGreeks {
    double delta = 0.0;
    double gamma = 0.0;
    double theta = 0.0;
    double vega = 0.0;
    double rho = 0.0;
    
    double volga = 0.0;
    double vanna = 0.0;
    double charm = 0.0;
    double speed = 0.0;
    double zomma = 0.0;
    double color = 0.0;
    
    double dollarDelta = 0.0;
    double dollarGamma = 0.0;
    double pinRisk = 0.0;
    double price = 0.0;
};

class MathEngine {
private:
    HardwareCapabilities capabilities_;
    
public:
    MathEngine();
    ~MathEngine() = default;
    
    MathEngine(const MathEngine&) = delete;
    MathEngine& operator=(const MathEngine&) = delete;
    
    double normalCDF(double x) const;
    double normalPDF(double x) const;
    
    void normalCDF_vector(const std::vector<double>& input, 
                         std::vector<double>& output) const;
    
    common::Greeks blackScholes(double spot, double strike, double timeToExpiry,
                               double riskFreeRate, double volatility, bool isCall) const;
    
    void blackScholesBatch(const std::vector<double>& spots,
                          const std::vector<double>& strikes,
                          const std::vector<double>& timeToExpiries,
                          const std::vector<double>& riskFreeRates,
                          const std::vector<double>& volatilities,
                          const std::vector<bool>& isCall,
                          std::vector<common::Greeks>& results) const;
    
    double impliedVolatility(double marketPrice, double spot, double strike,
                           double timeToExpiry, double riskFreeRate, bool isCall) const;
    
    double realizedVolatility(const std::vector<double>& prices, int windowDays = 20) const;
    
    VolatilityMetrics calculateVolatilityMetrics(const std::vector<double>& prices,
                                               const std::vector<double>& impliedVols) const;
    
    double percentileRank(double value, const std::vector<double>& dataset) const;
    
    common::Greeks aggregateGreeks(const std::vector<common::Greeks>& positions,
                                 const std::vector<double>& quantities) const;
    
    HardwareCapabilities getCapabilities() const { return capabilities_; }
    
private:
    void initializeOptimizations();
    
    double normalCDF_scalar(double x) const;
    void normalCDF_vector_scalar(const std::vector<double>& input, 
                               std::vector<double>& output) const;
    
#ifdef __SSE2__
    double normalCDF_sse2(double x) const;
    void normalCDF_vector_sse2(const std::vector<double>& input, 
                             std::vector<double>& output) const;
#endif

#ifdef __AVX2__
    double normalCDF_avx2(double x) const;
    void normalCDF_vector_avx2(const std::vector<double>& input, 
                             std::vector<double>& output) const;
#endif
};

inline double fastNormalCDF(double x) {
    static constexpr double a1 =  0.254829592;
    static constexpr double a2 = -0.284496736;
    static constexpr double a3 =  1.421413741;
    static constexpr double a4 = -1.453152027;
    static constexpr double a5 =  1.061405429;
    static constexpr double p  =  0.3275911;
    static constexpr double sqrt2 = 1.4142135623730950488;
    
    const double sign = (x >= 0) ? 1.0 : -1.0;
    x = std::abs(x) / sqrt2;
    
    const double t = 1.0 / (1.0 + p * x);
    const double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);
    
    return 0.5 * (1.0 + sign * y);
}

inline double fastNormalPDF(double x) {
    static constexpr double invSqrt2Pi = 0.3989422804014326779;
    return invSqrt2Pi * std::exp(-0.5 * x * x);
}

}