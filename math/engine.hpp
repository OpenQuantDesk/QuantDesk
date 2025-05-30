#pragma once

#include "common/greeks.hpp"
#include "math/hardware_capabilities.hpp"
#include <vector>
#include <memory>
#include <thread>

namespace math {

class MathEngine {
private:
    HardwareCapabilities capabilities_;
    std::unique_ptr<class HardwareOptimizer> optimizer_;
    std::unique_ptr<class MemoryPool> memoryPool_;
    
public:
    MathEngine();
    ~MathEngine();
    
    MathEngine(const MathEngine&) = delete;
    MathEngine& operator=(const MathEngine&) = delete;
    
    double normalCDF(double x) const;
    double normalPDF(double x) const;
    void normalCDF_vector(const std::vector<double>& input, std::vector<double>& output) const;
    
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
    
    double realizedVolatility(const std::vector<double>& prices, int windowDays = 252) const;
    
    common::Greeks aggregateGreeks(const std::vector<common::Greeks>& positions,
                                 const std::vector<double>& quantities) const;
    
    const HardwareCapabilities& getCapabilities() const { return capabilities_; }
    
private:
    double normalCDF_scalar(double x) const;
    void normalCDF_vector_scalar(const std::vector<double>& input, std::vector<double>& output) const;
    
#ifdef __SSE2__
    double normalCDF_sse2(double x) const;
    void normalCDF_vector_sse2(const std::vector<double>& input, std::vector<double>& output) const;
#endif

#ifdef __AVX2__
    double normalCDF_avx2(double x) const;
    void normalCDF_vector_avx2(const std::vector<double>& input, std::vector<double>& output) const;
#endif

    void initializeOptimizations();
    double percentileRank(double value, const std::vector<double>& dataset) const;
};

}