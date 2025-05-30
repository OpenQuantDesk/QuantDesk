#pragma once

#include "types.hpp"
#include "memory.hpp"
#include "threading.hpp"
#include <immintrin.h>
#include <cmath>
#include <algorithm>

namespace math::core {

class MathEngine {
private:
    HardwareInfo hardware_;
    std::unique_ptr<MemoryPool> memoryPool_;
    std::unique_ptr<GlobalThreadPool> threadPool_;
    ThreadLocalStorage<std::mt19937> rngStorage_;
    
public:
    MathEngine();
    ~MathEngine() = default;
    
    MathEngine(const MathEngine&) = delete;
    MathEngine& operator=(const MathEngine&) = delete;
    
    Real normalCDF(Real x) const;
    Real normalPDF(Real x) const;
    Real normalInverse(Real p) const;
    
    void normalCDF_vectorized(const Real* input, Real* output, Size count) const;
    void normalPDF_vectorized(const Real* input, Real* output, Size count) const;
    
    Greeks blackScholes(const MarketData& market) const;
    
    void blackScholesBatch(const MarketDataVector& markets, GreeksVector& results) const;
    
    Real impliedVolatility(Real marketPrice, const MarketData& market, 
                          Real tolerance = TOLERANCE, Integer maxIter = MAX_ITERATIONS) const;
    
    void impliedVolatilityBatch(const Vector& marketPrices, const MarketDataVector& markets,
                               Vector& impliedVols) const;
    
    Real realizedVolatility(const Vector& prices, Integer windowDays = 252) const;
    
    VolatilityMetrics calculateVolatilityMetrics(const Vector& prices, 
                                               const Vector& impliedVols) const;
    
    Greeks aggregateGreeks(const GreeksVector& positions, const Vector& quantities) const;
    
    template<typename T>
    T* allocate(Size count) {
        return static_cast<T*>(memoryPool_->allocate(count * sizeof(T)));
    }
    
    template<typename T>
    void deallocate(T* ptr, Size count) {
        memoryPool_->deallocate(ptr, count * sizeof(T));
    }
    
    GlobalThreadPool& getThreadPool() { return *threadPool_; }
    const HardwareInfo& getHardwareInfo() const { return hardware_; }
    
    std::mt19937& getRNG() { return rngStorage_.get(); }
    
private:
    void detectHardware();
    void initializeMemoryPool();
    void initializeThreadPool();
    
    Real normalCDF_scalar(Real x) const;
    Real normalCDF_avx512(Real x) const;
    Real normalCDF_avx2(Real x) const;
    
    void normalCDF_scalar_batch(const Real* input, Real* output, Size count) const;
    void normalCDF_avx512_batch(const Real* input, Real* output, Size count) const;
    void normalCDF_avx2_batch(const Real* input, Real* output, Size count) const;
    
    Real blackScholesPrice(Real S, Real K, Real T, Real r, Real v, bool isCall) const;
    Greeks blackScholesGreeks(Real S, Real K, Real T, Real r, Real v, bool isCall) const;
    
    Real percentileRank(Real value, const Vector& dataset) const;
};

inline MathEngine::MathEngine() 
    : rngStorage_([]() { return std::make_unique<std::mt19937>(std::random_device{}()); }) {
    detectHardware();
    initializeMemoryPool();
    initializeThreadPool();
}

inline void MathEngine::detectHardware() {
    hardware_.numCores = std::thread::hardware_concurrency();
    hardware_.numThreads = hardware_.numCores;
    
#ifdef __x86_64__
    unsigned int eax, ebx, ecx, edx;
    
    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        hardware_.hasSSE42 = (ecx & bit_SSE4_2) != 0;
        hardware_.hasAVX = (ecx & bit_AVX) != 0;
        hardware_.hasFMA = (ecx & bit_FMA) != 0;
    }
    
    if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
        hardware_.hasAVX2 = (ebx & bit_AVX2) != 0;
        hardware_.hasAVX512 = (ebx & bit_AVX512F) != 0;
    }
#endif

#ifdef __ARM_NEON
    hardware_.hasNEON = true;
#endif
}

inline void MathEngine::initializeMemoryPool() {
    memoryPool_ = std::make_unique<MemoryPool>();
}

inline void MathEngine::initializeThreadPool() {
    threadPool_ = std::make_unique<GlobalThreadPool>(hardware_.numThreads);
}

inline Real MathEngine::normalCDF(Real x) const {
    if (hardware_.hasAVX512) return normalCDF_avx512(x);
    if (hardware_.hasAVX2) return normalCDF_avx2(x);
    return normalCDF_scalar(x);
}

inline Real MathEngine::normalCDF_scalar(Real x) const {
    static constexpr Real a1 = 0.254829592;
    static constexpr Real a2 = -0.284496736;
    static constexpr Real a3 = 1.421413741;
    static constexpr Real a4 = -1.453152027;
    static constexpr Real a5 = 1.061405429;
    static constexpr Real p = 0.3275911;
    
    const Real sign = (x >= 0.0) ? 1.0 : -1.0;
    x = std::abs(x) / SQRT_2;
    
    const Real t = 1.0 / (1.0 + p * x);
    const Real y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);
    
    return 0.5 * (1.0 + sign * y);
}

inline Real MathEngine::normalPDF(Real x) const {
    return INV_SQRT_2PI * std::exp(-0.5 * x * x);
}

inline Real MathEngine::normalInverse(Real p) const {
    static constexpr Real a0 = -3.969683028665376e+01;
    static constexpr Real a1 =  2.209460984245205e+02;
    static constexpr Real a2 = -2.759285104469687e+02;
    static constexpr Real a3 =  1.383577518672690e+02;
    static constexpr Real a4 = -3.066479806614716e+01;
    static constexpr Real a5 =  2.506628277459239e+00;
    
    static constexpr Real b1 = -5.447609879822406e+01;
    static constexpr Real b2 =  1.615858368580409e+02;
    static constexpr Real b3 = -1.556989798598866e+02;
    static constexpr Real b4 =  6.680131188771972e+01;
    static constexpr Real b5 = -1.328068155288572e+01;
    
    if (p <= 0.0 || p >= 1.0) return 0.0;
    
    const Real q = p - 0.5;
    
    if (std::abs(q) <= 0.425) {
        const Real r = 0.180625 - q * q;
        return q * (((((a5 * r + a4) * r + a3) * r + a2) * r + a1) * r + a0) /
                   (((((b5 * r + b4) * r + b3) * r + b2) * r + b1) * r + 1.0);
    }
    
    const Real r = (q < 0.0) ? p : 1.0 - p;
    const Real t = std::sqrt(-2.0 * std::log(r));
    const Real sign = (q < 0.0) ? -1.0 : 1.0;
    
    return sign * (t - ((2.515517 + 0.802853 * t + 0.010328 * t * t) /
                       (1.0 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t)));
}

inline void MathEngine::normalCDF_vectorized(const Real* input, Real* output, Size count) const {
    if (hardware_.hasAVX512 && count >= SIMD_WIDTH_AVX512) {
        normalCDF_avx512_batch(input, output, count);
    } else if (hardware_.hasAVX2 && count >= SIMD_WIDTH_AVX2) {
        normalCDF_avx2_batch(input, output, count);
    } else {
        normalCDF_scalar_batch(input, output, count);
    }
}

inline void MathEngine::normalCDF_scalar_batch(const Real* input, Real* output, Size count) const {
    for (Size i = 0; i < count; ++i) {
        output[i] = normalCDF_scalar(input[i]);
    }
}

#ifdef __AVX512F__
inline Real MathEngine::normalCDF_avx512(Real x) const {
    return normalCDF_scalar(x);
}

inline void MathEngine::normalCDF_avx512_batch(const Real* input, Real* output, Size count) const {
    const Size simdCount = (count / SIMD_WIDTH_AVX512) * SIMD_WIDTH_AVX512;
    
    const __m512d a1_vec = _mm512_set1_pd(0.254829592);
    const __m512d a2_vec = _mm512_set1_pd(-0.284496736);
    const __m512d a3_vec = _mm512_set1_pd(1.421413741);
    const __m512d a4_vec = _mm512_set1_pd(-1.453152027);
    const __m512d a5_vec = _mm512_set1_pd(1.061405429);
    const __m512d p_vec = _mm512_set1_pd(0.3275911);
    const __m512d sqrt2_vec = _mm512_set1_pd(SQRT_2);
    const __m512d half_vec = _mm512_set1_pd(0.5);
    const __m512d one_vec = _mm512_set1_pd(1.0);
    
    for (Size i = 0; i < simdCount; i += SIMD_WIDTH_AVX512) {
        __m512d x_vec = _mm512_load_pd(&input[i]);
        __m512d sign_vec = _mm512_cmp_pd_mask(x_vec, _mm512_setzero_pd(), _CMP_GE_OQ) ? one_vec : _mm512_set1_pd(-1.0);
        x_vec = _mm512_div_pd(_mm512_abs_pd(x_vec), sqrt2_vec);
        
        __m512d t_vec = _mm512_div_pd(one_vec, _mm512_fmadd_pd(p_vec, x_vec, one_vec));
        __m512d x_sq = _mm512_mul_pd(x_vec, x_vec);
        __m512d exp_neg_x_sq = _mm512_exp_pd(_mm512_sub_pd(_mm512_setzero_pd(), x_sq));
        
        __m512d poly = _mm512_fmadd_pd(a5_vec, t_vec, a4_vec);
        poly = _mm512_fmadd_pd(poly, t_vec, a3_vec);
        poly = _mm512_fmadd_pd(poly, t_vec, a2_vec);
        poly = _mm512_fmadd_pd(poly, t_vec, a1_vec);
        poly = _mm512_mul_pd(_mm512_mul_pd(poly, t_vec), exp_neg_x_sq);
        
        __m512d y_vec = _mm512_sub_pd(one_vec, poly);
        __m512d result = _mm512_mul_pd(half_vec, _mm512_fmadd_pd(sign_vec, y_vec, one_vec));
        
        _mm512_store_pd(&output[i], result);
    }
    
    for (Size i = simdCount; i < count; ++i) {
        output[i] = normalCDF_scalar(input[i]);
    }
}
#else
inline Real MathEngine::normalCDF_avx512(Real x) const {
    return normalCDF_scalar(x);
}

inline void MathEngine::normalCDF_avx512_batch(const Real* input, Real* output, Size count) const {
    normalCDF_scalar_batch(input, output, count);
}
#endif

#ifdef __AVX2__
inline Real MathEngine::normalCDF_avx2(Real x) const {
    return normalCDF_scalar(x);
}

inline void MathEngine::normalCDF_avx2_batch(const Real* input, Real* output, Size count) const {
    const Size simdCount = (count / SIMD_WIDTH_AVX2) * SIMD_WIDTH_AVX2;
    
    for (Size i = 0; i < simdCount; i += SIMD_WIDTH_AVX2) {
        for (Size j = 0; j < SIMD_WIDTH_AVX2; ++j) {
            output[i + j] = normalCDF_scalar(input[i + j]);
        }
    }
    
    for (Size i = simdCount; i < count; ++i) {
        output[i] = normalCDF_scalar(input[i]);
    }
}
#else
inline Real MathEngine::normalCDF_avx2(Real x) const {
    return normalCDF_scalar(x);
}

inline void MathEngine::normalCDF_avx2_batch(const Real* input, Real* output, Size count) const {
    normalCDF_scalar_batch(input, output, count);
}
#endif

inline Greeks MathEngine::blackScholes(const MarketData& market) const {
    return blackScholesGreeks(market.spot, market.strike, market.timeToExpiry,
                             market.riskFreeRate, market.volatility, market.isCall);
}

inline Greeks MathEngine::blackScholesGreeks(Real S, Real K, Real T, Real r, Real v, bool isCall) const {
    Greeks greeks;
    
    if (T <= 0.0 || v <= 0.0 || S <= 0.0 || K <= 0.0) {
        return greeks;
    }
    
    const Real sqrtT = std::sqrt(T);
    const Real d1 = (std::log(S / K) + (r + 0.5 * v * v) * T) / (v * sqrtT);
    const Real d2 = d1 - v * sqrtT;
    
    const Real nd1 = normalCDF(d1);
    const Real nd2 = normalCDF(d2);
    const Real nPrimeD1 = normalPDF(d1);
    const Real discountFactor = std::exp(-r * T);
    
    if (isCall) {
        greeks.price = S * nd1 - K * discountFactor * nd2;
        greeks.delta = nd1;
    } else {
        greeks.price = K * discountFactor * normalCDF(-d2) - S * normalCDF(-d1);
        greeks.delta = nd1 - 1.0;
    }
    
    greeks.gamma = nPrimeD1 / (S * v * sqrtT);
    greeks.theta = (-S * nPrimeD1 * v / (2.0 * sqrtT) - 
                   (isCall ? 1.0 : -1.0) * r * K * discountFactor * 
                   (isCall ? nd2 : normalCDF(-d2))) / 365.0;
    greeks.vega = S * nPrimeD1 * sqrtT / 100.0;
    greeks.rho = (isCall ? 1.0 : -1.0) * K * T * discountFactor * 
                 (isCall ? nd2 : normalCDF(-d2)) / 100.0;
    greeks.impliedVol = v;
    
    return greeks;
}

inline void MathEngine::blackScholesBatch(const MarketDataVector& markets, GreeksVector& results) const {
    results.resize(markets.size());
    
    if (markets.size() < threadPool_->getNumThreads() * 4) {
        for (Size i = 0; i < markets.size(); ++i) {
            results[i] = blackScholes(markets[i]);
        }
        return;
    }
    
    threadPool_->parallelFor(markets.begin(), markets.end(),
        [this, &markets, &results](const MarketData& market) {
            Size index = &market - &markets[0];
            results[index] = blackScholes(market);
        });
}

inline Real MathEngine::impliedVolatility(Real marketPrice, const MarketData& market,
                                         Real tolerance, Integer maxIter) const {
    Real vol = 0.25;
    
    for (Integer i = 0; i < maxIter; ++i) {
        MarketData tempMarket = market;
        tempMarket.volatility = vol;
        
        const Greeks greeks = blackScholes(tempMarket);
        const Real price = greeks.price;
        const Real vega = greeks.vega;
        
        if (std::abs(vega) < tolerance) break;
        
        const Real diff = price - marketPrice;
        if (std::abs(diff) < tolerance) break;
        
        vol -= diff / (vega * 100.0);
        vol = std::clamp(vol, Real(0.001), Real(5.0));
    }
    
    return vol;
}

}