#include "math/engine.hpp"
#include <algorithm>
#include <immintrin.h>
#include <numeric>

namespace math {

    MathEngine::MathEngine(): capabilities_(HardwareCapabilities::detect())
    {
        initializeOptimizations();
    }

    double MathEngine::normalCDF(double x) const
    {
        if(capabilities_.avx2) { return normalCDF_avx2(x); }
        if(capabilities_.sse2) { return normalCDF_sse2(x); }
        return normalCDF_scalar(x);
    }

    void MathEngine::normalCDF_vector(const std::vector<double>& input,
                                      std::vector<double>& output) const
    {
        output.resize(input.size());

        if(capabilities_.avx2 && input.size() >= 8) {
            normalCDF_vector_avx2(input, output);
            return;
        }

        if(capabilities_.sse2 && input.size() >= 4) {
            normalCDF_vector_sse2(input, output);
            return;
        }

        normalCDF_vector_scalar(input, output);
    }

    Greeks MathEngine::blackScholes(double spot, double strike,
                                    double timeToExpiry, double riskFreeRate,
                                    double volatility, bool isCall) const
    {
        Greeks greeks;

        if(timeToExpiry <= 0 || volatility <= 0 || spot <= 0 || strike <= 0) {
            return greeks;
        }

        const double sqrtT = std::sqrt(timeToExpiry);
        const double d1
            = (std::log(spot / strike)
               + (riskFreeRate + 0.5 * volatility * volatility) * timeToExpiry)
              / (volatility * sqrtT);
        const double d2 = d1 - volatility * sqrtT;

        const double nd1 = normalCDF(d1);
        const double nd2 = normalCDF(d2);
        const double nPrimeD1 = normalPDF(d1);
        const double discountFactor = std::exp(-riskFreeRate * timeToExpiry);

        if(isCall) {
            greeks.price = spot * nd1 - strike * discountFactor * nd2;
            greeks.delta = nd1;
        }
        else {
            greeks.price = strike * discountFactor * normalCDF(-d2)
                           - spot * normalCDF(-d1);
            greeks.delta = nd1 - 1.0;
        }

        greeks.gamma = nPrimeD1 / (spot * volatility * sqrtT);
        greeks.theta = (-spot * nPrimeD1 * volatility / (2 * sqrtT)
                        - riskFreeRate * strike * discountFactor * nd2)
                       / 365.0;
        greeks.vega = spot * nPrimeD1 * sqrtT / 100.0;
        greeks.rho = strike * timeToExpiry * discountFactor * nd2 / 100.0;
        greeks.impliedVol = volatility;

        return greeks;
    }

    void MathEngine::blackScholesBatch(
        const std::vector<double>& spots, const std::vector<double>& strikes,
        const std::vector<double>& timeToExpiries,
        const std::vector<double>& riskFreeRates,
        const std::vector<double>& volatilities,
        const std::vector<bool>& isCall, std::vector<Greeks>& results) const
    {
        const size_t n = spots.size();
        results.resize(n);

        const size_t numThreads = std::thread::hardware_concurrency();
        const size_t chunkSize = (n + numThreads - 1) / numThreads;

        std::vector<std::thread> threads;
        threads.reserve(numThreads);

        for(size_t i = 0; i < numThreads; ++i) {
            const size_t start = i * chunkSize;
            const size_t end = std::min(start + chunkSize, n);

            if(start >= n) break;

            threads.emplace_back([this, &spots, &strikes, &timeToExpiries,
                                  &riskFreeRates, &volatilities, &isCall,
                                  &results, start, end]() {
                for(size_t j = start; j < end; ++j) {
                    results[j] = blackScholes(
                        spots[j], strikes[j], timeToExpiries[j],
                        riskFreeRates[j], volatilities[j], isCall[j]);
                }
            });
        }

        for(auto& thread: threads) { thread.join(); }
    }

    double MathEngine::impliedVolatility(double marketPrice, double spot,
                                         double strike, double timeToExpiry,
                                         double riskFreeRate, bool isCall) const
    {
        double vol = 0.25;
        const double tolerance = 1e-6;
        const int maxIterations = 50;

        for(int i = 0; i < maxIterations; ++i) {
            const auto greeks = blackScholes(spot, strike, timeToExpiry,
                                             riskFreeRate, vol, isCall);
            const double price = greeks.price;
            const double vega = greeks.vega;

            if(std::abs(vega) < tolerance) break;

            const double diff = price - marketPrice;
            if(std::abs(diff) < tolerance) break;

            vol -= diff / (vega * 100);
            vol = std::clamp(vol, 0.001, 5.0);
        }

        return vol;
    }

    double MathEngine::realizedVolatility(const std::vector<double>& prices,
                                          int windowDays) const
    {
        if(prices.size() < static_cast<size_t>(windowDays + 1)) { return 0.0; }

        std::vector<double> returns;
        returns.reserve(windowDays);

        const size_t start = prices.size() - windowDays - 1;
        for(size_t i = start + 1; i < prices.size(); ++i) {
            returns.push_back(std::log(prices[i] / prices[i - 1]));
        }

        const double mean = std::accumulate(returns.begin(), returns.end(), 0.0)
                            / returns.size();
        double variance = 0.0;

        for(double ret: returns) {
            const double diff = ret - mean;
            variance += diff * diff;
        }
        variance /= (returns.size() - 1);

        return std::sqrt(variance * 252);
    }

    double MathEngine::normalCDF_scalar(double x) const
    {
        static constexpr double a1 = 0.254829592;
        static constexpr double a2 = -0.284496736;
        static constexpr double a3 = 1.421413741;
        static constexpr double a4 = -1.453152027;
        static constexpr double a5 = 1.061405429;
        static constexpr double p = 0.3275911;
        static constexpr double sqrt2 = 1.4142135623730950488;

        const double sign = (x >= 0) ? 1.0 : -1.0;
        x = std::abs(x) / sqrt2;

        const double t = 1.0 / (1.0 + p * x);
        const double y = 1.0
                         - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t
                               * std::exp(-x * x);

        return 0.5 * (1.0 + sign * y);
    }

    double MathEngine::normalPDF(double x) const
    {
        static constexpr double invSqrt2Pi = 0.3989422804014326779;
        return invSqrt2Pi * std::exp(-0.5 * x * x);
    }

    void MathEngine::normalCDF_vector_scalar(const std::vector<double>& input,
                                             std::vector<double>& output) const
    {
        std::transform(input.begin(), input.end(), output.begin(),
                       [this](double x) { return normalCDF_scalar(x); });
    }

#ifdef __SSE2__
    double MathEngine::normalCDF_sse2(double x) const
    {
        return normalCDF_scalar(x);
    }

    void MathEngine::normalCDF_vector_sse2(const std::vector<double>& input,
                                           std::vector<double>& output) const
    {
        const size_t n = input.size();
        const size_t simdEnd = (n / 2) * 2;

        for(size_t i = 0; i < simdEnd; i += 2) {
            output[i] = normalCDF_scalar(input[i]);
            output[i + 1] = normalCDF_scalar(input[i + 1]);
        }

        for(size_t i = simdEnd; i < n; ++i) {
            output[i] = normalCDF_scalar(input[i]);
        }
    }
#endif

#ifdef __AVX2__
    double MathEngine::normalCDF_avx2(double x) const
    {
        return normalCDF_scalar(x);
    }

    void MathEngine::normalCDF_vector_avx2(const std::vector<double>& input,
                                           std::vector<double>& output) const
    {
        const size_t n = input.size();
        const size_t simdEnd = (n / 4) * 4;

        for(size_t i = 0; i < simdEnd; i += 4) {
            for(size_t j = 0; j < 4; ++j) {
                output[i + j] = normalCDF_scalar(input[i + j]);
            }
        }

        for(size_t i = simdEnd; i < n; ++i) {
            output[i] = normalCDF_scalar(input[i]);
        }
    }
#endif

    void MathEngine::initializeOptimizations() {}

    Greeks
    MathEngine::aggregateGreeks(const std::vector<Greeks>& positions,
                                const std::vector<double>& quantities) const
    {
        Greeks total;

        for(size_t i = 0; i < positions.size() && i < quantities.size(); ++i) {
            total.delta += positions[i].delta * quantities[i];
            total.gamma += positions[i].gamma * quantities[i];
            total.theta += positions[i].theta * quantities[i];
            total.vega += positions[i].vega * quantities[i];
            total.rho += positions[i].rho * quantities[i];
        }

        return total;
    }

    VolatilityMetrics MathEngine::calculateVolatilityMetrics(
        const std::vector<double>& prices,
        const std::vector<double>& impliedVols) const
    {
        VolatilityMetrics metrics;

        if(prices.size() < 20) return metrics;

        metrics.realized = realizedVolatility(prices);

        if(!impliedVols.empty()) {
            metrics.implied
                = std::accumulate(impliedVols.begin(), impliedVols.end(), 0.0)
                  / impliedVols.size();

            std::vector<double> sortedIVs = impliedVols;
            std::sort(sortedIVs.begin(), sortedIVs.end());
            metrics.ivRank = percentileRank(metrics.implied, sortedIVs);
        }

        return metrics;
    }

    double MathEngine::percentileRank(double value,
                                      const std::vector<double>& dataset) const
    {
        if(dataset.empty()) return 0.0;

        auto it = std::lower_bound(dataset.begin(), dataset.end(), value);
        double rank = static_cast<double>(std::distance(dataset.begin(), it))
                      / dataset.size();
        return std::clamp(rank, 0.0, 1.0);
    }

} // namespace math