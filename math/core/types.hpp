#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace math::core {

    using Real = double;
    using Integer = std::int64_t;
    using Size = std::size_t;

    struct Greeks {
        Real price = 0.0;
        Real delta = 0.0;
        Real gamma = 0.0;
        Real theta = 0.0;
        Real vega = 0.0;
        Real rho = 0.0;
        Real impliedVol = 0.0;
    };

    struct ExtendedGreeks : Greeks {
        Real volga = 0.0;
        Real vanna = 0.0;
        Real charm = 0.0;
        Real speed = 0.0;
        Real zomma = 0.0;
        Real color = 0.0;
        Real dollarDelta = 0.0;
        Real dollarGamma = 0.0;
        Real pinRisk = 0.0;
    };

    struct MarketData {
        Real spot = 0.0;
        Real strike = 0.0;
        Real timeToExpiry = 0.0;
        Real riskFreeRate = 0.0;
        Real volatility = 0.0;
        Real dividendYield = 0.0;
        bool isCall = true;
    };

    struct SimulationParams {
        Real spot = 100.0;
        Real volatility = 0.2;
        Real riskFreeRate = 0.05;
        Real timeToExpiry = 1.0;
        Real dividendYield = 0.0;
        Integer numSimulations = 100000;
        Integer timeSteps = 252;
        bool useAntitheticVariates = true;
        bool useControlVariates = false;
        uint32_t seed = 12345;
    };

    struct SimulationResult {
        Real expectedValue = 0.0;
        Real standardError = 0.0;
        Real probabilityProfit = 0.0;
        Real valueAtRisk95 = 0.0;
        Real valueAtRisk99 = 0.0;
        Real maxDrawdown = 0.0;
        std::vector<Real> pricePath;
        std::vector<Real> payoffs;
        Real convergenceRate = 0.0;
    };

    struct VolatilityMetrics {
        Real realized = 0.0;
        Real implied = 0.0;
        Real hvRank = 0.0;
        Real ivRank = 0.0;
        Real skew = 0.0;
        Real termStructure = 0.0;
        Real meanReversion = 0.0;
        Real persistence = 0.0;
    };

    struct ProbabilityMetrics {
        Real probabilityITM = 0.0;
        Real probabilityTouch = 0.0;
        Real expectedMove = 0.0;
        Real breakeven = 0.0;
        Real maxPain = 0.0;
        Real probabilityProfit = 0.0;
        Real probabilityMaxProfit = 0.0;
        Real probabilityMaxLoss = 0.0;
        Real expectedReturn = 0.0;
        Real sharpeRatio = 0.0;
        Real kelly = 0.0;
        Real winRate = 0.0;
        Real avgWin = 0.0;
        Real avgLoss = 0.0;
        Real profitFactor = 0.0;
    };

    struct RiskMetrics {
        Real var95 = 0.0;
        Real var99 = 0.0;
        Real cvar95 = 0.0;
        Real cvar99 = 0.0;
        Real maxDrawdown = 0.0;
        Real sharpeRatio = 0.0;
        Real sortinoRatio = 0.0;
        Real calmarRatio = 0.0;
        Real beta = 0.0;
        Real alpha = 0.0;
        Real trackingError = 0.0;
        Real informationRatio = 0.0;
    };

    enum class OptionType : uint8_t {
        Call = 0,
        Put = 1,
        Binary = 2,
        Barrier = 3,
        Asian = 4,
        Lookback = 5,
        Rainbow = 6
    };

    enum class ExerciseStyle : uint8_t {
        European = 0,
        American = 1,
        Bermudan = 2
    };

    enum class BarrierType : uint8_t {
        UpAndOut = 0,
        UpAndIn = 1,
        DownAndOut = 2,
        DownAndIn = 3
    };

    struct OptionSpec {
        OptionType type = OptionType::Call;
        ExerciseStyle exercise = ExerciseStyle::European;
        Real strike = 100.0;
        Real barrier = 0.0;
        Real rebate = 0.0;
        BarrierType barrierType = BarrierType::UpAndOut;
        std::vector<Real> exerciseDates;
    };

    template<typename T> using AlignedVector = std::vector<T>;

    using GreeksVector = AlignedVector<Greeks>;
    using ExtendedGreeksVector = AlignedVector<ExtendedGreeks>;
    using MarketDataVector = AlignedVector<MarketData>;

    template<Size N> using StaticArray = std::array<Real, N>;

    using Matrix = std::vector<std::vector<Real>>;
    using Vector = std::vector<Real>;

    struct HardwareInfo {
        bool hasAVX512 = false;
        bool hasAVX2 = false;
        bool hasAVX = false;
        bool hasSSE42 = false;
        bool hasFMA = false;
        bool hasNEON = false;
        Size numCores = 1;
        Size numThreads = 1;
        Size cacheLineSize = 64;
        Size l1CacheSize = 32768;
        Size l2CacheSize = 262144;
        Size l3CacheSize = 8388608;
    };

    constexpr Real PI = 3.14159265358979323846;
    constexpr Real SQRT_2PI = 2.50662827463100050241;
    constexpr Real INV_SQRT_2PI = 0.39894228040143267794;
    constexpr Real SQRT_2 = 1.41421356237309504880;
    constexpr Real LN_2 = 0.69314718055994530942;
    constexpr Real EULER = 2.71828182845904523536;

    constexpr Real EPSILON = 1e-15;
    constexpr Real TOLERANCE = 1e-8;
    constexpr Integer MAX_ITERATIONS = 100;

    constexpr Size CACHE_LINE_SIZE = 64;
    constexpr Size SIMD_WIDTH_AVX512 = 8;
    constexpr Size SIMD_WIDTH_AVX2 = 4;
    constexpr Size SIMD_WIDTH_SSE = 2;

} // namespace math::core