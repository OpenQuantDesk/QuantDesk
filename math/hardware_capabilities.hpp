#pragma once

#include <cstddef>

namespace math {

    class HardwareCapabilities {
    private:
        bool sse2_;
        bool sse4_1_;
        bool avx_;
        bool avx2_;
        bool fma_;
        bool neon_;
        bool opencl_;

        size_t cacheLineSize_;
        size_t l1CacheSize_;
        size_t l2CacheSize_;
        size_t l3CacheSize_;
        size_t numCores_;
        size_t numThreads_;
    public:
        HardwareCapabilities();
        ~HardwareCapabilities() = default;

        static HardwareCapabilities detect();

        bool hasSSE2() const { return sse2_; }
        bool hasSSE4_1() const { return sse4_1_; }
        bool hasAVX() const { return avx_; }
        bool hasAVX2() const { return avx2_; }
        bool hasFMA() const { return fma_; }
        bool hasNEON() const { return neon_; }
        bool hasOpenCL() const { return opencl_; }

        size_t getCacheLineSize() const { return cacheLineSize_; }
        size_t getL1CacheSize() const { return l1CacheSize_; }
        size_t getL2CacheSize() const { return l2CacheSize_; }
        size_t getL3CacheSize() const { return l3CacheSize_; }
        size_t getNumCores() const { return numCores_; }
        size_t getNumThreads() const { return numThreads_; }

        void print() const;
    private:
        void detectCpuFeatures();
        void detectCacheInfo();
        void detectCoreInfo();
        void detectOpenCL();
    };

} // namespace math