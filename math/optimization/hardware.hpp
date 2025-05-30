#pragma once
namespace math 
{
    struct HardwareCapabilities {
        bool sse2 = false;
        bool sse4_1 = false;
        bool avx = false;
        bool avx2 = false;
        bool fma = false;
        bool neon = false;
        bool opencl = false;
        
        static HardwareCapabilities detect();
    };
}