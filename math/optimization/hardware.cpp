#include "math/optimization/hardware.hpp"
#include <iostream>

#ifdef __x86_64__
#include <cpuid.h>
#endif

namespace math {

HardwareCapabilities HardwareCapabilities::detect() {
    HardwareCapabilities caps;
    
#ifdef __x86_64__
    unsigned int eax, ebx, ecx, edx;
    
    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        caps.sse2 = (edx & bit_SSE2) != 0;
        caps.sse4_1 = (ecx & bit_SSE4_1) != 0;
        caps.avx = (ecx & bit_AVX) != 0;
        caps.fma = (ecx & bit_FMA) != 0;
    }
    
    if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
        caps.avx2 = (ebx & bit_AVX2) != 0;
    }
#endif

#ifdef __ARM_NEON
    caps.neon = true;
#endif

#ifdef ENABLE_OPENCL
    try {
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        caps.opencl = !platforms.empty();
    } catch (...) {
        caps.opencl = false;
    }
#endif

    return caps;
}

void HardwareCapabilities::print() const {
    std::cout << "Hardware Capabilities:" << std::endl;
    std::cout << "  SSE2:   " << (sse2 ? "Yes" : "No") << std::endl;
    std::cout << "  SSE4.1: " << (sse4_1 ? "Yes" : "No") << std::endl;
    std::cout << "  AVX:    " << (avx ? "Yes" : "No") << std::endl;
    std::cout << "  AVX2:   " << (avx2 ? "Yes" : "No") << std::endl;
    std::cout << "  FMA:    " << (fma ? "Yes" : "No") << std::endl;
    std::cout << "  NEON:   " << (neon ? "Yes" : "No") << std::endl;
    std::cout << "  OpenCL: " << (opencl ? "Yes" : "No") << std::endl;
}

}