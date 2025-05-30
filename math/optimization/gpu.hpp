#pragma once

#include "../core/types.hpp"

namespace math::optimization {

class GPUAccelerator {
private:
    bool available_;
    
public:
    GPUAccelerator();
    
    bool isAvailable() const { return available_; }
    
    void vectorAdd(const core::Real* a, const core::Real* b, 
                   core::Real* result, core::Size count) const;
    
    void matrixMultiply(const core::Real* a, const core::Real* b, 
                       core::Real* c, core::Size m, core::Size n, core::Size k) const;
};

}