#pragma once

#include "../core/types.hpp"
#include <memory>

namespace math::optimization {

class CacheOptimizer {
private:
    core::Size cacheLineSize_;
    core::Size l1Size_;
    core::Size l2Size_;
    core::Size l3Size_;

public:
    explicit CacheOptimizer(core::Size cacheLineSize = 64);
    
    template<typename T>
    void optimizeMemoryLayout(std::vector<T>& data) const;
    
    core::Size getOptimalBlockSize(core::Size elementSize) const;
    core::Size getOptimalTileSize(core::Size matrixSize) const;
};

}
