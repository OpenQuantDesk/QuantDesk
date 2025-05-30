#include "memory.hpp"
#include <stdexcept>
#include <algorithm>

namespace math::core {

template<Size Alignment>
template<typename T>
T* AlignedAllocator<Alignment>::allocate(Size count) {
    Size size = count * sizeof(T);
    Size alignedSize = (size + Alignment - 1) & ~(Alignment - 1);
    void* ptr = std::aligned_alloc(Alignment, alignedSize);
    if (!ptr) throw std::bad_alloc();
    return static_cast<T*>(ptr);
}

template<Size Alignment>
template<typename T>
void AlignedAllocator<Alignment>::deallocate(T* ptr) {
    std::free(ptr);
}

template class AlignedAllocator<CACHE_LINE_SIZE>;
template class AlignedAllocator<32>;
template class AlignedAllocator<16>;

}