#pragma once

#include <vector>
#include <memory>

#ifdef ENABLE_OPENCL
#include <CL/cl.hpp>
#endif

namespace math {

struct HardwareCapabilities {
    bool sse2 = false;
    bool sse4_1 = false;
    bool avx = false;
    bool avx2 = false;
    bool fma = false;
    bool neon = false;
    bool opencl = false;
    
    size_t cacheLineSize = 64;
    size_t l1CacheSize = 32768;
    size_t l2CacheSize = 262144;
    size_t l3CacheSize = 8388608;
    size_t numCores = 1;
    size_t numThreads = 1;
    
    static HardwareCapabilities detect();
    void print() const;
    void optimize();
};

class HardwareOptimizer {
private:
    HardwareCapabilities capabilities_;
    
#ifdef ENABLE_OPENCL
    std::unique_ptr<cl::Context> clContext_;
    std::unique_ptr<cl::CommandQueue> clQueue_;
    std::unique_ptr<cl::Program> clProgram_;
#endif
    
public:
    explicit HardwareOptimizer(const HardwareCapabilities& caps);
    ~HardwareOptimizer() = default;
    
    HardwareOptimizer(const HardwareOptimizer&) = delete;
    HardwareOptimizer& operator=(const HardwareOptimizer&) = delete;
    
    void initializeOpenCL();
    bool isOpenCLAvailable() const;
    
    void vectorizedAdd(const std::vector<double>& a, const std::vector<double>& b, 
                      std::vector<double>& result) const;
    
    void vectorizedMultiply(const std::vector<double>& a, const std::vector<double>& b,
                           std::vector<double>& result) const;
    
    void matrixMultiply(const std::vector<std::vector<double>>& a,
                       const std::vector<std::vector<double>>& b,
                       std::vector<std::vector<double>>& result) const;
    
    size_t getOptimalChunkSize(size_t dataSize) const;
    size_t getOptimalThreadCount() const;
    
    HardwareCapabilities getCapabilities() const { return capabilities_; }
    
private:
    void vectorizedAdd_scalar(const std::vector<double>& a, const std::vector<double>& b,
                             std::vector<double>& result) const;
    
#ifdef __SSE2__
    void vectorizedAdd_sse2(const std::vector<double>& a, const std::vector<double>& b,
                           std::vector<double>& result) const;
#endif

#ifdef __AVX2__
    void vectorizedAdd_avx2(const std::vector<double>& a, const std::vector<double>& b,
                           std::vector<double>& result) const;
#endif

#ifdef ENABLE_OPENCL
    void vectorizedAdd_opencl(const std::vector<double>& a, const std::vector<double>& b,
                             std::vector<double>& result) const;
#endif
};

class MemoryPool {
private:
    struct Block {
        void* ptr;
        size_t size;
        bool inUse;
        Block* next;
    };
    
    Block* freeBlocks_;
    size_t totalSize_;
    size_t alignment_;
    void* basePtr_;
    
public:
    explicit MemoryPool(size_t totalSize, size_t alignment = 64);
    ~MemoryPool();
    
    MemoryPool(const MemoryPool&) = delete;
    MemoryPool& operator=(const MemoryPool&) = delete;
    
    void* allocate(size_t size);
    void deallocate(void* ptr);
    void reset();
    
    size_t getTotalSize() const { return totalSize_; }
    size_t getAvailableSize() const;
    double getFragmentation() const;
    
private:
    void* alignPointer(void* ptr) const;
    Block* findFreeBlock(size_t size);
    void splitBlock(Block* block, size_t size);
    void mergeBlocks();
};

template<typename T>
class AlignedVector {
private:
    T* data_;
    size_t size_;
    size_t capacity_;
    size_t alignment_;
    
public:
    explicit AlignedVector(size_t alignment = 64) 
        : data_(nullptr), size_(0), capacity_(0), alignment_(alignment) {}
    
    ~AlignedVector() {
        if (data_) {
            for (size_t i = 0; i < size_; ++i) {
                data_[i].~T();
            }
            std::aligned_free(data_);
        }
    }
    
    AlignedVector(const AlignedVector&) = delete;
    AlignedVector& operator=(const AlignedVector&) = delete;
    
    AlignedVector(AlignedVector&& other) noexcept
        : data_(other.data_), size_(other.size_), capacity_(other.capacity_), alignment_(other.alignment_) {
        other.data_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }
    
    AlignedVector& operator=(AlignedVector&& other) noexcept {
        if (this != &other) {
            this->~AlignedVector();
            data_ = other.data_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            alignment_ = other.alignment_;
            other.data_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }
        return *this;
    }
    
    void reserve(size_t newCapacity) {
        if (newCapacity > capacity_) {
            T* newData = static_cast<T*>(std::aligned_alloc(alignment_, newCapacity * sizeof(T)));
            if (!newData) throw std::bad_alloc();
            
            for (size_t i = 0; i < size_; ++i) {
                new(&newData[i]) T(std::move(data_[i]));
                data_[i].~T();
            }
            
            if (data_) std::aligned_free(data_);
            data_ = newData;
            capacity_ = newCapacity;
        }
    }
    
    void resize(size_t newSize) {
        if (newSize > capacity_) {
            reserve(newSize * 2);
        }
        
        for (size_t i = size_; i < newSize; ++i) {
            new(&data_[i]) T();
        }
        
        for (size_t i = newSize; i < size_; ++i) {
            data_[i].~T();
        }
        
        size_ = newSize;
    }
    
    void push_back(const T& value) {
        if (size_ >= capacity_) {
            reserve(capacity_ == 0 ? 1 : capacity_ * 2);
        }
        new(&data_[size_]) T(value);
        ++size_;
    }
    
    T& operator[](size_t index) { return data_[index]; }
    const T& operator[](size_t index) const { return data_[index]; }
    
    T* data() { return data_; }
    const T* data() const { return data_; }
    
    size_t size() const { return size_; }
    size_t capacity() const { return capacity_; }
    bool empty() const { return size_ == 0; }
    
    T* begin() { return data_; }
    T* end() { return data_ + size_; }
    const T* begin() const { return data_; }
    const T* end() const { return data_ + size_; }
};

class CacheOptimizer {
private:
    size_t cacheLineSize_;
    size_t l1CacheSize_;
    size_t l2CacheSize_;
    size_t l3CacheSize_;
    
public:
    explicit CacheOptimizer(const HardwareCapabilities& caps);
    
    template<typename T>
    void optimizeLayout(std::vector<T>& data) const {
        if (data.empty()) return;
        
        size_t elementsPerCacheLine = cacheLineSize_ / sizeof(T);
        if (elementsPerCacheLine == 0) elementsPerCacheLine = 1;
        
        if (data.size() % elementsPerCacheLine != 0) {
            size_t newSize = ((data.size() / elementsPerCacheLine) + 1) * elementsPerCacheLine;
            data.resize(newSize);
        }
    }
    
    size_t getOptimalBlockSize(size_t elementSize) const;
    size_t getOptimalTileSize(size_t matrixSize, size_t elementSize) const;
    
    template<typename Func>
    void blockwiseProcess(size_t dataSize, size_t elementSize, Func&& func) const {
        size_t blockSize = getOptimalBlockSize(elementSize);
        
        for (size_t i = 0; i < dataSize; i += blockSize) {
            size_t currentBlockSize = std::min(blockSize, dataSize - i);
            func(i, currentBlockSize);
        }
    }
};

class PrefetchOptimizer {
public:
    static void prefetchRead(const void* addr);
    static void prefetchWrite(void* addr);
    static void prefetchNTA(const void* addr);
    
    template<typename Iterator>
    static void prefetchRange(Iterator begin, Iterator end, size_t stride = 1) {
        auto it = begin;
        size_t count = 0;
        while (it != end) {
            if (count % stride == 0) {
                prefetchRead(&(*it));
            }
            ++it;
            ++count;
        }
    }
};

}