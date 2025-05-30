#pragma once

#include "types.hpp"
#include <memory>
#include <atomic>
#include <mutex>
#include <array>
#include <cstdlib>

namespace math::core {

template<Size Alignment = CACHE_LINE_SIZE>
class AlignedAllocator {
public:
    template<typename T>
    static T* allocate(Size count) {
        Size size = count * sizeof(T);
        Size alignedSize = (size + Alignment - 1) & ~(Alignment - 1);
        void* ptr = std::aligned_alloc(Alignment, alignedSize);
        if (!ptr) throw std::bad_alloc();
        return static_cast<T*>(ptr);
    }
    
    template<typename T>
    static void deallocate(T* ptr) {
        std::free(ptr);
    }
};

template<typename T, Size BlockSize = 4096>
class ObjectPool {
private:
    struct Block {
        alignas(CACHE_LINE_SIZE) std::array<std::byte, sizeof(T)> data;
        std::atomic<Block*> next{nullptr};
    };
    
    std::atomic<Block*> freeHead_{nullptr};
    std::vector<std::unique_ptr<Block[]>> chunks_;
    std::mutex chunkMutex_;
    static constexpr Size objectsPerChunk_ = BlockSize / sizeof(Block);
    
public:
    ObjectPool() {
        allocateChunk();
    }
    
    ~ObjectPool() = default;
    
    ObjectPool(const ObjectPool&) = delete;
    ObjectPool& operator=(const ObjectPool&) = delete;
    
    template<typename... Args>
    T* acquire(Args&&... args) {
        Block* block = freeHead_.load(std::memory_order_acquire);
        
        while (block) {
            if (freeHead_.compare_exchange_weak(block, block->next.load(std::memory_order_relaxed),
                                               std::memory_order_release, std::memory_order_acquire)) {
                return new(block->data.data()) T(std::forward<Args>(args)...);
            }
        }
        
        allocateChunk();
        return acquire(std::forward<Args>(args)...);
    }
    
    void release(T* obj) {
        if (!obj) return;
        
        obj->~T();
        Block* block = reinterpret_cast<Block*>(obj);
        Block* head = freeHead_.load(std::memory_order_relaxed);
        
        do {
            block->next.store(head, std::memory_order_relaxed);
        } while (!freeHead_.compare_exchange_weak(head, block,
                                                 std::memory_order_release,
                                                 std::memory_order_relaxed));
    }
    
private:
    void allocateChunk() {
        std::lock_guard<std::mutex> lock(chunkMutex_);
        
        auto chunk = std::make_unique<Block[]>(objectsPerChunk_);
        Block* blocks = chunk.get();
        
        for (Size i = 0; i < objectsPerChunk_ - 1; ++i) {
            blocks[i].next.store(&blocks[i + 1], std::memory_order_relaxed);
        }
        blocks[objectsPerChunk_ - 1].next.store(nullptr, std::memory_order_relaxed);
        
        Block* head = freeHead_.load(std::memory_order_relaxed);
        do {
            blocks[objectsPerChunk_ - 1].next.store(head, std::memory_order_relaxed);
        } while (!freeHead_.compare_exchange_weak(head, &blocks[0],
                                                 std::memory_order_release,
                                                 std::memory_order_relaxed));
        
        chunks_.push_back(std::move(chunk));
    }
};

template<Size NumSizeClasses = 16>
class SlabAllocator {
private:
    static constexpr std::array<Size, NumSizeClasses> sizeClasses_ = {
        16, 32, 64, 128, 256, 512, 1024, 2048,
        4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288
    };
    
    struct Slab {
        void* memory;
        Size objectSize;
        Size objectCount;
        std::atomic<void*> freeHead{nullptr};
        std::atomic<Size> freeCount{0};
        Slab* next{nullptr};
    };
    
    std::array<std::atomic<Slab*>, NumSizeClasses> slabs_;
    std::mutex slabMutex_;
    
public:
    SlabAllocator() {
        for (auto& slab : slabs_) {
            slab.store(nullptr, std::memory_order_relaxed);
        }
    }
    
    ~SlabAllocator() {
        for (Size i = 0; i < NumSizeClasses; ++i) {
            Slab* slab = slabs_[i].load(std::memory_order_relaxed);
            while (slab) {
                Slab* next = slab->next;
                std::free(slab->memory);
                delete slab;
                slab = next;
            }
        }
    }
    
    SlabAllocator(const SlabAllocator&) = delete;
    SlabAllocator& operator=(const SlabAllocator&) = delete;
    
    void* allocate(Size size) {
        Size sizeClass = getSizeClass(size);
        if (sizeClass == NumSizeClasses) {
            return std::aligned_alloc(CACHE_LINE_SIZE, size);
        }
        
        Slab* slab = slabs_[sizeClass].load(std::memory_order_acquire);
        
        while (slab) {
            void* obj = slab->freeHead.load(std::memory_order_acquire);
            if (obj) {
                void* next = *static_cast<void**>(obj);
                if (slab->freeHead.compare_exchange_weak(obj, next,
                                                        std::memory_order_release,
                                                        std::memory_order_acquire)) {
                    slab->freeCount.fetch_sub(1, std::memory_order_relaxed);
                    return obj;
                }
            }
            slab = slab->next;
        }
        
        return allocateFromNewSlab(sizeClass);
    }
    
    void deallocate(void* ptr, Size size) {
        if (!ptr) return;
        
        Size sizeClass = getSizeClass(size);
        if (sizeClass == NumSizeClasses) {
            std::free(ptr);
            return;
        }
        
        Slab* slab = findSlab(ptr, sizeClass);
        if (slab) {
            void* head = slab->freeHead.load(std::memory_order_relaxed);
            do {
                *static_cast<void**>(ptr) = head;
            } while (!slab->freeHead.compare_exchange_weak(head, ptr,
                                                          std::memory_order_release,
                                                          std::memory_order_relaxed));
            slab->freeCount.fetch_add(1, std::memory_order_relaxed);
        }
    }
    
private:
    Size getSizeClass(Size size) const {
        for (Size i = 0; i < NumSizeClasses; ++i) {
            if (size <= sizeClasses_[i]) {
                return i;
            }
        }
        return NumSizeClasses;
    }
    
    void* allocateFromNewSlab(Size sizeClass) {
        std::lock_guard<std::mutex> lock(slabMutex_);
        
        Size objectSize = sizeClasses_[sizeClass];
        Size objectCount = 4096 / objectSize;
        if (objectCount == 0) objectCount = 1;
        
        Size totalSize = objectCount * objectSize;
        void* memory = std::aligned_alloc(CACHE_LINE_SIZE, totalSize);
        if (!memory) throw std::bad_alloc();
        
        Slab* newSlab = new Slab{memory, objectSize, objectCount, nullptr, objectCount - 1, nullptr};
        
        char* obj = static_cast<char*>(memory);
        for (Size i = 0; i < objectCount - 1; ++i) {
            *reinterpret_cast<void**>(obj + i * objectSize) = obj + (i + 1) * objectSize;
        }
        *reinterpret_cast<void**>(obj + (objectCount - 1) * objectSize) = nullptr;
        newSlab->freeHead.store(obj + objectSize, std::memory_order_relaxed);
        
        Slab* head = slabs_[sizeClass].load(std::memory_order_relaxed);
        newSlab->next = head;
        slabs_[sizeClass].store(newSlab, std::memory_order_release);
        
        return obj;
    }
    
    Slab* findSlab(void* ptr, Size sizeClass) const {
        Slab* slab = slabs_[sizeClass].load(std::memory_order_acquire);
        while (slab) {
            char* slabStart = static_cast<char*>(slab->memory);
            char* slabEnd = slabStart + slab->objectCount * slab->objectSize;
            if (ptr >= slabStart && ptr < slabEnd) {
                return slab;
            }
            slab = slab->next;
        }
        return nullptr;
    }
};

template<typename T>
class LockFreeQueue {
private:
    struct Node {
        std::atomic<T*> data{nullptr};
        std::atomic<Node*> next{nullptr};
    };
    
    std::atomic<Node*> head_{nullptr};
    std::atomic<Node*> tail_{nullptr};
    
public:
    LockFreeQueue() {
        Node* dummy = new Node;
        head_.store(dummy, std::memory_order_relaxed);
        tail_.store(dummy, std::memory_order_relaxed);
    }
    
    ~LockFreeQueue() {
        while (Node* head = head_.load(std::memory_order_relaxed)) {
            head_.store(head->next.load(std::memory_order_relaxed), std::memory_order_relaxed);
            delete head;
        }
    }
    
    LockFreeQueue(const LockFreeQueue&) = delete;
    LockFreeQueue& operator=(const LockFreeQueue&) = delete;
    
    void enqueue(T item) {
        Node* newNode = new Node;
        T* data = new T(std::move(item));
        newNode->data.store(data, std::memory_order_relaxed);
        
        while (true) {
            Node* tail = tail_.load(std::memory_order_acquire);
            Node* next = tail->next.load(std::memory_order_acquire);
            
            if (tail == tail_.load(std::memory_order_acquire)) {
                if (next == nullptr) {
                    if (tail->next.compare_exchange_weak(next, newNode,
                                                        std::memory_order_release,
                                                        std::memory_order_relaxed)) {
                        break;
                    }
                } else {
                    tail_.compare_exchange_weak(tail, next,
                                              std::memory_order_release,
                                              std::memory_order_relaxed);
                }
            }
        }
        
        tail_.compare_exchange_weak(tail_.load(std::memory_order_acquire), newNode,
                                   std::memory_order_release,
                                   std::memory_order_relaxed);
    }
    
    bool dequeue(T& result) {
        while (true) {
            Node* head = head_.load(std::memory_order_acquire);
            Node* tail = tail_.load(std::memory_order_acquire);
            Node* next = head->next.load(std::memory_order_acquire);
            
            if (head == head_.load(std::memory_order_acquire)) {
                if (head == tail) {
                    if (next == nullptr) {
                        return false;
                    }
                    tail_.compare_exchange_weak(tail, next,
                                              std::memory_order_release,
                                              std::memory_order_relaxed);
                } else {
                    if (next == nullptr) {
                        continue;
                    }
                    T* data = next->data.load(std::memory_order_acquire);
                    if (data == nullptr) {
                        continue;
                    }
                    if (head_.compare_exchange_weak(head, next,
                                                   std::memory_order_release,
                                                   std::memory_order_relaxed)) {
                        result = *data;
                        delete data;
                        delete head;
                        return true;
                    }
                }
            }
        }
    }
    
    bool empty() const {
        Node* head = head_.load(std::memory_order_acquire);
        Node* tail = tail_.load(std::memory_order_acquire);
        return head == tail && head->next.load(std::memory_order_acquire) == nullptr;
    }
};

using MemoryPool = SlabAllocator<>;
using GreeksPool = ObjectPool<Greeks>;
using MarketDataPool = ObjectPool<MarketData>;

}