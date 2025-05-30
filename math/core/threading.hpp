#pragma once

#include "types.hpp"
#include "memory.hpp"
#include <thread>
#include <atomic>
#include <functional>
#include <future>
#include <condition_variable>
#include <mutex>
#include <queue>

namespace math::core {

template<typename T>
class WorkStealingQueue {
private:
    std::deque<T> queue_;
    mutable std::mutex mutex_;
    
public:
    WorkStealingQueue() = default;
    
    WorkStealingQueue(const WorkStealingQueue&) = delete;
    WorkStealingQueue& operator=(const WorkStealingQueue&) = delete;
    
    void push(T item) {
        std::lock_guard<std::mutex> lock(mutex_);
        queue_.push_front(std::move(item));
    }
    
    bool tryPop(T& item) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (queue_.empty()) {
            return false;
        }
        item = std::move(queue_.front());
        queue_.pop_front();
        return true;
    }
    
    bool trySteal(T& item) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (queue_.empty()) {
            return false;
        }
        item = std::move(queue_.back());
        queue_.pop_back();
        return true;
    }
    
    bool empty() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.empty();
    }
    
    Size size() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.size();
    }
};

class ThreadPool {
private:
    using Task = std::function<void()>;
    
    std::vector<std::thread> workers_;
    std::vector<std::unique_ptr<WorkStealingQueue<Task>>> localQueues_;
    LockFreeQueue<Task> globalQueue_;
    
    std::atomic<bool> done_{false};
    std::atomic<Size> activeThreads_{0};
    
    static thread_local Size workerIndex_;
    static thread_local WorkStealingQueue<Task>* localQueue_;
    
public:
    explicit ThreadPool(Size numThreads = std::thread::hardware_concurrency()) 
        : workers_(numThreads), localQueues_(numThreads) {
        
        try {
            for (Size i = 0; i < numThreads; ++i) {
                localQueues_[i] = std::make_unique<WorkStealingQueue<Task>>();
                workers_[i] = std::thread(&ThreadPool::workerThread, this, i);
            }
        } catch (...) {
            done_.store(true, std::memory_order_release);
            throw;
        }
    }
    
    ~ThreadPool() {
        done_.store(true, std::memory_order_release);
        
        for (auto& worker : workers_) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }
    
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    
    template<typename F, typename... Args>
    auto submit(F&& f, Args&&... args) -> std::future<std::invoke_result_t<F, Args...>> {
        using ReturnType = std::invoke_result_t<F, Args...>;
        
        auto task = std::make_shared<std::packaged_task<ReturnType()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        
        auto future = task->get_future();
        
        if (localQueue_) {
            localQueue_->push([task]() { (*task)(); });
        } else {
            globalQueue_.enqueue([task]() { (*task)(); });
        }
        
        return future;
    }
    
    template<typename Iterator, typename Function>
    void parallelFor(Iterator first, Iterator last, Function func) {
        const auto length = std::distance(first, last);
        if (length == 0) return;
        
        const Size numThreads = workers_.size();
        const Size blockSize = std::max(Size(1), Size(length) / numThreads);
        
        std::vector<std::future<void>> futures;
        futures.reserve(numThreads);
        
        auto current = first;
        while (current != last) {
            auto blockEnd = current;
            std::advance(blockEnd, std::min(blockSize, Size(std::distance(current, last))));
            
            futures.push_back(submit([current, blockEnd, func]() {
                for (auto it = current; it != blockEnd; ++it) {
                    func(*it);
                }
            }));
            
            current = blockEnd;
        }
        
        for (auto& future : futures) {
            future.wait();
        }
    }
    
    template<typename Iterator, typename Function, typename Reducer>
    auto parallelReduce(Iterator first, Iterator last, Function func, Reducer reducer) 
        -> std::invoke_result_t<Function, typename std::iterator_traits<Iterator>::value_type> {
        
        using ValueType = std::invoke_result_t<Function, typename std::iterator_traits<Iterator>::value_type>;
        
        const auto length = std::distance(first, last);
        if (length == 0) return ValueType{};
        if (length == 1) return func(*first);
        
        const Size numThreads = workers_.size();
        const Size blockSize = std::max(Size(1), Size(length) / numThreads);
        
        std::vector<std::future<ValueType>> futures;
        futures.reserve(numThreads);
        
        auto current = first;
        while (current != last) {
            auto blockEnd = current;
            std::advance(blockEnd, std::min(blockSize, Size(std::distance(current, last))));
            
            futures.push_back(submit([current, blockEnd, func, reducer]() -> ValueType {
                if (current == blockEnd) return ValueType{};
                
                ValueType result = func(*current);
                ++current;
                
                for (auto it = current; it != blockEnd; ++it) {
                    result = reducer(result, func(*it));
                }
                
                return result;
            }));
            
            current = blockEnd;
        }
        
        ValueType finalResult = futures[0].get();
        for (Size i = 1; i < futures.size(); ++i) {
            finalResult = reducer(finalResult, futures[i].get());
        }
        
        return finalResult;
    }
    
    Size getNumThreads() const {
        return workers_.size();
    }
    
    Size getActiveThreads() const {
        return activeThreads_.load(std::memory_order_acquire);
    }
    
    void waitForAll() {
        while (activeThreads_.load(std::memory_order_acquire) > 0) {
            std::this_thread::yield();
        }
    }
    
private:
    void workerThread(Size index) {
        workerIndex_ = index;
        localQueue_ = localQueues_[index].get();
        
        while (!done_.load(std::memory_order_acquire)) {
            runPendingTask();
        }
    }
    
    void runPendingTask() {
        Task task;
        
        if (popTaskFromLocalQueue(task) || 
            popTaskFromOtherQueue(task) || 
            popTaskFromGlobalQueue(task)) {
            
            activeThreads_.fetch_add(1, std::memory_order_acquire);
            task();
            activeThreads_.fetch_sub(1, std::memory_order_release);
        } else {
            std::this_thread::yield();
        }
    }
    
    bool popTaskFromLocalQueue(Task& task) {
        return localQueue_ && localQueue_->tryPop(task);
    }
    
    bool popTaskFromOtherQueue(Task& task) {
        for (Size i = 0; i < localQueues_.size(); ++i) {
            const Size index = (workerIndex_ + i + 1) % localQueues_.size();
            if (localQueues_[index]->trySteal(task)) {
                return true;
            }
        }
        return false;
    }
    
    bool popTaskFromGlobalQueue(Task& task) {
        return globalQueue_.dequeue(task);
    }
};

thread_local Size ThreadPool::workerIndex_;
thread_local WorkStealingQueue<ThreadPool::Task>* ThreadPool::localQueue_;

class TaskGroup {
private:
    std::atomic<Size> pendingTasks_{0};
    std::mutex completionMutex_;
    std::condition_variable completionCV_;
    
public:
    TaskGroup() = default;
    
    TaskGroup(const TaskGroup&) = delete;
    TaskGroup& operator=(const TaskGroup&) = delete;
    
    template<typename F, typename... Args>
    void run(ThreadPool& pool, F&& f, Args&&... args) {
        pendingTasks_.fetch_add(1, std::memory_order_acquire);
        
        pool.submit([this, f = std::forward<F>(f), args...]() mutable {
            try {
                f(args...);
            } catch (...) {
                
            }
            
            Size remaining = pendingTasks_.fetch_sub(1, std::memory_order_release) - 1;
            if (remaining == 0) {
                std::lock_guard<std::mutex> lock(completionMutex_);
                completionCV_.notify_all();
            }
        });
    }
    
    void wait() {
        std::unique_lock<std::mutex> lock(completionMutex_);
        completionCV_.wait(lock, [this] {
            return pendingTasks_.load(std::memory_order_acquire) == 0;
        });
    }
    
    bool waitFor(const std::chrono::milliseconds& timeout) {
        std::unique_lock<std::mutex> lock(completionMutex_);
        return completionCV_.wait_for(lock, timeout, [this] {
            return pendingTasks_.load(std::memory_order_acquire) == 0;
        });
    }
    
    Size getPendingTasks() const {
        return pendingTasks_.load(std::memory_order_acquire);
    }
};

template<typename T>
class ThreadLocalStorage {
private:
    thread_local static std::unique_ptr<T> instance_;
    std::function<std::unique_ptr<T>()> factory_;
    
public:
    template<typename F>
    explicit ThreadLocalStorage(F&& factory) 
        : factory_(std::forward<F>(factory)) {}
    
    T& get() {
        if (!instance_) {
            instance_ = factory_();
        }
        return *instance_;
    }
    
    T* operator->() {
        return &get();
    }
    
    T& operator*() {
        return get();
    }
};

template<typename T>
thread_local std::unique_ptr<T> ThreadLocalStorage<T>::instance_;

class SpinLock {
private:
    std::atomic_flag flag_ = ATOMIC_FLAG_INIT;
    
public:
    void lock() {
        while (flag_.test_and_set(std::memory_order_acquire)) {
            while (flag_.test(std::memory_order_relaxed)) {
                std::this_thread::yield();
            }
        }
    }
    
    bool try_lock() {
        return !flag_.test_and_set(std::memory_order_acquire);
    }
    
    void unlock() {
        flag_.clear(std::memory_order_release);
    }
};

class ReadWriteLock {
private:
    std::atomic<int> readers_{0};
    std::atomic<bool> writer_{false};
    
public:
    void readLock() {
        while (true) {
            while (writer_.load(std::memory_order_acquire)) {
                std::this_thread::yield();
            }
            
            readers_.fetch_add(1, std::memory_order_acquire);
            
            if (!writer_.load(std::memory_order_acquire)) {
                break;
            }
            
            readers_.fetch_sub(1, std::memory_order_release);
        }
    }
    
    void readUnlock() {
        readers_.fetch_sub(1, std::memory_order_release);
    }
    
    void writeLock() {
        bool expected = false;
        while (!writer_.compare_exchange_weak(expected, true, 
                                             std::memory_order_acquire,
                                             std::memory_order_relaxed)) {
            expected = false;
            std::this_thread::yield();
        }
        
        while (readers_.load(std::memory_order_acquire) > 0) {
            std::this_thread::yield();
        }
    }
    
    void writeUnlock() {
        writer_.store(false, std::memory_order_release);
    }
};

template<typename Mutex>
class SharedLockGuard {
private:
    Mutex& mutex_;
    
public:
    explicit SharedLockGuard(Mutex& m) : mutex_(m) {
        mutex_.readLock();
    }
    
    ~SharedLockGuard() {
        mutex_.readUnlock();
    }
    
    SharedLockGuard(const SharedLockGuard&) = delete;
    SharedLockGuard& operator=(const SharedLockGuard&) = delete;
};

template<typename Mutex>
class UniqueLockGuard {
private:
    Mutex& mutex_;
    
public:
    explicit UniqueLockGuard(Mutex& m) : mutex_(m) {
        mutex_.writeLock();
    }
    
    ~UniqueLockGuard() {
        mutex_.writeUnlock();
    }
    
    UniqueLockGuard(const UniqueLockGuard&) = delete;
    UniqueLockGuard& operator=(const UniqueLockGuard&) = delete;
};

class Barrier {
private:
    const Size count_;
    std::atomic<Size> waiting_{0};
    std::atomic<Size> generation_{0};
    
public:
    explicit Barrier(Size count) : count_(count) {}
    
    void wait() {
        const Size gen = generation_.load(std::memory_order_acquire);
        const Size waitingCount = waiting_.fetch_add(1, std::memory_order_acquire) + 1;
        
        if (waitingCount == count_) {
            waiting_.store(0, std::memory_order_release);
            generation_.fetch_add(1, std::memory_order_release);
        } else {
            while (generation_.load(std::memory_order_acquire) == gen) {
                std::this_thread::yield();
            }
        }
    }
};

using GlobalThreadPool = ThreadPool;

}