#include "threading.hpp"

namespace math::core {

thread_local Size ThreadPool::workerIndex_ = 0;
thread_local WorkStealingQueue<ThreadPool::Task>* ThreadPool::localQueue_ = nullptr;

template<typename T>
thread_local std::unique_ptr<T> ThreadLocalStorage<T>::instance_;

template class ThreadLocalStorage<std::mt19937>;
template class ThreadPool;
template class WorkStealingQueue<std::function<void()>>;
template class LockFreeQueue<std::function<void()>>;
template class ObjectPool<Greeks>;
template class ObjectPool<MarketData>;

}