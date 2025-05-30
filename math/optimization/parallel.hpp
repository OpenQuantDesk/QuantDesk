#pragma once

#include "../core/types.hpp"
#include "../core/threading.hpp"
#include <functional>

namespace math::optimization {

class ParallelProcessor {
private:
    std::shared_ptr<core::GlobalThreadPool> threadPool_;
    
public:
    explicit ParallelProcessor(core::Size numThreads = 0);
    
    template<typename Iterator, typename Function>
    void parallelFor(Iterator first, Iterator last, Function func) const;
    
    template<typename Iterator, typename Function, typename Reducer>
    auto parallelReduce(Iterator first, Iterator last, 
                       Function func, Reducer reducer) const;
};

}