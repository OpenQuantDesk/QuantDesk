#pragma once

#include "math/simulation/params.hpp"
#include "math/simulation/results.hpp"
#include <vector>
#include <functional>
#include <random>
#include <future>
#include <thread>

namespace math {

class MonteCarloEngine {
private:
    mutable std::mt19937 rng_;
    
public:
    MonteCarloEngine();
    ~MonteCarloEngine() = default;
    
    MonteCarloEngine(const MonteCarloEngine&) = delete;
    MonteCarloEngine& operator=(const MonteCarloEngine&) = delete;
    
    SimulationResult simulate(const SimulationParams& params, 
                             const std::function<double(double)>& payoffFunc);
    
    std::vector<double> generatePricePath(const SimulationParams& params);
    
    std::vector<SimulationResult> runBatch(
        const std::vector<SimulationParams>& paramsBatch,
        const std::vector<std::function<double(double)>>& payoffFuncs);
    
private:
    std::vector<double> runSimulationBatch(
        const SimulationParams& params,
        const std::function<double(double)>& payoffFunc,
        size_t numSims,
        size_t threadId) const;
};

class PathGenerator {
private:
    mutable std::mt19937 rng_;
    
public:
    PathGenerator();
    ~PathGenerator() = default;
    
    PathGenerator(const PathGenerator&) = delete;
    PathGenerator& operator=(const PathGenerator&) = delete;
    
    std::vector<double> generatePath(const SimulationParams& params);
    
    std::vector<std::vector<double>> generatePaths(const SimulationParams& params, 
                                                  int numPaths);
    
    std::vector<double> generateCorrelatedPath(const SimulationParams& params,
                                              const std::vector<double>& correlations);
};

}