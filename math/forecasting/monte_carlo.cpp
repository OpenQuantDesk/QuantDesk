#include "math/forecasting/monte_carlo.hpp"
#include <algorithm>
#include <numeric>
#include <future>
#include <thread>

namespace math {

MonteCarloEngine::MonteCarloEngine() : rng_(std::random_device{}()) {}

SimulationResult MonteCarloEngine::simulate(const SimulationParams& params, 
                                          const std::function<double(double)>& payoffFunc) {
    SimulationResult result;
    std::vector<double> payoffs;
    payoffs.reserve(params.numSimulations);
    
    const double dt = params.timeToExpiry / params.timeSteps;
    const double drift = (params.riskFreeRate - 0.5 * params.volatility * params.volatility) * dt;
    const double diffusion = params.volatility * std::sqrt(dt);
    
    std::normal_distribution<double> normal(0.0, 1.0);
    
    const size_t numThreads = std::thread::hardware_concurrency();
    const size_t simsPerThread = params.numSimulations / numThreads;
    
    std::vector<std::future<std::vector<double>>> futures;
    futures.reserve(numThreads);
    
    for (size_t i = 0; i < numThreads; ++i) {
        const size_t start = i * simsPerThread;
        const size_t end = (i == numThreads - 1) ? params.numSimulations : (i + 1) * simsPerThread;
        
        futures.emplace_back(std::async(std::launch::async, [=, &payoffFunc]() {
            std::mt19937 localRng(rng_() + i);
            std::normal_distribution<double> localNormal(0.0, 1.0);
            std::vector<double> localPayoffs;
            localPayoffs.reserve(end - start);
            
            for (size_t sim = start; sim < end; ++sim) {
                double spot = params.spot;
                
                for (int step = 0; step < params.timeSteps; ++step) {
                    double z = localNormal(localRng);
                    spot *= std::exp(drift + diffusion * z);
                    
                    if (params.useAntitheticVariates && sim % 2 == 1) {
                        spot = params.spot * std::exp(drift - diffusion * z);
                    }
                }
                
                localPayoffs.push_back(payoffFunc(spot));
            }
            
            return localPayoffs;
        }));
    }
    
    for (auto& future : futures) {
        auto localPayoffs = future.get();
        payoffs.insert(payoffs.end(), localPayoffs.begin(), localPayoffs.end());
    }
    
    result.payoffs = payoffs;
    result.expectedValue = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / payoffs.size();
    
    double variance = 0.0;
    for (double p : payoffs) {
        variance += (p - result.expectedValue) * (p - result.expectedValue);
    }
    variance /= (payoffs.size() - 1);
    result.standardError = std::sqrt(variance / payoffs.size());
    
    int profitable = std::count_if(payoffs.begin(), payoffs.end(), 
                                  [](double p) { return p > 0; });
    result.probabilityProfit = static_cast<double>(profitable) / payoffs.size();
    
    std::sort(payoffs.begin(), payoffs.end());
    result.valueAtRisk95 = payoffs[static_cast<size_t>(payoffs.size() * 0.05)];
    result.valueAtRisk99 = payoffs[static_cast<size_t>(payoffs.size() * 0.01)];
    
    return result;
}

PathGenerator::PathGenerator() : rng_(std::random_device{}()) {}

std::vector<double> PathGenerator::generatePath(const SimulationParams& params) {
    std::vector<double> path;
    path.reserve(params.timeSteps + 1);
    path.push_back(params.spot);
    
    const double dt = params.timeToExpiry / params.timeSteps;
    const double drift = (params.riskFreeRate - 0.5 * params.volatility * params.volatility) * dt;
    const double diffusion = params.volatility * std::sqrt(dt);
    
    std::normal_distribution<double> normal(0.0, 1.0);
    
    double spot = params.spot;
    for (int step = 0; step < params.timeSteps; ++step) {
        double z = normal(rng_);
        spot *= std::exp(drift + diffusion * z);
        path.push_back(spot);
    }
    
    return path;
}

std::vector<std::vector<double>> PathGenerator::generatePaths(const SimulationParams& params, 
                                                             int numPaths) {
    std::vector<std::vector<double>> paths;
    paths.reserve(numPaths);
    
    const size_t numThreads = std::thread::hardware_concurrency();
    const size_t pathsPerThread = numPaths / numThreads;
    
    std::vector<std::future<std::vector<std::vector<double>>>> futures;
    futures.reserve(numThreads);
    
    for (size_t i = 0; i < numThreads; ++i) {
        const size_t start = i * pathsPerThread;
        const size_t end = (i == numThreads - 1) ? numPaths : (i + 1) * pathsPerThread;
        
        futures.emplace_back(std::async(std::launch::async, [=]() {
            PathGenerator localGenerator;
            std::vector<std::vector<double>> localPaths;
            localPaths.reserve(end - start);
            
            for (size_t j = start; j < end; ++j) {
                localPaths.push_back(localGenerator.generatePath(params));
            }
            
            return localPaths;
        }));
    }
    
    for (auto& future : futures) {
        auto localPaths = future.get();
        paths.insert(paths.end(), localPaths.begin(), localPaths.end());
    }
    
    return paths;
}

}