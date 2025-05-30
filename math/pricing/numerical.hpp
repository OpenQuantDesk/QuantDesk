#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include "../core/threading.hpp"
#include <memory>
#include <vector>
#include <array>

namespace math::pricing {

template<core::Size MaxDepth = 1000>
class BinomialEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    mutable core::ObjectPool<core::Vector> vectorPool_;
    
public:
    explicit BinomialEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real american(const core::MarketData& market, core::Size steps = 100) const;
    core::Real european(const core::MarketData& market, core::Size steps = 100) const;
    core::Real barrier(const core::MarketData& market, core::Real barrierLevel, 
                      bool isKnockOut, bool isUp, core::Size steps = 100) const;
    
    core::Greeks americanGreeks(const core::MarketData& market, core::Size steps = 100) const;
    
    void priceBatch(const core::MarketDataVector& markets, core::Vector& prices,
                   core::Size steps = 100, bool isAmerican = false) const;
    
private:
    struct TreeParams {
        core::Real dt;
        core::Real u;
        core::Real d;
        core::Real p;
        core::Real discount;
    };
    
    TreeParams calculateParams(const core::MarketData& market, core::Size steps) const;
    core::Real calculatePayoff(core::Real spot, core::Real strike, bool isCall) const;
    core::Real earlyExerciseValue(core::Real spot, core::Real strike, bool isCall) const;
    
    void buildTree(const TreeParams& params, const core::MarketData& market,
                  core::Vector& tree, core::Size steps, bool isAmerican) const;
};

template<core::Size MaxDepth = 1000>
class TrinomialEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    mutable core::ObjectPool<core::Matrix> matrixPool_;
    
public:
    explicit TrinomialEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real american(const core::MarketData& market, core::Size steps = 100) const;
    core::Real european(const core::MarketData& market, core::Size steps = 100) const;
    core::Real bermudan(const core::MarketData& market, const core::Vector& exerciseTimes,
                       core::Size steps = 100) const;
    
    core::Greeks americanGreeks(const core::MarketData& market, core::Size steps = 100) const;
    
private:
    struct TrinomialParams {
        core::Real dt;
        core::Real u;
        core::Real d;
        core::Real pu;
        core::Real pm;
        core::Real pd;
        core::Real discount;
        core::Real lambda;
    };
    
    TrinomialParams calculateParams(const core::MarketData& market, core::Size steps) const;
    void buildTrinomialTree(const TrinomialParams& params, const core::MarketData& market,
                           core::Matrix& tree, core::Size steps, bool isAmerican) const;
};

class FiniteDifferenceEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit FiniteDifferenceEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real explicitFD(const core::MarketData& market, core::Size spotSteps = 100,
                         core::Size timeSteps = 100, core::Real spotMultiplier = 3.0) const;
    
    core::Real implicitFD(const core::MarketData& market, core::Size spotSteps = 100,
                         core::Size timeSteps = 100, core::Real spotMultiplier = 3.0) const;
    
    core::Real crankNicolsonFD(const core::MarketData& market, core::Size spotSteps = 100,
                              core::Size timeSteps = 100, core::Real spotMultiplier = 3.0) const;
    
    core::Greeks finiteDifferenceGreeks(const core::MarketData& market, 
                                       core::Size spotSteps = 100,
                                       core::Size timeSteps = 100) const;
    
    core::Real americanFD(const core::MarketData& market, core::Size spotSteps = 100,
                         core::Size timeSteps = 100) const;
    
private:
    struct FDGrid {
        core::Matrix values;
        core::Vector spotGrid;
        core::Vector timeGrid;
        core::Real dS;
        core::Real dt;
        core::Size M;
        core::Size N;
    };
    
    FDGrid setupGrid(const core::MarketData& market, core::Size spotSteps,
                    core::Size timeSteps, core::Real spotMultiplier) const;
    
    void applyBoundaryConditions(FDGrid& grid, const core::MarketData& market) const;
    void applyInitialConditions(FDGrid& grid, const core::MarketData& market) const;
    
    void explicitStep(FDGrid& grid, const core::MarketData& market, core::Size timeIndex) const;
    void implicitStep(FDGrid& grid, const core::MarketData& market, core::Size timeIndex) const;
    void crankNicolsonStep(FDGrid& grid, const core::MarketData& market, core::Size timeIndex) const;
    
    void solveTridiagonal(const core::Vector& a, const core::Vector& b, 
                         const core::Vector& c, core::Vector& d) const;
    
    core::Real interpolatePrice(const FDGrid& grid, core::Real spot) const;
    core::Real calculateEarlyExercise(const FDGrid& grid, const core::MarketData& market,
                                     core::Size timeIndex, core::Size spotIndex) const;
};

class AdaptiveEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    std::unique_ptr<BinomialEngine<>> binomialEngine_;
    std::unique_ptr<TrinomialEngine<>> trinomialEngine_;
    std::unique_ptr<FiniteDifferenceEngine> fdEngine_;
    
public:
    explicit AdaptiveEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real price(const core::MarketData& market, core::Real tolerance = 1e-6) const;
    core::Greeks greeks(const core::MarketData& market, core::Real tolerance = 1e-6) const;
    
    enum class Method { Auto, Binomial, Trinomial, FiniteDifference };
    
    core::Real priceWithMethod(const core::MarketData& market, Method method,
                              core::Size steps = 100) const;
    
private:
    Method selectOptimalMethod(const core::MarketData& market) const;
    core::Real refinePrice(const core::MarketData& market, Method method,
                          core::Real tolerance) const;
    bool hasConverged(const core::Vector& prices, core::Real tolerance) const;
};

template<core::Size MaxDepth>
inline BinomialEngine<MaxDepth>::BinomialEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

template<core::Size MaxDepth>
inline typename BinomialEngine<MaxDepth>::TreeParams 
BinomialEngine<MaxDepth>::calculateParams(const core::MarketData& market, core::Size steps) const {
    TreeParams params;
    params.dt = market.timeToExpiry / steps;
    params.u = std::exp(market.volatility * std::sqrt(params.dt));
    params.d = 1.0 / params.u;
    params.p = (std::exp(market.riskFreeRate * params.dt) - params.d) / (params.u - params.d);
    params.discount = std::exp(-market.riskFreeRate * params.dt);
    return params;
}

template<core::Size MaxDepth>
inline core::Real BinomialEngine<MaxDepth>::calculatePayoff(core::Real spot, core::Real strike, bool isCall) const {
    return isCall ? std::max(0.0, spot - strike) : std::max(0.0, strike - spot);
}

template<core::Size MaxDepth>
inline core::Real BinomialEngine<MaxDepth>::earlyExerciseValue(core::Real spot, core::Real strike, bool isCall) const {
    return calculatePayoff(spot, strike, isCall);
}

template<core::Size MaxDepth>
inline core::Real BinomialEngine<MaxDepth>::american(const core::MarketData& market, core::Size steps) const {
    if (steps > MaxDepth) steps = MaxDepth;
    
    const TreeParams params = calculateParams(market, steps);
    auto* tree = vectorPool_.acquire();
    tree->resize(steps + 1);
    
    for (core::Size i = 0; i <= steps; ++i) {
        const core::Real spot = market.spot * std::pow(params.u, static_cast<core::Real>(steps - i)) 
                                            * std::pow(params.d, static_cast<core::Real>(i));
        (*tree)[i] = calculatePayoff(spot, market.strike, market.isCall);
    }
    
    for (core::Integer step = static_cast<core::Integer>(steps) - 1; step >= 0; --step) {
        for (core::Size i = 0; i <= static_cast<core::Size>(step); ++i) {
            const core::Real spot = market.spot * std::pow(params.u, static_cast<core::Real>(step - i)) 
                                                * std::pow(params.d, static_cast<core::Real>(i));
            
            const core::Real holdValue = params.discount * (params.p * (*tree)[i] + (1.0 - params.p) * (*tree)[i + 1]);
            const core::Real exerciseValue = earlyExerciseValue(spot, market.strike, market.isCall);
            
            (*tree)[i] = std::max(holdValue, exerciseValue);
        }
    }
    
    const core::Real result = (*tree)[0];
    vectorPool_.release(tree);
    return result;
}

template<core::Size MaxDepth>
inline core::Real BinomialEngine<MaxDepth>::european(const core::MarketData& market, core::Size steps) const {
    if (steps > MaxDepth) steps = MaxDepth;
    
    const TreeParams params = calculateParams(market, steps);
    auto* tree = vectorPool_.acquire();
    tree->resize(steps + 1);
    
    for (core::Size i = 0; i <= steps; ++i) {
        const core::Real spot = market.spot * std::pow(params.u, static_cast<core::Real>(steps - i)) 
                                            * std::pow(params.d, static_cast<core::Real>(i));
        (*tree)[i] = calculatePayoff(spot, market.strike, market.isCall);
    }
    
    for (core::Integer step = static_cast<core::Integer>(steps) - 1; step >= 0; --step) {
        for (core::Size i = 0; i <= static_cast<core::Size>(step); ++i) {
            (*tree)[i] = params.discount * (params.p * (*tree)[i] + (1.0 - params.p) * (*tree)[i + 1]);
        }
    }
    
    const core::Real result = (*tree)[0];
    vectorPool_.release(tree);
    return result;
}

template<core::Size MaxDepth>
inline core::Greeks BinomialEngine<MaxDepth>::americanGreeks(const core::MarketData& market, core::Size steps) const {
    core::Greeks greeks;
    
    const core::Real basePrice = american(market, steps);
    const core::Real deltaShift = market.spot * 0.01;
    const core::Real vegaShift = 0.01;
    const core::Real thetaShift = 1.0 / 365.0;
    const core::Real rhoShift = 0.01;
    
    core::MarketData upMarket = market;
    upMarket.spot += deltaShift;
    const core::Real upPrice = american(upMarket, steps);
    
    core::MarketData downMarket = market;
    downMarket.spot -= deltaShift;
    const core::Real downPrice = american(downMarket, steps);
    
    greeks.price = basePrice;
    greeks.delta = (upPrice - downPrice) / (2.0 * deltaShift);
    greeks.gamma = (upPrice - 2.0 * basePrice + downPrice) / (deltaShift * deltaShift);
    
    if (market.timeToExpiry > thetaShift) {
        core::MarketData thetaMarket = market;
        thetaMarket.timeToExpiry -= thetaShift;
        greeks.theta = american(thetaMarket, steps) - basePrice;
    }
    
    core::MarketData vegaMarket = market;
    vegaMarket.volatility += vegaShift;
    greeks.vega = american(vegaMarket, steps) - basePrice;
    
    core::MarketData rhoMarket = market;
    rhoMarket.riskFreeRate += rhoShift;
    greeks.rho = american(rhoMarket, steps) - basePrice;
    
    greeks.impliedVol = market.volatility;
    
    return greeks;
}

inline FiniteDifferenceEngine::FiniteDifferenceEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

inline FiniteDifferenceEngine::FDGrid 
FiniteDifferenceEngine::setupGrid(const core::MarketData& market, core::Size spotSteps,
                                 core::Size timeSteps, core::Real spotMultiplier) const {
    FDGrid grid;
    grid.M = spotSteps;
    grid.N = timeSteps;
    
    const core::Real Smax = market.spot * spotMultiplier;
    grid.dS = Smax / grid.M;
    grid.dt = market.timeToExpiry / grid.N;
    
    grid.values.resize(grid.M + 1, core::Vector(grid.N + 1, 0.0));
    grid.spotGrid.resize(grid.M + 1);
    grid.timeGrid.resize(grid.N + 1);
    
    for (core::Size i = 0; i <= grid.M; ++i) {
        grid.spotGrid[i] = i * grid.dS;
    }
    
    for (core::Size j = 0; j <= grid.N; ++j) {
        grid.timeGrid[j] = j * grid.dt;
    }
    
    return grid;
}

inline void FiniteDifferenceEngine::applyInitialConditions(FDGrid& grid, const core::MarketData& market) const {
    for (core::Size i = 0; i <= grid.M; ++i) {
        const core::Real spot = grid.spotGrid[i];
        if (market.isCall) {
            grid.values[i][grid.N] = std::max(0.0, spot - market.strike);
        } else {
            grid.values[i][grid.N] = std::max(0.0, market.strike - spot);
        }
    }
}

inline void FiniteDifferenceEngine::applyBoundaryConditions(FDGrid& grid, const core::MarketData& market) const {
    for (core::Size j = 0; j <= grid.N; ++j) {
        const core::Real tau = market.timeToExpiry - grid.timeGrid[j];
        
        if (market.isCall) {
            grid.values[0][j] = 0.0;
            grid.values[grid.M][j] = grid.spotGrid[grid.M] - market.strike * std::exp(-market.riskFreeRate * tau);
        } else {
            grid.values[0][j] = market.strike * std::exp(-market.riskFreeRate * tau);
            grid.values[grid.M][j] = 0.0;
        }
    }
}

inline core::Real FiniteDifferenceEngine::explicitFD(const core::MarketData& market, 
                                                     core::Size spotSteps, core::Size timeSteps,
                                                     core::Real spotMultiplier) const {
    FDGrid grid = setupGrid(market, spotSteps, timeSteps, spotMultiplier);
    applyInitialConditions(grid, market);
    applyBoundaryConditions(grid, market);
    
    for (core::Size j = grid.N; j > 0; --j) {
        explicitStep(grid, market, j - 1);
    }
    
    return interpolatePrice(grid, market.spot);
}

inline void FiniteDifferenceEngine::explicitStep(FDGrid& grid, const core::MarketData& market, 
                                                 core::Size timeIndex) const {
    const core::Real r = market.riskFreeRate;
    const core::Real sigma = market.volatility;
    const core::Real dt = grid.dt;
    const core::Real dS = grid.dS;
    
    for (core::Size i = 1; i < grid.M; ++i) {
        const core::Real S = grid.spotGrid[i];
        const core::Real a = 0.5 * dt * (sigma * sigma * S * S / (dS * dS) - r * S / dS);
        const core::Real b = 1.0 - dt * (sigma * sigma * S * S / (dS * dS) + r);
        const core::Real c = 0.5 * dt * (sigma * sigma * S * S / (dS * dS) + r * S / dS);
        
        grid.values[i][timeIndex] = a * grid.values[i - 1][timeIndex + 1] +
                                   b * grid.values[i][timeIndex + 1] +
                                   c * grid.values[i + 1][timeIndex + 1];
    }
}

inline core::Real FiniteDifferenceEngine::interpolatePrice(const FDGrid& grid, core::Real spot) const {
    const core::Real index = spot / grid.dS;
    const core::Size i = static_cast<core::Size>(index);
    
    if (i >= grid.M) return grid.values[grid.M][0];
    if (i == 0) return grid.values[0][0];
    
    const core::Real weight = index - i;
    return (1.0 - weight) * grid.values[i][0] + weight * grid.values[i + 1][0];
}

}