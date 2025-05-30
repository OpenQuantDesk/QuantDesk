#pragma once

#include "../core/types.hpp"
#include "params.hpp"
#include "results.hpp"

namespace math::simulation {

class VarianceReducer {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit VarianceReducer(std::shared_ptr<core::MathEngine> engine);
    
    core::SimulationResult antitheticVariates(
        const core::SimulationParams& params,
        const std::function<core::Real(core::Real)>& payoffFunc) const;
    
    core::SimulationResult controlVariates(
        const core::SimulationParams& params,
        const std::function<core::Real(core::Real)>& payoffFunc,
        const std::function<core::Real(core::Real)>& controlFunc) const;
    
    core::SimulationResult importanceSampling(
        const core::SimulationParams& params,
        const std::function<core::Real(core::Real)>& payoffFunc,
        core::Real optimalDrift) const;
    
    struct VarianceMetrics {
        core::Real standardMC = 0.0;
        core::Real antitheticMC = 0.0;
        core::Real controlVariateMC = 0.0;
        core::Real varianceReduction = 0.0;
        core::Real efficiency = 0.0;
    };
    
    VarianceMetrics compareVarianceReduction(
        const core::SimulationParams& params,
        const std::function<core::Real(core::Real)>& payoffFunc) const;

private:
    core::SimulationResult standardMonteCarlo(
        const core::SimulationParams& params,
        const std::function<core::Real(core::Real)>& payoffFunc) const;
};

}