#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include <memory>
#include <functional>

namespace math::pricing {

class RainbowEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit RainbowEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real maxOption(const core::Vector& spots, const core::Vector& strikes,
                        core::Real timeToExpiry, core::Real riskFreeRate,
                        const core::Matrix& correlations, const core::Vector& vols,
                        const core::Vector& dividends = {}) const;
    
    core::Real minOption(const core::Vector& spots, const core::Vector& strikes,
                        core::Real timeToExpiry, core::Real riskFreeRate,
                        const core::Matrix& correlations, const core::Vector& vols,
                        const core::Vector& dividends = {}) const;
    
    core::Real spreadOption(core::Real S1, core::Real S2, core::Real K,
                           core::Real timeToExpiry, core::Real riskFreeRate,
                           core::Real vol1, core::Real vol2, core::Real correlation,
                           core::Real dividend1 = 0.0, core::Real dividend2 = 0.0) const;
    
    core::Real basketCall(const core::Vector& spots, const core::Vector& weights,
                         core::Real strike, core::Real timeToExpiry, core::Real riskFreeRate,
                         const core::Matrix& correlations, const core::Vector& vols,
                         const core::Vector& dividends = {}) const;
    
private:
    core::Real monteCarloRainbow(const core::Vector& spots, const core::Vector& strikes,
                                core::Real timeToExpiry, core::Real riskFreeRate,
                                const core::Matrix& correlations, const core::Vector& vols,
                                const core::Vector& dividends, bool isCall = true,
                                bool isMax = true, core::Integer numSimulations = 100000) const;
    
    core::Matrix choleskyDecomposition(const core::Matrix& correlation) const;
    core::Real bivariateCDF(core::Real x, core::Real y, core::Real rho) const;
};

class BermudanEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit BermudanEngine(std::shared_ptr<core::MathEngine> engine);
    
    struct ExerciseParams {
        core::Integer timeSteps = 252;
        core::Integer numPaths = 100000;
        bool useControlVariates = false;
        bool useAntitheticVariates = true;
    };
    
    core::Real price(const core::MarketData& market, const core::Vector& exerciseTimes,
                    const ExerciseParams& params = ExerciseParams{}) const;
    
    core::Real lsmPrice(const core::MarketData& market, const core::Vector& exerciseTimes,
                       const ExerciseParams& params = ExerciseParams{}) const;
    
    core::Greeks greeks(const core::MarketData& market, const core::Vector& exerciseTimes,
                       const ExerciseParams& params = ExerciseParams{}) const;
    
private:
    core::Real calculateContinuationValue(const core::Matrix& spotPaths,
                                         const core::Vector& payoffs,
                                         core::Size timeIndex) const;
};

class QuantoEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit QuantoEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real quantoOption(const core::MarketData& market, core::Real fxVol,
                           core::Real correlation, core::Real foreignRate) const;
    
    core::Greeks quantoGreeks(const core::MarketData& market, core::Real fxVol,
                             core::Real correlation, core::Real foreignRate) const;
    
    core::Real currencyTranslatedOption(const core::MarketData& market,
                                       core::Real fxSpot, core::Real fxVol,
                                       core::Real correlation) const;
};

class CompositeEngine {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit CompositeEngine(std::shared_ptr<core::MathEngine> engine);
    
    core::Real powerOption(const core::MarketData& market, core::Real power) const;
    
    core::Real logOption(const core::MarketData& market) const;
    
    core::Real chooseOption(const core::MarketData& market, core::Real chooseTimes) const;
    
    core::Real cliqueOption(const core::Vector& strikes, const core::Vector& spots,
                           const core::Vector& times, core::Real riskFreeRate,
                           const core::Matrix& correlations, const core::Vector& vols) const;
    
    core::Real shoutOption(const core::MarketData& market, 
                          const core::Vector& shoutTimes) const;
    
private:
    core::Real calculateOptimalExercise(const core::Vector& path,
                                       const core::Vector& exerciseTimes) const;
};

inline RainbowEngine::RainbowEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

inline BermudanEngine::BermudanEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

inline QuantoEngine::QuantoEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

inline CompositeEngine::CompositeEngine(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

}