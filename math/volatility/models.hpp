#pragma once

#include "../core/types.hpp"
#include "../core/engine.hpp"
#include "surface.hpp"
#include <complex>
#include <functional>
#include <memory>

namespace math::volatility {

class BlackScholesModel {
private:
    std::shared_ptr<core::MathEngine> engine_;
    core::Matrix impliedVols_;
    core::Vector strikes_;
    core::Vector maturities_;
    core::Real spot_;
    core::Real riskFreeRate_;
    
public:
    explicit BlackScholesModel(std::shared_ptr<core::MathEngine> engine);
    
    core::Real impliedVolatility(const core::MarketData& market, core::Real marketPrice) const;
    core::Real theoreticalPrice(const core::MarketData& market) const;
    
    void calibrate(const core::Vector& strikes, const core::Vector& maturities,
                  const core::Matrix& marketPrices, core::Real spot, core::Real riskFreeRate);
    
    core::Real getVolatility(core::Real strike, core::Real maturity) const;
    
    const core::Matrix& getImpliedVolatilities() const { return impliedVols_; }
    const core::Vector& getStrikes() const { return strikes_; }
    const core::Vector& getMaturities() const { return maturities_; }
};

class HestonModel {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit HestonModel(std::shared_ptr<core::MathEngine> engine);
    
    struct HestonParams {
        core::Real kappa = 2.0;     // Mean reversion speed
        core::Real theta = 0.04;    // Long-term variance
        core::Real sigma = 0.3;     // Volatility of volatility
        core::Real rho = -0.5;      // Correlation
        core::Real v0 = 0.04;       // Initial variance
    };
    
    std::complex<core::Real> characteristicFunction(const std::complex<core::Real>& u,
                                                   core::Real tau,
                                                   const HestonParams& params) const;
    
    core::Real optionPrice(core::Real spot, core::Real strike, core::Real timeToExpiry,
                          core::Real riskFreeRate, const HestonParams& params,
                          bool isCall = true) const;
    
    core::Real impliedVolatility(core::Real spot, core::Real strike, core::Real timeToExpiry,
                                core::Real riskFreeRate, const HestonParams& params,
                                bool isCall = true) const;
    
    HestonParams calibrate(const VolatilityPoints& marketData, core::Real spot,
                          core::Real riskFreeRate) const;
    
    bool validateParams(const HestonParams& params) const;
    HestonParams enforceConstraints(const HestonParams& params) const;
    
    core::Vector simulatePath(const HestonParams& params, core::Real spot,
                             core::Real timeToExpiry, core::Integer timeSteps,
                             std::mt19937& rng) const;
};

class StochasticVolatilityModel {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit StochasticVolatilityModel(std::shared_ptr<core::MathEngine> engine);
    
    struct SABRParams {
        core::Real alpha = 0.2;     // ATM volatility
        core::Real beta = 0.5;      // CEV exponent
        core::Real rho = 0.0;       // Correlation
        core::Real nu = 0.3;        // Volatility of volatility
    };
    
    core::Real sabr(core::Real forward, core::Real strike, core::Real timeToExpiry,
                   const SABRParams& params) const;
    
    SABRParams calibrateSABR(const core::Vector& strikes, const core::Vector& impliedVols,
                            core::Real forward, core::Real timeToExpiry) const;
    
    bool validateSABRParams(const SABRParams& params) const;
    SABRParams enforceSABRConstraints(const SABRParams& params) const;
    
    struct BergomiParams {
        core::Real H = 0.1;         // Hurst parameter
        core::Real eta = 0.3;       // Volatility of volatility
        core::Real rho = -0.5;      // Correlation
        core::Real xi = 0.04;       // Initial log-variance
    };
    
    core::Real bergomiVolatility(core::Real timeToExpiry, const BergomiParams& params) const;
    
private:
    core::Real sabrImpliedVolatility(core::Real forward, core::Real strike, core::Real timeToExpiry,
                                    const SABRParams& params) const;
};

class JumpDiffusionModel {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit JumpDiffusionModel(std::shared_ptr<core::MathEngine> engine);
    
    struct JumpParams {
        core::Real intensity = 0.1;          // Jump intensity (lambda)
        core::Real meanJumpSize = -0.05;     // Mean jump size (mu_J)
        core::Real jumpVolatility = 0.15;    // Jump volatility (sigma_J)
    };
    
    core::Real mertonPrice(core::Real spot, core::Real strike, core::Real timeToExpiry,
                          core::Real riskFreeRate, core::Real sigma,
                          const JumpParams& jumpParams, bool isCall = true) const;
    
    JumpParams calibrateJumps(const VolatilityPoints& marketData, core::Real spot,
                             core::Real riskFreeRate, core::Real sigma) const;
    
    struct KouParams {
        core::Real intensity = 0.1;
        core::Real p = 0.5;              // Probability of upward jump
        core::Real lambda1 = 10.0;       // Upward jump decay
        core::Real lambda2 = 5.0;        // Downward jump decay
    };
    
    core::Real kouPrice(core::Real spot, core::Real strike, core::Real timeToExpiry,
                       core::Real riskFreeRate, core::Real sigma,
                       const KouParams& kouParams, bool isCall = true) const;
    
private:
    core::Real calculateJumpCompensation(const JumpParams& params) const;
    core::Real kouCharacteristicFunction(const std::complex<core::Real>& u,
                                        const KouParams& params) const;
};

class LocalVolatilityModel {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit LocalVolatilityModel(std::shared_ptr<core::MathEngine> engine);
    
    core::Real dupireLocalVol(core::Real spot, core::Real timeToExpiry,
                             const std::function<core::Real(core::Real, core::Real)>& impliedVol) const;
    
    core::Matrix buildLocalVolSurface(const core::Vector& spots, const core::Vector& times,
                                     const std::function<core::Real(core::Real, core::Real)>& impliedVol) const;
    
    core::Real localVolatility(core::Real spot, core::Real timeToExpiry) const;
    
    void calibrateFromImpliedSurface(std::shared_ptr<VolatilitySurface> impliedSurface);
    
    core::Vector simulatePath(core::Real spot, core::Real timeToExpiry,
                             core::Integer timeSteps, std::mt19937& rng) const;
    
private:
    core::Matrix localVolSurface_;
    core::Vector spotGrid_;
    core::Vector timeGrid_;
    
    core::Real interpolateLocalVol(core::Real spot, core::Real time) const;
    core::Real numericalDerivative(const std::function<core::Real(core::Real)>& func,
                                  core::Real x, core::Real h = 1e-6) const;
};

class VarianceGammaModel {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit VarianceGammaModel(std::shared_ptr<core::MathEngine> engine);
    
    struct VGParams {
        core::Real sigma = 0.2;     // Volatility
        core::Real nu = 0.1;        // Variance rate
        core::Real theta = -0.1;    // Drift
    };
    
    std::complex<core::Real> characteristicFunction(const std::complex<core::Real>& u,
                                                   core::Real timeToExpiry,
                                                   const VGParams& params) const;
    
    core::Real optionPrice(core::Real spot, core::Real strike, core::Real timeToExpiry,
                          core::Real riskFreeRate, const VGParams& params,
                          bool isCall = true) const;
    
    VGParams calibrate(const VolatilityPoints& marketData, core::Real spot,
                      core::Real riskFreeRate) const;
    
    core::Vector simulatePath(core::Real spot, core::Real timeToExpiry,
                             core::Integer timeSteps, const VGParams& params,
                             std::mt19937& rng) const;
    
private:
    core::Real gammaVariate(core::Real shape, core::Real scale, std::mt19937& rng) const;
    core::Real calculateMartingaleCorrection(const VGParams& params) const;
};

class StochasticLocalVolModel {
private:
    std::shared_ptr<core::MathEngine> engine_;
    std::shared_ptr<HestonModel> hestonModel_;
    std::shared_ptr<LocalVolatilityModel> localVolModel_;
    
public:
    explicit StochasticLocalVolModel(std::shared_ptr<core::MathEngine> engine);
    
    struct SLVParams {
        HestonModel::HestonParams hestonParams;
        core::Real correlation = 0.0;    // Correlation between SV and LV
        core::Real leverage = 1.0;       // Leverage function parameter
    };
    
    core::Real optionPrice(core::Real spot, core::Real strike, core::Real timeToExpiry,
                          core::Real riskFreeRate, const SLVParams& params,
                          bool isCall = true) const;
    
    SLVParams calibrate(const VolatilityPoints& marketData, core::Real spot,
                       core::Real riskFreeRate) const;
    
    core::Vector simulatePath(core::Real spot, core::Real timeToExpiry,
                             core::Integer timeSteps, const SLVParams& params,
                             std::mt19937& rng) const;
    
private:
    core::Real leverageFunction(core::Real spot, core::Real variance, const SLVParams& params) const;
    void calibrateLeverageFunction(const VolatilityPoints& marketData);
};

class VolatilityInterpolator {
private:
    std::shared_ptr<core::MathEngine> engine_;
    
public:
    explicit VolatilityInterpolator(std::shared_ptr<core::MathEngine> engine);
    
    enum class InterpolationType {
        Linear,
        Cubic,
        Variance,
        SVI,
        SSVI
    };
    
    core::Real interpolate(const core::Vector& strikes, const core::Vector& vols,
                          core::Real targetStrike, InterpolationType type = InterpolationType::Cubic) const;
    
    core::Vector interpolateTermStructure(const core::Vector& times, const core::Vector& atmVols,
                                         core::Real targetTime, InterpolationType type = InterpolationType::Variance) const;
    
    core::Matrix interpolateSurface(const core::Vector& strikes, const core::Vector& times,
                                   const core::Matrix& vols, const core::Vector& targetStrikes,
                                   const core::Vector& targetTimes) const;
    
private:
    core::Real varianceInterpolation(const core::Vector& strikes, const core::Vector& vols,
                                    core::Real targetStrike, core::Real timeToExpiry) const;
    
    core::Real sviInterpolation(const core::Vector& strikes, const core::Vector& vols,
                               core::Real targetStrike, core::Real forward) const;
};

inline BlackScholesModel::BlackScholesModel(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)), spot_(0.0), riskFreeRate_(0.0) {}

inline bool StochasticVolatilityModel::validateSABRParams(const SABRParams& params) const {
    return params.alpha > 0.0 && params.beta >= 0.0 && params.beta <= 1.0 &&
           params.rho >= -1.0 && params.rho <= 1.0 && params.nu >= 0.0;
}

inline StochasticVolatilityModel::SABRParams 
StochasticVolatilityModel::enforceSABRConstraints(const SABRParams& params) const {
    SABRParams constrained = params;
    constrained.alpha = std::max(0.001, constrained.alpha);
    constrained.beta = std::clamp(constrained.beta, 0.0, 1.0);
    constrained.rho = std::clamp(constrained.rho, -0.99, 0.99);
    constrained.nu = std::max(0.0, constrained.nu);
    return constrained;
}

inline core::Real StochasticVolatilityModel::bergomiVolatility(core::Real timeToExpiry,
                                                              const BergomiParams& params) const {
    const core::Real variance = std::exp(params.xi + params.eta * std::sqrt(2.0 * params.H) * 
                                        std::pow(timeToExpiry, params.H));
    return std::sqrt(variance);
}

inline core::Real JumpDiffusionModel::calculateJumpCompensation(const JumpParams& params) const {
    return params.intensity * (std::exp(params.meanJumpSize + 0.5 * params.jumpVolatility * params.jumpVolatility) - 1.0);
}

inline VolatilityInterpolator::VolatilityInterpolator(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

}