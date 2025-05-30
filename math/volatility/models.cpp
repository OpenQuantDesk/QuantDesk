#include "models.hpp"
#include "../core/utils.hpp"
#include <complex>
#include <algorithm>
#include <numeric>

namespace math::volatility {

BlackScholesModel::BlackScholesModel(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real BlackScholesModel::impliedVolatility(const core::MarketData& market, 
                                               core::Real marketPrice) const {
    return engine_->impliedVolatility(marketPrice, market);
}

core::Real BlackScholesModel::theoreticalPrice(const core::MarketData& market) const {
    return engine_->blackScholes(market).price;
}

void BlackScholesModel::calibrate(const core::Vector& strikes, 
                                 const core::Vector& maturities,
                                 const core::Matrix& marketPrices,
                                 core::Real spot, core::Real riskFreeRate) {
    spot_ = spot;
    riskFreeRate_ = riskFreeRate;
    
    const core::Size numStrikes = strikes.size();
    const core::Size numMaturities = maturities.size();
    
    impliedVols_.resize(numStrikes, core::Vector(numMaturities));
    
    for (core::Size i = 0; i < numStrikes; ++i) {
        for (core::Size j = 0; j < numMaturities; ++j) {
            core::MarketData market;
            market.spot = spot;
            market.strike = strikes[i];
            market.timeToExpiry = maturities[j];
            market.riskFreeRate = riskFreeRate;
            market.isCall = true;
            
            impliedVols_[i][j] = impliedVolatility(market, marketPrices[i][j]);
        }
    }
}

core::Real BlackScholesModel::getVolatility(core::Real strike, core::Real maturity) const {
    return core::Interpolation::bilinear(strike, maturity, 
                                        strikes_[0], strikes_.back(),
                                        maturities_[0], maturities_.back(),
                                        impliedVols_[0][0], impliedVols_[0].back(),
                                        impliedVols_.back()[0], impliedVols_.back().back());
}

HestonModel::HestonModel(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

std::complex<core::Real> HestonModel::characteristicFunction(const std::complex<core::Real>& u,
                                                           core::Real tau,
                                                           const HestonParams& params) const {
    const auto i = std::complex<core::Real>(0.0, 1.0);
    
    const auto d = std::sqrt((params.rho * params.sigma * i * u - params.kappa) * 
                           (params.rho * params.sigma * i * u - params.kappa) - 
                           params.sigma * params.sigma * (i * u + u * u));
    
    const auto g = (params.kappa - params.rho * params.sigma * i * u - d) /
                   (params.kappa - params.rho * params.sigma * i * u + d);
    
    const auto C = params.kappa * params.theta / (params.sigma * params.sigma) *
                   ((params.kappa - params.rho * params.sigma * i * u - d) * tau -
                    2.0 * std::log((1.0 - g * std::exp(-d * tau)) / (1.0 - g)));
    
    const auto D = (params.kappa - params.rho * params.sigma * i * u - d) /
                   (params.sigma * params.sigma) * 
                   (1.0 - std::exp(-d * tau)) / (1.0 - g * std::exp(-d * tau));
    
    return std::exp(C + D * params.v0);
}

core::Real HestonModel::optionPrice(core::Real spot, core::Real strike, core::Real timeToExpiry,
                                   core::Real riskFreeRate, const HestonParams& params,
                                   bool isCall) const {
    const core::Real forward = spot * std::exp(riskFreeRate * timeToExpiry);
    const core::Real logMoneyness = std::log(forward / strike);
    
    auto integrand1 = [&](core::Real u) -> core::Real {
        const auto phi = characteristicFunction(std::complex<core::Real>(u - 1.0, -1.0), 
                                              timeToExpiry, params);
        return std::real(std::exp(-std::complex<core::Real>(0.0, 1.0) * u * logMoneyness) * 
                        phi / (std::complex<core::Real>(0.0, 1.0) * u));
    };
    
    auto integrand2 = [&](core::Real u) -> core::Real {
        const auto phi = characteristicFunction(std::complex<core::Real>(u, -1.0), 
                                              timeToExpiry, params);
        return std::real(std::exp(-std::complex<core::Real>(0.0, 1.0) * u * logMoneyness) * 
                        phi / (std::complex<core::Real>(0.0, 1.0) * u));
    };
    
    const core::Real P1 = 0.5 + core::NumericalMethods::trapezoidalRule(integrand1, 0.0, 100.0) / core::PI;
    const core::Real P2 = 0.5 + core::NumericalMethods::trapezoidalRule(integrand2, 0.0, 100.0) / core::PI;
    
    const core::Real discountFactor = std::exp(-riskFreeRate * timeToExpiry);
    
    if (isCall) {
        return spot * P1 - strike * discountFactor * P2;
    } else {
        return strike * discountFactor * (1.0 - P2) - spot * (1.0 - P1);
    }
}

HestonModel::HestonParams HestonModel::calibrate(const VolatilityPoints& marketData,
                                                core::Real spot, core::Real riskFreeRate) const {
    HestonParams params;
    params.kappa = 2.0;
    params.theta = 0.04;
    params.sigma = 0.3;
    params.rho = -0.5;
    params.v0 = 0.04;
    
    const core::Integer maxIterations = 100;
    const core::Real tolerance = 1e-6;
    
    for (core::Integer iter = 0; iter < maxIterations; ++iter) {
        core::Real totalError = 0.0;
        
        for (const auto& point : marketData) {
            const core::Real modelPrice = optionPrice(spot, point.strike, point.timeToExpiry,
                                                     riskFreeRate, params, point.isCall);
            const core::Real error = modelPrice - point.marketPrice;
            totalError += error * error;
        }
        
        if (std::sqrt(totalError / marketData.size()) < tolerance) {
            break;
        }
        
        const core::Real stepSize = 0.01;
        params.kappa += stepSize * (totalError > 0 ? -1.0 : 1.0);
        params = enforceConstraints(params);
    }
    
    return params;
}

HestonModel::HestonParams HestonModel::enforceConstraints(const HestonParams& params) const {
    HestonParams constrained = params;
    
    constrained.kappa = std::clamp(constrained.kappa, 0.1, 10.0);
    constrained.theta = std::clamp(constrained.theta, 0.01, 1.0);
    constrained.sigma = std::clamp(constrained.sigma, 0.01, 2.0);
    constrained.rho = std::clamp(constrained.rho, -0.99, 0.99);
    constrained.v0 = std::clamp(constrained.v0, 0.001, 1.0);
    
    if (2.0 * constrained.kappa * constrained.theta < constrained.sigma * constrained.sigma) {
        constrained.theta = constrained.sigma * constrained.sigma / (2.0 * constrained.kappa) + 0.001;
    }
    
    return constrained;
}

core::Vector HestonModel::simulatePath(const HestonParams& params, core::Real spot,
                                      core::Real timeToExpiry, core::Integer timeSteps,
                                      std::mt19937& rng) const {
    core::Vector spotPath(timeSteps + 1);
    core::Vector varPath(timeSteps + 1);
    
    spotPath[0] = spot;
    varPath[0] = params.v0;
    
    const core::Real dt = timeToExpiry / timeSteps;
    const core::Real sqrtDt = std::sqrt(dt);
    const core::Real rhoComp = std::sqrt(1.0 - params.rho * params.rho);
    
    std::normal_distribution<core::Real> normal(0.0, 1.0);
    
    for (core::Integer step = 1; step <= timeSteps; ++step) {
        const core::Real z1 = normal(rng);
        const core::Real z2 = normal(rng);
        const core::Real zCorr = params.rho * z1 + rhoComp * z2;
        
        const core::Real currentVar = std::max(varPath[step - 1], 0.0);
        const core::Real sqrtVar = std::sqrt(currentVar);
        
        varPath[step] = currentVar + params.kappa * (params.theta - currentVar) * dt + 
                       params.sigma * sqrtVar * sqrtDt * z1;
        varPath[step] = std::max(varPath[step], 0.0);
        
        const core::Real drift = -0.5 * currentVar * dt;
        spotPath[step] = spotPath[step - 1] * std::exp(drift + sqrtVar * sqrtDt * zCorr);
    }
    
    return spotPath;
}

StochasticVolatilityModel::StochasticVolatilityModel(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real StochasticVolatilityModel::sabr(core::Real forward, core::Real strike,
                                          core::Real timeToExpiry, const SABRParams& params) const {
    if (std::abs(params.alpha) < core::EPSILON) return 0.0;
    
    const core::Real logMoneyness = std::log(forward / strike);
    const core::Real fk = forward * strike;
    const core::Real fkAvg = std::sqrt(fk);
    
    core::Real vol = params.alpha;
    
    if (std::abs(logMoneyness) > core::EPSILON) {
        const core::Real z = params.nu / params.alpha * std::pow(fkAvg, 1.0 - params.beta) * logMoneyness;
        const core::Real x = std::log((std::sqrt(1.0 - 2.0 * params.rho * z + z * z) + z - params.rho) / (1.0 - params.rho));
        
        vol *= logMoneyness / x;
    }
    
    vol *= std::pow(fkAvg, params.beta - 1.0);
    
    const core::Real correction = 1.0 + timeToExpiry * (
        std::pow(1.0 - params.beta, 2) * params.alpha * params.alpha / (24.0 * std::pow(fkAvg, 2.0 - 2.0 * params.beta)) +
        0.25 * params.rho * params.beta * params.nu * params.alpha / std::pow(fkAvg, 1.0 - params.beta) +
        (2.0 - 3.0 * params.rho * params.rho) * params.nu * params.nu / 24.0
    );
    
    return vol * correction;
}

StochasticVolatilityModel::SABRParams StochasticVolatilityModel::calibrateSABR(
    const core::Vector& strikes, const core::Vector& impliedVols,
    core::Real forward, core::Real timeToExpiry) const {
    
    SABRParams params;
    params.alpha = 0.2;
    params.beta = 0.5;
    params.rho = 0.0;
    params.nu = 0.3;
    
    const core::Integer maxIterations = 100;
    const core::Real tolerance = 1e-6;
    
    for (core::Integer iter = 0; iter < maxIterations; ++iter) {
        core::Real totalError = 0.0;
        
        for (core::Size i = 0; i < strikes.size(); ++i) {
            const core::Real modelVol = sabr(forward, strikes[i], timeToExpiry, params);
            const core::Real error = modelVol - impliedVols[i];
            totalError += error * error;
        }
        
        if (std::sqrt(totalError / strikes.size()) < tolerance) {
            break;
        }
        
        const core::Real stepSize = 0.001;
        params.alpha += stepSize * (totalError > 0 ? -1.0 : 1.0);
        params.nu += stepSize * (totalError > 0 ? -1.0 : 1.0);
        
        params.alpha = std::clamp(params.alpha, 0.001, 2.0);
        params.nu = std::clamp(params.nu, 0.001, 2.0);
        params.rho = std::clamp(params.rho, -0.99, 0.99);
    }
    
    return params;
}

JumpDiffusionModel::JumpDiffusionModel(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real JumpDiffusionModel::mertonPrice(core::Real spot, core::Real strike,
                                          core::Real timeToExpiry, core::Real riskFreeRate,
                                          core::Real sigma, const JumpParams& jumpParams,
                                          bool isCall) const {
    core::Real price = 0.0;
    const core::Real lambda = jumpParams.intensity;
    const core::Real mu = jumpParams.meanJumpSize;
    const core::Real delta = jumpParams.jumpVolatility;
    
    const core::Integer maxTerms = 50;
    core::Real factorial = 1.0;
    
    for (core::Integer n = 0; n < maxTerms; ++n) {
        if (n > 0) factorial *= n;
        
        const core::Real poissonProb = std::exp(-lambda * timeToExpiry) * 
                                      std::pow(lambda * timeToExpiry, n) / factorial;
        
        const core::Real sigmaN = std::sqrt(sigma * sigma + n * delta * delta / timeToExpiry);
        const core::Real rN = riskFreeRate - lambda * (std::exp(mu + 0.5 * delta * delta) - 1.0) + 
                              n * std::log(1.0 + mu) / timeToExpiry;
        
        core::MarketData market;
        market.spot = spot;
        market.strike = strike;
        market.timeToExpiry = timeToExpiry;
        market.riskFreeRate = rN;
        market.volatility = sigmaN;
        market.isCall = isCall;
        
        const core::Real bsPrice = engine_->blackScholes(market).price;
        price += poissonProb * bsPrice;
    }
    
    return price;
}

LocalVolatilityModel::LocalVolatilityModel(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real LocalVolatilityModel::dupireLocalVol(core::Real spot, core::Real timeToExpiry,
                                               const std::function<core::Real(core::Real, core::Real)>& impliedVol) const {
    const core::Real h = 0.01;
    const core::Real dt = 1.0 / 365.0;
    
    const core::Real vol = impliedVol(spot, timeToExpiry);
    const core::Real vol_dt = impliedVol(spot, timeToExpiry + dt);
    const core::Real vol_dK_up = impliedVol(spot + h, timeToExpiry);
    const core::Real vol_dK_down = impliedVol(spot - h, timeToExpiry);
    const core::Real vol_d2K = impliedVol(spot + 2.0 * h, timeToExpiry);
    
    const core::Real dVol_dt = (vol_dt - vol) / dt;
    const core::Real dVol_dK = (vol_dK_up - vol_dK_down) / (2.0 * h);
    const core::Real d2Vol_dK2 = (vol_dK_up - 2.0 * vol + vol_dK_down) / (h * h);
    
    const core::Real d1 = std::log(spot / spot) / (vol * std::sqrt(timeToExpiry));
    
    const core::Real numerator = vol * vol + 2.0 * vol * timeToExpiry * dVol_dt;
    const core::Real denominator = 1.0 + d1 * std::sqrt(timeToExpiry) * dVol_dK + 
                                  0.25 * timeToExpiry * vol * (d1 * d1 - 1.0) * dVol_dK * dVol_dK +
                                  timeToExpiry * vol * d2Vol_dK2;
    
    return std::sqrt(std::max(0.01, numerator / denominator));
}

VarianceGammaModel::VarianceGammaModel(std::shared_ptr<core::MathEngine> engine)
    : engine_(std::move(engine)) {}

core::Real VarianceGammaModel::characteristicFunction(const std::complex<core::Real>& u,
                                                     core::Real timeToExpiry,
                                                     const VGParams& params) const {
    const auto i = std::complex<core::Real>(0.0, 1.0);
    
    const auto omega = std::log(1.0 - params.theta * params.nu - 0.5 * params.sigma * params.sigma * params.nu) / params.nu;
    
    const auto phi = std::exp(i * u * omega * timeToExpiry) * 
                     std::pow(1.0 - i * u * params.theta * params.nu + 0.5 * params.sigma * params.sigma * u * u * params.nu,
                             -timeToExpiry / params.nu);
    
    return phi;
}

core::Real VarianceGammaModel::optionPrice(core::Real spot, core::Real strike,
                                          core::Real timeToExpiry, core::Real riskFreeRate,
                                          const VGParams& params, bool isCall) const {
    const core::Real forward = spot * std::exp(riskFreeRate * timeToExpiry);
    const core::Real logMoneyness = std::log(forward / strike);
    
    auto integrand = [&](core::Real u) -> core::Real {
        const auto phi = characteristicFunction(std::complex<core::Real>(u, -0.5), timeToExpiry, params);
        return std::real(std::exp(-std::complex<core::Real>(0.0, 1.0) * u * logMoneyness) * phi);
    };
    
    const core::Real integral = core::NumericalMethods::trapezoidalRule(integrand, -50.0, 50.0);
    const core::Real price = std::exp(-riskFreeRate * timeToExpiry) * integral / core::PI;
    
    return isCall ? price : price - forward + strike * std::exp(-riskFreeRate * timeToExpiry);
}

}