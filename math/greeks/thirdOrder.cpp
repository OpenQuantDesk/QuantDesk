core::Real ThirdOrderCalculator::calculateDvannaDtime(const core::MarketData& market) const {
    const core::Real timeBump = 1.0 / 365.0;
    const core::Real spotBump = market.spot * defaultBumpSize_;
    
    auto vannaFunc = [this, &market, spotBump](core::Real time) {
        core::MarketData tempMarket = market;
        tempMarket.timeToExpiry = time;
        
        if (tempMarket.timeToExpiry <= 0) return 0.0;
        
        // Calculate vanna numerically
        core::MarketData upMarket = tempMarket;
        upMarket.spot += spotBump;
        const core::Real deltaUp = engine_->blackScholes(upMarket).delta;
        
        core::MarketData downMarket = tempMarket;
        downMarket.spot -= spotBump;
        const core::Real deltaDown = engine_->blackScholes(downMarket).delta;
        
        return (deltaUp - deltaDown) / (2.0 * spotBump);
    };
    
    if (market.timeToExpiry <= timeBump) return 0.0;
    
    const core::Real currentVanna = vannaFunc(market.timeToExpiry);
    const core::Real futureVanna = vannaFunc(market.timeToExpiry - timeBump);
    
    return (currentVanna - futureVanna) / timeBump;
}

core::Real ThirdOrderCalculator::calculateDcharmDvol(const core::MarketData& market) const {
    const core::Real volBump = 0.01;
    
    auto charmFunc = [this, &market](core::Real vol) {
        core::MarketData tempMarket = market;
        tempMarket.volatility = vol;
        
        const core::Real timeBump = 1.0 / 365.0;
        if (tempMarket.timeToExpiry <= timeBump) return 0.0;
        
        const core::Real currentDelta = engine_->blackScholes(tempMarket).delta;
        
        core::MarketData timeMarket = tempMarket;
        timeMarket.timeToExpiry -= timeBump;
        const core::Real futureDelta = engine_->blackScholes(timeMarket).delta;
        
        return (currentDelta - futureDelta) / timeBump;
    };
    
    const core::Real charmUp = charmFunc(market.volatility + volBump);
    const core::Real charmDown = charmFunc(market.volatility - volBump);
    
    return (charmUp - charmDown) / (2.0 * volBump);
}

ThirdOrderCalculator::ThirdOrderRiskMetrics ThirdOrderCalculator::analyzePortfolioRisk(
    const std::vector<ThirdOrderGreeks>& positions,
    const core::Vector& quantities) const {
    ThirdOrderRiskMetrics metrics;
    
    // Aggregate portfolio-level risks
    core::Real totalConvexity = 0.0;
    core::Real totalVolOfVol = 0.0;
    core::Real totalTimeDecay = 0.0;
    core::Real totalCrossRisk = 0.0;
    
    metrics.riskContributions.resize(positions.size());
    
    for (core::Size i = 0; i < positions.size(); ++i) {
        const core::Real qty = i < quantities.size() ? quantities[i] : 1.0;
        const ThirdOrderGreeks& greeks = positions[i];
        
        // Convexity risk (speed, zomma, color)
        totalConvexity += qty * (std::abs(greeks.speed) + 
                               std::abs(greeks.zomma) + 
                               std::abs(greeks.color));
        
        // Volatility of volatility risk (ultima, volga)
        totalVolOfVol += qty * (std::abs(greeks.ultima) + 
                               std::abs(greeks.volga));
        
        // Time decay acceleration
        totalTimeDecay += qty * (std::abs(greeks.dspeedDtime) + 
                                std::abs(greeks.dzommaDtime) + 
                                std::abs(greeks.dvannaDtime));
        
        // Cross-derivative risk
        totalCrossRisk += qty * (std::abs(greeks.dspeedDvol) + 
                                std::abs(greeks.dcolorDvol) + 
                                std::abs(greeks.dcharmDvol));
        
        // Per-position contribution
        metrics.riskContributions[i] = qty * (
            std::abs(greeks.speed) + 
            std::abs(greeks.zomma) + 
            std::abs(greeks.ultima) + 
            std::abs(greeks.totto)
        );
    }
    
    metrics.convexityRisk = totalConvexity;
    metrics.volOfVolRisk = totalVolOfVol;
    metrics.timeDecayAcceleration = totalTimeDecay;
    metrics.crossRisk = totalCrossRisk;
    
    // Estimate tail risk as weighted sum of third-order exposures
    metrics.tailRisk = 0.3 * totalConvexity + 
                      0.4 * totalVolOfVol + 
                      0.2 * totalTimeDecay + 
                      0.1 * totalCrossRisk;
    
    return metrics;
}

ThirdOrderCalculator::ConvexityProfile ThirdOrderCalculator::buildConvexityProfile(
    const core::MarketData& market,
    const core::Vector& spotRange,
    const core::Vector& volRange,
    const core::Vector& timeRange) const {
    ConvexityProfile profile;
    
    profile.spotLadder = spotRange;
    profile.volLadder = volRange;
    profile.timeLadder = timeRange;
    
    const core::Size spotSize = spotRange.size();
    const core::Size volSize = volRange.size();
    const core::Size timeSize = timeRange.size();
    
    // Initialize matrices
    profile.speedProfile.resize(spotSize, volSize);
    profile.zommaProfile.resize(spotSize, volSize);
    profile.colorProfile.resize(spotSize, timeSize);
    profile.ultimaProfile.resize(volSize, timeSize);
    
    // Calculate profiles
    core::MarketData tempMarket = market;
    
    // Speed and Zomma profiles (spot vs vol)
    for (core::Size i = 0; i < spotSize; ++i) {
        tempMarket.spot = spotRange[i];
        for (core::Size j = 0; j < volSize; ++j) {
            tempMarket.volatility = volRange[j];
            profile.speedProfile(i, j) = calculateSpeed(tempMarket);
            profile.zommaProfile(i, j) = calculateZomma(tempMarket);
        }
    }
    
    // Color profile (spot vs time)
    for (core::Size i = 0; i < spotSize; ++i) {
        tempMarket.spot = spotRange[i];
        for (core::Size j = 0; j < timeSize; ++j) {
            tempMarket.timeToExpiry = timeRange[j];
            profile.colorProfile(i, j) = calculateColor(tempMarket);
        }
    }
    
    // Ultima profile (vol vs time)
    for (core::Size i = 0; i < volSize; ++i) {
        tempMarket.volatility = volRange[i];
        for (core::Size j = 0; j < timeSize; ++j) {
            tempMarket.timeToExpiry = timeRange[j];
            profile.ultimaProfile(i, j) = calculateUltima(tempMarket);
        }
    }
    
    return profile;
}

ThirdOrderCalculator::HedgingRecommendation ThirdOrderCalculator::recommendThirdOrderHedge(
    const std::vector<ThirdOrderGreeks>& portfolio,
    const std::vector<ThirdOrderGreeks>& hedgeInstruments,
    const core::Vector& hedgeCosts) const {
    HedgingRecommendation rec;
    
    const core::Size n = portfolio.size();
    const core::Size m = hedgeInstruments.size();
    
    rec.thirdOrderExposures.resize(n);
    rec.hedgeInstruments.resize(m);
    rec.hedgeQuantities.resize(m);
    rec.hedgeEffectiveness.resize(n, m);
    
    // Calculate portfolio exposures
    for (core::Size i = 0; i < n; ++i) {
        const ThirdOrderGreeks& greeks = portfolio[i];
        rec.thirdOrderExposures[i] = std::abs(greeks.speed) + 
                                    std::abs(greeks.zomma) + 
                                    std::abs(greeks.ultima) + 
                                    std::abs(greeks.totto);
    }
    
    // Simple optimization: minimize exposure and cost
    for (core::Size j = 0; j < m; ++j) {
        const ThirdOrderGreeks& instr = hedgeInstruments[j];
        core::Real hedgeEffect = std::abs(instr.speed) + 
                                std::abs(instr.zomma) + 
                                std::abs(instr.ultima) + 
                                std::abs(instr.totto);
        
        // Calculate optimal quantity (simplified: inverse proportionality to cost)
        core::Real cost = j < hedgeCosts.size() ? hedgeCosts[j] : 1.0;
        rec.hedgeQuantities[j] = hedgeEffect / (cost + 1e-10);
        
        // Update effectiveness matrix
        for (core::Size i = 0; i < n; ++i) {
            rec.hedgeEffectiveness(i, j) = std::min(
                rec.thirdOrderExposures[i] / (hedgeEffect + 1e-10), 1.0);
        }
        
        rec.hedgeInstruments[j] = j;
    }
    
    // Calculate residual risk and total hedging cost
    rec.residualRisk = 0.0;
    for (core::Size i = 0; i < n; ++i) {
        core::Real hedgedExposure = rec.thirdOrderExposures[i];
        for (core::Size j = 0; j < m; ++j) {
            hedgedExposure -= rec.hedgeEffectiveness(i, j) * rec.hedgeQuantities[j];
        }
        rec.residualRisk += std::max(0.0, hedgedExposure);
    }
    
    rec.hedgingCost = 0.0;
    for (core::Size j = 0; j < m; ++j) {
        rec.hedgingCost += rec.hedgeQuantities[j] * 
                          (j < hedgeCosts.size() ? hedgeCosts[j] : 1.0);
    }
    
    return rec;
}

core::Real ThirdOrderCalculator::analyticalSpeed(const core::MarketData& market) const {
    const core::Real d1 = calculateD1(market);
    const core::Real sqrtT = std::sqrt(market.timeToExpiry);
    const core::Real sigma = market.volatility;
    const core::Real S = market.spot;
    
    const core::Real pdf = normalPDF(d1);
    const core::Real term1 = -pdf * (d1 * d1 + d1 * sigma * sqrtT + 1.0);
    const core::Real term2 = S * S * S * sigma * sigma * sqrtT;
    
    return term1 / term2;
}

core::Real ThirdOrderCalculator::analyticalZomma(const core::MarketData& market) const {
    const core::Real d1 = calculateD1(market);
    const core::Real d2 = calculateD2(market);
    const core::Real S = market.spot;
    const core::Real sigma = market.volatility;
    const core::Real sqrtT = std::sqrt(market.timeToExpiry);
    
    const core::Real pdf = normalPDF(d1);
    const core::Real gamma = pdf / (S * sigma * sqrtT);
    
    return gamma * (d1 * d2 - 1.0) / sigma;
}

core::Real ThirdOrderCalculator::analyticalColor(const core::MarketData& market) const {
    const core::Real d1 = calculateD1(market);
    const core::Real d2 = calculateD2(market);
    const core::Real S = market.spot;
    const core::Real sigma = market.volatility;
    const core::Real sqrtT = std::sqrt(market.timeToExpiry);
    const core::Real r = market.riskFreeRate;
    
    const core::Real pdf = normalPDF(d1);
    const core::Real gamma = pdf / (S * sigma * sqrtT);
    
    const core::Real term1 = -gamma * (0.5 + (r * d1 - d2 * sigma * sqrtT) / 
                                      (2.0 * sigma * sqrtT));
    
    return term1 / 365.0; // Convert to daily decay
}

core::Real ThirdOrderCalculator::analyticalUltima(const core::MarketData& market) const {
    const core::Real d1 = calculateD1(market);
    const core::Real d2 = calculateD2(market);
    const core::Real S = market.spot;
    const core::Real sigma = market.volatility;
    const core::Real T = market.timeToExpiry;
    
    const core::Real pdf = normalPDF(d1);
    const core::Real vega = S * pdf * std::sqrt(T);
    
    const core::Real term1 = -vega / (sigma * sigma);
    const core::Real term2 = d1 * d2 * (1.0 - d1 * d2) + d1 * d1 + d2 * d2;
    
    return term1 * term2;
}

core::Real ThirdOrderCalculator::numericalThirdDerivative(
    const std::function<core::Real(core::Real)>& func,
    core::Real x, core::Real h) const {
    const core::Real f1 = func(x + 2.0 * h);
    const core::Real f2 = func(x + h);
    const core::Real f3 = func(x - h);
    const core::Real f4 = func(x - 2.0 * h);
    
    return (f1 - 3.0 * f2 + 3.0 * f3 - f4) / (2.0 * h * h * h);
}

core::Real ThirdOrderCalculator::numericalMixedDerivative(
    const std::function<core::Real(core::Real, core::Real)>& func,
    core::Real x, core::Real y, 
    core::Real hx, core::Real hy,
    bool isThirdOrder) const {
    if (isThirdOrder) {
        // Third-order mixed derivative
        const core::Real f1 = func(x + hx, y + hy);
        const core::Real f2 = func(x + hx, y - hy);
        const core::Real f3 = func(x - hx, y + hy);
        const core::Real f4 = func(x - hx, y - hy);
        
        return (f1 - f2 - f3 + f4) / (4.0 * hx * hy * hy);
    } else {
        // Second-order mixed derivative
        const core::Real f1 = func(x + hx, y + hy);
        const core::Real f2 = func(x + hx, y - hy);
        const core::Real f3 = func(x - hx, y + hy);
        const core::Real f4 = func(x - hx, y - hy);
        
        return (f1 - f2 - f3 + f4) / (4.0 * hx * hy);
    }
}

core::Real ThirdOrderCalculator::calculatePortfolioConvexity(
    const std::vector<ThirdOrderGreeks>& positions,
    const core::Vector& quantities) const {
    core::Real totalConvexity = 0.0;
    
    for (core::Size i = 0; i < positions.size(); ++i) {
        const core::Real qty = i < quantities.size() ? quantities[i] : 1.0;
        totalConvexity += qty * (std::abs(positions[i].speed) + 
                               std::abs(positions[i].zomma) + 
                               std::abs(positions[i].color));
    }
    
    return totalConvexity;
}

core::Vector ThirdOrderCalculator::calculateRiskContributions(
    const std::vector<ThirdOrderGreeks>& positions,
    const core::Vector& quantities) const {
    core::Vector contributions(positions.size());
    
    for (core::Size i = 0; i < positions.size(); ++i) {
        const core::Real qty = i < quantities.size() ? quantities[i] : 1.0;
        contributions[i] = qty * (
            std::abs(positions[i].speed) + 
            std::abs(positions[i].zomma) + 
            std::abs(positions[i].ultima) + 
            std::abs(positions[i].totto)
        );
    }
    
    return contributions;
}

} // namespace math::greeks