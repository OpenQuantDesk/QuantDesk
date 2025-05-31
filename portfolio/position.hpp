/*
 * Filename: position.hpp
 * Developer: Benjamin Cance
 * Date: 5/31/2025
 * 
 * Copyright 2025 Open Quant Desk, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include "common/greeks.hpp"
#include "OptionsQuantLib/analytics/probability.hpp"
#include <string>
#include <chrono>

namespace portfolio {

class EnhancedPosition {
private:
    std::string symbol_;
    std::string optionSymbol_;
    std::string underlying_;
    double strike_;
    std::string expiration_;
    bool isCall_;
    double quantity_;
    double avgPrice_;
    double currentPrice_;
    double underlyingPrice_;
    
    common::Greeks greeks_;
    common::ExtendedGreeks extendedGreeks_;
    math::core::ProbabilityMetrics probMetrics_;
    
    double theoreticalValue_;
    double impliedVol_;
    double unrealizedPnL_;
    double dayChange_;
    double deltaAdjustedExposure_;
    double gammaRisk_;
    double vegaRisk_;
    double thetaDecay_;
    
    std::chrono::system_clock::time_point entryTime_;
    std::chrono::system_clock::time_point lastUpdate_;

public:
    EnhancedPosition();
    EnhancedPosition(const std::string& symbol, const std::string& underlying);
    ~EnhancedPosition() = default;
    
    const std::string& getSymbol() const { return symbol_; }
    const std::string& getOptionSymbol() const { return optionSymbol_; }
    const std::string& getUnderlying() const { return underlying_; }
    double getStrike() const { return strike_; }
    const std::string& getExpiration() const { return expiration_; }
    bool getIsCall() const { return isCall_; }
    double getQuantity() const { return quantity_; }
    double getAvgPrice() const { return avgPrice_; }
    double getCurrentPrice() const { return currentPrice_; }
    double getUnderlyingPrice() const { return underlyingPrice_; }
    const common::Greeks& getGreeks() const { return greeks_; }
    const common::ExtendedGreeks& getExtendedGreeks() const { return extendedGreeks_; }
    const math::core::ProbabilityMetrics& getProbMetrics() const { return probMetrics_; }
    double getTheoreticalValue() const { return theoreticalValue_; }
    double getImpliedVol() const { return impliedVol_; }
    double getUnrealizedPnL() const { return unrealizedPnL_; }
    double getDayChange() const { return dayChange_; }
    double getDeltaAdjustedExposure() const { return deltaAdjustedExposure_; }
    double getGammaRisk() const { return gammaRisk_; }
    double getVegaRisk() const { return vegaRisk_; }
    double getThetaDecay() const { return thetaDecay_; }
    const std::chrono::system_clock::time_point& getEntryTime() const { return entryTime_; }
    const std::chrono::system_clock::time_point& getLastUpdate() const { return lastUpdate_; }
    
    void setSymbol(const std::string& symbol) { symbol_ = symbol; }
    void setOptionSymbol(const std::string& optionSymbol) { optionSymbol_ = optionSymbol; }
    void setUnderlying(const std::string& underlying) { underlying_ = underlying; }
    void setStrike(double strike) { strike_ = strike; }
    void setExpiration(const std::string& expiration) { expiration_ = expiration; }
    void setIsCall(bool isCall) { isCall_ = isCall; }
    void setQuantity(double quantity) { quantity_ = quantity; }
    void setAvgPrice(double avgPrice) { avgPrice_ = avgPrice; }
    void setCurrentPrice(double currentPrice) { currentPrice_ = currentPrice; }
    void setUnderlyingPrice(double underlyingPrice) { underlyingPrice_ = underlyingPrice; }
    void setGreeks(const common::Greeks& greeks) { greeks_ = greeks; }
    void setExtendedGreeks(const common::ExtendedGreeks& extendedGreeks) { extendedGreeks_ = extendedGreeks; }
    void setProbMetrics(const math::core::ProbabilityMetrics& probMetrics) { probMetrics_ = probMetrics; }
    
    void updateAnalytics(double underlyingPx, double vol, double riskFreeRate,
                        const class math::core::MathEngine& engine);
    
    double getDaysToExpiry() const;
    double getNotionalValue() const;
    double getImpliedProbability() const;
    double getMoneyness() const;
    bool isExpiring(int days = 7) const;
    bool isOption() const;
    bool isEquity() const;
    
private:
    std::chrono::system_clock::time_point parseExpirationDate(const std::string& expiration) const;
    void updatePnL();
    void updateRiskMetrics();
};

}