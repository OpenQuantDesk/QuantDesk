/*
 * Filename: metrics.hpp
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

#include <map>
#include <string>
#include <chrono>

namespace portfolio {

class PortfolioMetrics {
private:
    double totalValue_;
    double totalPnL_;
    double dayPnL_;
    double cash_;
    
    double netDelta_;
    double netGamma_;
    double netTheta_;
    double netVega_;
    double netRho_;
    
    double netVolga_;
    double netVanna_;
    double dollarDelta_;
    double dollarGamma_;
    
    double var95_;
    double var99_;
    double expectedShortfall_;
    double maxDrawdown_;
    double leverageRatio_;
    
    double portfolioProbProfit_;
    double sharpeRatio_;
    double kellyOptimal_;
    double correlationRisk_;
    
    std::map<std::string, double> underlyingExposures_;
    std::map<std::string, double> sectorExposures_;
    std::map<std::string, double> strategyAllocations_;
    
    std::chrono::system_clock::time_point timestamp_;

public:
    PortfolioMetrics();
    ~PortfolioMetrics() = default;
    
    double getTotalValue() const { return totalValue_; }
    double getTotalPnL() const { return totalPnL_; }
    double getDayPnL() const { return dayPnL_; }
    double getCash() const { return cash_; }
    double getNetDelta() const { return netDelta_; }
    double getNetGamma() const { return netGamma_; }
    double getNetTheta() const { return netTheta_; }
    double getNetVega() const { return netVega_; }
    double getNetRho() const { return netRho_; }
    double getNetVolga() const { return netVolga_; }
    double getNetVanna() const { return netVanna_; }
    double getDollarDelta() const { return dollarDelta_; }
    double getDollarGamma() const { return dollarGamma_; }
    double getVar95() const { return var95_; }
    double getVar99() const { return var99_; }
    double getExpectedShortfall() const { return expectedShortfall_; }
    double getMaxDrawdown() const { return maxDrawdown_; }
    double getLeverageRatio() const { return leverageRatio_; }
    double getPortfolioProbProfit() const { return portfolioProbProfit_; }
    double getSharpeRatio() const { return sharpeRatio_; }
    double getKellyOptimal() const { return kellyOptimal_; }
    double getCorrelationRisk() const { return correlationRisk_; }
    const std::map<std::string, double>& getUnderlyingExposures() const { return underlyingExposures_; }
    const std::map<std::string, double>& getSectorExposures() const { return sectorExposures_; }
    const std::map<std::string, double>& getStrategyAllocations() const { return strategyAllocations_; }
    const std::chrono::system_clock::time_point& getTimestamp() const { return timestamp_; }
    
    void setTotalValue(double totalValue) { totalValue_ = totalValue; }
    void setTotalPnL(double totalPnL) { totalPnL_ = totalPnL; }
    void setDayPnL(double dayPnL) { dayPnL_ = dayPnL; }
    void setCash(double cash) { cash_ = cash; }
    void setNetDelta(double netDelta) { netDelta_ = netDelta; }
    void setNetGamma(double netGamma) { netGamma_ = netGamma; }
    void setNetTheta(double netTheta) { netTheta_ = netTheta; }
    void setNetVega(double netVega) { netVega_ = netVega; }
    void setNetRho(double netRho) { netRho_ = netRho; }
    void setNetVolga(double netVolga) { netVolga_ = netVolga; }
    void setNetVanna(double netVanna) { netVanna_ = netVanna; }
    void setDollarDelta(double dollarDelta) { dollarDelta_ = dollarDelta; }
    void setDollarGamma(double dollarGamma) { dollarGamma_ = dollarGamma; }
    void setVar95(double var95) { var95_ = var95; }
    void setVar99(double var99) { var99_ = var99; }
    void setExpectedShortfall(double expectedShortfall) { expectedShortfall_ = expectedShortfall; }
    void setMaxDrawdown(double maxDrawdown) { maxDrawdown_ = maxDrawdown; }
    void setLeverageRatio(double leverageRatio) { leverageRatio_ = leverageRatio; }
    void setPortfolioProbProfit(double portfolioProbProfit) { portfolioProbProfit_ = portfolioProbProfit; }
    void setSharpeRatio(double sharpeRatio) { sharpeRatio_ = sharpeRatio; }
    void setKellyOptimal(double kellyOptimal) { kellyOptimal_ = kellyOptimal; }
    void setCorrelationRisk(double correlationRisk) { correlationRisk_ = correlationRisk; }
    void setUnderlyingExposures(const std::map<std::string, double>& underlyingExposures) { underlyingExposures_ = underlyingExposures; }
    void setSectorExposures(const std::map<std::string, double>& sectorExposures) { sectorExposures_ = sectorExposures; }
    void setStrategyAllocations(const std::map<std::string, double>& strategyAllocations) { strategyAllocations_ = strategyAllocations; }
    void setTimestamp(const std::chrono::system_clock::time_point& timestamp) { timestamp_ = timestamp; }
    
    void addUnderlyingExposure(const std::string& underlying, double exposure);
    void addSectorExposure(const std::string& sector, double exposure);
    void addStrategyAllocation(const std::string& strategy, double allocation);
    
    double getTotalExposure() const;
    double getMaxSingleExposure() const;
    double getConcentrationRatio() const;
    bool isWithinRiskLimits() const;
    
    void updateTimestamp();
    void reset();
};

}