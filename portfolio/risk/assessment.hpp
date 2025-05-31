/*
 * Filename: assessment.hpp
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

#include "common/order_request.hpp"
#include <vector>
#include <string>
#include <map>
#include <chrono>

namespace portfolio {

class RiskAssessment {
private:
    bool withinLimits_;
    std::vector<std::string> violations_;
    double riskScore_;
    std::map<std::string, double> recommendations_;
    std::vector<common::OrderRequest> hedgeOrders_;
    double capitalEfficiency_;
    double portfolioHeat_;
    double stressTestScore_;
    std::chrono::system_clock::time_point timestamp_;

public:
    RiskAssessment();
    ~RiskAssessment() = default;
    
    bool isWithinLimits() const { return withinLimits_; }
    const std::vector<std::string>& getViolations() const { return violations_; }
    double getRiskScore() const { return riskScore_; }
    const std::map<std::string, double>& getRecommendations() const { return recommendations_; }
    const std::vector<common::OrderRequest>& getHedgeOrders() const { return hedgeOrders_; }
    double getCapitalEfficiency() const { return capitalEfficiency_; }
    double getPortfolioHeat() const { return portfolioHeat_; }
    double getStressTestScore() const { return stressTestScore_; }
    const std::chrono::system_clock::time_point& getTimestamp() const { return timestamp_; }
    
    void setWithinLimits(bool withinLimits) { withinLimits_ = withinLimits; }
    void setViolations(const std::vector<std::string>& violations) { violations_ = violations; }
    void setRiskScore(double riskScore) { riskScore_ = riskScore; }
    void setRecommendations(const std::map<std::string, double>& recommendations) { recommendations_ = recommendations; }
    void setHedgeOrders(const std::vector<common::OrderRequest>& hedgeOrders) { hedgeOrders_ = hedgeOrders; }
    void setCapitalEfficiency(double capitalEfficiency) { capitalEfficiency_ = capitalEfficiency; }
    void setPortfolioHeat(double portfolioHeat) { portfolioHeat_ = portfolioHeat; }
    void setStressTestScore(double stressTestScore) { stressTestScore_ = stressTestScore; }
    void setTimestamp(const std::chrono::system_clock::time_point& timestamp) { timestamp_ = timestamp; }
    
    void addViolation(const std::string& violation);
    void addRecommendation(const std::string& metric, double value);
    void addHedgeOrder(const common::OrderRequest& order);
    
    bool hasViolations() const { return !violations_.empty(); }
    bool isHighRisk() const { return riskScore_ > 0.8; }
    bool needsHedging() const { return !hedgeOrders_.empty(); }
    
    void updateTimestamp();
    void clear();
    
    std::string getSummary() const;
    std::string getDetailedReport() const;
};

}