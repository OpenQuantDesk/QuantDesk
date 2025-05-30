#pragma once

#include "portfolio/manager.hpp"
#include "limits.hpp"
#include "assessment.hpp"
#include "math/engine.hpp"
#include "common/order_request.hpp"
#include <memory>
#include <vector>
#include <functional>
#include <future>
#include <thread>
#include <atomic>
#include <shared_mutex>

namespace portfolio {

class RiskManager {
private:
    std::shared_ptr<PortfolioManager> portfolio_;
    std::shared_ptr<math::MathEngine> mathEngine_;
    RiskLimits limits_;
    
    std::vector<RiskAssessment> assessmentHistory_;
    mutable std::shared_mutex assessmentLock_;
    
    std::thread monitoringThread_;
    std::atomic<bool> monitoring_;
    
    std::vector<std::function<void(const RiskAssessment&)>> alertCallbacks_;

public:
    RiskManager(std::shared_ptr<PortfolioManager> portfolio,
               std::shared_ptr<math::MathEngine> engine);
    ~RiskManager();
    
    RiskManager(const RiskManager&) = delete;
    RiskManager& operator=(const RiskManager&) = delete;
    
    void startMonitoring();
    void stopMonitoring();
    
    std::future<RiskAssessment> assessRisk() const;
    std::future<bool> validateOrder(const common::OrderRequest& order) const;
    std::future<std::vector<common::OrderRequest>> generateHedgeOrders() const;
    
    void setRiskLimits(const RiskLimits& limits);
    const RiskLimits& getRiskLimits() const { return limits_; }
    
    std::future<double> calculatePortfolioVaR(double confidence = 0.95, int horizon = 1) const;
    std::future<double> calculateExpectedShortfall(double confidence = 0.95) const;
    std::future<double> calculateMaxDrawdown() const;
    
    std::future<std::map<std::string, double>> calculateComponentVaR() const;
    std::future<std::map<std::string, double>> calculateMarginalVaR() const;
    
    void addAlertCallback(std::function<void(const RiskAssessment&)> callback);
    std::vector<RiskAssessment> getAssessmentHistory(int days = 7) const;
    
    bool isMonitoring() const { return monitoring_.load(); }
    
private:
    void monitoringLoop();
    double calculateOrderImpact(const common::OrderRequest& order) const;
    double calculateRiskScore(const PortfolioMetrics& metrics) const;
    std::vector<common::OrderRequest> generateDeltaHedge(double targetDelta = 0.0) const;
    std::vector<common::OrderRequest> generateGammaHedge() const;
    std::vector<common::OrderRequest> generateVegaHedge() const;
    
    void notifyRiskAlert(const RiskAssessment& assessment);
    bool isViolation(const std::string& metric, double value, double limit) const;
    void checkPositionLimits(RiskAssessment& assessment) const;
    void checkPortfolioLimits(RiskAssessment& assessment) const;
    void checkConcentrationLimits(RiskAssessment& assessment) const;
};

}