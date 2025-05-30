#pragma once

namespace portfolio {

class RiskLimits {
private:
    double maxDelta_;
    double maxGamma_;
    double maxTheta_;
    double maxVega_;
    double maxPositionSize_;
    double maxPortfolioRisk_;
    double maxLeverage_;
    double maxConcentration_;
    double maxVaR_;
    double maxDrawdown_;
    double maxDailyLoss_;
    double maxWeeklyLoss_;
    double maxMonthlyLoss_;

public:
    RiskLimits();
    ~RiskLimits() = default;
    
    double getMaxDelta() const { return maxDelta_; }
    double getMaxGamma() const { return maxGamma_; }
    double getMaxTheta() const { return maxTheta_; }
    double getMaxVega() const { return maxVega_; }
    double getMaxPositionSize() const { return maxPositionSize_; }
    double getMaxPortfolioRisk() const { return maxPortfolioRisk_; }
    double getMaxLeverage() const { return maxLeverage_; }
    double getMaxConcentration() const { return maxConcentration_; }
    double getMaxVaR() const { return maxVaR_; }
    double getMaxDrawdown() const { return maxDrawdown_; }
    double getMaxDailyLoss() const { return maxDailyLoss_; }
    double getMaxWeeklyLoss() const { return maxWeeklyLoss_; }
    double getMaxMonthlyLoss() const { return maxMonthlyLoss_; }
    
    void setMaxDelta(double maxDelta) { maxDelta_ = maxDelta; }
    void setMaxGamma(double maxGamma) { maxGamma_ = maxGamma; }
    void setMaxTheta(double maxTheta) { maxTheta_ = maxTheta; }
    void setMaxVega(double maxVega) { maxVega_ = maxVega; }
    void setMaxPositionSize(double maxPositionSize) { maxPositionSize_ = maxPositionSize; }
    void setMaxPortfolioRisk(double maxPortfolioRisk) { maxPortfolioRisk_ = maxPortfolioRisk; }
    void setMaxLeverage(double maxLeverage) { maxLeverage_ = maxLeverage; }
    void setMaxConcentration(double maxConcentration) { maxConcentration_ = maxConcentration; }
    void setMaxVaR(double maxVaR) { maxVaR_ = maxVaR; }
    void setMaxDrawdown(double maxDrawdown) { maxDrawdown_ = maxDrawdown; }
    void setMaxDailyLoss(double maxDailyLoss) { maxDailyLoss_ = maxDailyLoss; }
    void setMaxWeeklyLoss(double maxWeeklyLoss) { maxWeeklyLoss_ = maxWeeklyLoss; }
    void setMaxMonthlyLoss(double maxMonthlyLoss) { maxMonthlyLoss_ = maxMonthlyLoss; }
    
    void setConservativeLimits();
    void setAggressiveLimits();
    void setDefaultLimits();
    
    bool isValid() const;
    void validate();
};

}