#pragma once

#include <vector>
#include <memory>

namespace math {

struct ProbabilityMetrics {
    double probabilityITM = 0.0;
    double probabilityTouch = 0.0;
    double expectedMove = 0.0;
    double breakeven = 0.0;
    double maxPain = 0.0;
    
    double probabilityProfit = 0.0;
    double probabilityMaxProfit = 0.0;
    double probabilityMaxLoss = 0.0;
    double expectedReturn = 0.0;
    double sharpeRatio = 0.0;
    
    double kelly = 0.0;
    double winRate = 0.0;
    double avgWin = 0.0;
    double avgLoss = 0.0;
    double profitFactor = 0.0;
};

class MathEngine;

class ProbabilityAnalyzer {
private:
    std::shared_ptr<MathEngine> mathEngine_;
    
public:
    explicit ProbabilityAnalyzer(std::shared_ptr<MathEngine> engine = nullptr);
    ~ProbabilityAnalyzer() = default;
    
    ProbabilityAnalyzer(const ProbabilityAnalyzer&) = delete;
    ProbabilityAnalyzer& operator=(const ProbabilityAnalyzer&) = delete;
    
    ProbabilityMetrics analyzeOption(double spot, double strike, double vol, 
                                   double timeToExpiry, bool isCall, double premium);
    
    ProbabilityMetrics analyzeStrategy(const std::vector<double>& strikes,
                                     const std::vector<bool>& isCall,
                                     const std::vector<double>& quantities,
                                     const std::vector<double>& premiums,
                                     double spot, double vol, double timeToExpiry);
    
    ProbabilityMetrics calculateHistoricalProbabilities(const std::vector<double>& prices,
                                                       double strike, int daysToExpiry);
    
private:
    double calculateTouchProbability(double spot, double barrier, double vol, double time);
    double calculateMaxPain(const std::vector<double>& strikes, 
                          const std::vector<double>& openInterest);
    double calculateKellyCriterion(double winRate, double avgWin, double avgLoss);
};

}