#pragma once
namespace math 
{
    struct ProbabilityMetrics 
    {
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
}