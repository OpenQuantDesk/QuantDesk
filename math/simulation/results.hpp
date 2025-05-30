#pragma once

namespace math 
{
    struct SimulationResult 
    {
    double expectedValue = 0.0;
    double standardError = 0.0;
    double probabilityProfit = 0.0;
    double valueAtRisk95 = 0.0;
    double valueAtRisk99 = 0.0;
    double maxDrawdown = 0.0;
    std::vector<double> pricePath;
    std::vector<double> payoffs;
    };
}