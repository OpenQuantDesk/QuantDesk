#pragma once
namespace math {
    struct SimulationParams {
        double spot;
        double volatility;
        double riskFreeRate;
        double timeToExpiry;
        int numSimulations = 100000;
        int timeSteps = 252;
        bool useAntitheticVariates = true;
    };
} // namespace math