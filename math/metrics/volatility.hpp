#pragma once
namespace math
{
    struct VolatilityMetrics 
    {
        double realized = 0.0;
        double implied = 0.0;
        double hvRank = 0.0;
        double ivRank = 0.0;
        double skew = 0.0;
        double termStructure = 0.0;
    };
}

