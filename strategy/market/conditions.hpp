#pragma once

#include <string>

namespace strategy {

class MarketConditions {
private:
    double vixLevel_;
    double ivRank_;
    double hvRank_;
    std::string regime_;
    bool highVolEnvironment_;
    double correlations_;
    double trendStrength_;
    double marketStress_;
    double liquidityCondition_;

public:
    MarketConditions();
    ~MarketConditions() = default;
    
    double getVixLevel() const { return vixLevel_; }
    double getIvRank() const { return ivRank_; }
    double getHvRank() const { return hvRank_; }
    const std::string& getRegime() const { return regime_; }
    bool isHighVolEnvironment() const { return highVolEnvironment_; }
    double getCorrelations() const { return correlations_; }
    double getTrendStrength() const { return trendStrength_; }
    double getMarketStress() const { return marketStress_; }
    double getLiquidityCondition() const { return liquidityCondition_; }
    
    void setVixLevel(double vixLevel) { vixLevel_ = vixLevel; }
    void setIvRank(double ivRank) { ivRank_ = ivRank; }
    void setHvRank(double hvRank) { hvRank_ = hvRank; }
    void setRegime(const std::string& regime) { regime_ = regime; }
    void setHighVolEnvironment(bool highVolEnvironment) { highVolEnvironment_ = highVolEnvironment; }
    void setCorrelations(double correlations) { correlations_ = correlations; }
    void setTrendStrength(double trendStrength) { trendStrength_ = trendStrength; }
    void setMarketStress(double marketStress) { marketStress_ = marketStress; }
    void setLiquidityCondition(double liquidityCondition) { liquidityCondition_ = liquidityCondition; }
    
    bool isFearRegime() const { return regime_ == "Fear"; }
    bool isGreedRegime() const { return regime_ == "Greed"; }
    bool isNeutralRegime() const { return regime_ == "Neutral"; }
    bool isTrendingMarket() const { return trendStrength_ > 0.7; }
    bool isStressedMarket() const { return marketStress_ > 0.8; }
    bool isHighLiquidity() const { return liquidityCondition_ > 0.7; }
    
    std::string getMarketSummary() const;
    double getOverallRiskLevel() const;
    bool isSuitableForPremiumSelling() const;
    bool isSuitableForDirectional() const;
    bool isSuitableForVolatilityPlays() const;
};

}