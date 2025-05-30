#pragma once

#include "../core/types.hpp"
#include "calculator.hpp"
#include <memory>

namespace math::greeks {

class SensitivityEngine {
private:
    std::shared_ptr<GreeksCalculator> calculator_;
    
public:
    explicit SensitivityEngine(std::shared_ptr<GreeksCalculator> calculator);
    
    struct SensitivityParams {
        core::Real spotBumpSize = 0.01;
        core::Real volBumpSize = 0.01;
        core::Real rateBumpSize = 0.0001;
        core::Real timeBumpSize = 1.0 / 365.0;
        bool useRelativeBumps = true;
        core::Integer numScenarios = 1000;
    };
    
    struct SensitivityProfile {
        core::Vector spotSensitivities;
        core::Vector volSensitivities;
        core::Vector rateSensitivities;
        core::Vector timeSensitivities;
        core::Matrix crossSensitivities;
        core::Real totalSpotSensitivity = 0.0;
        core::Real totalVolSensitivity = 0.0;
        core::Real totalRateSensitivity = 0.0;
        core::Real totalTimeSensitivity = 0.0;
        core::Vector riskFactors;
        core::Real totalRisk = 0.0;
    };
    
    SensitivityProfile calculateSensitivities(const core::MarketDataVector& markets,
                                             const core::Vector& quantities,
                                             const SensitivityParams& params = SensitivityParams{}) const;
    
    struct StressTestParams {
        core::Vector spotShocks = {-0.30, -0.20, -0.10, -0.05, 0.05, 0.10, 0.20, 0.30};
        core::Vector volShocks = {-0.10, -0.05, -0.02, 0.02, 0.05, 0.10};
        core::Vector rateShocks = {-0.02, -0.01, -0.005, 0.005, 0.01, 0.02};
        core::Vector timeShocks = {1.0/365.0, 7.0/365.0, 30.0/365.0, 90.0/365.0};
        bool includeCorrelationShocks = false;
        core::Real correlationShock = 0.2;
    };
    
    struct StressTestResult {
        core::Vector spotStresses;
        core::Vector volStresses;
        core::Vector rateStresses;
        core::Vector timeStresses;
        core::Real worstCase = 0.0;
        core::Real bestCase = 0.0;
        core::Real var95 = 0.0;
        core::Real var99 = 0.0;
        core::Real expectedShortfall95 = 0.0;
        core::Real expectedShortfall99 = 0.0;
        core::Vector stressContributions;
    };
    
    StressTestResult runStressTest(const core::MarketDataVector& markets,
                                  const core::Vector& quantities,
                                  const StressTestParams& params = StressTestParams{}) const;
    
    struct ScenarioAnalysis {
        core::Vector scenarioValues;
        core::Vector pnlValues;
        core::Real expectedPnl = 0.0;
        core::Real pnlVolatility = 0.0;
        core::Real worstCase = 0.0;
        core::Real bestCase = 0.0;
        core::Real probabilityOfLoss = 0.0;
        core::Real sharpeRatio = 0.0;
        core::Real skewness = 0.0;
        core::Real kurtosis = 0.0;
    };
    
    ScenarioAnalysis runScenarioAnalysis(const core::MarketDataVector& markets,
                                        const core::Vector& quantities,
                                        const core::Matrix& scenarios,
                                        const core::Vector& probabilities = {}) const;
    
    struct CorrelationAnalysis {
        core::Matrix correlationMatrix;
        core::Real averageCorrelation = 0.0;
        core::Real maxCorrelation = 0.0;
        core::Real minCorrelation = 0.0;
        core::Vector eigenvalues;
        core::Matrix eigenvectors;
        core::Real conditionalVaR = 0.0;
        core::Vector marginalContributions;
    };
    
    CorrelationAnalysis analyzeCorrelations(const core::MarketDataVector& markets,
                                          const core::Matrix& historicalReturns) const;
    
    struct HedgingSensitivity {
        core::Vector deltaHedgeRatios;
        core::Vector gammaHedgeRatios;
        core::Vector vegaHedgeRatios;
        core::Matrix hedgeEffectiveness;
        core::Real optimalHedgeRatio = 0.0;
        core::Real minimumVarianceRatio = 0.0;
        core::Real hedgingCost = 0.0;
    };
    
    HedgingSensitivity analyzeHedgingSensitivity(const core::MarketDataVector& portfolioMarkets,
                                               const core::Vector& portfolioQuantities,
                                               const core::MarketDataVector& hedgeInstruments) const;
    
    struct TailRiskAnalysis {
        core::Real tailRisk95 = 0.0;
        core::Real tailRisk99 = 0.0;
        core::Real expectedTailLoss95 = 0.0;
        core::Real expectedTailLoss99 = 0.0;
        core::Real extremeValueIndex = 0.0;
        core::Real tailDependence = 0.0;
        core::Vector tailContributions;
    };
    
    TailRiskAnalysis analyzeTailRisk(const core::MarketDataVector& markets,
                                    const core::Vector& quantities,
                                    const core::Matrix& historicalScenarios) const;
    
    struct TimeVaryingSensitivity {
        core::Matrix timeEvolution;
        core::Vector volatilityRegimes;
        core::Matrix regimeTransitions;
        core::Vector averageSensitivities;
        core::Vector sensitivityVolatilities;
        core::Real persistenceMetric = 0.0;
    };
    
    TimeVaryingSensitivity analyzeTimeVaryingSensitivity(const std::vector<core::MarketDataVector>& timeSeriesMarkets,
                                                        const std::vector<core::Vector>& timeSeriesQuantities) const;
    
    struct NonLinearSensitivity {
        core::Matrix hessianMatrix;
        core::Vector thirdOrderTerms;
        core::Real convexityAdjustment = 0.0;
        core::Real nonLinearityIndex = 0.0;
        core::Vector taylorExpansionTerms;
        core::Real approximationError = 0.0;
    };
    
    NonLinearSensitivity analyzeNonLinearSensitivity(const core::MarketDataVector& markets,
                                                    const core::Vector& quantities,
                                                    const SensitivityParams& params = SensitivityParams{}) const;
    
private:
    core::Real calculateSpotSensitivity(const core::MarketData& market, core::Real quantity,
                                       const SensitivityParams& params) const;
    
    core::Real calculateVolSensitivity(const core::MarketData& market, core::Real quantity,
                                      const SensitivityParams& params) const;
    
    core::Real calculateRateSensitivity(const core::MarketData& market, core::Real quantity,
                                       const SensitivityParams& params) const;
    
    core::Real calculateTimeSensitivity(const core::MarketData& market, core::Real quantity,
                                       const SensitivityParams& params) const;
    
    core::Matrix calculateCrossSensitivities(const core::MarketDataVector& markets,
                                            const core::Vector& quantities,
                                            const SensitivityParams& params) const;
    
    core::Real calculateSecondOrderSensitivity(const core::MarketData& market, core::Real quantity,
                                              const SensitivityParams& params) const;
    
    core::Real calculateCrossDerivative(const core::MarketDataVector& markets,
                                       const core::Vector& quantities,
                                       core::Size i, core::Size j,
                                       const SensitivityParams& params,
                                       core::Real baseValue) const;
    
    core::Real calculatePortfolioValue(const core::MarketDataVector& markets,
                                      const core::Vector& quantities) const;
    
    core::Matrix calculateCovarianceMatrix(const core::Matrix& returns) const;
    core::Vector calculateEigenvalues(const core::Matrix& matrix) const;
    core::Matrix calculateEigenvectors(const core::Matrix& matrix) const;
    
    core::Real calculateConditionalVaR(const core::Vector& returns, const core::Matrix& correlations,
                                      core::Real confidence = 0.95) const;
    
    core::Vector calculateMarginalRiskContributions(const core::MarketDataVector& markets,
                                                   const core::Vector& quantities,
                                                   const core::Matrix& correlations) const;
};

class SensitivityReporting {
private:
    std::shared_ptr<SensitivityEngine> engine_;
    
public:
    explicit SensitivityReporting(std::shared_ptr<SensitivityEngine> engine);
    
    struct RiskReport {
        std::string reportDate;
        std::string portfolioId;
        SensitivityEngine::SensitivityProfile sensitivities;
        SensitivityEngine::StressTestResult stressResults;
        SensitivityEngine::CorrelationAnalysis correlations;
        SensitivityEngine::TailRiskAnalysis tailRisk;
        std::vector<std::string> riskAlerts;
        core::Real overallRiskScore = 0.0;
    };
    
    RiskReport generateRiskReport(const core::MarketDataVector& markets,
                                 const core::Vector& quantities,
                                 const std::string& portfolioId = "DEFAULT") const;
    
    struct ComplianceReport {
        bool withinLimits = true;
        std::vector<std::string> violations;
        core::Real utilizationPercentage = 0.0;
        core::Vector limitUsage;
        std::vector<std::string> recommendations;
    };
    
    ComplianceReport checkCompliance(const SensitivityEngine::SensitivityProfile& profile,
                                    const core::Vector& limits) const;
    
    void exportToHTML(const RiskReport& report, const std::string& filename) const;
    void exportToJSON(const RiskReport& report, const std::string& filename) const;
    void exportToPDF(const RiskReport& report, const std::string& filename) const;
    
    struct AlertThresholds {
        core::Real deltaThreshold = 1000.0;
        core::Real gammaThreshold = 500.0;
        core::Real vegaThreshold = 5000.0;
        core::Real var95Threshold = 100000.0;
        core::Real concentrationThreshold = 0.25;
        core::Real correlationThreshold = 0.8;
    };
    
    std::vector<std::string> generateAlerts(const SensitivityEngine::SensitivityProfile& profile,
                                           const SensitivityEngine::StressTestResult& stressResults,
                                           const AlertThresholds& thresholds = AlertThresholds{}) const;
    
private:
    core::Real calculateRiskScore(const SensitivityEngine::SensitivityProfile& profile,
                                 const SensitivityEngine::StressTestResult& stressResults) const;
    
    std::string formatSensitivityValue(core::Real value, const std::string& type) const;
    std::string generateExecutiveSummary(const RiskReport& report) const;
};

class RealTimeSensitivityMonitor {
private:
    std::shared_ptr<SensitivityEngine> engine_;
    std::shared_ptr<SensitivityReporting> reporting_;
    mutable core::ReadWriteLock dataLock_;
    
    struct MonitoringData {
        core::MarketDataVector currentMarkets;
        core::Vector currentQuantities;
        SensitivityEngine::SensitivityProfile lastProfile;
        std::vector<SensitivityReporting::RiskReport> historicalReports;
        core::Integer lastUpdateTime = 0;
    };
    
    mutable MonitoringData data_;
    bool isMonitoring_ = false;
    core::Integer updateIntervalSeconds_ = 60;
    
public:
    explicit RealTimeSensitivityMonitor(std::shared_ptr<SensitivityEngine> engine,
                                       std::shared_ptr<SensitivityReporting> reporting);
    
    void startMonitoring(const core::MarketDataVector& markets, const core::Vector& quantities);
    void stopMonitoring();
    void updatePositions(const core::MarketDataVector& markets, const core::Vector& quantities);
    
    void setUpdateInterval(core::Integer seconds) { updateIntervalSeconds_ = seconds; }
    core::Integer getUpdateInterval() const { return updateIntervalSeconds_; }
    
    SensitivityEngine::SensitivityProfile getCurrentSensitivities() const;
    SensitivityReporting::RiskReport getLatestReport() const;
    
    std::vector<std::string> getActiveAlerts() const;
    
    struct MonitoringSettings {
        bool enableAlerts = true;
        bool enableReporting = true;
        bool enableHistoricalTracking = true;
        core::Integer maxHistoricalReports = 100;
        SensitivityReporting::AlertThresholds alertThresholds;
    };
    
    void setMonitoringSettings(const MonitoringSettings& settings) { settings_ = settings; }
    const MonitoringSettings& getMonitoringSettings() const { return settings_; }
    
private:
    MonitoringSettings settings_;
    
    void performUpdate();
    void checkForAlerts(const SensitivityEngine::SensitivityProfile& profile);
    void storeHistoricalData(const SensitivityReporting::RiskReport& report);
};

inline SensitivityEngine::SensitivityEngine(std::shared_ptr<GreeksCalculator> calculator)
    : calculator_(std::move(calculator)) {}

inline SensitivityReporting::SensitivityReporting(std::shared_ptr<SensitivityEngine> engine)
    : engine_(std::move(engine)) {}

inline RealTimeSensitivityMonitor::RealTimeSensitivityMonitor(
    std::shared_ptr<SensitivityEngine> engine,
    std::shared_ptr<SensitivityReporting> reporting)
    : engine_(std::move(engine)), reporting_(std::move(reporting)) {}

}