#pragma once

#include "../core/engine.hpp"
#include "../core/types.hpp"
#include "models.hpp"
#include "surface.hpp"
#include <memory>
#include <string>

namespace math::volatility {

    class VolatilityCalibrator {
    private:
        std::shared_ptr<core::MathEngine> engine_;
    public:
        explicit VolatilityCalibrator(std::shared_ptr<core::MathEngine> engine);

        struct CalibrationParams {
            core::Real tolerance = 1e-4;
            core::Integer maxIterations = 100;
            core::Real stepSize = 0.01;
            bool enforceArbitrage = true;
            core::Real smoothing = 0.1;
            bool useWeights = false;
            core::Vector weights;
        };

        struct CalibrationResult {
            bool success = false;
            core::Real calibrationError = 0.0;
            std::string errorMessage;
            core::Vector strikes;
            core::Vector maturities;
            core::Matrix impliedVols;
            core::Vector residuals;
            core::Integer iterations = 0;
            core::Real finalRMSE = 0.0;
        };

        struct ModelCalibrationResult {
            bool success = false;
            core::Real calibrationError = 0.0;
            std::string errorMessage;
            HestonModel::HestonParams hestonParams;
            StochasticVolatilityModel::SABRParams sabrParams;
            JumpDiffusionModel::JumpParams jumpParams;
            core::Integer iterations = 0;
            core::Vector modelPrices;
            core::Vector marketPrices;
        };

        struct ValidationResult {
            bool isValid = true;
            core::Real averageError = 0.0;
            core::Real maxError = 0.0;
            core::Size largeErrorCount = 0;
            core::Vector priceErrors;
            std::vector<std::string> issues;
        };

        CalibrationResult
        calibrateImpliedSurface(const VolatilityPoints& marketData,
                                core::Real spot, core::Real riskFreeRate,
                                const CalibrationParams& params
                                = CalibrationParams{}) const;

        ModelCalibrationResult
        calibrateHeston(const VolatilityPoints& marketData, core::Real spot,
                        core::Real riskFreeRate,
                        const CalibrationParams& params
                        = CalibrationParams{}) const;

        ModelCalibrationResult calibrateSABR(const VolatilityPoints& marketData,
                                             core::Real forward,
                                             core::Real timeToExpiry,
                                             const CalibrationParams& params
                                             = CalibrationParams{}) const;

        ModelCalibrationResult calibrateJumpDiffusion(
            const VolatilityPoints& marketData, core::Real spot,
            core::Real riskFreeRate, core::Real baseVol,
            const CalibrationParams& params = CalibrationParams{}) const;

        ValidationResult
        validateCalibration(const VolatilityPoints& marketData,
                            const CalibrationResult& calibrationResult,
                            core::Real spot, core::Real riskFreeRate) const;

        struct OptimizationResult {
            bool converged = false;
            core::Integer iterations = 0;
            core::Real finalObjective = 0.0;
            core::Vector finalParameters;
            core::Matrix jacobian;
            core::Matrix hessian;
        };

        OptimizationResult optimizeParameters(
            const core::Vector& initialParams,
            const std::function<core::Real(const core::Vector&)>& objective,
            const std::function<core::Vector(const core::Vector&)>& gradient,
            const CalibrationParams& params = CalibrationParams{}) const;

        struct CrossValidationResult {
            core::Real averageError = 0.0;
            core::Real standardDeviation = 0.0;
            core::Vector foldErrors;
            bool isStable = true;
        };

        CrossValidationResult crossValidate(const VolatilityPoints& marketData,
                                            core::Real spot,
                                            core::Real riskFreeRate,
                                            core::Integer numFolds = 5) const;

        struct BenchmarkResult {
            core::Real blackScholesRMSE = 0.0;
            core::Real hestonRMSE = 0.0;
            core::Real sabrRMSE = 0.0;
            core::Real jumpDiffusionRMSE = 0.0;
            std::string bestModel;
            core::Real improvementRatio = 0.0;
        };

        BenchmarkResult benchmarkModels(const VolatilityPoints& marketData,
                                        core::Real spot,
                                        core::Real riskFreeRate) const;
    private:
        void interpolateMissingPoints(core::Matrix& impliedVols,
                                      const core::Vector& strikes,
                                      const core::Vector& maturities) const;

        void enforceArbitrageConstraints(core::Matrix& impliedVols,
                                         const core::Vector& strikes,
                                         const core::Vector& maturities,
                                         core::Real spot,
                                         core::Real riskFreeRate) const;

        void applySmoothening(core::Matrix& impliedVols,
                              core::Real factor) const;

        core::Real calculateCalibrationError(const VolatilityPoints& marketData,
                                             const core::Matrix& impliedVols,
                                             const core::Vector& strikes,
                                             const core::Vector& maturities,
                                             core::Real spot,
                                             core::Real riskFreeRate) const;

        void optimizeHestonParams(HestonModel::HestonParams& params,
                                  const VolatilityPoints& marketData,
                                  core::Real spot, core::Real riskFreeRate,
                                  core::Real stepSize) const;

        core::Vector calculateGradient(
            const core::Vector& params,
            const std::function<core::Real(const core::Vector&)>& objective,
            core::Real epsilon = 1e-8) const;

        core::Matrix calculateHessian(
            const core::Vector& params,
            const std::function<core::Real(const core::Vector&)>& objective,
            core::Real epsilon = 1e-6) const;

        bool checkConvergence(const core::Vector& oldParams,
                              const core::Vector& newParams,
                              core::Real tolerance) const;

        core::Real lineSearch(
            const core::Vector& params, const core::Vector& direction,
            const std::function<core::Real(const core::Vector&)>& objective,
            core::Real maxStep = 1.0) const;
    };

    class AdvancedCalibrator {
    private:
        std::shared_ptr<core::MathEngine> engine_;
        std::shared_ptr<VolatilityCalibrator> baseCalibrator_;
    public:
        explicit AdvancedCalibrator(std::shared_ptr<core::MathEngine> engine);

        struct MultiObjectiveParams {
            core::Real priceWeight = 1.0;
            core::Real greeksWeight = 0.5;
            core::Real smoothnessWeight = 0.1;
            core::Real arbitrageWeight = 2.0;
            bool penalizeExtrapolation = true;
        };

        struct RobustCalibrationResult {
            VolatilityCalibrator::CalibrationResult baseResult;
            core::Matrix confidenceIntervals;
            core::Vector parameterStability;
            core::Real robustnessScore = 0.0;
            bool isRobust = true;
        };

        RobustCalibrationResult
        robustCalibration(const VolatilityPoints& marketData, core::Real spot,
                          core::Real riskFreeRate,
                          const MultiObjectiveParams& params
                          = MultiObjectiveParams{}) const;

        struct TimeSeriesCalibrationResult {
            std::vector<VolatilityCalibrator::CalibrationResult>
                timeSeriesResults;
            core::Matrix parameterEvolution;
            core::Vector stabilityMetrics;
            core::Real averageCalibrationTime = 0.0;
        };

        TimeSeriesCalibrationResult calibrateTimeSeries(
            const std::vector<VolatilityPoints>& timeSeriesData,
            const core::Vector& spots, const core::Vector& riskFreeRates,
            const MultiObjectiveParams& params = MultiObjectiveParams{}) const;

        struct RegimeDetectionResult {
            core::Vector regimeIndicators;
            core::Integer numRegimes = 0;
            core::Vector regimeVolatilities;
            core::Vector transitionProbabilities;
        };

        RegimeDetectionResult detectVolatilityRegimes(
            const std::vector<VolatilityPoints>& timeSeriesData,
            core::Integer maxRegimes = 3) const;

        struct UncertaintyQuantification {
            core::Matrix parameterCovariance;
            core::Vector parameterStdErrors;
            core::Matrix priceConfidenceIntervals;
            core::Real modelUncertainty = 0.0;
        };

        UncertaintyQuantification quantifyUncertainty(
            const VolatilityPoints& marketData,
            const VolatilityCalibrator::CalibrationResult& calibResult,
            core::Real spot, core::Real riskFreeRate,
            core::Integer numBootstraps = 1000) const;
    private:
        core::Real multiObjectiveFunction(const core::Vector& params,
                                          const VolatilityPoints& marketData,
                                          const MultiObjectiveParams& objParams,
                                          core::Real spot,
                                          core::Real riskFreeRate) const;

        core::Real
        calculateSmoothnessPenalty(const core::Matrix& impliedVols) const;
        core::Real calculateArbitragePenalty(const core::Matrix& impliedVols,
                                             const core::Vector& strikes,
                                             const core::Vector& maturities,
                                             core::Real spot,
                                             core::Real riskFreeRate) const;

        core::Vector bootstrapSample(const VolatilityPoints& marketData,
                                     std::mt19937& rng) const;

        core::Real calculateRobustnessScore(
            const std::vector<VolatilityCalibrator::CalibrationResult>&
                bootstrapResults) const;
    };

    class CalibrationDiagnostics {
    private:
        std::shared_ptr<core::MathEngine> engine_;
    public:
        explicit CalibrationDiagnostics(
            std::shared_ptr<core::MathEngine> engine);

        struct DiagnosticResult {
            core::Vector residuals;
            core::Vector standardizedResiduals;
            core::Vector leverageValues;
            core::Vector cookDistances;
            core::Vector outlierIndices;
            bool hasOutliers = false;
            bool passesNormalityTest = true;
            bool passesHomoscedasticityTest = true;
        };

        DiagnosticResult analyzeResiduals(
            const VolatilityPoints& marketData,
            const VolatilityCalibrator::CalibrationResult& calibResult,
            core::Real spot, core::Real riskFreeRate) const;

        struct SensitivityAnalysis {
            core::Matrix parameterSensitivities;
            core::Vector elasticities;
            core::Real conditionNumber = 0.0;
            bool isWellConditioned = true;
        };

        SensitivityAnalysis analyzeSensitivity(
            const VolatilityPoints& marketData,
            const VolatilityCalibrator::CalibrationResult& calibResult,
            core::Real spot, core::Real riskFreeRate) const;

        struct GoodnessOfFit {
            core::Real rSquared = 0.0;
            core::Real adjustedRSquared = 0.0;
            core::Real aic = 0.0;
            core::Real bic = 0.0;
            core::Real logLikelihood = 0.0;
            bool isGoodFit = true;
        };

        GoodnessOfFit assessGoodnessOfFit(
            const VolatilityPoints& marketData,
            const VolatilityCalibrator::CalibrationResult& calibResult,
            core::Real spot, core::Real riskFreeRate) const;
    private:
        core::Real calculateLeverageValue(core::Size index,
                                          const core::Matrix& jacobian) const;
        core::Real calculateCookDistance(core::Size index,
                                         const core::Vector& residuals,
                                         const core::Matrix& jacobian,
                                         core::Real mse) const;

        bool jarqueBeraTest(const core::Vector& residuals,
                            core::Real alpha = 0.05) const;
        bool breuschPaganTest(const core::Vector& residuals,
                              const core::Matrix& regressors,
                              core::Real alpha = 0.05) const;
    };

    inline VolatilityCalibrator::VolatilityCalibrator(
        std::shared_ptr<core::MathEngine> engine)
        : engine_(std::move(engine))
    {}

    inline AdvancedCalibrator::AdvancedCalibrator(
        std::shared_ptr<core::MathEngine> engine)
        : engine_(engine),
          baseCalibrator_(std::make_shared<VolatilityCalibrator>(engine))
    {}

    inline CalibrationDiagnostics::CalibrationDiagnostics(
        std::shared_ptr<core::MathEngine> engine)
        : engine_(std::move(engine))
    {}

} // namespace math::volatility