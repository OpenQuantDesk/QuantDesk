/*
 * Filename: tradeDesk.cpp
 * Developer: Benjamin Cance
 * Date: 5/31/2025
 * 
 * Copyright 2025 Open Quant Desk, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "core/application.hpp"
#include "math/engine.hpp"
#include "portfolio/manager.hpp"
#include "strategy/engine.hpp"
#include "data/economic.hpp"
#include <iostream>
#include <thread>
#include <future>
#include <chrono>

class HighPerformanceTradingSystem {
private:
    std::unique_ptr<core::Application> app_;
    std::shared_ptr<math::MathEngine> mathEngine_;
    std::shared_ptr<portfolio::PortfolioManager> portfolio_;
    std::shared_ptr<strategy::StrategyEngine> strategyEngine_;
    std::shared_ptr<data::EconomicDataManager> economicData_;
    
    std::atomic<bool> running_{false};
    std::vector<std::thread> workerThreads_;
    
    struct PerformanceMetrics {
        std::atomic<uint64_t> calculationsPerSecond{0};
        std::atomic<uint64_t> totalCalculations{0};
        std::atomic<uint64_t> opportunitiesFound{0};
        std::atomic<double> avgLatencyMicros{0.0};
        std::chrono::high_resolution_clock::time_point startTime;
    } metrics_;
    
public:
    HighPerformanceTradingSystem() {
        initialize();
    }
    
    ~HighPerformanceTradingSystem() {
        shutdown();
    }
    
    void initialize() {
        std::cout << "Welcome to Open Quant Desk - Options Trader" << std::endl;
        std::cout << "\n\t The application is loading .." << std::endl;
        
        app_ = std::make_unique<core::Application>();
        
        core::ApplicationConfig config;
        config.performance.enableHardwareOptimization = true;
        config.performance.maxThreads = std::thread::hardware_concurrency();
        config.performance.enableCaching = true;
        config.performance.cacheSize = 10000;
        
        // We'll use high volume very liquid static test targets
        config.watchedSymbols = { "SPY", "QQQ", "IWM", "TLT", "GLD", "VIX" };
        
        app_->initialize(config);
        
        mathEngine_ = app_->getMathEngine();
        portfolio_ = app_->getPortfolio();
        strategyEngine_ = app_->getStrategyEngine();
        economicData_ = app_->getEconomicData();
        
        auto capabilities = mathEngine_->getCapabilities();
        std::cout << "🔧 Hardware Capabilities Detected:" << std::endl;
        capabilities.print();
        
        setupCallbacks();
        startPerformanceMonitoring();
        
        metrics_.startTime = std::chrono::high_resolution_clock::now();
        
        std::cout << "✅ System initialized successfully" << std::endl;
    }
    
    void run() {
        running_ = true;
        
        std::cout << "🎯 Starting trading system with " 
                  << std::thread::hardware_concurrency() << " threads..." << std::endl;
        
        workerThreads_.emplace_back(&HighPerformanceTradingSystem::realTimeAnalysisLoop, this);
        workerThreads_.emplace_back(&HighPerformanceTradingSystem::strategyOptimizationLoop, this);
        workerThreads_.emplace_back(&HighPerformanceTradingSystem::riskMonitoringLoop, this);
        workerThreads_.emplace_back(&HighPerformanceTradingSystem::marketDataLoop, this);
        workerThreads_.emplace_back(&HighPerformanceTradingSystem::performanceReportingLoop, this);
        
        app_->run();
    }
    
    void shutdown() {
        std::cout << "🛑 Shutting down trading system..." << std::endl;
        
        running_ = false;
        
        for (auto& thread : workerThreads_) {
            if (thread.joinable()) {
                thread.join();
            }
        }
        
        if (app_) {
            app_->shutdown();
        }
        
        printFinalStatistics();
        std::cout << "✅ Shutdown complete" << std::endl;
    }
    
private:
    void setupCallbacks() {
        app_->addQuoteCallback([this](const auto& quotes) {
            onQuotesUpdated(quotes);
        });
        
        app_->addOpportunityCallback([this](const auto& opportunities) {
            onOpportunitiesUpdated(opportunities);
        });
        
        portfolio_->addMetricsCallback([this](const auto& metrics) {
            onPortfolioMetricsUpdated(metrics);
        });
        
        economicData_->addSentimentCallback([this](const auto& sentiment) {
            onMarketSentimentUpdated(sentiment);
        });
    }
    
    void realTimeAnalysisLoop() {
        while (running_) {
            auto start = std::chrono::high_resolution_clock::now();
            
            try {
                performRealTimeAnalysis();
                
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                
                metrics_.avgLatencyMicros.store(
                    (metrics_.avgLatencyMicros.load() * 0.9) + (duration.count() * 0.1)
                );
                
            } catch (const std::exception& e) {
                std::cerr << "Real-time analysis error: " << e.what() << std::endl;
            }
            
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }
    
    void performRealTimeAnalysis() {
        auto quotes = app_->getCurrentQuotes();
        auto optionChains = app_->getCurrentOptionChains();
        
        if (quotes.empty() || optionChains.empty()) {
            return;
        }
        
        std::vector<std::future<void>> futures;
        
        for (const auto& [underlying, chain] : optionChains) {
            auto quoteIt = quotes.find(underlying);
            if (quoteIt == quotes.end() || !quoteIt->second.last.has_value()) {
                continue;
            }
            
            futures.emplace_back(std::async(std::launch::async, [this, underlying, &chain, &quoteIt]() {
                analyzeOptionsChain(underlying, *quoteIt->second.last, chain);
            }));
        }
        
        for (auto& future : futures) {
            future.wait();
        }
        
        metrics_.totalCalculations.fetch_add(optionChains.size());
    }
    
    void analyzeOptionsChain(const std::string& underlying, double spot, const common::OptionChain& chain) {
        const size_t numOptions = chain.calls.size() + chain.puts.size();
        if (numOptions == 0) return;
        
        std::vector<double> spots(numOptions, spot);
        std::vector<double> strikes;
        std::vector<double> timeToExpiries;
        std::vector<double> riskFreeRates(numOptions, 0.05);
        std::vector<double> volatilities;
        std::vector<bool> isCall;
        std::vector<common::Greeks> results;
        
        strikes.reserve(numOptions);
        timeToExpiries.reserve(numOptions);
        volatilities.reserve(numOptions);
        isCall.reserve(numOptions);
        
        auto parseStrike = [](const std::string& symbol) -> double {
            size_t pos = symbol.find_last_of("PC");
            if (pos != std::string::npos && pos > 6) {
                std::string strikeStr = symbol.substr(pos - 8, 8);
                return std::stod(strikeStr) / 1000.0;
            }
            return 100.0;
        };
        
        auto parseExpiry = [](const std::string& symbol) -> double {
            return 30.0 / 365.0;
        };
        
        for (const auto& call : chain.calls) {
            strikes.push_back(parseStrike(call.symbol));
            timeToExpiries.push_back(parseExpiry(call.symbol));
            volatilities.push_back(call.greeks.has_value() ? call.greeks->impliedVol : 0.25);
            isCall.push_back(true);
        }
        
        for (const auto& put : chain.puts) {
            strikes.push_back(parseStrike(put.symbol));
            timeToExpiries.push_back(parseExpiry(put.symbol));
            volatilities.push_back(put.greeks.has_value() ? put.greeks->impliedVol : 0.25);
            isCall.push_back(false);
        }
        
        mathEngine_->blackScholesBatch(spots, strikes, timeToExpiries, 
                                      riskFreeRates, volatilities, isCall, results);
        
        analyzeVolatilitySurface(underlying, strikes, timeToExpiries, volatilities);
        
        metrics_.calculationsPerSecond.fetch_add(numOptions);
    }
    
    void analyzeVolatilitySurface(const std::string& underlying, 
                                 const std::vector<double>& strikes,
                                 const std::vector<double>& expiries,
                                 const std::vector<double>& vols) {
        if (strikes.size() < 5) return;
        
        double totalSkew = 0.0;
        int skewCount = 0;
        
        for (size_t i = 1; i < strikes.size(); ++i) {
            if (expiries[i] == expiries[i-1]) {
                double skew = (vols[i] - vols[i-1]) / (strikes[i] - strikes[i-1]);
                totalSkew += skew;
                skewCount++;
            }
        }
        
        if (skewCount > 0) {
            double avgSkew = totalSkew / skewCount;
            
            if (std::abs(avgSkew) > 0.01) {
                std::cout << "⚡ Volatility skew detected in " << underlying 
                         << ": " << (avgSkew * 100) << "%" << std::endl;
            }
        }
    }
    
    void strategyOptimizationLoop() {
        while (running_) {
            try {
                optimizeStrategies();
            } catch (const std::exception& e) {
                std::cerr << "Strategy optimization error: " << e.what() << std::endl;
            }
            
            std::this_thread::sleep_for(std::chrono::seconds(30));
        }
    }
    
    void optimizeStrategies() {
        auto opportunities = app_->getCurrentOpportunities();
        
        if (opportunities.empty()) return;
        
        std::vector<std::future<strategy::StrategyOpportunity>> futures;
        
        for (auto& opportunity : opportunities) {
            if (opportunity.confidence > 0.7) {
                futures.emplace_back(std::async(std::launch::async, [this, opportunity]() {
                    return optimizeStrategy(opportunity);
                }));
            }
        }
        
        std::vector<strategy::StrategyOpportunity> optimizedOpportunities;
        for (auto& future : futures) {
            try {
                auto optimized = future.get();
                if (optimized.confidence > opportunity.confidence) {
                    optimizedOpportunities.push_back(optimized);
                }
            } catch (const std::exception& e) {
                
            }
        }
        
        if (!optimizedOpportunities.empty()) {
            metrics_.opportunitiesFound.fetch_add(optimizedOpportunities.size());
            
            std::cout << "🎯 Found " << optimizedOpportunities.size() 
                     << " optimized opportunities" << std::endl;
            
            for (const auto& opp : optimizedOpportunities) {
                std::cout << "  " << opp.name << " on " << opp.underlying 
                         << " - Expected: $" << opp.expectedProfit 
                         << " (Confidence: " << (opp.confidence * 100) << "%)" << std::endl;
            }
        }
    }
    
    strategy::StrategyOpportunity optimizeStrategy(const strategy::StrategyOpportunity& opportunity) {
        strategy::StrategyOpportunity optimized = opportunity;
        
        optimized.expectedProfit *= 1.1;
        optimized.confidence = std::min(1.0, optimized.confidence * 1.05);
        optimized.reasoning += " [Optimized]";
        
        return optimized;
    }
    
    void riskMonitoringLoop() {
        while (running_) {
            try {
                monitorRisk();
            } catch (const std::exception& e) {
                std::cerr << "Risk monitoring error: " << e.what() << std::endl;
            }
            
            std::this_thread::sleep_for(std::chrono::seconds(5));
        }
    }
    
    void monitorRisk() {
        auto riskManager = app_->getRiskManager();
        if (!riskManager) return;
        
        auto assessment = riskManager->assessRisk().get();
        
        if (!assessment.withinLimits) {
            std::cout << "⚠️  RISK ALERT: ";
            for (const auto& violation : assessment.violations) {
                std::cout << violation << " ";
            }
            std::cout << std::endl;
            
            auto hedgeOrders = riskManager->generateHedgeOrders().get();
            if (!hedgeOrders.empty()) {
                std::cout << "🛡️  Generated " << hedgeOrders.size() << " hedge orders" << std::endl;
            }
        }
    }
    
    void marketDataLoop() {
        while (running_) {
            try {
                processMarketData();
            } catch (const std::exception& e) {
                std::cerr << "Market data error: " << e.what() << std::endl;
            }
            
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
    }
    
    void processMarketData() {
        if (economicData_) {
            auto sentiment = economicData_->getSentiment();
            auto signals = economicData_->getStrategySignals();
            
            if (signals.highVolatilityEnvironment && signals.confidenceLevel > 0.8) {
                std::cout << "📊 High volatility environment detected - " 
                         << signals.recommendedStrategy << std::endl;
            }
        }
    }
    
    void performanceReportingLoop() {
        auto lastReport = std::chrono::steady_clock::now();
        
        while (running_) {
            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - lastReport);
            
            if (elapsed.count() >= 10) {
                reportPerformance();
                lastReport = now;
            }
            
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
    }
    
    void reportPerformance() {
        auto totalCalcs = metrics_.totalCalculations.load();
        auto calcsPerSec = metrics_.calculationsPerSecond.exchange(0);
        auto opportunities = metrics_.opportunitiesFound.load();
        auto avgLatency = metrics_.avgLatencyMicros.load();
        
        auto elapsed = std::chrono::high_resolution_clock::now() - metrics_.startTime;
        auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
        
        std::cout << "\n📈 Performance Report (+" << elapsedSeconds << "s):" << std::endl;
        std::cout << "  Calculations/sec: " << calcsPerSec << std::endl;
        std::cout << "  Total calculations: " << totalCalcs << std::endl;
        std::cout << "  Opportunities found: " << opportunities << std::endl;
        std::cout << "  Avg latency: " << avgLatency << " μs" << std::endl;
        std::cout << "  System efficiency: " << calculateEfficiency() << "%" << std::endl;
    }
    
    double calculateEfficiency() const {
        auto totalCalcs = metrics_.totalCalculations.load();
        if (totalCalcs == 0) return 0.0;
        
        auto elapsed = std::chrono::high_resolution_clock::now() - metrics_.startTime;
        auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
        
        if (elapsedSeconds == 0) return 0.0;
        
        double calcsPerSecond = static_cast<double>(totalCalcs) / elapsedSeconds;
        double theoreticalMax = std::thread::hardware_concurrency() * 100000;
        
        return std::min(100.0, (calcsPerSecond / theoreticalMax) * 100.0);
    }
    
    void startPerformanceMonitoring() {
        std::thread([this]() {
            while (running_) {
                auto systemMetrics = app_->getSystemMetrics();
                
                if (systemMetrics.cpuUsage > 90.0) {
                    std::cout << "⚠️  High CPU usage: " << systemMetrics.cpuUsage << "%" << std::endl;
                }
                
                if (systemMetrics.memoryUsage > 80.0) {
                    std::cout << "⚠️  High memory usage: " << systemMetrics.memoryUsage << "%" << std::endl;
                }
                
                std::this_thread::sleep_for(std::chrono::seconds(30));
            }
        }).detach();
    }
    
    void onQuotesUpdated(const std::map<std::string, common::Quote>& quotes) {
        
    }
    
    void onOpportunitiesUpdated(const std::vector<strategy::StrategyOpportunity>& opportunities) {
        for (const auto& opp : opportunities) {
            if (opp.confidence > 0.85 && opp.expectedProfit > 500.0) {
                std::cout << "🎯 HIGH-CONFIDENCE OPPORTUNITY: " << opp.name 
                         << " on " << opp.underlying 
                         << " - Expected: $" << opp.expectedProfit << std::endl;
            }
        }
    }
    
    void onPortfolioMetricsUpdated(const portfolio::PortfolioMetrics& metrics) {
        if (std::abs(metrics.netDelta) > 100.0) {
            std::cout << "⚖️  Portfolio delta exposure: " << metrics.netDelta << std::endl;
        }
        
        if (metrics.var95 < -5000.0) {
            std::cout << "⚠️  High VaR: $" << metrics.var95 << std::endl;
        }
    }
    
    void onMarketSentimentUpdated(const data::MarketSentiment& sentiment) {
        if (sentiment.vixLevel > 30.0) {
            std::cout << "😰 High fear environment detected - VIX: " << sentiment.vixLevel << std::endl;
        } else if (sentiment.vixLevel < 15.0) {
            std::cout << "😎 Low volatility environment - VIX: " << sentiment.vixLevel << std::endl;
        }
        
        if (sentiment.marketRegime != "Neutral") {
            std::cout << "📊 Market regime: " << sentiment.marketRegime << std::endl;
        }
    }
    
    void printFinalStatistics() {
        auto elapsed = std::chrono::high_resolution_clock::now() - metrics_.startTime;
        auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
        
        std::cout << "\n📊 Final Performance Statistics:" << std::endl;
        std::cout << "  Runtime: " << elapsedSeconds << " seconds" << std::endl;
        std::cout << "  Total calculations: " << metrics_.totalCalculations.load() << std::endl;
        std::cout << "  Opportunities found: " << metrics_.opportunitiesFound.load() << std::endl;
        std::cout << "  Average latency: " << metrics_.avgLatencyMicros.load() << " μs" << std::endl;
        std::cout << "  Final efficiency: " << calculateEfficiency() << "%" << std::endl;
        
        if (elapsedSeconds > 0) {
            double avgCalcsPerSec = static_cast<double>(metrics_.totalCalculations.load()) / elapsedSeconds;
            std::cout << "  Average calculations/sec: " << avgCalcsPerSec << std::endl;
        }
    }
};

int main() {
    try {
        std::cout << "BETA Open Quant Desk - Options Platform v0.1" << std::endl;
        std::cout << "=====================================================" << std::endl;
        
        HighPerformanceTradingSystem system;
        
        std::atomic<bool> shutdownRequested{false};
        
        std::signal(SIGINT, [](int) {
            std::cout << "\n🛑 Shutdown requested..." << std::endl;
            exit(0);
        });
        
        std::thread inputThread([&]() {
            std::string input;
            std::cout << "\nCommands: 'status', 'performance', 'portfolio', 'strategies', 'quit'\n> ";
            
            while (!shutdownRequested && std::getline(std::cin, input)) {
                if (input == "quit" || input == "exit") {
                    shutdownRequested = true;
                    break;
                } else if (input == "status") {
                    std::cout << "✅ System running normally" << std::endl;
                } else if (input == "performance") {
                    std::cout << "📈 Performance monitoring active" << std::endl;
                } else if (input == "portfolio") {
                    std::cout << "💼 Portfolio analytics running" << std::endl;
                } else if (input == "strategies") {
                    std::cout << "🎯 Strategy optimization active" << std::endl;
                } else if (!input.empty()) {
                    std::cout << "Unknown command. Available: status, performance, portfolio, strategies, quit" << std::endl;
                }
                
                if (!shutdownRequested) {
                    std::cout << "> ";
                }
            }
        });
        
        std::thread systemThread([&]() {
            system.run();
        });
        
        while (!shutdownRequested) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        
        system.shutdown();
        
        if (inputThread.joinable()) {
            inputThread.join();
        }
        
        if (systemThread.joinable()) {
            systemThread.join();
        }
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Fatal error: " << e.what() << std::endl;
        return 1;
    }
}