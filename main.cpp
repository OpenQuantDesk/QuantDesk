#include "core/application.hpp"
#include "gui/main_window.hpp"
#include <QApplication>
#include <iostream>
#include <memory>
#include <csignal>
#include <atomic>
#include <thread>
#include <chrono>

namespace {
    std::atomic<bool> shutdownRequested{false};
    std::unique_ptr<core::Application> globalApp;
}

void signalHandler(int signal) {
    std::cout << "\nReceived signal " << signal << ", initiating graceful shutdown..." << std::endl;
    shutdownRequested = true;
    if (globalApp) {
        globalApp->shutdown();
    }
}

void setupSignalHandlers() {
    std::signal(SIGINT, signalHandler);
    std::signal(SIGTERM, signalHandler);
#ifndef _WIN32
    std::signal(SIGHUP, signalHandler);
    std::signal(SIGQUIT, signalHandler);
#endif
}

void printBanner() {
    std::cout << R"(
  ╔═══════════════════════════════════════════════════════════╗
  ║                                                           ║
  ║     🎯 QUANTITATIVE OPTIONS TRADING PLATFORM v4.0        ║
  ║                                                           ║
  ║  • Advanced Mathematical Engine with Hardware Optimization║
  ║  • Real-time Portfolio Analytics & Risk Management       ║
  ║  • Multi-threaded Strategy Analysis & Optimization       ║
  ║  • Economic Data Integration & Market Regime Detection   ║
  ║  • High-Performance Monte Carlo Simulation Engine        ║
  ║                                                           ║
  ╚═══════════════════════════════════════════════════════════╝
)" << std::endl;
}

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [OPTIONS]\n\n"
              << "Options:\n"
              << "  --gui              Launch with graphical interface\n"
              << "  --config FILE      Use specific configuration file\n"
              << "  --watchlist FILE   Load watchlist from file\n"
              << "  --threads N        Set maximum number of threads\n"
              << "  --log-level LEVEL  Set logging level (DEBUG, INFO, WARN, ERROR)\n"
              << "  --benchmark        Run performance benchmark\n"
              << "  --help             Show this help message\n"
              << "  --version          Show version information\n\n"
              << "Environment Variables:\n"
              << "  BROKER_API_KEY     API key for broker connection\n"
              << "  FRED_API_KEY       Federal Reserve economic data API key\n"
              << "  ALPHA_VANTAGE_KEY  Alpha Vantage API key\n"
              << "  LOG_LEVEL          Logging level override\n"
              << std::endl;
}

void runBenchmark() {
    std::cout << "🔬 Running Performance Benchmark..." << std::endl;
    
    auto app = std::make_unique<core::Application>();
    app->initialize();
    
    auto mathEngine = app->getMathEngine();
    if (!mathEngine) {
        std::cerr << "Failed to initialize math engine for benchmark" << std::endl;
        return;
    }
    
    auto capabilities = mathEngine->getCapabilities();
    capabilities.print();
    
    const size_t numOptions = 10000;
    const size_t iterations = 100;
    
    std::vector<double> spots(numOptions, 100.0);
    std::vector<double> strikes(numOptions, 100.0);
    std::vector<double> timeToExpiries(numOptions, 30.0/365.0);
    std::vector<double> riskFreeRates(numOptions, 0.05);
    std::vector<double> volatilities(numOptions, 0.25);
    std::vector<bool> isCall(numOptions, true);
    std::vector<common::Greeks> results(numOptions);
    
    for (size_t i = 0; i < numOptions; ++i) {
        strikes[i] = 90.0 + (i % 20);
        volatilities[i] = 0.15 + (i % 10) * 0.01;
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (size_t iter = 0; iter < iterations; ++iter) {
        mathEngine->blackScholesBatch(spots, strikes, timeToExpiries, 
                                     riskFreeRates, volatilities, isCall, results);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    double totalCalculations = static_cast<double>(numOptions * iterations);
    double calculationsPerSecond = totalCalculations / (duration.count() / 1e6);
    
    std::cout << "\n📊 Benchmark Results:" << std::endl;
    std::cout << "  Options calculated: " << totalCalculations << std::endl;
    std::cout << "  Total time: " << duration.count() / 1000.0 << " ms" << std::endl;
    std::cout << "  Calculations/second: " << std::fixed << std::setprecision(0) 
              << calculationsPerSecond << std::endl;
    std::cout << "  Average per option: " << std::setprecision(2) 
              << (duration.count() / totalCalculations) << " μs" << std::endl;
    
    if (calculationsPerSecond > 1000000) {
        std::cout << "  Performance: 🚀 Excellent" << std::endl;
    } else if (calculationsPerSecond > 500000) {
        std::cout << "  Performance: ✅ Good" << std::endl;
    } else if (calculationsPerSecond > 100000) {
        std::cout << "  Performance: ⚠️  Moderate" << std::endl;
    } else {
        std::cout << "  Performance: 🐌 Needs optimization" << std::endl;
    }
}

int runConsoleMode(const core::ApplicationConfig& config) {
    try {
        globalApp = std::make_unique<core::Application>();
        globalApp->initialize(config);
        
        std::cout << "🚀 Starting Options Trading Platform..." << std::endl;
        std::cout << "📡 Hardware capabilities detected:" << std::endl;
        
        auto mathEngine = globalApp->getMathEngine();
        if (mathEngine) {
            mathEngine->getCapabilities().print();
        }
        
        std::thread appThread([&]() {
            globalApp->run();
        });
        
        std::thread statusThread([&]() {
            while (!shutdownRequested && globalApp->isRunning()) {
                auto metrics = globalApp->getSystemMetrics();
                
                std::cout << "\r📊 CPU: " << std::fixed << std::setprecision(1) 
                         << metrics.cpuUsage << "% | "
                         << "RAM: " << metrics.memoryUsage << "% | "
                         << "Threads: " << metrics.threadCount << " | "
                         << "Quotes: " << globalApp->getCurrentQuotes().size() << " | "
                         << "Opportunities: " << globalApp->getCurrentOpportunities().size()
                         << std::flush;
                
                std::this_thread::sleep_for(std::chrono::seconds(2));
            }
        });
        
        std::cout << "\n💡 Commands: 'help', 'status', 'quotes', 'portfolio', 'strategies', 'quit'" << std::endl;
        std::cout << "> ";
        
        std::string input;
        while (!shutdownRequested && std::getline(std::cin, input)) {
            if (input == "quit" || input == "exit") {
                break;
            } else if (input == "help") {
                std::cout << "\nAvailable commands:" << std::endl;
                std::cout << "  help       - Show this help" << std::endl;
                std::cout << "  status     - Show system status" << std::endl;
                std::cout << "  quotes     - Show current quotes" << std::endl;
                std::cout << "  portfolio  - Show portfolio summary" << std::endl;
                std::cout << "  strategies - Show strategy opportunities" << std::endl;
                std::cout << "  reload     - Reload watchlist" << std::endl;
                std::cout << "  benchmark  - Run performance benchmark" << std::endl;
                std::cout << "  quit/exit  - Shutdown application" << std::endl;
            } else if (input == "status") {
                auto metrics = globalApp->getSystemMetrics();
                std::cout << "\n📊 System Status:" << std::endl;
                std::cout << "  CPU Usage: " << metrics.cpuUsage << "%" << std::endl;
                std::cout << "  Memory Usage: " << metrics.memoryUsage << "%" << std::endl;
                std::cout << "  Active Threads: " << metrics.threadCount << std::endl;
                std::cout << "  Last Update: " << std::chrono::duration_cast<std::chrono::seconds>(
                    std::chrono::system_clock::now() - metrics.lastUpdate).count() << "s ago" << std::endl;
            } else if (input == "quotes") {
                auto quotes = globalApp->getCurrentQuotes();
                std::cout << "\n📈 Current Quotes (" << quotes.size() << "):" << std::endl;
                for (const auto& [symbol, quote] : quotes) {
                    if (quote.last.has_value()) {
                        std::cout << "  " << symbol << ": $" << std::fixed << std::setprecision(2) 
                                 << *quote.last;
                        if (quote.bid.has_value() && quote.ask.has_value()) {
                            std::cout << " (bid: " << *quote.bid << ", ask: " << *quote.ask << ")";
                        }
                        std::cout << std::endl;
                    }
                }
            } else if (input == "portfolio") {
                auto portfolio = globalApp->getPortfolio();
                if (portfolio) {
                    auto metrics = portfolio->calculateMetrics().get();
                    std::cout << "\n💼 Portfolio Summary:" << std::endl;
                    std::cout << "  Total Value: $" << std::fixed << std::setprecision(2) << metrics.totalValue << std::endl;
                    std::cout << "  Cash: $" << metrics.cash << std::endl;
                    std::cout << "  P&L: $" << metrics.totalPnL << std::endl;
                    std::cout << "  Net Delta: " << std::setprecision(1) << metrics.netDelta << std::endl;
                    std::cout << "  Net Theta: $" << std::setprecision(0) << metrics.netTheta << "/day" << std::endl;
                    std::cout << "  VaR (95%): $" << metrics.var95 << std::endl;
                }
            } else if (input == "strategies") {
                auto opportunities = globalApp->getCurrentOpportunities();
                std::cout << "\n🎯 Strategy Opportunities (" << opportunities.size() << "):" << std::endl;
                for (size_t i = 0; i < std::min(size_t(5), opportunities.size()); ++i) {
                    const auto& opp = opportunities[i];
                    std::cout << "  " << (i+1) << ". " << opp.name << " on " << opp.underlying << std::endl;
                    std::cout << "     Expected: $" << std::fixed << std::setprecision(0) << opp.expectedProfit
                             << " | Risk: $" << opp.maxRisk 
                             << " | Prob: " << std::setprecision(0) << (opp.probabilityProfit * 100) << "%"
                             << " | Conf: " << (opp.confidence * 100) << "%" << std::endl;
                }
            } else if (input == "reload") {
                globalApp->reloadWatchlist();
                std::cout << "✅ Watchlist reloaded" << std::endl;
            } else if (input == "benchmark") {
                runBenchmark();
            } else if (!input.empty()) {
                std::cout << "Unknown command: " << input << ". Type 'help' for available commands." << std::endl;
            }
            
            if (!shutdownRequested) {
                std::cout << "> ";
            }
        }
        
        shutdownRequested = true;
        statusThread.join();
        
        std::cout << "\n🛑 Shutting down..." << std::endl;
        globalApp->shutdown();
        
        if (appThread.joinable()) {
            appThread.join();
        }
        
        std::cout << "✅ Shutdown complete" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Fatal error: " << e.what() << std::endl;
        return 1;
    }
}

int runGuiMode(int argc, char* argv[], const core::ApplicationConfig& config) {
    QApplication app(argc, argv);
    app.setApplicationName("QuantTrader Pro");
    app.setApplicationVersion("4.0");
    app.setApplicationDisplayName("Quantitative Options Trading Platform");
    
    try {
        globalApp = std::make_unique<core::Application>();
        globalApp->initialize(config);
        
        auto mainWindow = std::make_unique<gui::MainWindow>();
        mainWindow->show();
        
        std::thread appThread([&]() {
            globalApp->run();
        });
        
        int result = app.exec();
        
        shutdownRequested = true;
        globalApp->shutdown();
        
        if (appThread.joinable()) {
            appThread.join();
        }
        
        return result;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ GUI Error: " << e.what() << std::endl;
        return 1;
    }
}

int main(int argc, char* argv[]) {
    setupSignalHandlers();
    printBanner();
    
    bool useGui = false;
    bool runBench = false;
    std::string configFile;
    std::string watchlistFile;
    std::string logLevel = "INFO";
    size_t maxThreads = 0;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--gui") {
            useGui = true;
        } else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "--version") {
            std::cout << "Quantitative Options Trading Platform v4.0" << std::endl;
            std::cout << "Built with advanced mathematical optimization" << std::endl;
            return 0;
        } else if (arg == "--benchmark") {
            runBench = true;
        } else if (arg == "--config" && i + 1 < argc) {
            configFile = argv[++i];
        } else if (arg == "--watchlist" && i + 1 < argc) {
            watchlistFile = argv[++i];
        } else if (arg == "--log-level" && i + 1 < argc) {
            logLevel = argv[++i];
        } else if (arg == "--threads" && i + 1 < argc) {
            maxThreads = std::stoul(argv[++i]);
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }
    
    if (runBench) {
        runBenchmark();
        return 0;
    }
    
    core::ApplicationConfig config;
    
    if (!configFile.empty()) {
        try {
            config = core::ConfigManager::loadFromFile(configFile);
            std::cout << "📄 Configuration loaded from " << configFile << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "⚠️  Failed to load config: " << e.what() << std::endl;
        }
    }
    
    config.logging.logLevel = logLevel;
    if (maxThreads > 0) {
        config.performance.maxThreads = maxThreads;
    }
    
    if (!watchlistFile.empty()) {
        std::ifstream file(watchlistFile);
        if (file.is_open()) {
            config.watchedSymbols.clear();
            std::string symbol;
            while (std::getline(file, symbol)) {
                if (!symbol.empty() && symbol[0] != '#') {
                    config.watchedSymbols.push_back(symbol);
                }
            }
            std::cout << "📋 Watchlist loaded: " << config.watchedSymbols.size() << " symbols" << std::endl;
        }
    }
    
    std::cout << "🔧 Configuration:" << std::endl;
    std::cout << "  Mode: " << (useGui ? "GUI" : "Console") << std::endl;
    std::cout << "  Log Level: " << config.logging.logLevel << std::endl;
    std::cout << "  Max Threads: " << (config.performance.maxThreads == 0 ? 
                                      std::thread::hardware_concurrency() : 
                                      config.performance.maxThreads) << std::endl;
    std::cout << "  Watched Symbols: " << config.watchedSymbols.size() << std::endl;
    std::cout << "  Hardware Optimization: " << (config.performance.enableHardwareOptimization ? "Enabled" : "Disabled") << std::endl;
    
    if (useGui) {
        return runGuiMode(argc, argv, config);
    } else {
        return runConsoleMode(config);
    }
}