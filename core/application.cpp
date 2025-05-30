#include "core/application.hpp"
#include "broker/registry.hpp"
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

namespace core {

Application::Application() {
    config_ = loadDefaultConfig();
}

Application::~Application() {
    shutdown();
}

void Application::notifyStatusUpdate(const std::string& status) {
    std::shared_lock<std::shared_mutex> lock(dataLock_);
    for (const auto& callback : statusCallbacks_) {
        try {
            callback(status);
        } catch (const std::exception& e) {
            
        }
    }
}

ApplicationConfig Application::loadDefaultConfig() const {
    ApplicationConfig config;
    
    config.watchedSymbols = {"SPY", "QQQ", "IWM"};
    config.updateIntervalSeconds = 10;
    config.enableBackgroundAnalysis = true;
    config.enableRiskManagement = true;
    config.configFile = "config.json";
    
    return config;
}

std::vector<std::string> Application::loadWatchlistFromFile(const std::string& filename) const {
    std::vector<std::string> symbols;
    
    try {
        std::ifstream file(filename);
        if (file.is_open()) {
            std::string line;
            while (std::getline(file, line)) {
                line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
                if (!line.empty() && line[0] != '#') {
                    std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                    symbols.push_back(line);
                }
            }
        }
    } catch (const std::exception& e) {
        
    }
    
    if (symbols.empty()) {
        symbols = {"SPY", "QQQ", "IWM"};
    }
    
    return symbols;
}

void Application::saveWatchlistToFile(const std::string& filename) const {
    try {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << "# Options Trading Watchlist\n";
            file << "# Lines starting with # are comments\n";
            
            std::shared_lock<std::shared_mutex> lock(dataLock_);
            for (const auto& symbol : config_.watchedSymbols) {
                file << symbol << "\n";
            }
        }
    } catch (const std::exception& e) {
        
    }
}

ApplicationConfig ConfigManager::loadFromFile(const std::string& filename) {
    ApplicationConfig config;
    
    try {
        std::ifstream file(filename);
        if (file.is_open()) {
            nlohmann::json j;
            file >> j;
            
            if (j.contains("watchedSymbols")) {
                config.watchedSymbols = j["watchedSymbols"].get<std::vector<std::string>>();
            }
            
            if (j.contains("updateIntervalSeconds")) {
                config.updateIntervalSeconds = j["updateIntervalSeconds"];
            }
            
            if (j.contains("enableBackgroundAnalysis")) {
                config.enableBackgroundAnalysis = j["enableBackgroundAnalysis"];
            }
            
            if (j.contains("enableRiskManagement")) {
                config.enableRiskManagement = j["enableRiskManagement"];
            }
            
            if (j.contains("brokers")) {
                for (const auto& brokerJson : j["brokers"]) {
                    common::BrokerConfig brokerConfig;
                    brokerConfig.brokerName = brokerJson["name"];
                    brokerConfig.sandboxMode = brokerJson.value("sandboxMode", true);
                    
                    if (brokerJson.contains("parameters")) {
                        brokerConfig.parameters = brokerJson["parameters"].get<std::map<std::string, std::string>>();
                    }
                    
                    config.brokers.push_back(brokerConfig);
                }
            }
            
            if (j.contains("economicDataKeys")) {
                config.economicDataKeys = j["economicDataKeys"].get<std::map<std::string, std::string>>();
            }
        }
    } catch (const std::exception& e) {
        
    }
    
    return config;
}

void ConfigManager::saveToFile(const ApplicationConfig& config, const std::string& filename) {
    try {
        nlohmann::json j;
        
        j["watchedSymbols"] = config.watchedSymbols;
        j["updateIntervalSeconds"] = config.updateIntervalSeconds;
        j["enableBackgroundAnalysis"] = config.enableBackgroundAnalysis;
        j["enableRiskManagement"] = config.enableRiskManagement;
        j["economicDataKeys"] = config.economicDataKeys;
        
        nlohmann::json brokersJson = nlohmann::json::array();
        for (const auto& broker : config.brokers) {
            nlohmann::json brokerJson;
            brokerJson["name"] = broker.brokerName;
            brokerJson["sandboxMode"] = broker.sandboxMode;
            brokerJson["parameters"] = broker.parameters;
            brokersJson.push_back(brokerJson);
        }
        j["brokers"] = brokersJson;
        
        std::ofstream file(filename);
        if (file.is_open()) {
            file << j.dump(4);
        }
    } catch (const std::exception& e) {
        
    }
}

ApplicationConfig ConfigManager::loadFromEnvironment() {
    ApplicationConfig config;
    
    if (const char* symbols = std::getenv("WATCHED_SYMBOLS")) {
        std::stringstream ss(symbols);
        std::string symbol;
        while (std::getline(ss, symbol, ',')) {
            symbol.erase(std::remove_if(symbol.begin(), symbol.end(), ::isspace), symbol.end());
            if (!symbol.empty()) {
                config.watchedSymbols.push_back(symbol);
            }
        }
    }
    
    if (const char* interval = std::getenv("UPDATE_INTERVAL")) {
        try {
            config.updateIntervalSeconds = std::stoi(interval);
        } catch (const std::exception& e) {
            
        }
    }
    
    if (const char* analysis = std::getenv("ENABLE_ANALYSIS")) {
        config.enableBackgroundAnalysis = (std::string(analysis) == "true");
    }
    
    if (const char* risk = std::getenv("ENABLE_RISK")) {
        config.enableRiskManagement = (std::string(risk) == "true");
    }
    
    return config;
}

common::BrokerConfig ConfigManager::parseBrokerConfig(const std::map<std::string, std::string>& params) {
    common::BrokerConfig config;
    
    auto nameIt = params.find("name");
    if (nameIt != params.end()) {
        config.brokerName = nameIt->second;
    }
    
    auto sandboxIt = params.find("sandbox");
    if (sandboxIt != params.end()) {
        config.sandboxMode = (sandboxIt->second == "true");
    }
    
    config.parameters = params;
    
    return config;
}

std::map<std::string, std::string> ConfigManager::parseEconomicDataConfig(const std::map<std::string, std::string>& params) {
    return params;
}

}initialize() {
    initialize(config_);
}

void Application::initialize(const ApplicationConfig& config) {
    config_ = config;
    
    loadConfiguration();
    initializeComponents();
    loadPlugins();
    initializeBrokers();
}

void Application::run() {
    if (running_.load()) return;
    
    running_ = true;
    startBackgroundThreads();
    
    while (running_.load()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

void Application::shutdown() {
    running_ = false;
    stopBackgroundThreads();
    saveConfiguration();
}

void Application::addQuoteCallback(std::function<void(const std::map<std::string, common::Quote>&)> callback) {
    std::unique_lock<std::shared_mutex> lock(dataLock_);
    quoteCallbacks_.push_back(callback);
}

void Application::addOpportunityCallback(std::function<void(const std::vector<strategy::StrategyOpportunity>&)> callback) {
    std::unique_lock<std::shared_mutex> lock(dataLock_);
    opportunityCallbacks_.push_back(callback);
}

void Application::addStatusCallback(std::function<void(const std::string&)> callback) {
    std::unique_lock<std::shared_mutex> lock(dataLock_);
    statusCallbacks_.push_back(callback);
}

std::map<std::string, common::Quote> Application::getCurrentQuotes() const {
    std::shared_lock<std::shared_mutex> lock(dataLock_);
    return currentQuotes_;
}

std::map<std::string, common::OptionChain> Application::getCurrentOptionChains() const {
    std::shared_lock<std::shared_mutex> lock(dataLock_);
    return currentOptionChains_;
}

std::vector<strategy::StrategyOpportunity> Application::getCurrentOpportunities() const {
    std::shared_lock<std::shared_mutex> lock(dataLock_);
    return currentOpportunities_;
}

std::future<broker::Result<common::OrderResponse>> Application::placeOrder(const common::OrderRequest& request) {
    return std::async(std::launch::async, [this, request]() -> broker::Result<common::OrderResponse> {
        if (!riskManager_ || !brokerManager_) {
            return std::string("System not initialized");
        }
        
        auto validationResult = riskManager_->validateOrder(request);
        bool isValid = validationResult.get();
        
        if (!isValid) {
            return std::string("Order rejected by risk management");
        }
        
        return brokerManager_->placeOrder(request).get();
    });
}

std::future<broker::Result<common::OrderResponse>> Application::cancelOrder(const std::string& accountId, 
                                                                           const std::string& orderId) {
    return std::async(std::launch::async, [this, accountId, orderId]() -> broker::Result<common::OrderResponse> {
        if (!brokerManager_) {
            return std::string("Broker manager not initialized");
        }
        
        return brokerManager_->cancelOrder(accountId, orderId).get();
    });
}

void Application::reloadWatchlist() {
    auto symbols = loadWatchlistFromFile("watchlist.txt");
    
    std::unique_lock<std::shared_mutex> lock(dataLock_);
    config_.watchedSymbols = symbols;
    currentQuotes_.clear();
    currentOptionChains_.clear();
    
    notifyStatusUpdate("Watchlist reloaded");
}

void Application::addSymbol(const std::string& symbol) {
    std::unique_lock<std::shared_mutex> lock(dataLock_);
    
    auto& symbols = config_.watchedSymbols;
    if (std::find(symbols.begin(), symbols.end(), symbol) == symbols.end()) {
        symbols.push_back(symbol);
        notifyStatusUpdate("Added symbol: " + symbol);
    }
}

void Application::removeSymbol(const std::string& symbol) {
    std::unique_lock<std::shared_mutex> lock(dataLock_);
    
    auto& symbols = config_.watchedSymbols;
    symbols.erase(std::remove(symbols.begin(), symbols.end(), symbol), symbols.end());
    
    currentQuotes_.erase(symbol);
    currentOptionChains_.erase(symbol);
    
    notifyStatusUpdate("Removed symbol: " + symbol);
}

std::vector<std::string> Application::getWatchedSymbols() const {
    std::shared_lock<std::shared_mutex> lock(dataLock_);
    return config_.watchedSymbols;
}

void Application::initializeComponents() {
    mathEngine_ = std::make_shared<math::MathEngine>();
    
    economicData_ = std::make_shared<data::EconomicDataManager>();
    
    portfolio_ = std::make_shared<portfolio::PortfolioManager>(mathEngine_);
    
    riskManager_ = std::make_shared<portfolio::RiskManager>(portfolio_, mathEngine_);
    
    strategyEngine_ = std::make_shared<strategy::StrategyEngine>(mathEngine_, economicData_);
    
    brokerManager_ = std::make_unique<broker::BrokerManager>();
    
    pluginLoader_ = std::make_unique<broker::PluginLoader>();
    
    notifyStatusUpdate("Components initialized");
}

void Application::loadConfiguration() {
    try {
        if (std::filesystem::exists(config_.configFile)) {
            config_ = ConfigManager::loadFromFile(config_.configFile);
        }
    } catch (const std::exception& e) {
        notifyStatusUpdate("Failed to load configuration: " + std::string(e.what()));
    }
}

void Application::saveConfiguration() const {
    try {
        ConfigManager::saveToFile(config_, config_.configFile);
    } catch (const std::exception& e) {
        
    }
}

void Application::loadPlugins() {
    try {
        auto& registry = broker::BrokerRegistry::getInstance();
        
        std::filesystem::path pluginDir("plugins");
        if (std::filesystem::exists(pluginDir)) {
            for (const auto& entry : std::filesystem::directory_iterator(pluginDir)) {
                if (entry.path().extension() == ".so" || entry.path().extension() == ".dll") {
                    try {
                        auto factory = pluginLoader_->loadPlugin(entry.path().string());
                        if (factory) {
                            registry.registerBroker(std::move(factory));
                        }
                    } catch (const std::exception& e) {
                        notifyStatusUpdate("Failed to load plugin: " + entry.path().string());
                    }
                }
            }
        }
        
        notifyStatusUpdate("Plugins loaded");
    } catch (const std::exception& e) {
        notifyStatusUpdate("Plugin loading failed: " + std::string(e.what()));
    }
}

void Application::initializeBrokers() {
    if (!brokerManager_) return;
    
    std::vector<std::future<void>> brokerFutures;
    
    for (const auto& brokerConfig : config_.brokers) {
        brokerFutures.emplace_back(std::async(std::launch::async, [this, brokerConfig]() {
            auto result = brokerManager_->addBroker(brokerConfig.brokerName, brokerConfig);
            auto success = result.get();
            
            if (broker::isSuccess(success)) {
                notifyStatusUpdate("Connected to broker: " + brokerConfig.brokerName);
            } else {
                notifyStatusUpdate("Failed to connect to broker: " + brokerConfig.brokerName + 
                                 " - " + broker::getError(success));
            }
        }));
    }
    
    for (auto& future : brokerFutures) {
        future.wait();
    }
}

void Application::startBackgroundThreads() {
    if (config_.enableBackgroundAnalysis) {
        threads_.emplace_back(&Application::dataCollectionLoop, this);
        threads_.emplace_back(&Application::analysisLoop, this);
    }
    
    if (config_.enableRiskManagement && riskManager_) {
        threads_.emplace_back(&Application::riskMonitoringLoop, this);
        riskManager_->startMonitoring();
    }
    
    if (portfolio_) {
        portfolio_->start();
    }
    
    notifyStatusUpdate("Background threads started");
}

void Application::stopBackgroundThreads() {
    for (auto& thread : threads_) {
        if (thread.joinable()) {
            thread.join();
        }
    }
    threads_.clear();
    
    if (riskManager_) {
        riskManager_->stopMonitoring();
    }
    
    if (portfolio_) {
        portfolio_->stop();
    }
}

void Application::dataCollectionLoop() {
    while (running_.load()) {
        try {
            updateQuotes();
            updateOptionChains();
            updatePortfolioAnalytics();
        } catch (const std::exception& e) {
            notifyStatusUpdate("Data collection error: " + std::string(e.what()));
        }
        
        std::this_thread::sleep_for(std::chrono::seconds(config_.updateIntervalSeconds));
    }
}

void Application::analysisLoop() {
    while (running_.load()) {
        try {
            runStrategyAnalysis();
        } catch (const std::exception& e) {
            notifyStatusUpdate("Analysis error: " + std::string(e.what()));
        }
        
        std::this_thread::sleep_for(std::chrono::seconds(30));
    }
}

void Application::riskMonitoringLoop() {
    while (running_.load()) {
        try {
            performRiskCheck();
        } catch (const std::exception& e) {
            notifyStatusUpdate("Risk monitoring error: " + std::string(e.what()));
        }
        
        std::this_thread::sleep_for(std::chrono::seconds(5));
    }
}

void Application::updateQuotes() {
    if (!brokerManager_) return;
    
    auto quotesResult = brokerManager_->getQuotes(config_.watchedSymbols);
    auto quotes = quotesResult.get();
    
    if (broker::isSuccess(quotes)) {
        std::map<std::string, common::Quote> quoteMap;
        for (const auto& quote : broker::getValue(quotes)) {
            quoteMap[quote.symbol] = quote;
        }
        
        {
            std::unique_lock<std::shared_mutex> lock(dataLock_);
            currentQuotes_ = quoteMap;
        }
        
        notifyQuoteUpdate(quoteMap);
    }
}

void Application::updateOptionChains() {
    if (!brokerManager_) return;
    
    std::map<std::string, common::OptionChain> chains;
    
    for (const auto& symbol : config_.watchedSymbols) {
        try {
            auto chainResult = brokerManager_->getOptionChain(symbol, "");
            auto chain = chainResult.get();
            
            if (broker::isSuccess(chain)) {
                chains[symbol] = broker::getValue(chain);
            }
        } catch (const std::exception& e) {
            
        }
    }
    
    {
        std::unique_lock<std::shared_mutex> lock(dataLock_);
        currentOptionChains_ = chains;
    }
}

void Application::updatePortfolioAnalytics() {
    if (!portfolio_ || !brokerManager_) return;
    
    try {
        std::shared_lock<std::shared_mutex> lock(dataLock_);
        portfolio_->updateAllPositions(currentQuotes_, currentOptionChains_);
    } catch (const std::exception& e) {
        
    }
}

void Application::runStrategyAnalysis() {
    if (!strategyEngine_) return;
    
    std::shared_lock<std::shared_mutex> lock(dataLock_);
    auto quotesSnapshot = currentQuotes_;
    auto chainsSnapshot = currentOptionChains_;
    lock.unlock();
    
    auto opportunitiesResult = strategyEngine_->scanOpportunities(chainsSnapshot, quotesSnapshot);
    auto opportunities = opportunitiesResult.get();
    
    {
        std::unique_lock<std::shared_mutex> updateLock(dataLock_);
        currentOpportunities_ = opportunities;
    }
    
    notifyOpportunityUpdate(opportunities);
}

void Application::performRiskCheck() {
    if (!riskManager_) return;
    
    try {
        auto assessmentResult = riskManager_->assessRisk();
        auto assessment = assessmentResult.get();
        
        if (!assessment.withinLimits) {
            std::string alertMessage = "Risk limits violated: ";
            for (const auto& violation : assessment.violations) {
                alertMessage += violation + "; ";
            }
            notifyStatusUpdate(alertMessage);
        }
    } catch (const std::exception& e) {
        
    }
}

void Application::notifyQuoteUpdate(const std::map<std::string, common::Quote>& quotes) {
    std::shared_lock<std::shared_mutex> lock(dataLock_);
    for (const auto& callback : quoteCallbacks_) {
        try {
            callback(quotes);
        } catch (const std::exception& e) {
            
        }
    }
}

void Application::notifyOpportunityUpdate(const std::vector<strategy::StrategyOpportunity>& opportunities) {
    std::shared_lock<std::shared_mutex> lock(dataLock_);
    for (const auto& callback : opportunityCallbacks_) {
        try {
            callback(opportunities);
        } catch (const std::exception& e) {
            
        }
    }
}