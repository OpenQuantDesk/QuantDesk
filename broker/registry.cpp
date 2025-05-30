#include "broker/registry.hpp"
#include <mutex>

namespace broker {

std::unique_ptr<BrokerRegistry> BrokerRegistry::instance_;
std::once_flag BrokerRegistry::instanceFlag_;

BrokerRegistry& BrokerRegistry::getInstance() {
    std::call_once(instanceFlag_, []() {
        instance_ = std::unique_ptr<BrokerRegistry>(new BrokerRegistry());
    });
    return *instance_;
}

void BrokerRegistry::registerBroker(std::unique_ptr<IBrokerFactory> factory) {
    std::unique_lock<std::shared_mutex> lock(factoriesLock_);
    std::string name = factory->getName();
    factories_[name] = std::move(factory);
}

void BrokerRegistry::unregisterBroker(const std::string& name) {
    std::unique_lock<std::shared_mutex> lock(factoriesLock_);
    factories_.erase(name);
}

std::vector<std::string> BrokerRegistry::getAvailableBrokers() const {
    std::shared_lock<std::shared_mutex> lock(factoriesLock_);
    std::vector<std::string> brokers;
    brokers.reserve(factories_.size());
    
    for (const auto& [name, factory] : factories_) {
        brokers.push_back(name);
    }
    
    return brokers;
}

std::unique_ptr<IBroker> BrokerRegistry::createBroker(const common::BrokerConfig& config) const {
    std::shared_lock<std::shared_mutex> lock(factoriesLock_);
    
    auto it = factories_.find(config.brokerName);
    if (it == factories_.end()) {
        return nullptr;
    }
    
    return it->second->createBroker(config);
}

IBrokerFactory* BrokerRegistry::getFactory(const std::string& brokerName) const {
    std::shared_lock<std::shared_mutex> lock(factoriesLock_);
    
    auto it = factories_.find(brokerName);
    return (it != factories_.end()) ? it->second.get() : nullptr;
}

bool BrokerRegistry::isRegistered(const std::string& brokerName) const {
    std::shared_lock<std::shared_mutex> lock(factoriesLock_);
    return factories_.find(brokerName) != factories_.end();
}

std::vector<std::string> BrokerRegistry::getRequiredParameters(const std::string& brokerName) const {
    std::shared_lock<std::shared_mutex> lock(factoriesLock_);
    
    auto it = factories_.find(brokerName);
    if (it != factories_.end()) {
        return it->second->getRequiredParameters();
    }
    
    return {};
}

BrokerManager::~BrokerManager() {
    std::unique_lock<std::shared_mutex> lock(brokersLock_);
    brokers_.clear();
}

std::future<Result<bool>> BrokerManager::addBroker(const std::string& name, const common::BrokerConfig& config) {
    return std::async(std::launch::async, [this, name, config]() -> Result<bool> {
        try {
            auto& registry = BrokerRegistry::getInstance();
            auto broker = registry.createBroker(config);
            
            if (!broker) {
                return std::string("Failed to create broker: " + config.brokerName);
            }
            
            auto connectResult = broker->connect().get();
            if (!isSuccess(connectResult)) {
                return getError(connectResult);
            }
            
            {
                std::unique_lock<std::shared_mutex> lock(brokersLock_);
                brokers_[name] = std::move(broker);
                
                if (primaryBroker_.empty()) {
                    primaryBroker_ = name;
                }
            }
            
            notifyConnectionChange(name, true);
            return true;
            
        } catch (const std::exception& e) {
            return std::string("Exception adding broker: " + std::string(e.what()));
        }
    });
}

void BrokerManager::removeBroker(const std::string& name) {
    std::unique_lock<std::shared_mutex> lock(brokersLock_);
    
    auto it = brokers_.find(name);
    if (it != brokers_.end()) {
        it->second->disconnect();
        brokers_.erase(it);
        
        if (primaryBroker_ == name) {
            primaryBroker_ = brokers_.empty() ? "" : brokers_.begin()->first;
        }
        
        notifyConnectionChange(name, false);
    }
}

IBroker* BrokerManager::getBroker(const std::string& name) const {
    std::shared_lock<std::shared_mutex> lock(brokersLock_);
    
    if (name.empty()) {
        if (!primaryBroker_.empty()) {
            auto it = brokers_.find(primaryBroker_);
            return (it != brokers_.end()) ? it->second.get() : nullptr;
        }
        return brokers_.empty() ? nullptr : brokers_.begin()->second.get();
    }
    
    auto it = brokers_.find(name);
    return (it != brokers_.end()) ? it->second.get() : nullptr;
}

std::vector<std::string> BrokerManager::getActiveBrokers() const {
    std::shared_lock<std::shared_mutex> lock(brokersLock_);
    std::vector<std::string> brokers;
    brokers.reserve(brokers_.size());
    
    for (const auto& [name, broker] : brokers_) {
        brokers.push_back(name);
    }
    
    return brokers;
}

std::vector<std::string> BrokerManager::getConnectedBrokers() const {
    std::shared_lock<std::shared_mutex> lock(brokersLock_);
    std::vector<std::string> connected;
    
    for (const auto& [name, broker] : brokers_) {
        if (broker->isConnected()) {
            connected.push_back(name);
        }
    }
    
    return connected;
}

void BrokerManager::setPrimaryBroker(const std::string& name) {
    std::unique_lock<std::shared_mutex> lock(brokersLock_);
    
    if (brokers_.find(name) != brokers_.end()) {
        primaryBroker_ = name;
    }
}

IBroker* BrokerManager::getOptimalBroker(const common::Instrument& instrument) const {
    std::shared_lock<std::shared_mutex> lock(brokersLock_);
    
    for (const auto& [name, broker] : brokers_) {
        if (broker->isConnected() && broker->supportsInstrumentType(instrument.type)) {
            auto exchanges = broker->getSupportedExchanges();
            if (std::find(exchanges.begin(), exchanges.end(), instrument.exchange) != exchanges.end()) {
                return broker.get();
            }
        }
    }
    
    return getBroker();
}

std::future<Result<std::vector<common::Quote>>> BrokerManager::getQuotes(const std::vector<std::string>& symbols) const {
    return std::async(std::launch::async, [this, symbols]() -> Result<std::vector<common::Quote>> {
        auto broker = getBroker();
        if (!broker) {
            return std::string("No broker available");
        }
        
        if (!broker->isConnected()) {
            return std::string("Broker not connected");
        }
        
        return broker->getQuotes(symbols).get();
    });
}

std::future<Result<common::Quote>> BrokerManager::getQuote(const std::string& symbol) const {
    return std::async(std::launch::async, [this, symbol]() -> Result<common::Quote> {
        auto broker = getBroker();
        if (!broker) {
            return std::string("No broker available");
        }
        
        if (!broker->isConnected()) {
            return std::string("Broker not connected");
        }
        
        return broker->getQuote(symbol).get();
    });
}

std::future<Result<common::OptionChain>> BrokerManager::getOptionChain(const std::string& underlying, 
                                                                      const std::string& expiration) const {
    return std::async(std::launch::async, [this, underlying, expiration]() -> Result<common::OptionChain> {
        auto broker = getBroker();
        if (!broker) {
            return std::string("No broker available");
        }
        
        if (!broker->isConnected()) {
            return std::string("Broker not connected");
        }
        
        return broker->getOptionChain(underlying, expiration).get();
    });
}

std::future<Result<common::OrderResponse>> BrokerManager::placeOrder(const common::OrderRequest& request) const {
    return std::async(std::launch::async, [this, request]() -> Result<common::OrderResponse> {
        auto broker = getOptimalBroker(request.instrument);
        if (!broker) {
            return std::string("No suitable broker available");
        }
        
        if (!broker->isConnected()) {
            return std::string("Broker not connected");
        }
        
        return broker->placeOrder(request).get();
    });
}

std::future<Result<common::OrderResponse>> BrokerManager::cancelOrder(const std::string& accountId, 
                                                                     const std::string& orderId) const {
    return std::async(std::launch::async, [this, accountId, orderId]() -> Result<common::OrderResponse> {
        auto broker = getBroker();
        if (!broker) {
            return std::string("No broker available");
        }
        
        if (!broker->isConnected()) {
            return std::string("Broker not connected");
        }
        
        return broker->cancelOrder(accountId, orderId).get();
    });
}

std::future<Result<std::vector<common::Account>>> BrokerManager::getAccounts() const {
    return std::async(std::launch::async, [this]() -> Result<std::vector<common::Account>> {
        auto broker = getBroker();
        if (!broker) {
            return std::string("No broker available");
        }
        
        if (!broker->isConnected()) {
            return std::string("Broker not connected");
        }
        
        return broker->getAccounts().get();
    });
}

std::future<Result<std::vector<common::Position>>> BrokerManager::getPositions(const std::string& accountId) const {
    return std::async(std::launch::async, [this, accountId]() -> Result<std::vector<common::Position>> {
        auto broker = getBroker();
        if (!broker) {
            return std::string("No broker available");
        }
        
        if (!broker->isConnected()) {
            return std::string("Broker not connected");
        }
        
        return broker->getPositions(accountId).get();
    });
}

void BrokerManager::subscribeQuotes(const std::vector<std::string>& symbols,
                                   std::function<void(const common::Quote&)> callback) {
    std::shared_lock<std::shared_mutex> lock(brokersLock_);
    
    for (const auto& [name, broker] : brokers_) {
        if (broker->isConnected() && broker->supportsStreaming()) {
            broker->subscribeQuotes(symbols, callback);
            break;
        }
    }
}

void BrokerManager::subscribeExecutions(std::function<void(const common::ExecutionReport&)> callback) {
    std::shared_lock<std::shared_mutex> lock(brokersLock_);
    
    for (const auto& [name, broker] : brokers_) {
        if (broker->isConnected() && broker->supportsStreaming()) {
            broker->subscribeExecutions(callback);
        }
    }
}

void BrokerManager::addConnectionCallback(std::function<void(const std::string&, bool)> callback) {
    std::unique_lock<std::shared_mutex> lock(brokersLock_);
    connectionCallbacks_.push_back(callback);
}

void BrokerManager::notifyConnectionChange(const std::string& brokerName, bool connected) {
    for (const auto& callback : connectionCallbacks_) {
        try {
            callback(brokerName, connected);
        } catch (const std::exception& e) {
            
        }
    }
}

}