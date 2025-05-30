#pragma once

#include "common/types.hpp"
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <future>
#include <variant>
#include <shared_mutex>
#include <atomic>

namespace broker {

template<typename T>
using Result = std::variant<T, std::string>;

template<typename T>
bool isSuccess(const Result<T>& result) {
    return std::holds_alternative<T>(result);
}

template<typename T>
const T& getValue(const Result<T>& result) {
    return std::get<T>(result);
}

template<typename T>
const std::string& getError(const Result<T>& result) {
    return std::get<std::string>(result);
}

class IBroker {
public:
    virtual ~IBroker() = default;
    
    virtual std::future<Result<bool>> connect() = 0;
    virtual void disconnect() = 0;
    virtual bool isConnected() const = 0;
    virtual std::string getName() const = 0;
    
    virtual std::future<Result<std::vector<common::Account>>> getAccounts() = 0;
    virtual std::future<Result<std::vector<common::Position>>> getPositions(const std::string& accountId) = 0;
    virtual std::future<Result<common::Quote>> getQuote(const std::string& symbol) = 0;
    virtual std::future<Result<std::vector<common::Quote>>> getQuotes(const std::vector<std::string>& symbols) = 0;
    virtual std::future<Result<common::OptionChain>> getOptionChain(const std::string& underlying, const std::string& expiration) = 0;
    virtual std::future<Result<common::OrderResponse>> placeOrder(const common::OrderRequest& request) = 0;
    virtual std::future<Result<common::OrderResponse>> cancelOrder(const std::string& accountId, const std::string& orderId) = 0;
    
    virtual bool supportsStreaming() const { return false; }
    virtual void subscribeQuotes(const std::vector<std::string>& symbols, 
                               std::function<void(const common::Quote&)> callback) {}
    virtual void subscribeExecutions(std::function<void(const common::ExecutionReport&)> callback) {}
    
    virtual bool supportsInstrumentType(common::InstrumentType type) const = 0;
    virtual std::vector<std::string> getSupportedExchanges() const = 0;
};

class IBrokerFactory {
public:
    virtual ~IBrokerFactory() = default;
    virtual std::unique_ptr<IBroker> createBroker(const common::BrokerConfig& config) = 0;
    virtual std::string getName() const = 0;
    virtual std::vector<std::string> getRequiredParameters() const = 0;
};

class BrokerRegistry {
private:
    std::map<std::string, std::unique_ptr<IBrokerFactory>> factories_;
    mutable std::shared_mutex factoriesLock_;
    
    static std::unique_ptr<BrokerRegistry> instance_;
    static std::once_flag instanceFlag_;
    
    BrokerRegistry() = default;
    
public:
    ~BrokerRegistry() = default;
    
    BrokerRegistry(const BrokerRegistry&) = delete;
    BrokerRegistry& operator=(const BrokerRegistry&) = delete;
    
    static BrokerRegistry& getInstance();
    
    void registerBroker(std::unique_ptr<IBrokerFactory> factory);
    void unregisterBroker(const std::string& name);
    
    std::vector<std::string> getAvailableBrokers() const;
    std::unique_ptr<IBroker> createBroker(const common::BrokerConfig& config) const;
    IBrokerFactory* getFactory(const std::string& brokerName) const;
    
    bool isRegistered(const std::string& brokerName) const;
    std::vector<std::string> getRequiredParameters(const std::string& brokerName) const;
};

class BrokerManager {
private:
    std::map<std::string, std::unique_ptr<IBroker>> brokers_;
    std::string primaryBroker_;
    mutable std::shared_mutex brokersLock_;
    
    std::vector<std::function<void(const std::string&, bool)>> connectionCallbacks_;
    
public:
    BrokerManager() = default;
    ~BrokerManager();
    
    BrokerManager(const BrokerManager&) = delete;
    BrokerManager& operator=(const BrokerManager&) = delete;
    
    std::future<Result<bool>> addBroker(const std::string& name, const common::BrokerConfig& config);
    void removeBroker(const std::string& name);
    
    IBroker* getBroker(const std::string& name = "") const;
    std::vector<std::string> getActiveBrokers() const;
    std::vector<std::string> getConnectedBrokers() const;
    
    void setPrimaryBroker(const std::string& name);
    std::string getPrimaryBroker() const { return primaryBroker_; }
    
    IBroker* getOptimalBroker(const common::Instrument& instrument) const;
    
    std::future<Result<std::vector<common::Quote>>> getQuotes(const std::vector<std::string>& symbols) const;
    std::future<Result<common::Quote>> getQuote(const std::string& symbol) const;
    std::future<Result<common::OptionChain>> getOptionChain(const std::string& underlying, 
                                                           const std::string& expiration) const;
    std::future<Result<common::OrderResponse>> placeOrder(const common::OrderRequest& request) const;
    std::future<Result<common::OrderResponse>> cancelOrder(const std::string& accountId, 
                                                          const std::string& orderId) const;
    std::future<Result<std::vector<common::Account>>> getAccounts() const;
    std::future<Result<std::vector<common::Position>>> getPositions(const std::string& accountId) const;
    
    void subscribeQuotes(const std::vector<std::string>& symbols,
                        std::function<void(const common::Quote&)> callback);
    void subscribeExecutions(std::function<void(const common::ExecutionReport&)> callback);
    
    void addConnectionCallback(std::function<void(const std::string&, bool)> callback);
    
private:
    void notifyConnectionChange(const std::string& brokerName, bool connected);
};

class PluginLoader {
private:
    std::vector<void*> loadedLibraries_;
    mutable std::shared_mutex librariesLock_;
    
public:
    ~PluginLoader();
    
    std::unique_ptr<IBrokerFactory> loadPlugin(const std::string& libraryPath);
    void unloadPlugin(const std::string& libraryPath);
    void unloadAllPlugins();
};

}

extern "C" {
    broker::IBrokerFactory* createBrokerFactory();
}