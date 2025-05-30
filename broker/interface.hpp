#pragma once

#include "common/types.hpp"
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <future>
#include <variant>

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

}

extern "C" {
    broker::IBrokerFactory* createBrokerFactory();
}