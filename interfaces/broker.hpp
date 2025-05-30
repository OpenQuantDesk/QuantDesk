#pragma once

#include "common/quote.hpp"
#include "common/option_chain.hpp"
#include "common/order_request.hpp"
#include "common/order_response.hpp"
#include "common/instrument.hpp"
#include "broker/result.hpp"
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <future>

namespace broker {

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
    
    virtual std::string getApiVersion() const { return "1.0"; }
    virtual bool isMarketDataEnabled() const { return true; }
    virtual bool isTradingEnabled() const { return true; }
};

}