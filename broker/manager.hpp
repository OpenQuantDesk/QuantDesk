#pragma once

#include "broker/result.hpp"
#include "common/instrument.hpp"
#include "interfaces/broker.hpp"
#include <functional>
#include <future>
#include <map>
#include <memory>
#include <shared_mutex>
#include <string>
#include <vector>

namespace broker {

    class BrokerManager {
    private:
        std::map<std::string, std::unique_ptr<IBroker>> brokers_;
        std::string primaryBroker_;
        mutable std::shared_mutex brokersLock_;

        std::vector<std::function<void(const std::string&, bool)>>
            connectionCallbacks_;
    public:
        BrokerManager() = default;
        ~BrokerManager();

        BrokerManager(const BrokerManager&) = delete;
        BrokerManager& operator=(const BrokerManager&) = delete;

        std::future<Result<bool>> addBroker(const std::string& name,
                                            const common::BrokerConfig& config);
        void removeBroker(const std::string& name);

        IBroker* getBroker(const std::string& name = "") const;
        std::vector<std::string> getActiveBrokers() const;
        std::vector<std::string> getConnectedBrokers() const;

        void setPrimaryBroker(const std::string& name);
        std::string getPrimaryBroker() const { return primaryBroker_; }

        IBroker* getOptimalBroker(const common::Instrument& instrument) const;

        std::future<Result<std::vector<common::Quote>>>
        getQuotes(const std::vector<std::string>& symbols) const;
        std::future<Result<common::Quote>>
        getQuote(const std::string& symbol) const;
        std::future<Result<common::OptionChain>>
        getOptionChain(const std::string& underlying,
                       const std::string& expiration) const;
        std::future<Result<common::OrderResponse>>
        placeOrder(const common::OrderRequest& request) const;
        std::future<Result<common::OrderResponse>>
        cancelOrder(const std::string& accountId,
                    const std::string& orderId) const;
        std::future<Result<std::vector<common::Account>>> getAccounts() const;
        std::future<Result<std::vector<common::Position>>>
        getPositions(const std::string& accountId) const;

        void
        subscribeQuotes(const std::vector<std::string>& symbols,
                        std::function<void(const common::Quote&)> callback);
        void subscribeExecutions(
            std::function<void(const common::ExecutionReport&)> callback);

        void addConnectionCallback(
            std::function<void(const std::string&, bool)> callback);

        size_t getBrokerCount() const;
        size_t getConnectedBrokerCount() const;
        bool hasConnectedBrokers() const;
    private:
        void notifyConnectionChange(const std::string& brokerName,
                                    bool connected);
    };

} // namespace broker