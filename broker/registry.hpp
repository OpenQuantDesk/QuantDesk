/*
 * Filename: registry.hpp
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

#pragma once

#include "interfaces/broker.hpp"
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <shared_mutex>
#include <functional>

namespace broker {

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
    
    std::future<Result<std::vector<common::Quote>>> getQu