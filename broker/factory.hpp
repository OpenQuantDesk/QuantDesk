/*
 * Filename: factory.hpp
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
#include "common/instrument.hpp"
#include <string>
#include <vector>
#include <memory>

namespace broker {

class IBrokerFactory {
public:
    virtual ~IBrokerFactory() = default;
    
    virtual std::unique_ptr<IBroker> createBroker(const common::BrokerConfig& config) = 0;
    virtual std::string getName() const = 0;
    virtual std::vector<std::string> getRequiredParameters() const = 0;
    virtual std::string getDescription() const { return ""; }
    virtual std::string getVersion() const { return "1.0"; }
    virtual bool isAvailable() const { return true; }
    virtual std::vector<common::InstrumentType> getSupportedInstruments() const = 0;
    virtual std::vector<std::string> getSupportedExchanges() const = 0;
};

}

extern "C" {
    broker::IBrokerFactory* createBrokerFactory();
}