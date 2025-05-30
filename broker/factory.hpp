#pragma once

#include "broker/broker_interface.hpp"
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