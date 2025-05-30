#pragma once

#include <string>
#include <optional>

namespace common {

enum class InstrumentType {
    EQUITY, 
    OPTION, 
    FUTURE, 
    FOREX, 
    CRYPTO
};

class Instrument {
private:
    std::string symbol_;
    std::string name_;
    InstrumentType type_;
    std::string exchange_;
    std::string currency_;
    std::optional<std::string> underlying_;
    std::optional<double> strike_;
    std::optional<std::string> expiration_;
    std::optional<bool> isCall_;
    double multiplier_;
    double tickSize_;

public:
    Instrument() : type_(InstrumentType::EQUITY), currency_("USD"), multiplier_(1.0), tickSize_(0.01) {}
    
    Instrument(const std::string& symbol, InstrumentType type, const std::string& exchange)
        : symbol_(symbol), type_(type), exchange_(exchange), currency_("USD"), 
          multiplier_(1.0), tickSize_(0.01) {}

    const std::string& getSymbol() const { return symbol_; }
    const std::string& getName() const { return name_; }
    InstrumentType getType() const { return type_; }
    const std::string& getExchange() const { return exchange_; }
    const std::string& getCurrency() const { return currency_; }
    const std::optional<std::string>& getUnderlying() const { return underlying_; }
    const std::optional<double>& getStrike() const { return strike_; }
    const std::optional<std::string>& getExpiration() const { return expiration_; }
    const std::optional<bool>& getIsCall() const { return isCall_; }
    double getMultiplier() const { return multiplier_; }
    double getTickSize() const { return tickSize_; }

    void setSymbol(const std::string& symbol) { symbol_ = symbol; }
    void setName(const std::string& name) { name_ = name; }
    void setType(InstrumentType type) { type_ = type; }
    void setExchange(const std::string& exchange) { exchange_ = exchange; }
    void setCurrency(const std::string& currency) { currency_ = currency; }
    void setUnderlying(const std::string& underlying) { underlying_ = underlying; }
    void setStrike(double strike) { strike_ = strike; }
    void setExpiration(const std::string& expiration) { expiration_ = expiration; }
    void setIsCall(bool isCall) { isCall_ = isCall; }
    void setMultiplier(double multiplier) { multiplier_ = multiplier; }
    void setTickSize(double tickSize) { tickSize_ = tickSize; }

    bool isOption() const { return type_ == InstrumentType::OPTION; }
    bool isEquity() const { return type_ == InstrumentType::EQUITY; }
    bool isFuture() const { return type_ == InstrumentType::FUTURE; }
};

}