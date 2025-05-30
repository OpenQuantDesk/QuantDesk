#pragma once

#include "common/instrument.hpp"
#include <string>
#include <vector>
#include <optional>
#include <map>

namespace common {

enum class OrderSide {
    BUY, 
    SELL, 
    BUY_TO_OPEN, 
    BUY_TO_CLOSE, 
    SELL_TO_OPEN, 
    SELL_TO_CLOSE
};

enum class OrderType {
    MARKET, 
    LIMIT, 
    STOP, 
    STOP_LIMIT
};

enum class TimeInForce {
    DAY, 
    GTC, 
    IOC, 
    FOK
};

class OrderRequest {
private:
    std::string accountId_;
    Instrument instrument_;
    OrderSide side_;
    OrderType type_;
    double quantity_;
    TimeInForce timeInForce_;
    std::optional<double> limitPrice_;
    std::optional<double> stopPrice_;
    std::vector<OrderRequest> legs_;
    std::map<std::string, std::string> metadata_;

public:
    OrderRequest() : side_(OrderSide::BUY), type_(OrderType::MARKET), 
                    quantity_(0.0), timeInForce_(TimeInForce::DAY) {}

    OrderRequest(const std::string& accountId, const Instrument& instrument, 
                OrderSide side, OrderType type, double quantity)
        : accountId_(accountId), instrument_(instrument), side_(side), 
          type_(type), quantity_(quantity), timeInForce_(TimeInForce::DAY) {}

    const std::string& getAccountId() const { return accountId_; }
    const Instrument& getInstrument() const { return instrument_; }
    OrderSide getSide() const { return side_; }
    OrderType getType() const { return type_; }
    double getQuantity() const { return quantity_; }
    TimeInForce getTimeInForce() const { return timeInForce_; }
    const std::optional<double>& getLimitPrice() const { return limitPrice_; }
    const std::optional<double>& getStopPrice() const { return stopPrice_; }
    const std::vector<OrderRequest>& getLegs() const { return legs_; }
    const std::map<std::string, std::string>& getMetadata() const { return metadata_; }

    void setAccountId(const std::string& accountId) { accountId_ = accountId; }
    void setInstrument(const Instrument& instrument) { instrument_ = instrument; }
    void setSide(OrderSide side) { side_ = side; }
    void setType(OrderType type) { type_ = type; }
    void setQuantity(double quantity) { quantity_ = quantity; }
    void setTimeInForce(TimeInForce timeInForce) { timeInForce_ = timeInForce; }
    void setLimitPrice(double limitPrice) { limitPrice_ = limitPrice; }
    void setStopPrice(double stopPrice) { stopPrice_ = stopPrice; }
    void setLegs(const std::vector<OrderRequest>& legs) { legs_ = legs; }
    void setMetadata(const std::map<std::string, std::string>& metadata) { metadata_ = metadata; }

    void addLeg(const OrderRequest& leg) { legs_.push_back(leg); }
    void addMetadata(const std::string& key, const std::string& value) { metadata_[key] = value; }
    
    bool isMultiLeg() const { return !legs_.empty(); }
    bool isMarketOrder() const { return type_ == OrderType::MARKET; }
    bool isLimitOrder() const { return type_ == OrderType::LIMIT; }
    bool isStopOrder() const { return type_ == OrderType::STOP || type_ == OrderType::STOP_LIMIT; }
    
    double getEstimatedNotional() const;
    bool isValid() const;
};

}