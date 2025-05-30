#pragma once

#include <string>
#include <vector>
#include <optional>
#include <chrono>
#include <map>
#include <atomic>

namespace common {

enum class InstrumentType {
    EQUITY, OPTION, FUTURE, FOREX, CRYPTO
};

enum class OrderSide {
    BUY, SELL, BUY_TO_OPEN, BUY_TO_CLOSE, SELL_TO_OPEN, SELL_TO_CLOSE
};

enum class OrderType {
    MARKET, LIMIT, STOP, STOP_LIMIT
};

enum class OrderStatus {
    PENDING, OPEN, FILLED, CANCELED, REJECTED
};

enum class TimeInForce {
    DAY, GTC, IOC, FOK
};

struct Money {
    double amount = 0.0;
    std::string currency = "USD";
    
    Money() = default;
    Money(double amt, const std::string& curr = "USD") : amount(amt), currency(curr) {}
    
    Money operator+(const Money& other) const;
    Money operator-(const Money& other) const;
    Money operator*(double multiplier) const;
};

struct Greeks {
    double delta = 0.0;
    double gamma = 0.0;
    double theta = 0.0;
    double vega = 0.0;
    double rho = 0.0;
    double impliedVol = 0.0;
    double price = 0.0;
};

struct ExtendedGreeks {
    double price = 0.0;
    double delta = 0.0;
    double gamma = 0.0;
    double theta = 0.0;
    double vega = 0.0;
    double rho = 0.0;
    
    double volga = 0.0;
    double vanna = 0.0;
    double charm = 0.0;
    double speed = 0.0;
    double zomma = 0.0;
    double color = 0.0;
    
    double dollarDelta = 0.0;
    double dollarGamma = 0.0;
    double pinRisk = 0.0;
};

struct Instrument {
    std::string symbol;
    std::string name;
    InstrumentType type;
    std::string exchange;
    std::string currency = "USD";
    
    std::optional<std::string> underlying;
    std::optional<double> strike;
    std::optional<std::string> expiration;
    std::optional<bool> isCall;
    
    double multiplier = 1.0;
    double tickSize = 0.01;
};

struct Quote {
    std::string symbol;
    std::optional<double> bid;
    std::optional<double> ask;
    std::optional<double> last;
    std::optional<double> open;
    std::optional<double> high;
    std::optional<double> low;
    std::optional<double> close;
    long volume = 0;
    std::chrono::system_clock::time_point timestamp;
    std::optional<Greeks> greeks;
    std::optional<long> openInterest;
};

struct OptionChain {
    std::string underlying;
    std::string expiration;
    std::vector<Quote> calls;
    std::vector<Quote> puts;
    std::chrono::system_clock::time_point timestamp;
};

struct Account {
    std::string id;
    std::string name;
    std::string type;
    Money totalEquity;
    Money buyingPower;
    Money cashBalance;
    bool dayTradingEnabled = false;
    int optionLevel = 0;
};

struct Position {
    std::string symbol;
    InstrumentType instrumentType;
    double quantity = 0.0;
    Money avgCost;
    Money marketValue;
    Money unrealizedPnL;
    std::chrono::system_clock::time_point acquired;
};

struct OrderRequest {
    std::string accountId;
    Instrument instrument;
    OrderSide side;
    OrderType type;
    double quantity;
    TimeInForce timeInForce = TimeInForce::DAY;
    
    std::optional<double> limitPrice;
    std::optional<double> stopPrice;
    
    std::vector<OrderRequest> legs;
    std::map<std::string, std::string> metadata;
};

struct OrderResponse {
    std::string orderId;
    std::string accountId;
    OrderStatus status;
    std::string statusMessage;
    std::chrono::system_clock::time_point timestamp;
    
    double filledQuantity = 0.0;
    double remainingQuantity = 0.0;
    std::optional<double> avgFillPrice;
    Money commission;
    Money fees;
};

struct ExecutionReport {
    std::string orderId;
    std::string executionId;
    double quantity;
    double price;
    std::chrono::system_clock::time_point timestamp;
    std::string venue;
};

struct BrokerConfig {
    std::string brokerName;
    std::map<std::string, std::string> parameters;
    bool sandboxMode = true;
};

}