#pragma once

#include "greeks.hpp"
#include <string>
#include <optional>
#include <chrono>

namespace common {

class Quote {
private:
    std::string symbol_;
    std::optional<double> bid_;
    std::optional<double> ask_;
    std::optional<double> last_;
    std::optional<double> open_;
    std::optional<double> high_;
    std::optional<double> low_;
    std::optional<double> close_;
    long volume_;
    std::chrono::system_clock::time_point timestamp_;
    std::optional<Greeks> greeks_;
    std::optional<long> openInterest_;

public:
    Quote() : volume_(0), timestamp_(std::chrono::system_clock::now()) {}
    
    explicit Quote(const std::string& symbol) 
        : symbol_(symbol), volume_(0), timestamp_(std::chrono::system_clock::now()) {}

    const std::string& getSymbol() const { return symbol_; }
    const std::optional<double>& getBid() const { return bid_; }
    const std::optional<double>& getAsk() const { return ask_; }
    const std::optional<double>& getLast() const { return last_; }
    const std::optional<double>& getOpen() const { return open_; }
    const std::optional<double>& getHigh() const { return high_; }
    const std::optional<double>& getLow() const { return low_; }
    const std::optional<double>& getClose() const { return close_; }
    long getVolume() const { return volume_; }
    const std::chrono::system_clock::time_point& getTimestamp() const { return timestamp_; }
    const std::optional<Greeks>& getGreeks() const { return greeks_; }
    const std::optional<long>& getOpenInterest() const { return openInterest_; }

    void setSymbol(const std::string& symbol) { symbol_ = symbol; }
    void setBid(double bid) { bid_ = bid; }
    void setAsk(double ask) { ask_ = ask; }
    void setLast(double last) { last_ = last; }
    void setOpen(double open) { open_ = open; }
    void setHigh(double high) { high_ = high; }
    void setLow(double low) { low_ = low; }
    void setClose(double close) { close_ = close; }
    void setVolume(long volume) { volume_ = volume; }
    void setTimestamp(const std::chrono::system_clock::time_point& timestamp) { timestamp_ = timestamp; }
    void setGreeks(const Greeks& greeks) { greeks_ = greeks; }
    void setOpenInterest(long openInterest) { openInterest_ = openInterest; }

    double getMidPrice() const;
    double getSpread() const;
    bool isValid() const;
    bool isStale(std::chrono::seconds threshold = std::chrono::seconds(60)) const;
};

}