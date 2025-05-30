#pragma once

#include "money.hpp"
#include <string>
#include <optional>
#include <chrono>

namespace common {

enum class OrderStatus {
    PENDING, 
    OPEN, 
    PARTIALLY_FILLED,
    FILLED, 
    CANCELED, 
    REJECTED,
    EXPIRED
};

class OrderResponse {
private:
    std::string orderId_;
    std::string accountId_;
    OrderStatus status_;
    std::string statusMessage_;
    std::chrono::system_clock::time_point timestamp_;
    double filledQuantity_;
    double remainingQuantity_;
    std::optional<double> avgFillPrice_;
    Money commission_;
    Money fees_;

public:
    OrderResponse() : status_(OrderStatus::PENDING), 
                     timestamp_(std::chrono::system_clock::now()),
                     filledQuantity_(0.0), remainingQuantity_(0.0) {}

    OrderResponse(const std::string& orderId, const std::string& accountId, OrderStatus status)
        : orderId_(orderId), accountId_(accountId), status_(status),
          timestamp_(std::chrono::system_clock::now()),
          filledQuantity_(0.0), remainingQuantity_(0.0) {}

    const std::string& getOrderId() const { return orderId_; }
    const std::string& getAccountId() const { return accountId_; }
    OrderStatus getStatus() const { return status_; }
    const std::string& getStatusMessage() const { return statusMessage_; }
    const std::chrono::system_clock::time_point& getTimestamp() const { return timestamp_; }
    double getFilledQuantity() const { return filledQuantity_; }
    double getRemainingQuantity() const { return remainingQuantity_; }
    const std::optional<double>& getAvgFillPrice() const { return avgFillPrice_; }
    const Money& getCommission() const { return commission_; }
    const Money& getFees() const { return fees_; }

    void setOrderId(const std::string& orderId) { orderId_ = orderId; }
    void setAccountId(const std::string& accountId) { accountId_ = accountId; }
    void setStatus(OrderStatus status) { status_ = status; }
    void setStatusMessage(const std::string& statusMessage) { statusMessage_ = statusMessage; }
    void setTimestamp(const std::chrono::system_clock::time_point& timestamp) { timestamp_ = timestamp; }
    void setFilledQuantity(double filledQuantity) { filledQuantity_ = filledQuantity; }
    void setRemainingQuantity(double remainingQuantity) { remainingQuantity_ = remainingQuantity; }
    void setAvgFillPrice(double avgFillPrice) { avgFillPrice_ = avgFillPrice; }
    void setCommission(const Money& commission) { commission_ = commission; }
    void setFees(const Money& fees) { fees_ = fees; }

    bool isFilled() const { return status_ == OrderStatus::FILLED; }
    bool isPartiallyFilled() const { return status_ == OrderStatus::PARTIALLY_FILLED; }
    bool isActive() const { return status_ == OrderStatus::OPEN || status_ == OrderStatus::PARTIALLY_FILLED; }
    bool isCanceled() const { return status_ == OrderStatus::CANCELED; }
    bool isRejected() const { return status_ == OrderStatus::REJECTED; }
    bool isFinal() const { return isFilled() || isCanceled() || isRejected(); }
    
    double getFillPercentage() const;
    Money getTotalCost() const;
};

}