#pragma once

#include <string>

namespace common {

class Money {
private:
    double amount_;
    std::string currency_;

public:
    Money() : amount_(0.0), currency_("USD") {}
    Money(double amount, const std::string& currency = "USD") 
        : amount_(amount), currency_(currency) {}

    double getAmount() const { return amount_; }
    const std::string& getCurrency() const { return currency_; }
    
    void setAmount(double amount) { amount_ = amount; }
    void setCurrency(const std::string& currency) { currency_ = currency; }

    Money operator+(const Money& other) const;
    Money operator-(const Money& other) const;
    Money operator*(double multiplier) const;
    Money& operator+=(const Money& other);
    Money& operator-=(const Money& other);
    Money& operator*=(double multiplier);
    
    bool operator==(const Money& other) const;
    bool operator!=(const Money& other) const;
    bool operator<(const Money& other) const;
    bool operator>(const Money& other) const;
};

}