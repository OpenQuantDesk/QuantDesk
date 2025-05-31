/*
 * Filename: money.hpp
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