/*
 * Filename: types.cpp
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

#include "types.hpp"

namespace common {

Money Money::operator+(const Money& other) const {
    if (currency != other.currency) {
        // For simplicity, assume USD conversion rate of 1.0
        return Money(amount + other.amount, currency);
    }
    return Money(amount + other.amount, currency);
}

Money Money::operator-(const Money& other) const {
    if (currency != other.currency) {
        // For simplicity, assume USD conversion rate of 1.0
        return Money(amount - other.amount, currency);
    }
    return Money(amount - other.amount, currency);
}

Money Money::operator*(double multiplier) const {
    return Money(amount * multiplier, currency);
}

}