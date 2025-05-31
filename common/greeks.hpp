/*
 * Filename: greeks.hpp
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

namespace common {

struct Greeks {
    double delta = 0.0;
    double gamma = 0.0;
    double theta = 0.0;
    double vega = 0.0;
    double rho = 0.0;
    double impliedVol = 0.0;
    double price = 0.0;
    
    Greeks() = default;
    Greeks(double d, double g, double t, double v, double r, double iv, double p)
        : delta(d), gamma(g), theta(t), vega(v), rho(r), impliedVol(iv), price(p) {}
    
    Greeks operator+(const Greeks& other) const;
    Greeks operator-(const Greeks& other) const;
    Greeks operator*(double multiplier) const;
    Greeks& operator+=(const Greeks& other);
    Greeks& operator-=(const Greeks& other);
    Greeks& operator*=(double multiplier);
};

struct ExtendedGreeks : public Greeks {
    double volga = 0.0;
    double vanna = 0.0;
    double charm = 0.0;
    double speed = 0.0;
    double zomma = 0.0;
    double color = 0.0;
    double dollarDelta = 0.0;
    double dollarGamma = 0.0;
    double pinRisk = 0.0;
    
    ExtendedGreeks() = default;
    ExtendedGreeks(const Greeks& base) : Greeks(base) {}
};

}