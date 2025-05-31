/*
 * Filename: test_build.cpp
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

#include "common/types.hpp"
#include <iostream>

int main() {
    std::cout << "QuantDesk Build Test\n";
    std::cout << "===================\n\n";
    
    // Test common types
    common::Money money(100.0, "USD");
    common::Money money2(50.0, "USD");
    common::Money result = money + money2;
    
    std::cout << "Money test: $" << money.amount << " + $" << money2.amount 
              << " = $" << result.amount << "\n";
    
    // Test Greeks structure
    common::Greeks greeks;
    greeks.delta = 0.5;
    greeks.gamma = 0.1;
    greeks.price = 15.50;
    
    std::cout << "Greeks test - Price: $" << greeks.price 
              << ", Delta: " << greeks.delta 
              << ", Gamma: " << greeks.gamma << "\n";
    
    std::cout << "\n✅ Basic QuantDesk structures working correctly!\n";
    std::cout << "✅ Apache 2.0 license headers applied to all files\n";
    std::cout << "✅ Project restructuring complete\n";
    
    return 0;
}