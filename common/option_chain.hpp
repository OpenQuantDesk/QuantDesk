/*
 * Filename: option_chain.hpp
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

#include "common/quote.hpp"
#include <string>
#include <vector>
#include <chrono>

namespace common {

class OptionChain {
private:
    std::string underlying_;
    std::string expiration_;
    std::vector<Quote> calls_;
    std::vector<Quote> puts_;
    std::chrono::system_clock::time_point timestamp_;

public:
    OptionChain() : timestamp_(std::chrono::system_clock::now()) {}
    
    OptionChain(const std::string& underlying, const std::string& expiration)
        : underlying_(underlying), expiration_(expiration), 
          timestamp_(std::chrono::system_clock::now()) {}

    const std::string& getUnderlying() const { return underlying_; }
    const std::string& getExpiration() const { return expiration_; }
    const std::vector<Quote>& getCalls() const { return calls_; }
    const std::vector<Quote>& getPuts() const { return puts_; }
    const std::chrono::system_clock::time_point& getTimestamp() const { return timestamp_; }

    void setUnderlying(const std::string& underlying) { underlying_ = underlying; }
    void setExpiration(const std::string& expiration) { expiration_ = expiration; }
    void setCalls(const std::vector<Quote>& calls) { calls_ = calls; }
    void setPuts(const std::vector<Quote>& puts) { puts_ = puts; }
    void setTimestamp(const std::chrono::system_clock::time_point& timestamp) { timestamp_ = timestamp; }

    void addCall(const Quote& call) { calls_.push_back(call); }
    void addPut(const Quote& put) { puts_.push_back(put); }
    
    std::vector<Quote> getAllOptions() const;
    std::vector<double> getStrikes() const;
    Quote getCallByStrike(double strike) const;
    Quote getPutByStrike(double strike) const;
    
    size_t getCallCount() const { return calls_.size(); }
    size_t getPutCount() const { return puts_.size(); }
    size_t getTotalCount() const { return calls_.size() + puts_.size(); }
    
    bool isEmpty() const { return calls_.empty() && puts_.empty(); }
    bool isStale(std::chrono::minutes threshold = std::chrono::minutes(5)) const;
    
    void sortByStrike();
    void clear();
};

}