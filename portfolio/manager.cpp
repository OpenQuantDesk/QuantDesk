/*
 * Filename: manager.cpp
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

#include "portfolio/manager.hpp"
#include <future>
#include <algorithm>

namespace portfolio {

PortfolioManager::PortfolioManager(std::shared_ptr<math::MathEngine> engine) 
    : mathEngine_(engine), running_(false) {}

PortfolioManager::~PortfolioManager() {
    stop();
}

void PortfolioManager::start() {
    std::unique_lock<std::shared_mutex> lock(positionsLock_);
    if (running_.load()) return;
    
    running_ = true;
    analyticsThread_ = std::thread(&PortfolioManager::analyticsLoop, this);
}

void PortfolioManager::stop() {
    running_ = false;
    if (analyticsThread_.joinable()) {
        analyticsThread_.join();
    }
}

void PortfolioManager::updatePosition(const common::Position& position) {
    std::unique_lock<std::shared_mutex> lock(positionsLock_);
    
    EnhancedPosition& pos = positions_[position.symbol];
    pos.symbol = position.symbol;
    pos.quantity = position.quantity;
    pos.avgPrice = position.avgCost.amount / (position.quantity * 100);
    pos.lastUpdate = std::chrono::system_clock::now();
    
    if (position.quantity == 0.0) {
        positions_.erase(position.symbol);
    }
    
    updateTotalValues();
    notifyPositionUpdated(pos);
}

void PortfolioManager::updatePosition(const std::string& symbol, double quantity, double price) {
    std::unique_lock<std::shared_mutex> lock(positionsLock_);
    
    EnhancedPosition& pos = positions_[symbol];
    pos.symbol = symbol;
    pos.quantity = quantity;
    pos.avgPrice = price;
    pos.lastUpdate = std::chrono::system_clock::now();
    
    if (quantity == 0.0) {
        positions_.erase(symbol);
    }
    
    updateTotalValues();
    notifyPositionUpdated(pos);
}

void PortfolioManager::updatePositionPrice(const std::string& symbol, double price, double underlyingPrice) {
    std::unique_lock<std::shared_mutex> lock(positionsLock_);
    
    auto it = positions_.find(symbol);
    if (it != positions_.end()) {
        double oldPrice = it->second.currentPrice;
        it->second.currentPrice = price;
        it->second.underlyingPrice = underlyingPrice;
        it->second.dayChange = price - oldPrice;
        it->second.lastUpdate = std::chrono::system_clock::now();
        
        updateTotalValues();
        notifyPositionUpdated(it->second);
    }
}

void PortfolioManager::updatePositionGreeks(const std::string& symbol, const common::Greeks& greeks) {
    std::unique_lock<std::shared_mutex> lock(positionsLock_);
    
    auto it = positions_.find(symbol);
    if (it != positions_.end()) {
        it->second.greeks.delta = greeks.delta;
        it->second.greeks.gamma = greeks.gamma;
        it->second.greeks.theta = greeks.theta;
        it->second.greeks.vega = greeks.vega;
        it->second.greeks.rho = greeks.rho;
        it->second.greeks.impliedVol = greeks.impliedVol;
        it->second.greeks.price = greeks.price;
        
        notifyPositionUpdated(it->second);
    }
}

void PortfolioManager::updateAllPositions(const std::map<std::string, common::Quote>& quotes,
                                         const std::map<std::string, common::OptionChain>& chains) {
    std::unique_lock<std::shared_mutex> lock(positionsLock_);
    
    for (auto& [symbol, position] : positions_) {
        auto quoteIt = quotes.find(position.underlying);
        if (quoteIt != quotes.end() && quoteIt->second.last.has_value()) {
            position.underlyingPrice = *quoteIt->second.last;
            
            auto chainIt = chains.find(position.underlying);
            if (chainIt != chains.end()) {
                for (const auto& option : chainIt->second.calls) {
                    if (option.symbol == symbol) {
                        if (option.last.has_value()) {
                            position.currentPrice = *option.last;
                        } else if (option.bid.has_value() && option.ask.has_value()) {
                            position.currentPrice = (*option.bid + *option.ask) / 2.0;
                        }
                        
                        if (option.greeks.has_value()) {
                            position.greeks = *option.greeks;
                        }
                        break;
                    }
                }
                
                for (const auto& option : chainIt->second.puts) {
                    if (option.symbol == symbol) {
                        if (option.last.has_value()) {
                            position.currentPrice = *option.last;
                        } else if (option.bid.has_value() && option.ask.has_value()) {
                            position.currentPrice = (*option.bid + *option.ask) / 2.0;
                        }
                        
                        if (option.greeks.has_value()) {
                            position.greeks = *option.greeks;
                        }
                        break;
                    }
                }
            }
            
            position.updateAnalytics(position.underlyingPrice, position.greeks.impliedVol, 0.05, *mathEngine_);
        }
    }
    
    updateTotalValues();
}

std::future<PortfolioMetrics> PortfolioManager::calculateMetrics() const {
    return std::async(std::launch::async, [this]() {
        std::shared_lock<std::shared_mutex> lock(positionsLock_);
        PortfolioMetrics metrics;
        
        metrics.cash = cash_;
        metrics.timestamp = std::chrono::system_clock::now();
        
        for (const auto& [symbol, pos] : positions_) {
            metrics.totalValue += pos.getNotionalValue();
            metrics.totalPnL += pos.unrealizedPnL;
            metrics.dayPnL += pos.dayChange * pos.quantity * 100;
            
            metrics.netDelta += pos.extendedGreeks.delta;
            metrics.netGamma += pos.extendedGreeks.gamma;
            metrics.netTheta += pos.extendedGreeks.theta;
            metrics.netVega += pos.extendedGreeks.vega;
            metrics.netRho += pos.extendedGreeks.rho;
            metrics.netVolga += pos.extendedGreeks.volga;
            metrics.netVanna += pos.extendedGreeks.vanna;
            metrics.dollarDelta += pos.extendedGreeks.dollarDelta;
            metrics.dollarGamma += pos.extendedGreeks.dollarGamma;
            
            metrics.underlyingExposures[pos.underlying] += pos.getNotionalValue();
        }
        
        metrics.totalValue += metrics.cash;
        
        calculatePortfolioRisk(metrics);
        calculateAdvancedMetrics(metrics);
        
        return metrics;
    });
}

std::future<std::vector<EnhancedPosition>> PortfolioManager::getPositions() const {
    return std::async(std::launch::async, [this]() {
        std::shared_lock<std::shared_mutex> lock(positionsLock_);
        std::vector<EnhancedPosition> result;
        result.reserve(positions_.size());
        
        for (const auto& [symbol, pos] : positions_) {
            result.push_back(pos);
        }
        
        return result;
    });
}

std::future<EnhancedPosition> PortfolioManager::getPosition(const std::string& symbol) const {
    return std::async(std::launch::async, [this, symbol]() {
        std::shared_lock<std::shared_mutex> lock(positionsLock_);
        auto it = positions_.find(symbol);
        if (it != positions_.end()) {
            return it->second;
        }
        return EnhancedPosition{};
    });
}

std::future<std::vector<EnhancedPosition>> PortfolioManager::getPositionsByUnderlying(const std::string& underlying) const {
    return std::async(std::launch::async, [this, underlying]() {
        std::shared_lock<std::shared_mutex> lock(positionsLock_);
        std::vector<EnhancedPosition> result;
        
        for (const auto& [symbol, pos] : positions_) {
            if (pos.underlying == underlying) {
                result.push_back(pos);
            }
        }
        
        return result;
    });
}

std::future<std::vector<EnhancedPosition>> PortfolioManager::getExpiringPositions(int days) const {
    return std::async(std::launch::async, [this, days]() {
        std::shared_lock<std::shared_mutex> lock(positionsLock_);
        std::vector<EnhancedPosition> result;
        
        for (const auto& [symbol, pos] : positions_) {
            if (pos.isExpiring(days)) {
                result.push_back(pos);
            }
        }
        
        return result;
    });
}

std::future<std::vector<EnhancedPosition>> PortfolioManager::getHighRiskPositions(double threshold) const {
    return std::async(std::launch::async, [this, threshold]() {
        std::shared_lock<std::shared_mutex> lock(positionsLock_);
        std::vector<EnhancedPosition> result;
        
        for (const auto& [symbol, pos] : positions_) {
            double riskScore = calculatePositionRisk(pos);
            if (riskScore > threshold) {
                result.push_back(pos);
            }
        }
        
        return result;
    });
}

double PortfolioManager::getCash() const {
    std::shared_lock<std::shared_mutex> lock(positionsLock_);
    return cash_;
}

void PortfolioManager::setCash(double cash) {
    std::unique_lock<std::shared_mutex> lock(positionsLock_);
    cash_ = cash;
    updateTotalValues();
}

void PortfolioManager::addCash(double amount) {
    std::unique_lock<std::shared_mutex> lock(positionsLock_);
    cash_ += amount;
    updateTotalValues();
}

common::Greeks PortfolioManager::getPortfolioGreeks() const {
    std::shared_lock<std::shared_mutex> lock(positionsLock_);
    common::Greeks total;
    
    for (const auto& [symbol, pos] : positions_) {
        total.delta += pos.greeks.delta * pos.quantity;
        total.gamma += pos.greeks.gamma * pos.quantity;
        total.theta += pos.greeks.theta * pos.quantity;
        total.vega += pos.greeks.vega * pos.quantity;
        total.rho += pos.greeks.rho * pos.quantity;
    }
    
    return total;
}

std::vector<std::string> PortfolioManager::getUnderlyings() const {
    std::shared_lock<std::shared_mutex> lock(positionsLock_);
    std::vector<std::string> underlyings;
    
    for (const auto& [symbol, pos] : positions_) {
        if (std::find(underlyings.begin(), underlyings.end(), pos.underlying) == underlyings.end()) {
            underlyings.push_back(pos.underlying);
        }
    }
    
    return underlyings;
}

std::map<std::string, double> PortfolioManager::getUnderlyingExposures() const {
    std::shared_lock<std::shared_mutex> lock(positionsLock_);
    std::map<std::string, double> exposures;
    
    for (const auto& [symbol, pos] : positions_) {
        exposures[pos.underlying] += pos.getNotionalValue();
    }
    
    return exposures;
}

void PortfolioManager::addPositionCallback(std::function<void(const EnhancedPosition&)> callback) {
    std::unique_lock<std::shared_mutex> lock(positionsLock_);
    positionCallbacks_.push_back(callback);
}

void PortfolioManager::addMetricsCallback(std::function<void(const PortfolioMetrics&)> callback) {
    std::unique_lock<std::shared_mutex> lock(metricsLock_);
    metricsCallbacks_.push_back(callback);
}

void PortfolioManager::analyticsLoop() {
    while (running_.load()) {
        try {
            auto metrics = calculateMetrics().get();
            updateMetricsHistory(metrics);
            updatePnLHistory();
            
            std::shared_lock<std::shared_mutex> lock(metricsLock_);
            for (const auto& callback : metricsCallbacks_) {
                callback(metrics);
            }
        } catch (const std::exception& e) {
            
        }
        
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }
}

void PortfolioManager::calculatePortfolioRisk(PortfolioMetrics& metrics) const {
    const int numSims = 10000;
    std::vector<double> pnlSims;
    pnlSims.reserve(numSims);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> marketDist(0.0, 0.02);
    
    for (int i = 0; i < numSims; ++i) {
        double portfolioPnL = 0.0;
        double marketMove = marketDist(gen);
        
        for (const auto& [symbol, pos] : positions_) {
            double pnl = pos.extendedGreeks.delta * pos.underlyingPrice * marketMove +
                        0.5 * pos.extendedGreeks.gamma * pos.underlyingPrice * pos.underlyingPrice * 
                        marketMove * marketMove;
            portfolioPnL += pnl;
        }
        pnlSims.push_back(portfolioPnL);
    }
    
    std::sort(pnlSims.begin(), pnlSims.end());
    metrics.var95 = pnlSims[static_cast<size_t>(numSims * 0.05)];
    metrics.var99 = pnlSims[static_cast<size_t>(numSims * 0.01)];
    
    double es = 0.0;
    for (int i = 0; i < numSims * 0.05; ++i) {
        es += pnlSims[i];
    }
    metrics.expectedShortfall = es / (numSims * 0.05);
}

void PortfolioManager::calculateAdvancedMetrics(PortfolioMetrics& metrics) const {
    metrics.leverageRatio = calculateLeverageRatio();
    metrics.correlationRisk = calculateConcentrationRisk();
    
    if (!metricsHistory_.empty()) {
        std::vector<double> returns;
        for (size_t i = 1; i < metricsHistory_.size(); ++i) {
            double ret = (metricsHistory_[i].totalValue - metricsHistory_[i-1].totalValue) / 
                        metricsHistory_[i-1].totalValue;
            returns.push_back(ret);
        }
        
        if (!returns.empty()) {
            double avgReturn = std::accumulate(returns.begin(), returns.end(), 0.0) / returns.size();
            double variance = 0.0;
            for (double ret : returns) {
                variance += (ret - avgReturn) * (ret - avgReturn);
            }
            variance /= returns.size();
            double vol = std::sqrt(variance * 252);
            
            if (vol > 0) {
                metrics.sharpeRatio = (avgReturn * 252 - 0.02) / vol;
            }
        }
    }
}

void PortfolioManager::updateMetricsHistory(const PortfolioMetrics& metrics) {
    std::unique_lock<std::shared_mutex> lock(metricsLock_);
    metricsHistory_.push_back(metrics);
    
    if (metricsHistory_.size() > 252) {
        metricsHistory_.erase(metricsHistory_.begin());
    }
}

void PortfolioManager::updatePnLHistory() {
    std::unique_lock<std::shared_mutex> lock(positionsLock_);
    
    for (const auto& [symbol, pos] : positions_) {
        pnlHistory_[symbol].push_back(pos.unrealizedPnL);
        
        if (pnlHistory_[symbol].size() > 252) {
            pnlHistory_[symbol].erase(pnlHistory_[symbol].begin());
        }
    }
}

void PortfolioManager::notifyPositionUpdated(const EnhancedPosition& position) {
    for (const auto& callback : positionCallbacks_) {
        callback(position);
    }
}

void PortfolioManager::updateTotalValues() {
    double totalValue = cash_;
    double totalPnL = 0.0;
    double riskScore = 0.0;
    
    for (const auto& [symbol, pos] : positions_) {
        totalValue += pos.getNotionalValue();
        totalPnL += pos.unrealizedPnL;
        riskScore += calculatePositionRisk(pos);
    }
    
    totalValue_.store(totalValue);
    totalPnL_.store(totalPnL);
    riskScore_.store(riskScore / positions_.size());
}

double PortfolioManager::calculateConcentrationRisk() const {
    auto exposures = getUnderlyingExposures();
    double totalExposure = 0.0;
    double maxExposure = 0.0;
    
    for (const auto& [underlying, exposure] : exposures) {
        totalExposure += std::abs(exposure);
        maxExposure = std::max(maxExposure, std::abs(exposure));
    }
    
    return (totalExposure > 0) ? maxExposure / totalExposure : 0.0;
}

double PortfolioManager::calculateLeverageRatio() const {
    double notionalValue = 0.0;
    for (const auto& [symbol, pos] : positions_) {
        notionalValue += pos.getNotionalValue();
    }
    
    double equity = totalValue_.load();
    return (equity > 0) ? notionalValue / equity : 0.0;
}

double PortfolioManager::calculatePositionRisk(const EnhancedPosition& position) const {
    double timeRisk = position.isExpiring(7) ? 0.5 : 0.0;
    double deltaRisk = std::min(1.0, std::abs(position.greeks.delta) / 0.5);
    double thetaRisk = std::min(1.0, std::abs(position.greeks.theta) / 50.0);
    
    return (timeRisk + deltaRisk + thetaRisk) / 3.0;
}

}