/*
 * Filename: main_window.hpp
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

#include "core/application.hpp"
#include "widgets/portfolio.hpp"
#include "widgets/market_data.hpp"
#include "widgets/strategy.hpp"
#include "widgets/risk.hpp"
#include "widgets/order_entry.hpp"
#include <QMainWindow>
#include <QTabWidget>
#include <QMenuBar>
#include <QStatusBar>
#include <QTimer>
#include <QLabel>
#include <memory>

namespace gui {

class MainWindow : public QMainWindow {
    Q_OBJECT

private:
    std::unique_ptr<core::Application> application_;
    
    QTabWidget* centralTabs_;
    widgets::PortfolioWidget* portfolioWidget_;
    widgets::MarketDataWidget* marketDataWidget_;
    widgets::StrategyWidget* strategyWidget_;
    widgets::RiskWidget* riskWidget_;
    widgets::OrderEntryWidget* orderEntryWidget_;
    
    QTimer* updateTimer_;
    QLabel* statusLabel_;
    QLabel* connectionLabel_;
    QLabel* timeLabel_;
    
public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow() override;

private slots:
    void onApplicationInitialized();
    void onQuotesUpdated(const std::map<std::string, common::Quote>& quotes);
    void onOpportunitiesUpdated(const std::vector<strategy::StrategyOpportunity>& opportunities);
    void onStatusUpdated(const QString& status);
    
    void updateDisplay();
    void updateTime();
    
    void openOrderDialog();
    void openSettingsDialog();
    void openAboutDialog();
    void openBacktestDialog();
    void openRiskDialog();
    
    void exportData();
    void importData();
    void saveLayout();
    void loadLayout();
    
    void connectToBroker();
    void disconnectFromBroker();
    void refreshData();
    
protected:
    void closeEvent(QCloseEvent* event) override;
    void showEvent(QShowEvent* event) override;

private:
    void setupUI();
    void setupMenuBar();
    void setupStatusBar();
    void setupConnections();
    void setupTimer();
    
    void initializeApplication();
    void connectSignals();
    
    void updateConnectionStatus();
    void updatePortfolioDisplay();
    void updateMarketDataDisplay();
    void updateStrategyDisplay();
    void updateRiskDisplay();
    
    void showErrorMessage(const QString& title, const QString& message);
    void showInfoMessage(const QString& title, const QString& message);
    bool showConfirmDialog(const QString& title, const QString& message);
};

}