#pragma once

#include "core/application.hpp"
#include "gui/widgets/portfolio_widget.hpp"
#include "gui/widgets/market_data_widget.hpp"
#include "gui/widgets/strategy_widget.hpp"
#include "gui/widgets/risk_widget.hpp"
#include "gui/widgets/order_entry_widget.hpp"
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