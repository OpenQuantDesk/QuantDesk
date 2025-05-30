#include "core/application.hpp"
#include "gui/main_window.hpp"
#include <QApplication>
#include <iostream>
#include <memory>

int main(int argc, char* argv[]) {
    bool useGui = argc > 1 && std::string(argv[1]) == "--gui";
    
    if (useGui) {
        QApplication app(argc, argv);
        app.setApplicationName("QuantTrader Pro");
        app.setApplicationVersion("1.0");
        
        auto mainWindow = std::make_unique<gui::MainWindow>();
        mainWindow->show();
        
        return app.exec();
    } else {
        try {
            auto app = std::make_unique<core::Application>();
            app->initialize();
            app->run();
            return 0;
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    }
}