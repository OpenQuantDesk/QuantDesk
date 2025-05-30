cmake_minimum_required(VERSION 3.20)
project(Open Quant Desk - Options VERSION 4.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_HARDWARE_OPTIMIZATION "Enable hardware-specific optimizations" ON)
option(ENABLE_OPENCL "Enable OpenCL GPU acceleration" OFF)
option(ENABLE_GUI "Build with Qt GUI" ON)
option(ENABLE_TESTING "Build unit tests" ON)
option(ENABLE_BENCHMARKS "Build performance benchmarks" ON)
option(ENABLE_SANITIZERS "Enable address/thread sanitizers" OFF)
option(ENABLE_LTO "Enable link-time optimization" ON)

include(CheckCXXCompilerFlag)
include(FindPackageHandleStandardArgs)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ${ENABLE_LTO})
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")

if(ENABLE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,thread,undefined")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address,thread,undefined")
endif()

find_package(Threads REQUIRED)

if(ENABLE_GUI)
    find_package(Qt6 REQUIRED COMPONENTS Core Widgets Charts)
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTORCC ON)
    set(CMAKE_AUTOUIC ON)
endif()

find_package(PkgConfig REQUIRED)
pkg_check_modules(CURL REQUIRED libcurl)

find_package(nlohmann_json REQUIRED)

if(ENABLE_OPENCL)
    find_package(OpenCL REQUIRED)
    add_definitions(-DENABLE_OPENCL)
endif()

check_cxx_compiler_flag("-msse2" COMPILER_SUPPORTS_SSE2)
check_cxx_compiler_flag("-msse4.1" COMPILER_SUPPORTS_SSE4_1)
check_cxx_compiler_flag("-mavx" COMPILER_SUPPORTS_AVX)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
check_cxx_compiler_flag("-mfma" COMPILER_SUPPORTS_FMA)

if(ENABLE_HARDWARE_OPTIMIZATION)
    if(COMPILER_SUPPORTS_SSE2)
        add_definitions(-D__SSE2__)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    endif()
    
    if(COMPILER_SUPPORTS_SSE4_1)
        add_definitions(-D__SSE4_1__)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
    endif()
    
    if(COMPILER_SUPPORTS_AVX)
        add_definitions(-D__AVX__)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    endif()
    
    if(COMPILER_SUPPORTS_AVX2)
        add_definitions(-D__AVX2__)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
    endif()
    
    if(COMPILER_SUPPORTS_FMA)
        add_definitions(-D__FMA__)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfma")
    endif()
endif()

add_library(common STATIC
    common/types.cpp
)

target_include_directories(common PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(math_engine STATIC
    math/engine.cpp
    math/optimization/hardware.cpp
    math/greeks/firstOrder.cpp
    math/greeks/secondOrder.cpp
    math/forecasting/monte_carlo.cpp
    math/metrics/volatility.cpp
    math/metrics/probability.cpp
    math/simulation/params.cpp
    math/simulation/results.cpp
)

target_link_libraries(math_engine PUBLIC
    common
    Threads::Threads
)

target_include_directories(math_engine PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if(ENABLE_OPENCL)
    target_link_libraries(math_engine PRIVATE OpenCL::OpenCL)
endif()

add_library(broker_interface STATIC
    broker/interface.cpp
    broker/registry.cpp
)

target_link_libraries(broker_interface PUBLIC
    common
    Threads::Threads
    ${CMAKE_DL_LIBS}
)

target_include_directories(broker_interface PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(portfolio_manager STATIC
    portfolio/manager.cpp
    portfolio/position.cpp
    portfolio/risk.cpp
)

target_link_libraries(portfolio_manager PUBLIC
    common
    math_engine
    Threads::Threads
)

target_include_directories(portfolio_manager PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(strategy_engine STATIC
    strategy/engine.cpp
    strategy/optimization.cpp
    strategy/opportunity.cpp
)

target_link_libraries(strategy_engine PUBLIC
    common
    math_engine
    portfolio_manager
    Threads::Threads
)

target_include_directories(strategy_engine PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(data_providers STATIC
    data/economic.cpp
    data/providers/fred/getFredData.cpp
    data/providers/yahoo/getYahooData.cpp
    data/providers/avantage/getAVantageData.cpp
)

target_link_libraries(data_providers PUBLIC
    common
    ${CURL_LIBRARIES}
    nlohmann_json::nlohmann_json
    Threads::Threads
)

target_include_directories(data_providers PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CURL_INCLUDE_DIRS}
)

target_compile_options(data_providers PRIVATE ${CURL_CFLAGS_OTHER})

add_library(core_application STATIC
    core/application.cpp
)

target_link_libraries(core_application PUBLIC
    common
    math_engine
    broker_interface
    portfolio_manager
    strategy_engine
    data_providers
    nlohmann_json::nlohmann_json
    Threads::Threads
)

target_include_directories(core_application PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if(ENABLE_GUI)
    add_library(gui_components STATIC
        gui/main_window.cpp
        gui/widgets/portfolio_widget.cpp
        gui/widgets/market_data_widget.cpp
        gui/widgets/strategy_widget.cpp
        gui/widgets/risk_widget.cpp
        gui/widgets/order_entry_widget.cpp
    )
    
    target_link_libraries(gui_components PUBLIC
        core_application
        Qt6::Core
        Qt6::Widgets
        Qt6::Charts
    )
    
    target_include_directories(gui_components PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

add_executable(options_desk
    main.cpp
)

target_link_libraries(options_desk PRIVATE
    core_application
)

if(ENABLE_GUI)
    target_link_libraries(options_desk PRIVATE gui_components)
endif()

set_target_properties(options_desk PROPERTIES
    OUTPUT_NAME "options_desk"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

add_subdirectory(broker/brokers)

if(ENABLE_TESTING)
    find_package(GTest REQUIRED)
    enable_testing()
    
    add_executable(unit_tests
        tests/math/test_engine.cpp
        tests/math/test_greeks.cpp
        tests/portfolio/test_manager.cpp
        tests/strategy/test_engine.cpp
        tests/core/test_application.cpp
        tests/tradeDesk.cpp
    )
    
    target_link_libraries(unit_tests PRIVATE
        core_application
        GTest::GTest
        GTest::Main
        Threads::Threads
    )
    
    target_include_directories(unit_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    
    add_test(NAME UnitTests COMMAND unit_tests)
    
    set_target_properties(unit_tests PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
    )
endif()

if(ENABLE_BENCHMARKS)
    find_package(benchmark REQUIRED)
    
    add_executable(performance_benchmarks
        benchmarks/math_benchmark.cpp
        benchmarks/portfolio_benchmark.cpp
        benchmarks/strategy_benchmark.cpp
    )
    
    target_link_libraries(performance_benchmarks PRIVATE
        core_application
        benchmark::benchmark
        Threads::Threads
    )
    
    target_include_directories(performance_benchmarks PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    
    set_target_properties(performance_benchmarks PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/benchmarks
    )
endif()

install(TARGETS options_desk
    RUNTIME DESTINATION bin
)

install(FILES
    config/default_config.json
    config/watchlist.txt
    DESTINATION etc/options_desk
)

install(DIRECTORY docs/
    DESTINATION share/doc/options_desk
)

set(CPACK_PACKAGE_NAME "Open Quant Desk - Options")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Quantitative Options Trading Platform")
set(CPACK_PACKAGE_VENDOR "OpenQuantDesk, Inc")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_CONTACT "developers@openquantdesk.com")

if(WIN32)
    set(CPACK_GENERATOR "NSIS")
elseif(APPLE)
    set(CPACK_GENERATOR "DragNDrop")
else()
    set(CPACK_GENERATOR "DEB;RPM;TGZ")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libcurl4, libqt6core6, libqt6widgets6")
    set(CPACK_RPM_PACKAGE_REQUIRES "libcurl, qt6-qtbase")
endif()

include(CPack)

function(add_broker_plugin PLUGIN_NAME)
    add_library(${PLUGIN_NAME} SHARED
        broker/brokers/${PLUGIN_NAME}/${PLUGIN_NAME}.cpp
    )
    
    target_link_libraries(${PLUGIN_NAME} PRIVATE
        broker_interface
        ${ARGN}
    )
    
    target_include_directories(${PLUGIN_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    
    set_target_properties(${PLUGIN_NAME} PROPERTIES
        PREFIX ""
        SUFFIX ".so"
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plugins
    )
    
    install(TARGETS ${PLUGIN_NAME}
        LIBRARY DESTINATION lib/options_desk/plugins
    )
endfunction()

message(STATUS "")
message(STATUS "Open Quant Desk - Options Build Configuration:")
message(STATUS "  Version: ${PROJECT_VERSION}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Hardware Optimization: ${ENABLE_HARDWARE_OPTIMIZATION}")
message(STATUS "  OpenCL Support: ${ENABLE_OPENCL}")
message(STATUS "  GUI Support: ${ENABLE_GUI}")
message(STATUS "  Testing: ${ENABLE_TESTING}")
message(STATUS "  Benchmarks: ${ENABLE_BENCHMARKS}")
message(STATUS "  Link-Time Optimization: ${ENABLE_LTO}")
message(STATUS "  Sanitizers: ${ENABLE_SANITIZERS}")
message(STATUS "")

if(ENABLE_HARDWARE_OPTIMIZATION)
    message(STATUS "Hardware Optimizations:")
    message(STATUS "  SSE2: ${COMPILER_SUPPORTS_SSE2}")
    message(STATUS "  SSE4.1: ${COMPILER_SUPPORTS_SSE4_1}")
    message(STATUS "  AVX: ${COMPILER_SUPPORTS_AVX}")
    message(STATUS "  AVX2: ${COMPILER_SUPPORTS_AVX2}")
    message(STATUS "  FMA: ${COMPILER_SUPPORTS_FMA}")
    message(STATUS "")
endif()

message(STATUS "Dependencies:")
message(STATUS "  Threads: ${CMAKE_THREAD_LIBS_INIT}")
message(STATUS "  CURL: ${CURL_VERSION}")
message(STATUS "  nlohmann_json: Found")
if(ENABLE_GUI)
    message(STATUS "  Qt6: ${Qt6_VERSION}")
endif()
if(ENABLE_OPENCL)
    message(STATUS "  OpenCL: ${OpenCL_VERSION}")
endif()
message(STATUS "")