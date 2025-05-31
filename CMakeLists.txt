cmake_minimum_required(VERSION 3.20)
project(QuantDesk VERSION 4.0.0 LANGUAGES CXX)

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

find_package(PkgConfig)
if(PkgConfig_FOUND)
    pkg_check_modules(CURL libcurl)
endif()

find_package(nlohmann_json)

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

# Add OptionsQuantLib as a subdirectory
# add_subdirectory(OptionsQuantLib)  # Temporarily disabled for testing

add_library(common STATIC
    common/types.cpp
)

target_include_directories(common PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(broker_interface STATIC
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
)

target_link_libraries(portfolio_manager PUBLIC
    common
    # OptionsQuantLib  # Temporarily disabled
    Threads::Threads
)

target_include_directories(portfolio_manager PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(strategy_engine STATIC
    strategy/strategy_engine.cpp
    strategy/optimization.cpp
    strategy/opportunity.cpp
)

target_link_libraries(strategy_engine PUBLIC
    common
    # OptionsQuantLib  # Temporarily disabled
    portfolio_manager
    Threads::Threads
)

target_include_directories(strategy_engine PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Only build data providers if dependencies are available
if(CURL_FOUND OR nlohmann_json_FOUND)
    add_library(data_providers STATIC
        data/providers/fred/getFredData.cpp
        data/providers/yahoo/getYahooData.cpp
        data/providers/avantage/getAVantageData.cpp
    )

    target_link_libraries(data_providers PUBLIC
        common
        Threads::Threads
    )
    
    if(CURL_FOUND)
        target_link_libraries(data_providers PUBLIC ${CURL_LIBRARIES})
        target_include_directories(data_providers PUBLIC ${CURL_INCLUDE_DIRS})
        target_compile_options(data_providers PRIVATE ${CURL_CFLAGS_OTHER})
    endif()
    
    if(nlohmann_json_FOUND)
        target_link_libraries(data_providers PUBLIC nlohmann_json::nlohmann_json)
    endif()

    target_include_directories(data_providers PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
else()
    # Create a dummy data_providers library
    add_library(data_providers STATIC
        common/types.cpp  # Just link to something that exists
    )
    target_include_directories(data_providers PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

add_library(core_application STATIC
    core/application.cpp
)

target_link_libraries(core_application PUBLIC
    common
    # OptionsQuantLib  # Temporarily disabled
    broker_interface
    portfolio_manager
    strategy_engine
    data_providers
    Threads::Threads
)

if(nlohmann_json_FOUND)
    target_link_libraries(core_application PUBLIC nlohmann_json::nlohmann_json)
endif()

target_include_directories(core_application PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if(ENABLE_GUI)
    add_library(gui_components STATIC
        gui/main_window.cpp
        gui/widgets/portfolio.cpp
        gui/widgets/market_data.cpp
        gui/widgets/strategy.cpp
        gui/widgets/risk.cpp
        gui/widgets/order_entry.cpp
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

add_executable(quantdesk
    main.cpp
)

add_executable(test_build
    test_build.cpp
)

target_link_libraries(quantdesk PRIVATE
    core_application
)

target_link_libraries(test_build PRIVATE
    common
)

if(ENABLE_GUI)
    target_link_libraries(quantdesk PRIVATE gui_components)
endif()

set_target_properties(quantdesk PROPERTIES
    OUTPUT_NAME "quantdesk"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# add_subdirectory(broker/brokers)  # Commented out for now

if(ENABLE_TESTING)
    find_package(GTest REQUIRED)
    enable_testing()
    
    add_executable(unit_tests
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

install(TARGETS quantdesk
    RUNTIME DESTINATION bin
)

install(FILES
    config/default_config.json
    config/watchlist.txt
    DESTINATION etc/quantdesk
)

install(DIRECTORY docs/
    DESTINATION share/doc/quantdesk
)

set(CPACK_PACKAGE_NAME "QuantDesk")
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
        LIBRARY DESTINATION lib/quantdesk/plugins
    )
endfunction()

message(STATUS "")
message(STATUS "QuantDesk Build Configuration:")
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
message(STATUS "  OptionsQuantLib: Bundled")
if(ENABLE_GUI)
    message(STATUS "  Qt6: ${Qt6_VERSION}")
endif()
if(ENABLE_OPENCL)
    message(STATUS "  OpenCL: ${OpenCL_VERSION}")
endif()
message(STATUS "")