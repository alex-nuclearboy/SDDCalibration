cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# Set the project name
project(SDDCalibration VERSION 1.0)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Gpad Tree Spectrum Rint Physics MathCore Thread)

# Add the include directory to the list of include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Include ROOT headers
include(${ROOT_USE_FILE})

# Add executable file to the project
add_executable(peak_finder 
    src/SDDCalibration.cxx
    peak_finder.cpp)

# Link ROOT libraries
target_link_libraries(peak_finder PRIVATE ${ROOT_LIBRARIES})

# Set output directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

# Additional compile options
add_compile_options(peak_finder PRIVATE -Wall -Wextra -pedantic)
