cmake_minimum_required(VERSION 3.22)

project(nBody VERSION 1.0 LANGUAGES CXX)

# Ensure cmake can find conan libraries
set(Eigen3_DIR ${CMAKE_BINARY_DIR})
set(Catch2_DIR ${CMAKE_BINARY_DIR})

# Make executables appear in build, not build/src
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Build application
add_subdirectory(app)

# Build library
add_subdirectory(src)

# Build tests
add_subdirectory(test)