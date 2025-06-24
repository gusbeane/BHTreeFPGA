cmake_minimum_required(VERSION 3.10)
project(BHTreeFPGA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Export compile commands for IDE support (generates compile_commands.json)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set output directories - executables go to the top level (sw directory)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR})

# Include directories
include_directories(include)

# Add main executable
add_executable(bhtree
    src/peano_hilbert.cpp
    src/pointcloud.cpp
    src/bhtree.cpp
)

# Add test executable
add_executable(test_peano_hilbert
    src/peano_hilbert.cpp
    tests/test_peano_hilbert.cpp
)

# Optional: Add compiler warnings
target_compile_options(bhtree PRIVATE -Wall -Wextra)

# Optional: Add debug flags for development
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(bhtree PRIVATE -g -O0)
else()
    target_compile_options(bhtree PRIVATE -O3)
endif() 