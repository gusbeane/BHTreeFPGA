cmake_minimum_required(VERSION 3.10)
project(BHTreeFPGA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)

# Add main executable
add_executable(bhtree
    src/peano_hilbert.cpp
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