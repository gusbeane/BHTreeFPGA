#!/bin/bash

# Compile script for FPGA host application
# This script compiles the host code that runs on the ARM processor

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Compiling FPGA Host Application ===${NC}"

# Check if we're in the correct directory
if [ ! -f "host/main.cpp" ]; then
    echo -e "${RED}Error: Please run this script from the hw/ directory${NC}"
    exit 1
fi

# Check if Vitis is sourced
if [ -z "$XILINX_VITIS" ]; then
    echo -e "${YELLOW}Warning: Vitis environment not detected. Sourcing setup script...${NC}"
    if [ -f "./setup_vitis.sh" ]; then
        source ./setup_vitis.sh
    else
        echo -e "${RED}Error: setup_vitis.sh not found. Please source Vitis manually.${NC}"
        exit 1
    fi
fi

# Create build directory for host code
mkdir -p build/host

# Compiler settings
CXX=g++
CXXFLAGS="-std=c++17 -O2 -Wall -Wextra"
INCLUDES="-I${XILINX_XRT}/include -I${XILINX_VITIS}/include/ -I./include"
LDFLAGS="-L${XILINX_XRT}/lib -L${XILINX_VITIS}/lib -lOpenCL -lxrt_coreutil -pthread"

# Source files
SRC_FILES="host/main.cpp"

# Output executable
OUTPUT="build/host/bhtree_host"

echo -e "${YELLOW}Compiling with:${NC}"
echo "  Compiler: $CXX"
echo "  Flags: $CXXFLAGS"
echo "  Includes: $INCLUDES"
echo "  Libraries: $LDFLAGS"
echo "  Source: $SRC_FILES"
echo "  Output: $OUTPUT"

# Compile
echo -e "${YELLOW}Compiling...${NC}"
$CXX $CXXFLAGS $INCLUDES $SRC_FILES $LDFLAGS -o $OUTPUT

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✅ Compilation successful!${NC}"
    echo -e "${GREEN}Executable created: $OUTPUT${NC}"
    echo ""
    echo -e "${YELLOW}To run the host application:${NC}"
    echo "  cd build/host"
    echo "  ./bhtree_host [path_to_xclbin]"
    echo ""
    echo -e "${YELLOW}Or use the run script:${NC}"
    echo "  ./scripts/run_fpga.sh"
else
    echo -e "${RED}❌ Compilation failed!${NC}"
    exit 1
fi 