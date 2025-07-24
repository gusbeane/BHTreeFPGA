#!/bin/bash

# Run script for FPGA host application
# This script loads the xclbin onto the FPGA and runs the host application

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Running FPGA Barnes-Hut Tree Application ===${NC}"

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

# Paths
XCLBIN_PATH="build/bhtree/create_bhtree_kernel.xclbin"
HOST_EXECUTABLE="build/host/bhtree_host"

# Check if xclbin exists
if [ ! -f "$XCLBIN_PATH" ]; then
    echo -e "${RED}Error: xclbin file not found at $XCLBIN_PATH${NC}"
    echo -e "${YELLOW}Please make sure the kernel has been built successfully.${NC}"
    exit 1
fi

# Check if host executable exists
if [ ! -f "$HOST_EXECUTABLE" ]; then
    echo -e "${YELLOW}Host executable not found. Compiling...${NC}"
    if [ -f "scripts/compile_host.sh" ]; then
        ./scripts/compile_host.sh
    else
        echo -e "${RED}Error: compile_host.sh not found${NC}"
        exit 1
    fi
fi

# Display FPGA information
echo -e "${BLUE}=== FPGA Information ===${NC}"
if command -v xbutil &> /dev/null; then
    echo -e "${YELLOW}Available FPGA devices:${NC}"
    xbutil examine
else
    echo -e "${YELLOW}xbutil not available. Using basic device detection.${NC}"
fi

# Check if we have XRT environment
if [ -z "$XILINX_XRT" ]; then
    echo -e "${RED}Error: XRT environment not found. Please source XRT setup script.${NC}"
    exit 1
fi

echo -e "${BLUE}=== Loading and Running Application ===${NC}"
echo -e "${YELLOW}XCLBIN file: $XCLBIN_PATH${NC}"
echo -e "${YELLOW}Host executable: $HOST_EXECUTABLE${NC}"

# Check file sizes
echo -e "${YELLOW}File information:${NC}"
echo "  XCLBIN size: $(du -h $XCLBIN_PATH | cut -f1)"
echo "  Host executable size: $(du -h $HOST_EXECUTABLE | cut -f1)"

# Run the application
echo -e "${BLUE}=== Executing on FPGA ===${NC}"
echo -e "${YELLOW}Starting execution...${NC}"

# Change to the directory containing the executable
cd build/host

# Run with timing
start_time=$(date +%s.%N)
./bhtree_host ../../$XCLBIN_PATH
exit_code=$?
end_time=$(date +%s.%N)

# Calculate execution time
execution_time=$(echo "$end_time - $start_time" | bc)

# Go back to original directory
cd ../../

echo -e "${BLUE}=== Execution Complete ===${NC}"
echo -e "${YELLOW}Total execution time: ${execution_time} seconds${NC}"

if [ $exit_code -eq 0 ]; then
    echo -e "${GREEN}✅ FPGA execution completed successfully!${NC}"
else
    echo -e "${RED}❌ FPGA execution failed with exit code $exit_code${NC}"
fi

exit $exit_code 