#!/bin/bash

# Build script for bhtree kernel
# Usage: ./build_bhtree.sh [csim|syn|cosim|package]

set -e

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HW_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="$HW_DIR/config/bhtree"

# Default operation
OPERATION=${1:-csim}

echo "Building bhtree kernel with operation: $OPERATION"
echo "HW Directory: $HW_DIR"
echo "Config Directory: $CONFIG_DIR"

# Change to config directory
cd "$CONFIG_DIR"

# Run the appropriate vitis-run command
case $OPERATION in
    csim)
        echo "Running C simulation..."
        vitis-run --mode hls --csim --config hls_config.cfg --work_dir ../../build/bhtree
        ;;
    syn)
        echo "Running synthesis..."
        v++ -c --mode hls --config hls_config.cfg --work_dir ../../build/bhtree
        ;;
    cosim)
        echo "Running C/RTL co-simulation..."
        vitis-run --mode hls --cosim --config hls_config.cfg --work_dir ../../build/bhtree
        ;;
    package)
        echo "Running package..."
        vitis-run --mode hls --package --config hls_config.cfg --work_dir ../../build/bhtree
        ;;
    impl)
        echo "Running implementation..."
        vitis-run --mode hls --impl --config hls_config.cfg --work_dir ../../build/bhtree
        ;;
    link)
        echo "Running link..."
        v++ -l ../../build/bhtree/create_bhtree_kernel.xo --platform xilinx_u200_gen3x16_xdma_2_202110_1 -o ../../build/bhtree/create_bhtree_kernel.xclbin
        ;;
    *)
        echo "Error: Unknown operation '$OPERATION'"
        echo "Usage: $0 [csim|syn|cosim|package|impl|link]"
        exit 1
        ;;
esac

echo "Build completed successfully!" 
