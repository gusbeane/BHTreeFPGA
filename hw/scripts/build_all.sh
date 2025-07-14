#!/bin/bash

# Build script for all kernels
# Usage: ./build_all.sh [csim|syn|cosim|package] [kernel_name]

set -e

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HW_DIR="$(dirname "$SCRIPT_DIR")"

# Default operation
OPERATION=${1:-csim}
KERNEL=${2:-all}

# Available kernels
KERNELS=("bhtree")

echo "=== Multi-Kernel Build Script ==="
echo "Operation: $OPERATION"
echo "Kernel: $KERNEL"
echo "HW Directory: $HW_DIR"

# Function to build a single kernel
build_kernel() {
    local kernel_name=$1
    local config_dir="$HW_DIR/config/$kernel_name"
    
    if [[ ! -d "$config_dir" ]]; then
        echo "Error: Config directory not found: $config_dir"
        return 1
    fi
    
    echo "Building kernel: $kernel_name"
    cd "$config_dir"
    
    case $OPERATION in
        csim)
            echo "  Running C simulation..."
            vitis-run --mode hls --csim --config hls_config.cfg
            ;;
        syn)
            echo "  Running synthesis..."
            vitis-run --mode hls --syn --config hls_config.cfg
            ;;
        cosim)
            echo "  Running C/RTL co-simulation..."
            vitis-run --mode hls --cosim --config hls_config.cfg
            ;;
        package)
            echo "  Running package..."
            vitis-run --mode hls --package --config hls_config.cfg
            ;;
        *)
            echo "Error: Unknown operation '$OPERATION'"
            echo "Usage: $0 [csim|syn|cosim|package] [kernel_name]"
            exit 1
            ;;
    esac
    
    echo "  ✅ $kernel_name completed successfully!"
}

# Build specified kernel(s)
if [[ "$KERNEL" == "all" ]]; then
    for kernel in "${KERNELS[@]}"; do
        build_kernel "$kernel"
    done
else
    # Check if specified kernel exists
    if [[ " ${KERNELS[*]} " =~ " ${KERNEL} " ]]; then
        build_kernel "$KERNEL"
    else
        echo "Error: Unknown kernel '$KERNEL'"
        echo "Available kernels: ${KERNELS[*]}"
        exit 1
    fi
fi

echo "=== Build completed successfully! ===" 