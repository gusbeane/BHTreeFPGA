#!/usr/bin/env bash

# Source the Vitis setup
source hw/setup_vitis.sh

# Launch Cursor with the Vitis environment
echo "Launching Cursor with Vitis environment..."
echo "XILINX_VITIS: $XILINX_VITIS"
echo "XILINX_HLS: $XILINX_HLS"

# Launch Cursor
cursor . 