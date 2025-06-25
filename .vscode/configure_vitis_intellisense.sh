#!/usr/bin/env bash

# Script to configure Cursor IntelliSense for Vitis HLS development
# This script sources the Vitis environment and updates the c_cpp_properties.json with actual paths

echo "Configuring Cursor IntelliSense for Vitis HLS..."

# Source the Vitis environment (temporarily for this script)
export PYTHONPATH="${PYTHONPATH:-}"
export MATLABPATH="${MATLABPATH:-}"
export XILINX_HOME="/opt/xilinx"
VITIS_VERSION="2025.1"
export VITIS_DIR="${XILINX_HOME}/${VITIS_VERSION}/Vitis"
export XILINX_VITIS="${VITIS_DIR}"
export XILINX_HLS="${VITIS_DIR}"

# Verify paths exist
if [[ ! -d "${VITIS_DIR}/include" ]]; then
    echo "Error: Vitis include directory not found at ${VITIS_DIR}/include"
    exit 1
fi

echo "Found Vitis installation at: ${VITIS_DIR}"
echo "Include directory: ${VITIS_DIR}/include"

# Update the c_cpp_properties.json with the actual paths
cat > .vscode/c_cpp_properties.json << EOF
{
    "configurations": [
        {
            "name": "Linux-Vitis",
            "includePath": [
                "\${workspaceFolder}/**",
                "\${workspaceFolder}/hw/include",
                "${VITIS_DIR}/include",
                "/usr/include",
                "/usr/include/c++/11",
                "/usr/include/x86_64-linux-gnu/c++/11"
            ],
            "defines": [
                "__SYNTHESIS__",
                "__VITIS_HLS__",
                "_GLIBCXX_USE_CXX11_ABI=1"
            ],
            "compilerPath": "/usr/bin/g++",
            "cStandard": "c17",
            "cppStandard": "c++11",
            "intelliSenseMode": "linux-gcc-x64"
        }
    ],
    "version": 4
}
EOF

echo "Updated .vscode/c_cpp_properties.json with Vitis paths"
echo "Please reload the C++ extension in Cursor for changes to take effect."
echo "You can do this by:"
echo "1. Open Command Palette (Ctrl+Shift+P)"
echo "2. Type 'C/C++: Reload IntelliSense Database'"
echo "3. Or restart Cursor" 