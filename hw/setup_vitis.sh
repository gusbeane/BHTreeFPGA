#!/usr/bin/env bash
#
# ~/.vitis/setup_vitis.sh

# Only bail out of this file—never the entire shell
# -----------------------------------------------
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
  _SOURCED=1
else
  echo "Error: please source this script, not execute it."
  exit 1
fi

# Helper: print an error and return from the script
_error() {
  echo "Error: $1" >&2
  return 1
}

# 1) catch undefined variables, but do NOT exit on every failed command
# Note: We don't use set -u because Vitis scripts may reference undefined variables
# set -u

# 2) ensure XILINX_HOME is set
export XILINX_HOME="/opt/xilinx"
[[ -n "${XILINX_HOME:-}" ]] || _error "XILINX_HOME is not set."

# 2.5) Set commonly undefined variables to avoid errors
export PYTHONPATH="${PYTHONPATH:-}"
export MATLABPATH="${MATLABPATH:-}"

# 3) pin version & compute paths
VITIS_VERSION="2025.1"
VITIS_DIR="${XILINX_HOME}/${VITIS_VERSION}/Vitis"
XRT_DIR="${XILINX_HOME}/xrt"

# 4) source Vitis settings
if [[ -r "${VITIS_DIR}/settings64.sh" ]]; then
  source "${VITIS_DIR}/settings64.sh" || _error "failed to source settings64.sh"
else
  _error "cannot find ${VITIS_DIR}/settings64.sh"
fi

# 5) source XRT setup
if [[ -r "${XRT_DIR}/setup.sh" ]]; then
  source "${XRT_DIR}/setup.sh" || _error "failed to source xrt/setup.sh"
else
  _error "cannot find ${XRT_DIR}/setup.sh"
fi

# 6) export your platform path
export PLATFORM_REPO_PATHS="${XILINX_HOME}/platforms/xilinx_u200_gen3x16_xdma_2_202110_1"
export GIT_EXEC_PATH=${XILINX_VITIS}/tps/lnx64/git-2.45.0
