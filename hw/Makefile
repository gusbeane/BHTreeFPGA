# Makefile for Vitis HLS multi-kernel build
# Place this Makefile at your project root directory
# Usage: make [target] KERNEL=<kernel_name>
# Example: make syn KERNEL=treecon
# Available targets: csim, syn, cosim, package, impl, link, clean

# Default kernel if none specified
# You can override this by running: make csim KERNEL=my_other_kernel
KERNEL ?= treecon

# Variables - these store values we'll reuse throughout the Makefile
# $(shell ...) runs a shell command and captures its output
SCRIPT_DIR := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
HW_DIR := $(shell dirname $(SCRIPT_DIR))/hw
CONFIG_DIR := $(HW_DIR)/config/$(KERNEL)
BUILD_DIR := $(HW_DIR)/build/$(KERNEL)

# Default target - this runs when you just type "make" with no arguments
# The first target in a Makefile is automatically the default
csim: 
	@echo "Building $(KERNEL) kernel - C simulation..."
	@echo "HW Directory: $(HW_DIR)"
	@echo "Config Directory: $(CONFIG_DIR)"
	@echo "Build Directory: $(BUILD_DIR)"
	@if [ ! -d "$(CONFIG_DIR)" ]; then \
		echo "Error: Config directory $(CONFIG_DIR) does not exist!"; \
		echo "Available kernels:"; \
		ls $(HW_DIR)/config/ 2>/dev/null || echo "No config directories found"; \
		exit 1; \
	fi
	cd $(CONFIG_DIR) && vitis-run --mode hls --csim --config hls_config.cfg --work_dir $(BUILD_DIR)
	@echo "$(KERNEL) C simulation completed successfully!"

# Synthesis target
syn:
	@echo "Building $(KERNEL) kernel - synthesis..."
	@if [ ! -d "$(CONFIG_DIR)" ]; then \
		echo "Error: Config directory $(CONFIG_DIR) does not exist!"; \
		exit 1; \
	fi
	cd $(CONFIG_DIR) && v++ -c --mode hls --config hls_config.cfg --work_dir $(BUILD_DIR)
	@echo "$(KERNEL) synthesis completed successfully!"

# C/RTL co-simulation target  
cosim:
	@echo "Building $(KERNEL) kernel - C/RTL co-simulation..."
	@if [ ! -d "$(CONFIG_DIR)" ]; then \
		echo "Error: Config directory $(CONFIG_DIR) does not exist!"; \
		exit 1; \
	fi
	cd $(CONFIG_DIR) && vitis-run --mode hls --cosim --config hls_config.cfg --work_dir $(BUILD_DIR)
	@echo "$(KERNEL) co-simulation completed successfully!"

# Package target
package:
	@echo "Building $(KERNEL) kernel - package..."
	@if [ ! -d "$(CONFIG_DIR)" ]; then \
		echo "Error: Config directory $(CONFIG_DIR) does not exist!"; \
		exit 1; \
	fi
	cd $(CONFIG_DIR) && vitis-run --mode hls --package --config hls_config.cfg --work_dir $(BUILD_DIR)
	@echo "$(KERNEL) package completed successfully!"

# Implementation target
impl:
	@echo "Building $(KERNEL) kernel - implementation..."
	@if [ ! -d "$(CONFIG_DIR)" ]; then \
		echo "Error: Config directory $(CONFIG_DIR) does not exist!"; \
		exit 1; \
	fi
	cd $(CONFIG_DIR) && vitis-run --mode hls --impl --config hls_config.cfg --work_dir $(BUILD_DIR)
	@echo "$(KERNEL) implementation completed successfully!"

# Link target
link:
	@echo "Building $(KERNEL) kernel - link..."
	v++ -l $(BUILD_DIR)/create_$(KERNEL)_kernel.xo \
		--platform xilinx_u200_gen3x16_xdma_2_202110_1 \
		-o $(BUILD_DIR)/create_$(KERNEL)_kernel.xclbin
	@echo "$(KERNEL) link completed successfully!"

# Clean target - removes build artifacts for specified kernel
clean:
	@echo "Cleaning $(KERNEL) build directory..."
	rm -rf $(BUILD_DIR)
	@echo "$(KERNEL) clean completed!"

# Clean all kernels
clean-all:
	@echo "Cleaning all build directories..."
	rm -rf $(HW_DIR)/build/
	@echo "All kernels cleaned!"

# List available kernels (based on config directories)
list-kernels:
	@echo "Available kernels:"
	@ls $(HW_DIR)/config/ 2>/dev/null || echo "No config directories found"

# Convenience targets for common kernels (add your kernel names here)
# These let you run "make treecon-syn" instead of "make syn KERNEL=treecon"
treecon-csim:
	@$(MAKE) csim KERNEL=treecon

treecon-syn:
	@$(MAKE) syn KERNEL=treecon

treecon-cosim:
	@$(MAKE) cosim KERNEL=treecon

# Add more convenience targets for other kernels as needed:
# example_kernel-csim:
#	@$(MAKE) csim KERNEL=example_kernel
#
# example_kernel-syn:
#	@$(MAKE) syn KERNEL=example_kernel

# Help target - shows available commands
help:
	@echo "Multi-kernel Vitis HLS Makefile"
	@echo ""
	@echo "Usage: make [target] KERNEL=<kernel_name>"
	@echo ""
	@echo "Available targets:"
	@echo "  csim    - Run C simulation (default)"
	@echo "  syn     - Run synthesis" 
	@echo "  cosim   - Run C/RTL co-simulation"
	@echo "  package - Run package"
	@echo "  impl    - Run implementation"
	@echo "  link    - Run link"
	@echo "  clean   - Remove build artifacts for specified kernel"
	@echo "  clean-all - Remove all build artifacts"
	@echo "  list-kernels - Show available kernels"
	@echo "  help    - Show this help message"
	@echo ""
	@echo "Examples:"
	@echo "  make syn KERNEL=treecon"
	@echo "  make csim KERNEL=my_other_kernel" 
	@echo "  make treecon-syn    (convenience target)"
	@echo ""
	@echo "Default kernel: $(KERNEL)"

# .PHONY tells Make these targets don't create files with the same name
.PHONY: csim syn cosim package impl link clean clean-all list-kernels help treecon-csim treecon-syn treecon-cosim