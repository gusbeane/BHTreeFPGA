# Makefile for Multi-Kernel HLS BHTree Project

# Paths
HW_ROOT := ..
SCRIPTS_DIR := $(HW_ROOT)/scripts
BUILD_DIR := $(HW_ROOT)/build

# Available kernels
KERNELS := bhtree

# Default target
.PHONY: all
all: csim

# C Simulation targets
.PHONY: csim csim-bhtree
csim: csim-bhtree

csim-bhtree:
	@echo "Running C simulation for bhtree..."
	$(SCRIPTS_DIR)/build_bhtree.sh csim

# Synthesis targets
.PHONY: syn syn-bhtree
syn: syn-bhtree

syn-bhtree:
	@echo "Running synthesis for bhtree..."
	$(SCRIPTS_DIR)/build_bhtree.sh syn

# C/RTL Co-simulation targets
.PHONY: cosim cosim-bhtree
cosim: cosim-bhtree

cosim-bhtree:
	@echo "Running C/RTL co-simulation for bhtree..."
	$(SCRIPTS_DIR)/build_bhtree.sh cosim

# Package targets
.PHONY: package package-bhtree
package: package-bhtree

package-bhtree:
	@echo "Running package for bhtree..."
	$(SCRIPTS_DIR)/build_bhtree.sh package

# Build all kernels with specific operation
.PHONY: all-csim all-syn all-cosim all-package
all-csim:
	$(SCRIPTS_DIR)/build_all.sh csim

all-syn:
	$(SCRIPTS_DIR)/build_all.sh syn

all-cosim:
	$(SCRIPTS_DIR)/build_all.sh cosim

all-package:
	$(SCRIPTS_DIR)/build_all.sh package

# Clean targets
.PHONY: clean clean-bhtree clean-all
clean: clean-all

clean-bhtree:
	@echo "Cleaning bhtree build artifacts..."
	rm -rf $(BUILD_DIR)/bhtree

clean-all:
	@echo "Cleaning all build artifacts..."
	rm -rf $(BUILD_DIR)/*

# Help target
.PHONY: help
help:
	@echo "Multi-Kernel HLS BHTree Build System"
	@echo "===================================="
	@echo ""
	@echo "Available targets:"
	@echo "  all          - Run C simulation for all kernels (default)"
	@echo "  csim         - Run C simulation for all kernels"
	@echo "  syn          - Run synthesis for all kernels"
	@echo "  cosim        - Run C/RTL co-simulation for all kernels"
	@echo "  package      - Run package for all kernels"
	@echo ""
	@echo "Per-kernel targets:"
	@echo "  csim-bhtree  - Run C simulation for bhtree kernel"
	@echo "  syn-bhtree   - Run synthesis for bhtree kernel"
	@echo "  cosim-bhtree - Run C/RTL co-simulation for bhtree kernel"
	@echo "  package-bhtree - Run package for bhtree kernel"
	@echo ""
	@echo "Batch targets:"
	@echo "  all-csim     - Run C simulation for all kernels"
	@echo "  all-syn      - Run synthesis for all kernels"
	@echo "  all-cosim    - Run C/RTL co-simulation for all kernels"
	@echo "  all-package  - Run package for all kernels"
	@echo ""
	@echo "Clean targets:"
	@echo "  clean        - Clean all build artifacts"
	@echo "  clean-bhtree - Clean bhtree build artifacts"
	@echo ""
	@echo "Available kernels: $(KERNELS)" 