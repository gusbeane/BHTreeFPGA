# Makefile for HLS BHTree Testbench

# Paths
HLS_ROOT := ..
VITIS_HLS := $(XILINX_HLS)/bin/vitis_hls

# Compiler settings
CXX := g++
CXXFLAGS := -std=c++11 -Wall -O2 -g

# Include paths for manual test (no sw dependencies)
INCLUDES := -I$(HLS_ROOT)/include \
           -I$(XILINX_HLS)/include

# Source files  
HLS_SOURCES := $(HLS_ROOT)/src/bhtree_hls.cpp
TEST_SOURCES := test_bhtree_manual.cpp

# Output
TARGET := test_bhtree_manual
LOG_FILE := test_results.log

# Default target
all: $(TARGET)

# Build the manual testbench (no sw dependencies)
$(TARGET): $(TEST_SOURCES) $(HLS_SOURCES)
	@echo "Building HLS BHTree manual testbench..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) \
		$(TEST_SOURCES) $(HLS_SOURCES) \
		-o $(TARGET)
	@echo "Build complete: $(TARGET)"

# Run the manual testbench
run: $(TARGET)
	@echo "Running HLS BHTree manual testbench..."
	./$(TARGET) | tee $(LOG_FILE)
	@echo "Results saved to $(LOG_FILE)"

# Clean up
clean:
	rm -f $(TARGET) $(LOG_FILE)
	rm -rf *.o

# Run with Vitis HLS C simulation
hls-csim: 
	@echo "Creating HLS project for C simulation..."
	@mkdir -p hls_csim
	@cd hls_csim && \
	echo "open_project -reset test_bhtree" > run_csim.tcl && \
	echo "set_top create_bhtree_kernel" >> run_csim.tcl && \
	echo "add_files $(HLS_ROOT)/src/bhtree_hls.cpp" >> run_csim.tcl && \
	echo "add_files -tb ../test_bhtree_manual.cpp" >> run_csim.tcl && \
	echo "open_solution solution1" >> run_csim.tcl && \
	echo "set_part xcu200-fsgd2104-2-e" >> run_csim.tcl && \
	echo "create_clock -period 10ns" >> run_csim.tcl && \
	echo "csim_design" >> run_csim.tcl && \
	echo "exit" >> run_csim.tcl && \
	$(VITIS_HLS) -f run_csim.tcl

# Help
help:
	@echo "Available targets:"
	@echo "  all       - Build the manual testbench"
	@echo "  run       - Build and run the manual testbench"
	@echo "  hls-csim  - Run HLS C simulation"
	@echo "  clean     - Clean up build artifacts"
	@echo "  help      - Show this help message"

.PHONY: all run clean hls-csim help 