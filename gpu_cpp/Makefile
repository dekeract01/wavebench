# HDF5 C++ compiler wrapper
CXX = h5c++

# Compiler flags for macOS
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra

# OpenCL framework (macOS)
LDFLAGS = -framework OpenCL

# For Linux, use instead:
# LDFLAGS = -lOpenCL

# Target executable
TARGET = wave_opencl

# Source file
SRC = wave_opencl.cpp
KERNEL = wave_opencl.cl

# Output file
OUTPUT = wave_solution_opencl.h5
LOG = wave_opencl.out.txt

# Default target
all: $(TARGET)

# Build executable
$(TARGET): $(SRC) $(KERNEL)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# Run the simulation
run: $(TARGET)
	@./$(TARGET) > $(LOG); status=$$?; cat $(LOG); exit $$status

# Clean build artifacts
clean:
	rm -f $(TARGET) $(OUTPUT) $(LOG) *.png *.pdf

# Build and run
test: all run

# Show OpenCL devices
info:
	@echo "Checking OpenCL availability..."
	@clinfo 2>/dev/null || echo "clinfo not installed (brew install clinfo)"

# Help message
help:
	@echo "Available targets:"
	@echo "  make       - Compile OpenCL program"
	@echo "  make run   - Compile, run, and save wave_opencl.out.txt"
	@echo "  make test  - Same as 'make run'"
	@echo "  make info  - Show OpenCL devices"
	@echo "  make clean - Remove executable and output"
	@echo ""
	@echo "Note: On macOS, OpenCL is built-in"
	@echo "      On Linux, install: sudo apt-get install opencl-headers ocl-icd-opencl-dev"

.PHONY: all run clean test info help