# Homebrew GCC with OpenMP support
CXX = /opt/homebrew/bin/g++-16

# Compiler flags with OpenMP
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra -fopenmp -I/opt/homebrew/include

# Homebrew HDF5 C library
LDFLAGS = -L/opt/homebrew/lib -Wl,-rpath,/opt/homebrew/lib -lhdf5

# Target executable
TARGET = wave_openmp

# Source file
SRC = wave_openmp.cpp

# Output file
OUTPUT = wave_solution_openmp.h5
LOG = wave_openmp.out.txt

# Default target
all: $(TARGET)

# Build executable
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# Run with default number of threads
run: $(TARGET)
	@./$(TARGET) > $(LOG); status=$$?; cat $(LOG); exit $$status

# Run with specific number of threads
run1: $(TARGET)
	@OMP_NUM_THREADS=1 ./$(TARGET) > wave_openmp_1.out.txt; status=$$?; cat wave_openmp_1.out.txt; exit $$status

run2: $(TARGET)
	@OMP_NUM_THREADS=2 ./$(TARGET) > wave_openmp_2.out.txt; status=$$?; cat wave_openmp_2.out.txt; exit $$status

run4: $(TARGET)
	@OMP_NUM_THREADS=4 ./$(TARGET) > wave_openmp_4.out.txt; status=$$?; cat wave_openmp_4.out.txt; exit $$status

run8: $(TARGET)
	@OMP_NUM_THREADS=8 ./$(TARGET) > wave_openmp_8.out.txt; status=$$?; cat wave_openmp_8.out.txt; exit $$status

run16: $(TARGET)
	@OMP_NUM_THREADS=16 ./$(TARGET) > wave_openmp_16.out.txt; status=$$?; cat wave_openmp_16.out.txt; exit $$status

# Scaling test
scaling: $(TARGET)
	@echo "Running scaling test..."
	@echo ""
	@echo "1 thread:"
	@OMP_NUM_THREADS=1 ./$(TARGET) | grep "total_time:"
	@echo ""
	@echo "2 threads:"
	@OMP_NUM_THREADS=2 ./$(TARGET) | grep "total_time:"
	@echo ""
	@echo "4 threads:"
	@OMP_NUM_THREADS=4 ./$(TARGET) | grep "total_time:"
	@echo ""
	@echo "8 threads:"
	@OMP_NUM_THREADS=8 ./$(TARGET) | grep "total_time:"

# Detailed scaling analysis
scaling-full: $(TARGET)
	@echo "Detailed scaling analysis..."
	@echo "Threads | Time (s) | Speedup | Efficiency"
	@echo "--------|----------|---------|------------"
	@T1=$$(OMP_NUM_THREADS=1 ./$(TARGET) 2>&1 | grep "total_time:" | tail -1 | awk '{print $$2}'); \
	for t in 1 2 4 8; do \
		T=$$(OMP_NUM_THREADS=$$t ./$(TARGET) 2>&1 | grep "total_time:" | tail -1 | awk '{print $$2}'); \
		S=$$(echo "scale=2; $$T1 / $$T" | bc); \
		E=$$(echo "scale=2; $$S / $$t * 100" | bc); \
		printf "%7d | %8.3f | %7.2f | %9.1f%%\n" $$t $$T $$S $$E; \
	done

# Clean build artifacts
clean:
	rm -f $(TARGET) $(OUTPUT) $(LOG) wave_openmp_*.out.txt *.png *.pdf

# Build and run
test: all run

# Benchmark with different thread counts
bench: $(TARGET)
	@echo "Running benchmark with different thread counts..."
	@for t in 1 2 4 8; do \
		echo ""; \
		echo "=== $$t threads ==="; \
		OMP_NUM_THREADS=$$t ./$(TARGET); \
	done

# Show compiler and CPU info
info:
	@echo "Checking for GCC installations..."
	@ls /opt/homebrew/bin/gcc-* 2>/dev/null || ls /usr/local/bin/gcc-* 2>/dev/null || echo "No GCC found - run: brew install gcc"
	@echo ""
	@echo "Available CPU cores:"
	@sysctl -n hw.ncpu

# Help message
help:
	@echo "Available targets:"
	@echo "  make             - Compile with OpenMP"
	@echo "  make run         - Run with default threads and save wave_openmp.out.txt"
	@echo "  make run1/2/4/8  - Run with specific thread count"
	@echo "  make scaling     - Quick scaling test"
	@echo "  make scaling-full- Detailed scaling analysis"
	@echo "  make bench       - Benchmark multiple thread counts"
	@echo "  make info        - Show compiler and CPU info"
	@echo "  make clean       - Remove executable and output"
	@echo ""
	@echo "Note: Make sure CXX is set to the correct GCC version in the Makefile"
	@echo "      Run 'make info' to see available GCC versions"

.PHONY: all run run1 run2 run4 run8 run16 scaling scaling-full clean test bench info help