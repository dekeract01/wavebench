# MPI C++ compiler wrapper
MPICXX = mpic++
# Alternative: MPICXX = mpicxx

# Compiler flags
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra -I/opt/homebrew/include

# Homebrew HDF5 C library
LDFLAGS = -L/opt/homebrew/lib -Wl,-rpath,/opt/homebrew/lib -lhdf5

# Target executable
TARGET = wave_mpi

# Source file
SRC = wave_mpi.cpp

# Output file
OUTPUT = wave_solution_mpi.h5
LOG = wave_mpi.out.txt

# Default target
all: $(TARGET)

# Build executable
$(TARGET): $(SRC)
	$(MPICXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# Run with different number of processes
run: $(TARGET)
	@mpirun -np 4 ./$(TARGET) > $(LOG); status=$$?; cat $(LOG); exit $$status

run1: $(TARGET)
	mpirun -np 1 ./$(TARGET)

run2: $(TARGET)
	mpirun -np 2 ./$(TARGET)

run4: $(TARGET)
	mpirun -np 4 ./$(TARGET)

run8: $(TARGET)
	mpirun -np 8 ./$(TARGET)

run16: $(TARGET)
	mpirun -np 16 ./$(TARGET)

# Scaling test
scaling: $(TARGET)
	@echo "Running scaling test..."
	@echo ""
	@echo "1 process:"
	@mpirun -np 1 ./$(TARGET) | grep "Total time"
	@echo ""
	@echo "2 processes:"
	@mpirun -np 2 ./$(TARGET) | grep "Total time"
	@echo ""
	@echo "4 processes:"
	@mpirun -np 4 ./$(TARGET) | grep "Total time"
	@echo ""
	@echo "8 processes:"
	@mpirun -np 8 ./$(TARGET) | grep "Total time"

# Detailed scaling analysis
scaling-full: $(TARGET)
	@echo "Detailed scaling analysis..."
	@echo "Processes | Time (s) | Speedup | Efficiency"
	@echo "----------|----------|---------|------------"
	@T1=$$(mpirun -np 1 ./$(TARGET) 2>&1 | grep "Total time" | awk '{print $$3}'); \
	for p in 1 2 4 8; do \
		T=$$(mpirun -np $$p ./$(TARGET) 2>&1 | grep "Total time" | awk '{print $$3}'); \
		S=$$(echo "scale=2; $$T1 / $$T" | bc); \
		E=$$(echo "scale=2; $$S / $$p * 100" | bc); \
		printf "%9d | %8.3f | %7.2f | %9.1f%%\n" $$p $$T $$S $$E; \
	done

# Clean build artifacts
clean:
	rm -f $(TARGET) $(OUTPUT) $(LOG) *.png *.pdf

# Build and run
test: all run

# Benchmark with different process counts
bench: $(TARGET)
	@echo "Running benchmark with different process counts..."
	@for p in 1 2 4 8; do \
		echo ""; \
		echo "=== $$p processes ==="; \
		mpirun -np $$p ./$(TARGET); \
	done

# Show MPI info
info:
	@echo "MPI compiler wrapper: $(MPICXX)"
	@which $(MPICXX) || echo "$(MPICXX) not found"
	@echo ""
	@echo "MPI version:"
	@mpirun --version || echo "mpirun not found"
	@echo ""
	@echo "Available CPU cores:"
	@sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo "Unknown"

# Help message
help:
	@echo "Available targets:"
	@echo "  make             - Compile with MPI"
	@echo "  make run         - Run with 4 processes and save wave_mpi.out.txt"
	@echo "  make run1/2/4/8  - Run with specific process count"
	@echo "  make scaling     - Quick scaling test"
	@echo "  make scaling-full- Detailed scaling analysis"
	@echo "  make bench       - Benchmark multiple process counts"
	@echo "  make info        - Show MPI and CPU info"
	@echo "  make clean       - Remove executable and output"
	@echo ""
	@echo "Examples:"
	@echo "  make run4                  # Run with 4 MPI processes"
	@echo "  mpirun -np 16 ./wave_mpi   # Custom process count"
	@echo "  make scaling               # Test parallel scaling"

.PHONY: all run run1 run2 run4 run8 run16 scaling scaling-full clean test bench info help