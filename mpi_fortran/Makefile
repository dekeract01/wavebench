# MPI Fortran compiler
MPIFC = /opt/homebrew/bin/mpif90
# Or use: mpifort

# Compiler flags
FFLAGS = -O3 -march=native -Wall -I/opt/homebrew/include

# Homebrew HDF5 Fortran libraries
LDFLAGS = -L/opt/homebrew/lib -Wl,-rpath,/opt/homebrew/lib -lhdf5_fortran -lhdf5

# Target executable
TARGET = wave_mpi

# Source file
SRC = wave_mpi.f90

# Default target
all: $(TARGET)

# Build executable
$(TARGET): $(SRC)
	$(MPIFC) $(FFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

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

# Benchmark scaling
scaling: $(TARGET)
	@echo "Running scaling test..."
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

# Clean build artifacts
clean:
	rm -f $(TARGET) $(OUTPUT) $(LOG) *.mod *.o

# Build and run
test: all run

# Help message
help:
	@echo "Available targets:"
	@echo "  make          - Compile the program"
	@echo "  make run      - Run with 4 processes (default)"
	@echo "  make run1     - Run with 1 process"
	@echo "  make run2     - Run with 2 processes"
	@echo "  make run4     - Run with 4 processes"
	@echo "  make run8     - Run with 8 processes"
	@echo "  make scaling  - Run scaling benchmark"
	@echo "  make clean    - Remove executable and output files"
	@echo ""
	@echo "Examples:"
	@echo "  make run4             # Run with 4 MPI processes"
	@echo "  make scaling          # Test parallel scaling"
	@echo "  mpirun -np 16 ./wave_mpi  # Custom process count"

.PHONY: all run run1 run2 run4 run8 scaling clean test help