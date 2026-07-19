# Compiler selection
FC = gfortran
# FC = ifort

# Compiler flags
ifeq ($(FC),gfortran)
	FFLAGS = -O3 -march=native -Wall
else ifeq ($(FC),ifort)
	FFLAGS = -O3 -xHost
endif

# Target executable
TARGET = wave

# Source file
SRC = wave.f90

# Generated output
OUTPUT = wave_solution_fortran.bin
LOG = wave_fortran.out.txt

# Default target
all: $(TARGET)

# Build executable
$(TARGET): $(SRC)
	$(FC) $(FFLAGS) -o $(TARGET) $(SRC)

# Run the simulation
run: $(TARGET)
	@./$(TARGET) > $(LOG); status=$$?; cat $(LOG); exit $$status

# Clean build artifacts
clean:
	rm -f $(TARGET) $(OUTPUT) $(LOG) *.mod *.o

# Build and run
test: all run

# Help message
help:
	@echo "Available targets:"
	@echo "  make           - Compile the program (default: gfortran)"
	@echo "  make run       - Compile, run, and save wave_fortran.out.txt"
	@echo "  make test      - Same as 'make run'"
	@echo "  make clean     - Remove executable and output files"
	@echo "  make FC=ifort  - Compile using Intel Fortran compiler"
	@echo ""
	@echo "Examples:"
	@echo "  make                    # Compile with gfortran"
	@echo "  make FC=ifort           # Compile with ifort"
	@echo "  make run                # Compile and run"
	@echo "  make clean              # Clean up"

.PHONY: all run clean test help