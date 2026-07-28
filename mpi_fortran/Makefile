# ================================================================
# Portable MPI Fortran Makefile
# macOS (Homebrew) + Linux (Fedora/OpenMPI)
# ================================================================

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)

    # ------------------------------------------------------------
    # macOS + Homebrew
    # ------------------------------------------------------------

    MPIFC = /opt/homebrew/bin/mpif90

    FFLAGS = -O3 -march=native -Wall \
             -I/opt/homebrew/include

    LDFLAGS = -L/opt/homebrew/lib \
              -Wl,-rpath,/opt/homebrew/lib \
              -lhdf5_fortran -lhdf5


else ifeq ($(UNAME_S),Linux)

    # ------------------------------------------------------------
    # Linux + OpenMPI + Fedora HDF5
    # ------------------------------------------------------------

    MPIFC = mpif90

    # Find HDF5 Fortran module directory automatically
    HDF5_MOD_DIR := $(shell find /usr/lib64/gfortran/modules -path "*openmpi*" -name hdf5.mod -printf '%h' 2>/dev/null | head -1)

    FFLAGS = -O3 -march=native -Wall \
             -I$(HDF5_MOD_DIR)

    LDFLAGS = -L/usr/lib64/openmpi/lib \
              -lhdf5_fortran -lhdf5


else

    $(error Unsupported operating system: $(UNAME_S))

endif


# ================================================================
# Program settings
# ================================================================

TARGET = wave_mpi
SRC = wave_mpi.f90

LOG = wave_mpi.out.txt
OUTPUT = wave_solution_mpi.h5


# ================================================================
# Build
# ================================================================

all: $(TARGET)


$(TARGET): $(SRC)
	$(MPIFC) $(FFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)


# ================================================================
# Run
# ================================================================

run: $(TARGET)
	@mpirun -np 4 ./$(TARGET) > $(LOG); \
	status=$$?; \
	cat $(LOG); \
	exit $$status


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


# ================================================================
# Scaling benchmark
# ================================================================

scaling: $(TARGET)
	@echo "Running MPI scaling test..."
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

	@echo ""
	@echo "16 processes:"
	@mpirun -np 16 ./$(TARGET) | grep "Total time"



scaling-full: $(TARGET)

	@echo "Processes | Time | Speedup"
	@echo "----------|------|--------"

	@T1=$$(mpirun -np 1 ./$(TARGET) 2>&1 | grep "Total time" | awk '{print $$3}'); \
	for p in 1 2 4 8 16; do \
		T=$$(mpirun -np $$p ./$(TARGET) 2>&1 | grep "Total time" | awk '{print $$3}'); \
		S=$$(echo "scale=3; $$T1/$$T" | bc); \
		printf "%9d | %4.3f | %7.3f\n" $$p $$T $$S; \
	done


# ================================================================
# Information
# ================================================================

info:
	@echo "Operating system: $(UNAME_S)"
	@echo "MPI compiler: $(MPIFC)"
	@echo "HDF5 module directory: $(HDF5_MOD_DIR)"
	@echo ""
	@echo "MPI version:"
	@mpirun --version | head -1
	@echo ""
	@echo "Fortran compiler:"
	@$(MPIFC) --version | head -1


# ================================================================
# Cleaning
# ================================================================

clean:
	rm -f $(TARGET) $(OUTPUT) $(LOG)
	rm -f *.mod *.o


# ================================================================
# Build and run
# ================================================================

test: all run


# ================================================================
# Help
# ================================================================

help:
	@echo "Targets:"
	@echo "  make              Build MPI Fortran program"
	@echo "  make run          Run with 4 MPI processes"
	@echo "  make run1         Run with 1 process"
	@echo "  make run2         Run with 2 processes"
	@echo "  make run4         Run with 4 processes"
	@echo "  make run8         Run with 8 processes"
	@echo "  make run16        Run with 16 processes"
	@echo "  make scaling      MPI scaling test"
	@echo "  make scaling-full Detailed scaling"
	@echo "  make info         Show compiler information"
	@echo "  make clean        Remove build files"


.PHONY: all run run1 run2 run4 run8 run16 \
        scaling scaling-full info clean test help