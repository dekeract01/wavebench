# Portable Makefile for macOS and Linux
# Supports Homebrew macOS and Fedora/Linux

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    # macOS (Homebrew)
    CXX = /opt/homebrew/bin/g++-16
    CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra -fopenmp \
               -I/opt/homebrew/include
    LDFLAGS = -L/opt/homebrew/lib -Wl,-rpath,/opt/homebrew/lib -lhdf5

else ifeq ($(UNAME_S),Linux)
    # Linux (Fedora)
    CXX = h5c++
    CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra -fopenmp
    LDFLAGS =

else
    $(error Unsupported operating system: $(UNAME_S))
endif


TARGET = wave_openmp
SRC = wave_openmp.cpp

OUTPUT = wave_solution_openmp.h5
LOG = wave_openmp.out.txt


all: $(TARGET)


$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)


run: $(TARGET)
	@./$(TARGET) > $(LOG); status=$$?; cat $(LOG); exit $$status


run1: $(TARGET)
	@OMP_NUM_THREADS=1 ./$(TARGET)

run2: $(TARGET)
	@OMP_NUM_THREADS=2 ./$(TARGET)

run4: $(TARGET)
	@OMP_NUM_THREADS=4 ./$(TARGET)

run8: $(TARGET)
	@OMP_NUM_THREADS=8 ./$(TARGET)

run16: $(TARGET)
	@OMP_NUM_THREADS=16 ./$(TARGET)


scaling: $(TARGET)
	@echo "Running scaling test..."
	@for t in 1 2 4 8 16; do \
		echo ""; \
		echo "=== $$t threads ==="; \
		OMP_NUM_THREADS=$$t ./$(TARGET) | grep "total_time:"; \
	done


scaling-full: $(TARGET)
	@echo "Threads | Time (s) | Speedup | Efficiency"
	@echo "--------|----------|---------|------------"
	@T1=$$(OMP_NUM_THREADS=1 ./$(TARGET) 2>&1 | grep "total_time:" | awk '{print $$2}'); \
	for t in 1 2 4 8 16; do \
		T=$$(OMP_NUM_THREADS=$$t ./$(TARGET) 2>&1 | grep "total_time:" | awk '{print $$2}'); \
		S=$$(echo "scale=3; $$T1/$$T" | bc); \
		E=$$(echo "scale=2; $$S/$$t*100" | bc); \
		printf "%7d | %8.3f | %7.2f | %8.1f%%\n" $$t $$T $$S $$E; \
	done


bench: $(TARGET)
	@for t in 1 2 4 8 16; do \
		echo ""; \
		echo "=== $$t threads ==="; \
		OMP_NUM_THREADS=$$t ./$(TARGET); \
	done


info:
	@echo "Operating system: $(UNAME_S)"
	@echo "Compiler: $(CXX)"
	@$(CXX) --version | head -1
	@echo "CPU cores:"
	@command -v nproc >/dev/null && nproc || sysctl -n hw.ncpu


clean:
	rm -f $(TARGET) $(OUTPUT) $(LOG) wave_openmp_*.out.txt *.png *.pdf


test: all run


help:
	@echo "Targets:"
	@echo "  make              Build"
	@echo "  make run          Run"
	@echo "  make run1/2/4/8/16 Run with threads"
	@echo "  make scaling      Scaling test"
	@echo "  make scaling-full Detailed scaling"
	@echo "  make clean        Remove files"
	@echo "  make info         System information"


.PHONY: all run run1 run2 run4 run8 run16 scaling scaling-full bench info clean test help