# C++ compiler
CXX = g++
# Or use: clang++

# Compiler flags
CXXFLAGS = -std=c++17 -O3 -march=native

# Target executable
TARGET = wave_cpp

# Source file
SRC = wave.cpp

# Output file
OUTPUT = wave_solution_cpp.dat

# Benchmark report
LOG = wave_cpp.out.txt

# Default target
all: $(TARGET)

# Build executable
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Run the simulation
run: $(TARGET)
	@./$(TARGET) > $(LOG); status=$$?; cat $(LOG); exit $$status

# Clean build artifacts
clean:
	rm -f $(TARGET) $(OUTPUT) $(LOG) *.png *.pdf

# Build and run
test: all run

# Benchmark - run multiple times
bench: $(TARGET)
	@echo "Running 5 iterations for benchmark..."
	@for i in 1 2 3 4 5; do \
		echo "Run $$i:"; \
		./$(TARGET); \
		echo ""; \
	done

# Compile with different compilers
gcc: 
	g++ $(CXXFLAGS) -o $(TARGET) $(SRC)

clang:
	clang++ $(CXXFLAGS) -o $(TARGET) $(SRC)

# Profile with optimization report (GCC)
profile-gcc:
	g++ -std=c++17 -O3 -march=native -fopt-info-vec-optimized -o $(TARGET) $(SRC)

# Profile with optimization report (Clang)
profile-clang:
	clang++ -std=c++17 -O3 -march=native -Rpass=loop-vectorize -o $(TARGET) $(SRC)

# Help message
help:
	@echo "Available targets:"
	@echo "  make           - Compile with default compiler"
	@echo "  make run       - Compile, run, and save wave_cpp.out.txt"
	@echo "  make test      - Same as 'make run'"
	@echo "  make bench     - Run 5 times for benchmarking"
	@echo "  make gcc       - Compile with GCC"
	@echo "  make clang     - Compile with Clang"
	@echo "  make clean     - Remove executable and generated output"
	@echo ""
	@echo "Examples:"
	@echo "  make           # Compile"
	@echo "  make run       # Compile and run"
	@echo "  make CXX=clang++ run  # Use Clang"

.PHONY: all run clean test bench gcc clang profile-gcc profile-clang help