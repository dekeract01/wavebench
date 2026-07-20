# CUDA compiler
NVCC = nvcc

# GPU architecture (RTX A4000 = Ampere, compute capability 8.6)
ARCH = sm_86

# Compiler flags
# -O3 for device code; -Xcompiler passes flags to the host compiler (g++)
NVCCFLAGS = -std=c++17 -O3 -arch=$(ARCH) -Xcompiler -O3,-march=native

# Target executable
TARGET = wave_cuda

# Source file
SRC = wave.cu

# Output file
OUTPUT = wave_solution_cuda.dat

# Benchmark report
LOG = wave_cuda.out.txt

# Default target
all: $(TARGET)

# Build executable
$(TARGET): $(SRC)
	$(NVCC) $(NVCCFLAGS) -o $(TARGET) $(SRC)

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

# Build for the GPU detected in this machine (CUDA >= 11.6)
native:
	$(NVCC) -std=c++17 -O3 -arch=native -Xcompiler -O3,-march=native -o $(TARGET) $(SRC)

# Show per-kernel register/memory usage from ptxas
profile-ptxas:
	$(NVCC) $(NVCCFLAGS) -Xptxas -v -o $(TARGET) $(SRC)

# Line-info build for use with Nsight Compute / Nsight Systems
profile-nsight:
	$(NVCC) $(NVCCFLAGS) -lineinfo -o $(TARGET) $(SRC)

# Profile a run with Nsight Systems (timeline) if installed
nsys: profile-nsight
	nsys profile --stats=true -o wave_cuda_nsys ./$(TARGET)

# Profile kernels with Nsight Compute if installed
ncu: profile-nsight
	ncu --set full -o wave_cuda_ncu ./$(TARGET)

# Check the toolchain and GPU are visible
check:
	@which $(NVCC) || echo "nvcc not on PATH - try: export PATH=/usr/local/cuda/bin:$$PATH"
	@$(NVCC) --version | tail -1
	@nvidia-smi -L

# Help message
help:
	@echo "Available targets:"
	@echo "  make                - Compile for $(ARCH) (RTX A4000)"
	@echo "  make run            - Compile, run, and save wave_cuda.out.txt"
	@echo "  make test           - Same as 'make run'"
	@echo "  make bench          - Run 5 times for benchmarking"
	@echo "  make native         - Compile for whatever GPU is in this machine"
	@echo "  make profile-ptxas  - Compile with register/shared-mem usage report"
	@echo "  make nsys           - Profile timeline with Nsight Systems"
	@echo "  make ncu            - Profile kernels with Nsight Compute"
	@echo "  make check          - Verify nvcc and GPU are visible"
	@echo "  make clean          - Remove executable and generated output"
	@echo ""
	@echo "Examples:"
	@echo "  make                # Compile"
	@echo "  make run            # Compile and run"
	@echo "  make ARCH=sm_80 run # Target a different GPU (e.g. A100)"

.PHONY: all run clean test bench native profile-ptxas profile-nsight nsys ncu check help