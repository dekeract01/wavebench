# wave_stdpar Makefile
# backend: serial (gfortran) | multicore (nvfortran CPU) | parallel (nvfortran GPU)
#   make                          build (default: multicore)
#   make backend=serial run       serial reference via gfortran
#   make backend=multicore run    CPU threads
#   make backend=parallel run     GPU (A4000), writes wave_stdpar.out.txt
# gfortran cannot target the GPU; stdpar multicore/parallel need nvfortran (NVHPC).

SHELL := /bin/bash

# gpu_arch:  RTX A4000 = cc86  (A100 = cc80, Hopper = cc90)
# cuda_home: external CUDA toolkit for nvfortran to use, overriding the newer
#   CUDA that NVHPC 26.5 bundles. Driver 545 supports 12.3, so point at that.
#   Passed via NVHPC_CUDA_HOME -- do NOT also add -gpu=cudaX.Y (that forces the
#   SDK-bundled toolkit and ignores this).
backend   ?= multicore
gpu_arch  ?= cc86
cuda_home ?= /usr/local/cuda-12.3

ifeq ($(backend),serial)
    FC     = gfortran
    FFLAGS = -O3 -march=native -Wall
    CUDA_ENV =
else ifeq ($(backend),multicore)
    FC     = nvfortran
    FFLAGS = -O3 -stdpar=multicore
    CUDA_ENV =
else ifeq ($(backend),parallel)
    FC     = nvfortran
    FFLAGS = -O3 -stdpar=gpu -gpu=$(strip $(gpu_arch)) -Minfo=accel
    CUDA_ENV = NVHPC_CUDA_HOME=$(strip $(cuda_home))
else
    $(error Unknown backend '$(backend)': use serial, multicore, or parallel)
endif

TARGET = wave_stdpar
SRC    = wave_stdpar.f90
OUTPUT = wave_solution_stdpar.bin
LOG    = wave_stdpar.out.txt

all: $(TARGET)

$(TARGET): $(SRC)
	$(CUDA_ENV) $(FC) $(FFLAGS) -o $(TARGET) $(SRC)

run: $(TARGET)
	set -o pipefail; ./$(TARGET) 2>&1 | tee $(LOG)

clean:
	rm -f $(TARGET) $(OUTPUT) $(LOG) *.mod *.o

test: all run

help:
	@echo "backends: serial | multicore | parallel   (default: multicore)"
	@echo "  make                          build (backend=multicore)"
	@echo "  make backend=serial run       serial reference (gfortran)"
	@echo "  make backend=multicore run    CPU threads (nvfortran)"
	@echo "  make backend=parallel run     GPU, writes $(LOG)"
	@echo "  make backend=parallel gpu_arch=cc90 cuda_home=/usr/local/cuda-12.6 run  # another box"
	@echo "  make clean                    remove build + output artifacts"

.PHONY: all run clean test help