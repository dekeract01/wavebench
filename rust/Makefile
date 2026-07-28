CARGO ?= $(HOME)/.cargo/bin/cargo

TARGET = target/release/wave
OUTPUT = wave_solution_rust.h5
LOG = wave_rust.out.txt

# Detect operating system
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    # macOS + Homebrew
    export HDF5_DIR = /opt/homebrew
    export PKG_CONFIG_PATH = /opt/homebrew/lib/pkgconfig

else ifeq ($(UNAME_S),Linux)
    # Linux/Fedora
    export HDF5_DIR = /usr
    export PKG_CONFIG_PATH = /usr/lib64/pkgconfig:/usr/lib/pkgconfig

else
    $(error Unsupported operating system: $(UNAME_S))
endif


all: $(TARGET)


$(TARGET): Cargo.toml wave.rs
	$(CARGO) build --release


run: $(TARGET)
	@./$(TARGET) > $(LOG); status=$$?; cat $(LOG); exit $$status


clean:
	$(CARGO) clean
	rm -f $(OUTPUT) $(LOG)


test: all run


bench: $(TARGET)
	@for i in 1 2 3 4 5; do \
		echo "run $$i:"; \
		./$(TARGET); \
		echo ""; \
	done


info:
	@echo "OS: $(UNAME_S)"
	@echo "HDF5_DIR: $(HDF5_DIR)"
	@echo "PKG_CONFIG_PATH: $(PKG_CONFIG_PATH)"
	@echo "Cargo: $(CARGO)"


help:
	@echo "Available targets:"
	@echo "  make       - Compile the Rust HDF5 solver"
	@echo "  make run   - Run solver"
	@echo "  make test  - Build and run"
	@echo "  make bench - Run benchmark"
	@echo "  make info  - Show configuration"
	@echo "  make clean - Remove build/output files"


.PHONY: all run clean test bench info help