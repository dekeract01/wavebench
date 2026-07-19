CARGO ?= $(HOME)/.cargo/bin/cargo
TARGET = target/release/wave
OUTPUT = wave_solution_rust.h5
LOG = wave_rust.out.txt

export HDF5_DIR = /opt/homebrew
export PKG_CONFIG_PATH = /opt/homebrew/lib/pkgconfig

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

help:
	@echo "Available targets:"
	@echo "  make       - Compile the Rust HDF5 solver"
	@echo "  make run   - Compile, run, and save wave_rust.out.txt"
	@echo "  make test  - Same as make run"
	@echo "  make bench - Run five benchmark iterations"
	@echo "  make clean - Remove build and generated output"

.PHONY: all run clean test help bench