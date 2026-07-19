#!/usr/bin/env python
"""
plot_wave.py

plots the numerical and analytical wave solutions, along with the
absolute numerical error, from a rust hdf5 output file.

input: wave_solution_rust.h5
output: result_wave_comparison_rust.png
author: dekeract01, 2025
university of southampton
"""

import h5py
import matplotlib.pyplot as plt

data_file = "wave_solution_rust.h5"
implementation = data_file.removeprefix("wave_solution_").removesuffix(".h5")
output_file = f"result_wave_comparison_{implementation}.png"


def read_solution(filename):
    """Return x, numerical solution, analytical solution, and absolute error."""
    with h5py.File(filename, "r") as input_file:
        required_datasets = ("x", "phi_numerical", "phi_exact", "error")
        missing_datasets = [name for name in required_datasets if name not in input_file]
        if missing_datasets:
            raise ValueError(
                f"{filename} is missing datasets: {', '.join(missing_datasets)}"
            )
        return tuple(input_file[name][:] for name in required_datasets)


def main():
    figure, (solution_axis, error_axis) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True, layout="constrained"
    )

    x, phi_numerical, phi_exact, error = read_solution(data_file)
    solution_axis.plot(x, phi_numerical, linewidth=2, label="numerical")
    solution_axis.plot(
        x,
        phi_exact,
        "--",
        color="black",
        linewidth=1.5,
        label="analytical",
    )
    error_axis.semilogy(x, error, linewidth=2, label="absolute error")

    solution_axis.set_ylabel("phi")
    solution_axis.set_title("one-dimensional wave equation")
    solution_axis.grid(True, alpha=0.3)
    solution_axis.legend()

    error_axis.set_xlabel("x")
    error_axis.set_ylabel("absolute error")
    error_axis.set_title("numerical error")
    error_axis.grid(True, alpha=0.3, which="both")
    error_axis.legend()

    figure.savefig(output_file, dpi=150)
    print(f"plot written to {output_file}")


if __name__ == "__main__":
    main()