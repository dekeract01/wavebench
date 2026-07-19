#!/usr/bin/env python
"""
plot_wavef90.py

plots the numerical and analytical wave solutions, along with the
absolute numerical error, from a fortran stream-binary output file.

input: wave_solution_fortran.bin
output: result_wave_comparison_fortran.png
author: dekeract01, 2022
university of southampton
"""

import matplotlib.pyplot as plt
import numpy as np

# now, this is just an attempt at reading the binary file, one can just output .dat file
# but I want to see if I can read the binary file as a practice
data_file = "wave_solution_fortran.bin"
implementation = data_file.removeprefix("wave_solution_").removesuffix(".bin")
output_file = f"result_wave_comparison_{implementation}.png"


def read_solution(filename):
    """Return x, numerical solution, analytical solution, and absolute error."""
    data = np.fromfile(filename, dtype=np.float64)
    if data.size == 0 or data.size % 4 != 0:
        raise ValueError(
            f"{filename} must contain groups of four float64 values: x, "
            "phi_numerical, phi_exact, error"
        )
    return data.reshape(-1, 4).T


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