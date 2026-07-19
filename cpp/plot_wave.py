#!/usr/bin/env python
"""
    plot_wave.py

    plots the numerical and analytical wave solutions, along with the
    absolute numerical error, from a wave solver output file.

    input: none
    output: result_wave_comparison_<implementation>.png
    author: dekeract01, 2021
    modified by: dekeract01, 2026
    university of southampton
"""

import matplotlib.pyplot as plt
import numpy as np

data_file = "wave_solution_cpp.dat"
implementation = data_file.removeprefix("wave_solution_").removesuffix(".dat")
output_file = f"result_wave_comparison_{implementation}.png"


def read_solution(filename):
    """Return x, numerical solution, analytical solution, and absolute error."""
    data = np.loadtxt(filename, comments="#")
    if data.ndim != 2 or data.shape[1] != 4:
        raise ValueError(
            f"{filename} must contain four columns: x, phi_numerical, "
            "phi_exact, error"
        )
    return data.T


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