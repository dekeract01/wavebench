#!/usr/bin/env python
"""
    plot_wave.py

    plots all available one-dimensional wave solver results.

    input: wave_solution_*.h5, wave_solution_*.dat, or wave_solution_*.bin
    output: result_wave_comparison_<implementation>.pdf beside each input file
    author: dekeract01, 2021
    modified by: dekeract01, 2026
    university of southampton
"""

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy

root_directory = Path(__file__).resolve().parent
required_datasets = ("x", "phi_numerical", "phi_exact", "error")
format_priority = {".h5": 0, ".dat": 1, ".bin": 2}


def read_solution(data_file):
    """Return x, numerical solution, analytical solution, and absolute error."""
    if data_file.suffix == ".h5":
        with h5py.File(data_file, "r") as input_file:
            missing_datasets = [name for name in required_datasets if name not in input_file]
            if missing_datasets:
                raise ValueError(f"{data_file} is missing datasets: {', '.join(missing_datasets)}")
            return tuple(input_file[name][:] for name in required_datasets)

    if data_file.suffix == ".dat":
        data = numpy.loadtxt(data_file, comments="#")  # Use the correct path to your data file
    elif data_file.suffix == ".bin":
        data = numpy.fromfile(data_file, dtype=numpy.float64)
        if data.size % 4 != 0:
            raise ValueError(f"{data_file} must contain groups of four float64 values")
        data = data.reshape(-1, 4)
    else:
        raise ValueError(f"unsupported solution format: {data_file.suffix}")

    if data.ndim != 2 or data.shape[1] != 4:
        raise ValueError(f"{data_file} must contain four columns: x, phi_numerical, phi_exact, error")
    return data.T


def find_solution_files():
    """Find one preferred solution file in each solver directory."""
    solution_files = []
    for directory in sorted(path for path in root_directory.iterdir() if path.is_dir()):
        candidates = [path for path in directory.glob("wave_solution_*") if path.suffix in format_priority]
        if candidates:
            solution_files.append(min(candidates, key=lambda path: (format_priority[path.suffix], path.name)))
    return solution_files


solution_files = find_solution_files()
if not solution_files:
    raise FileNotFoundError("no wave_solution_* result files found")

for data_file in solution_files:

    x, phi_num, phi_ex, error = read_solution(data_file)

    implementation = data_file.stem.removeprefix("wave_solution_")
    output_file = data_file.with_name(f"result_wave_comparison_{implementation}.png")

    # --------------------------------------------------------------
    fig01, axs01 = plt.subplots(2, 1, sharex=True)

    axs01[0].plot(x, phi_num, linewidth=2.0, label='Numerical')
    axs01[0].plot(x, phi_ex, markevery=3, markersize=5, linewidth=2.0,
                  linestyle='dashed', marker='o', color='k', label='Analytical')
    axs01[0].legend(frameon=False, loc="best")

    axs01[1].semilogy(x, numpy.maximum(error, numpy.finfo(float).tiny),
                      linewidth=2.0, color='k', label='Absolute error')
    axs01[1].legend(frameon=False, loc="best")

    axs01[0].set_ylabel(r'$\phi$')
    axs01[1].set_ylabel(r'$|\phi - \phi_{exact}|$')

    axs01[1].set_xlabel(r'$x \ (m)$')

    # Enable minor ticks
    axs01[0].minorticks_on()
    axs01[1].minorticks_on()

    # Set fontsize (equivalent to 'fontsize(gcf, fontsize, "points")')
    fontsize = 12  # Change as needed
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    fig01.set_size_inches(10, 8)  # Adjust the width and height as needed

    plt.savefig(output_file)
    plt.close(fig01)
    print(f"plot written to {output_file.relative_to(root_directory)}")