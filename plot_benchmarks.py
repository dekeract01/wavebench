#!/usr/bin/env python
"""
    plot_benchmarks.py

    finds wave benchmark reports below this directory and plots the measured
    time per iteration and iterations per second, and reports the max error
    of each implementation against the analytical solution.

    input: */wave_*.out.txt
    output: result_benchmark_comparison.png, result_benchmark_error.png
    author: dekeract01, 2025
    modified by: dekeract01, 2026
    university of southampton
"""

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy

project_directory = Path(__file__).resolve().parent
output_file = project_directory / "result_benchmark_comparison_macos.png"
error_output_file = project_directory / "result_benchmark_error_macos.png"

metric_patterns = {
    "time_per_iteration":    re.compile(r"^\s*time[ _]per[ _]iteration:\s*([0-9.eE+-]+)", re.IGNORECASE),
    "iterations_per_second": re.compile(r"^\s*iterations[ _]per[ _]second:\s*([0-9.eE+-]+)", re.IGNORECASE),
    "max_error":             re.compile(r"^\s*max[ _]error:\s*([0-9.eE+-]+)", re.IGNORECASE),
}
grid_pattern = re.compile(r"^\s*(?:total )?grid points:\s*(\d+)", re.IGNORECASE)

directory_names = {
    "cpp":          "c++",
    "fortran":      "fortran",
    "gpu_cpp":      "opencl",
    "julia":        "julia",
    "mp_cpp":       "openmp c++",
    "mpi_cpp":      "mpi c++",
    "mpi_fortran":  "mpi fortran",
    "python":       "python",
    "python_numba": "python numba",
    "rust":         "rust",
}


def implementation_name(log_file):
    """Return a readable name from the directory and benchmark filename."""
    name = directory_names.get(log_file.parent.name, log_file.parent.name)
    thread_match = re.fullmatch(r"wave_openmp_(\d+)\.out\.txt", log_file.name)
    if thread_match:
        name += f" ({thread_match.group(1)} threads)"
    return name


def read_benchmark(log_file):
    """Extract the final measured metrics and grid size from a report file."""
    metrics = {name: [] for name in metric_patterns}
    grid_points = None

    for line in log_file.read_text(encoding="utf-8", errors="replace").splitlines():
        grid_match = grid_pattern.match(line)
        if grid_match:
            grid_points = int(grid_match.group(1))

        for metric_name, pattern in metric_patterns.items():
            metric_match = pattern.match(line)
            if metric_match:
                metrics[metric_name].append(float(metric_match.group(1)))

    if any(not values for values in metrics.values()):
        return None

    return {
        "name":                  implementation_name(log_file),
        "grid_points":           grid_points,
        "time_per_iteration":    metrics["time_per_iteration"][-1],
        "iterations_per_second": metrics["iterations_per_second"][-1],
        "max_error":             metrics["max_error"][-1],
    }


log_files = sorted(project_directory.rglob("wave_*.out.txt"))
benchmarks = [benchmark for log_file in log_files if (benchmark := read_benchmark(log_file))]

if not benchmarks:
    raise RuntimeError("no complete wave_*.out.txt benchmark reports were found")

benchmarks.sort(key=lambda benchmark: benchmark["time_per_iteration"])

labels = []
for benchmark in benchmarks:
    label = benchmark["name"]

    labels.append(label)

positions = numpy.arange(len(benchmarks))
time_per_iteration    = [benchmark["time_per_iteration"]    for benchmark in benchmarks]
iterations_per_second = [benchmark["iterations_per_second"] for benchmark in benchmarks]

# --------------------------------------------------------------
fig01, axs01 = plt.subplots(2, 1, sharex=True)

bars0 = axs01[0].bar(positions, time_per_iteration, width=0.65, color='b')
bars1 = axs01[1].bar(positions, iterations_per_second, width=0.65, color='r')

# compact value labels above each bar
for axis, bars in ((axs01[0], bars0), (axs01[1], bars1)):
    for bar in bars:
        value = bar.get_height()
        axis.annotate(f"{value:.3g}",
                      (bar.get_x() + bar.get_width()/2, value),
                      ha="center", va="bottom",
                      xytext=(0, 3), textcoords="offset points", fontsize=8)

axs01[0].set_ylabel('time per iteration $(ms)$')
axs01[1].set_ylabel('iterations per second $(1/s)$')

axs01[1].set_xticks(positions, labels, rotation=45, ha='right', rotation_mode='anchor')

# Enable minor ticks
axs01[0].minorticks_on()
axs01[1].minorticks_on()

# Set fontsize (equivalent to 'fontsize(gcf, fontsize, "points")')
fontsize = 12  # Change as needed
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

fig01.set_size_inches(11, 8)  # Adjust the width and height as needed

axs01[0].set_ylim(0, 13)        # time per iteration axis
axs01[1].set_ylim(0, 1500)     # iterations per second axis


fig01.tight_layout()  # keep the angled labels from being clipped
plt.savefig(output_file)

# --------------------------------------------------------------
fig02, axs02 = plt.subplots(1, 1)

max_error = [benchmark["max_error"] for benchmark in benchmarks]

bars2 = axs02.bar(positions, max_error, width=0.65, color='k')
axs02.set_yscale('log')

axs02.set_ylim(1e-15, 8e-6)     # error axis (log scale — give positive values)


# compact value labels above each bar
for bar in bars2:
    value = bar.get_height()
    axs02.annotate(f"{value:.3g}",
                   (bar.get_x() + bar.get_width()/2, value),
                   ha="center", va="bottom",
                   xytext=(0, 3), textcoords="offset points", fontsize=8)

axs02.set_ylabel(r'$max \ |\phi - \phi_{exact}|$')

axs02.set_xticks(positions, labels, rotation=45, ha='right', rotation_mode='anchor')


fig02.set_size_inches(11, 5)  # Adjust the width and height as needed

fig02.tight_layout()  # keep the angled labels from being clipped
plt.savefig(error_output_file)


print(f"benchmark plot written to {output_file}")
print(f"error plot written to {error_output_file}")
print("included implementations:")
for benchmark in benchmarks:
    print(f"  {benchmark['name']}: "
          f"{benchmark['time_per_iteration']:.8f} ms/iteration, "
          f"{benchmark['iterations_per_second']:.8f} iterations/s, "
          f"max_error: {benchmark['max_error']:.6e}")