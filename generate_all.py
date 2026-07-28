#!/usr/bin/env python3
"""
    generate_all.py

    cleans, builds, and runs every wave-equation implementation.

    set clean_setting to 1 to run "make clean" in every project that has a
    Makefile and exit without compiling or running anything.
    set clean_setting to 0 to build and run all implementations.

    make commands are run with conda removed from the environment
    (equivalent to "conda deactivate"), because the conda venv3 environment
    interferes with the system compilers and libraries.

    author: dekeract01, 2026
    university of southampton
"""

import os
import subprocess
import sys
from pathlib import Path

clean_setting = 0
root_directory = Path(__file__).resolve().parent

# name, directory, run command
implementations = [
    ("c++",          root_directory / "cpp",          ("make", "run")),
    ("fortran",      root_directory / "fortran",      ("make", "run")),
    # ("opencl c++",   root_directory / "gpu_cpp",      ("make", "run")),
    ("julia",        root_directory / "julia",        ("julia", "wave.jl")),
    ("openmp c++",   root_directory / "mp_cpp",       ("make", "run")),
    ("mpi c++",      root_directory / "mpi_cpp",      ("make", "run")),
    ("cuda c++",      root_directory / "cuda_cpp",      ("make", "run")),
    ("mpi fortran",  root_directory / "mpi_fortran",  ("make", "run")),
    # ("python",       root_directory / "python",       (sys.executable, "wave.py")),
    ("numba", root_directory / "python_numba", (sys.executable, "wave.py")),
    ("rust",         root_directory / "rust",         ("make", "run")),
]

failures = []

# --------------------------------------------------------------
# build an environment with conda deactivated for the make commands:
# remove the conda environment's directories from PATH and drop all
# CONDA_* variables, so make picks up the system compilers instead
make_environment = os.environ.copy()
conda_prefix = make_environment.get("CONDA_PREFIX", "")
if conda_prefix:
    path_entries = make_environment.get("PATH", "").split(os.pathsep)
    path_entries = [entry for entry in path_entries if not entry.startswith(conda_prefix)]
    make_environment["PATH"] = os.pathsep.join(path_entries)
    for key in list(make_environment):
        if key.startswith("CONDA_"):
            del make_environment[key]

# --------------------------------------------------------------
# clean only, then stop (nothing is built or run)
if clean_setting == 1:
    for name, directory, command in implementations:
        if (directory / "Makefile").exists():
            print(f"\n{'-' * 72}\ncleaning {name}\n{'=' * 72}")
            try:
                subprocess.run(("make", "clean"), cwd=directory, check=True, env=make_environment)
            except FileNotFoundError:
                print("command not found: make", file=sys.stderr)
                failures.append(f"clean: {name}")
            except subprocess.CalledProcessError as error:
                print(f"command failed with exit code {error.returncode}", file=sys.stderr)
                failures.append(f"clean: {name}")

    if failures:
        print(f"\ncleaning finished with failures: {', '.join(failures)}", file=sys.stderr)
        sys.exit(1)

    print("\nall implementations cleaned successfully")
    sys.exit(0)

# --------------------------------------------------------------
# build and run
for name, directory, command in implementations:
    print(f"\n{'=' * 72}\nbuilding and running {name}\n{'=' * 72}")

    # make commands run without conda, everything else keeps the normal environment
    if command[0] == "make":
        environment = make_environment
    else:
        environment = os.environ

    try:
        subprocess.run(command, cwd=directory, check=True, env=environment)
    except FileNotFoundError:
        print(f"command not found: {command[0]}", file=sys.stderr)
        failures.append(name)
    except subprocess.CalledProcessError as error:
        print(f"command failed with exit code {error.returncode}", file=sys.stderr)
        failures.append(name)

# --------------------------------------------------------------
if failures:
    print(f"\ncompleted with failures: {', '.join(failures)}", file=sys.stderr)
    sys.exit(1)

print("\nall implementations completed successfully")