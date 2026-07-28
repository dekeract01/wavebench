# Wave Equation Benchmarks

I use this repository to solve the one-dimensional wave equation, more precisely the linear advection equation,

$$
\frac{\partial \phi}{\partial t} = -c_0 \frac{\partial \phi}{\partial x}.
$$

I have been updating this repository from time to time. I use the wave equation to learn how to code different things each time, including new languages, numerical libraries, parallel programming models, GPU programming, build systems, data output, and performance measurement.

The wave equation is a useful starting point because it is simple enough to understand completely, while still requiring the core parts of a realistic numerical solver: a spatial grid, boundary conditions, a stable time integrator, numerical accuracy checks, data output, and careful performance measurement. Its analytical sine-wave solution also lets me measure the numerical error directly instead of treating runtime as the only result.

## Implementations

I maintain versions in C++, CUDA C++, Fortran, OpenCL C++, Julia, OpenMP C++, MPI C++, MPI Fortran, Python, Python with Numba, and Rust. The serial versions provide reference implementations, while the OpenMP, MPI, OpenCL, CUDA, and Numba versions let me compare different ways of parallelising or accelerating the same problem.

All cases use the same numerical problem:

| Setting | Value |
| --- | --- |
| Domain | Periodic, $0 \leq x < 1$ |
| Initial condition | $\phi(x, 0) = \sin(2\pi x)$ |
| Wave speed | $c_0 = 0.5$ |
| Grid points | $n_x = 2{,}000{,}000$ |
| Time step | $\Delta t = 0.0000001$ |
| Iterations | 1,000 |
| Spatial discretisation | Fourth-order central finite difference |
| Time integration | Third-order Runge-Kutta |
| CFL number | $0.1$ |

At the final time, I compare the computed solution with

$$
\phi_{\mathrm{exact}}(x,t) = \sin\left(2\pi(x-c_0t)\right).
$$

## Running The Benchmarks

[`generate_all.py`](generate_all.py) is the root runner. With `clean_setting = 0`, it builds and runs the configured implementations in turn. With `clean_setting = 1`, it runs `make clean` in every directory that has a Makefile and then exits. This keeps cleaning separate from a new benchmark run. The CUDA implementation can be run independently with `make run` in [`cuda_cpp/`](cuda_cpp/).

```bash
python generate_all.py
python plot_wave.py
python plot_benchmarks.py
```

The solver runs write solution data in their own directories. The shared plotting scripts then create a solution comparison plot for each available result and two benchmark summary images in this directory.

## Benchmark Plots

### Fedora Linux Performance Results

![Fedora Linux performance comparison: time per iteration and iterations per second](https://dekeract01.github.io/images/result_benchmark_comparison_fedora.png)

The upper panel shows the time per iteration, where lower values are better. The lower panel shows iterations per second, where higher values are better. The bars are sorted by time per iteration. CUDA C++ is fastest in the current Fedora results, at approximately $0.767\,\mathrm{ms}$ per iteration and $1{,}300$ iterations per second on the RTX A4000. The CPU results reflect differences in compiler optimisation, memory access, runtime overhead, and parallel execution. These are results for this machine and configuration, not universal language rankings.

### Fedora Linux Accuracy Results

![Fedora Linux maximum error comparison](https://dekeract01.github.io/images/result_benchmark_error_fedora.png)

The error figure shows

$$
\max_x \left|\phi(x,t) - \phi_{\mathrm{exact}}(x,t)\right|
$$

on a logarithmic scale. The CUDA C++ implementation has a maximum error of $1.23 \times 10^{-14}$, consistent with the double-precision CPU implementations. The OpenCL implementation uses single precision, so its error is larger, around $10^{-5}$. The MPI Fortran result in this run has a larger error, around $6.17 \times 10^{-7}$, which demonstrates why accuracy needs to be assessed alongside timing.

## Output Format

Most implementations write HDF5 solution files containing `x`, `phi_numerical`, `phi_exact`, and `error`. The root [`plot_wave.py`](plot_wave.py) also supports the older text and binary output files, so I can plot results from every implementation through one script.

## Primary Test System

The Fedora benchmark results shown above were produced on the following machine:

| Component | Details |
| --- | --- |
| Machine | Dell Precision 3660 |
| CPU | 13th Gen Intel Core i9-13900K, 32 logical processors @ 5.50 GHz |
| GPU | NVIDIA RTX A4000 |
| Integrated GPU | Intel UHD Graphics 770 |
| Memory | 64 GB |
| OS | Fedora Linux 37 Workstation Edition (x86_64) |
| Kernel | 6.5.12-100.fc37.x86_64 |

The CUDA results run on the discrete RTX A4000. Their performance includes the implementation's host-device data handling, so results on other discrete GPUs, integrated GPUs, and CPU-only systems will differ.

## Earlier Reference System

Earlier results were collected on the following Apple Silicon system. They are useful as a separate reference point, but should not be compared as direct language or hardware rankings against the Fedora results above.

| Component | Details |
| --- | --- |
| Machine | MacBook Pro (14-inch, 2021) |
| CPU | Apple M1 Pro, 10 cores (8 performance + 2 efficiency) @ 3.23 GHz |
| GPU | Apple M1 Pro, 14-core integrated GPU |
| Memory | 16 GB unified |
| OS | macOS Sequoia 15.7.5 (arm64) |

On this system, OpenCL runs on the integrated M1 Pro GPU using the same unified memory as the CPU, so no host-device transfer over PCIe is involved. Results on discrete-GPU systems differ accordingly.