# Wave Equation Benchmarks

I use this repository to solve the one-dimensional wave equation, more precisely the linear advection equation,

$$
\frac{\partial \phi}{\partial t} = -c_0 \frac{\partial \phi}{\partial x}.
$$

I have been updating this repository from time to time. I use the wave equation to learn how to code different things each time, including new languages, numerical libraries, parallel programming models, GPU programming, build systems, data output, and performance measurement.

The wave equation is a useful starting point because it is simple enough to understand completely, while still requiring the core parts of a realistic numerical solver: a spatial grid, boundary conditions, a stable time integrator, numerical accuracy checks, data output, and careful performance measurement. Its analytical sine-wave solution also lets me measure the numerical error directly instead of treating runtime as the only result.

## Implementations

I maintain versions in C++, CUDA C++, Fortran, standard-parallelism Fortran (`do concurrent`), OpenCL C++, Julia, OpenMP C++, MPI C++, MPI Fortran, Python, and Python with Numba, and Rust. The serial versions provide reference implementations, while the OpenMP, MPI, OpenCL, CUDA, standard-parallelism, and Numba versions let me compare different ways of parallelising or accelerating the same problem.

The standard-parallelism Fortran version is the newest addition, and the one I understand least well so far. It is just ordinary Fortran with the loops written as `do concurrent`, and the same source compiles to serial, multi-core CPU, or GPU depending only on the compiler flag I pass to `nvfortran` (`-stdpar=multicore` or `-stdpar=gpu`). I did not write any CUDA or add any directives for it, so I mostly wanted to see how far the compiler could get on its own.

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

The upper panel shows the time per iteration, where lower values are better. The lower panel shows iterations per second, where higher values are better. The bars are sorted by time per iteration.

CUDA C++ is still the fastest here, at about $0.767\,\mathrm{ms}$ per iteration and roughly $1{,}300$ iterations per second on the RTX A4000. The result that surprised me is the standard-parallelism Fortran version, which came second at about $1.24\,\mathrm{ms}$ per iteration and around $808$ iterations per second on the same A4000. I did not expect plain Fortran with `do concurrent` to land that close to the hand-written CUDA C++, within roughly 1.6 times, for that little effort. Both of those run on the same GPU, so at least that part is a fair comparison.

I am not completely sure why the gap to CUDA C++ is there. When I built the standard-parallelism version, the compiler reported that it was copying the arrays between the CPU and the GPU around each step, so my guess is that a lot of the difference is data movement rather than the actual computation. I have not properly measured that yet, so I am treating it as a guess for now, and it is the next thing I want to look into.

The remaining CPU results reflect differences in compiler optimisation, memory access, runtime overhead, and parallel execution. These are results for this machine and configuration, not universal language rankings.

### Fedora Linux Accuracy Results

![Fedora Linux maximum error comparison](https://dekeract01.github.io/images/result_benchmark_error_fedora.png)

The error figure shows

$$
\max_x \left|\phi(x,t) - \phi_{\mathrm{exact}}(x,t)\right|
$$

on a logarithmic scale. In this run everything sits around $10^{-14}$, at the double-precision floor, including the standard-parallelism Fortran version at $1.23 \times 10^{-14}$, which is reassuring given how little I changed to get it onto the GPU.

Earlier I had a bug in the MPI Fortran version that is worth mentioning, because it is exactly why I keep the accuracy plot next to the timing plot. I was only exchanging one ghost point between the sub-domains, but the fourth-order stencil reaches two points either side, so the points near each internal boundary were using a stale value. That pushed the MPI Fortran error up to about $6.17 \times 10^{-7}$ while every other result stayed near $10^{-14}$. The timing looked perfectly reasonable, so without the accuracy plot I might not have noticed. I have since widened the halo to two points, and MPI Fortran is now back down with everything else at about $1.17 \times 10^{-14}$.

### macOS Performance Results

![macOS performance comparison: time per iteration and iterations per second](https://dekeract01.github.io/images/result_benchmark_comparison.png)

The OpenCL C++ implementation is fastest in these results, at about $1.26\,\mathrm{ms}$ per iteration, because its kernels run on the M1 Pro integrated GPU. This one is expected rather than surprising, and it is not really a fair comparison with the others for two reasons. First, the OpenCL version uses single precision, so it moves and computes half as much data as the double-precision implementations. Second, on the M1 the GPU shares unified memory with the CPU, so there is no separate host-device transfer. The accuracy plot below shows the trade-off for the single precision. The remaining CPU results reflect differences in compiler optimisation, memory access, runtime overhead, and parallel execution on this Apple Silicon system.

### macOS Accuracy Results

![macOS maximum error comparison](https://dekeract01.github.io/images/result_benchmark_error.png)

Most CPU implementations cluster around $10^{-14}$ maximum error. The OpenCL implementation uses single precision, so its error is much larger, around $10^{-5}$. That larger error is expected and is the direct trade-off for the speed it gets on the integrated GPU. These results belong to the M1 Pro configuration described below and should be interpreted separately from the Fedora measurements.

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

The CUDA and standard-parallelism GPU results run on the discrete RTX A4000. Their performance includes the implementation's host-device data handling, so results on other discrete GPUs, integrated GPUs, and CPU-only systems will differ.

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