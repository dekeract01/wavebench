#!/usr/bin/env python
"""
wave.py

solves the one-dimensional linear advection equation
    dphi/dt = -c0*dphi/dx
on a periodic domain using a fourth-order central spatial derivative
and third-order runge-kutta time integration.

this numba baseline reports timing and numerical error against the
analytical sine-wave solution. it is intended as a reference for
comparing parallel implementations and optimisation strategies.

output: wave_solution_python_numba.h5
author: dekeract01, 2025
university of southampton
"""

from time import perf_counter

import h5py
import numpy as np
from numba import njit, prange

c0 = 0.5
dt = 0.0000001
niter = 1000
nx = 2000000
dx = 1.0/nx
pi = np.pi


# @njit(parallel=True)
# def spatial_derivative_fourth(dphi_dx, phi, grid_spacing):
#     """Calculate a fourth-order periodic central spatial derivative."""
#     n = phi.size
#     inv_12dx = 1.0/(12.0*grid_spacing)

#     for index in prange(n):
#         im2 = (index - 2 + n) % n
#         im1 = (index - 1 + n) % n
#         ip1 = (index + 1) % n
#         ip2 = (index + 2) % n
#         dphi_dx[index] = (
#             -phi[ip2] + 8.0*phi[ip1] - 8.0*phi[im1] + phi[im2]
#         )*inv_12dx

@njit(parallel=True)
def spatial_derivative_fourth(dphi_dx, phi, grid_spacing):
    """Fourth-order periodic central derivative: vectorized interior + explicit wrap."""
    n = phi.size
    inv_12dx = 1.0/(12.0*grid_spacing)

    # interior: contiguous stencil, SIMD + parallel
    for i in prange(2, n - 2):
        dphi_dx[i] = (-phi[i+2] + 8.0*phi[i+1] - 8.0*phi[i-1] + phi[i-2])*inv_12dx

    # periodic wrap: four edge points
    dphi_dx[0]   = (-phi[2] + 8.0*phi[1]   - 8.0*phi[n-1] + phi[n-2])*inv_12dx
    dphi_dx[1]   = (-phi[3] + 8.0*phi[2]   - 8.0*phi[0]   + phi[n-1])*inv_12dx
    dphi_dx[n-2] = (-phi[0] + 8.0*phi[n-1] - 8.0*phi[n-3] + phi[n-4])*inv_12dx
    dphi_dx[n-1] = (-phi[1] + 8.0*phi[0]   - 8.0*phi[n-2] + phi[n-3])*inv_12dx


@njit(parallel=True)
def rhs(dphi, phi, dphi_dx, wave_speed, grid_spacing):
    """Calculate dphi/dt for the linear advection equation."""
    spatial_derivative_fourth(dphi_dx, phi, grid_spacing)
    for index in prange(phi.size):
        dphi[index] = -wave_speed*dphi_dx[index]


@njit(parallel=True)
def rk3_step(phi, phi1, phi2, k1, k2, k3, dphi_dx, wave_speed, grid_spacing, time_step):
    """Advance phi by one third-order runge-kutta step."""
    rhs(k1, phi, dphi_dx, wave_speed, grid_spacing)
    for index in prange(phi.size):
        phi1[index] = phi[index] + time_step*k1[index]

    rhs(k2, phi1, dphi_dx, wave_speed, grid_spacing)
    for index in prange(phi.size):
        phi2[index] = phi[index] + time_step*(0.25*k1[index] + 0.25*k2[index])

    rhs(k3, phi2, dphi_dx, wave_speed, grid_spacing)
    for index in prange(phi.size):
        phi[index] += time_step*(k1[index]/6.0 + k2[index]/6.0 + 2.0*k3[index]/3.0)


def main():
    print("1d wave equation simulation (python numba)")
    print("="*50)
    print(f"grid points: {nx}")
    print(f"dx: {dx:.8f}")
    print(f"dt: {dt:.8f}")
    print(f"wave_speed_c0: {c0:.8f}")
    print(f"cfl_number: {c0*dt/dx:.8f}")
    print(f"number_of_iterations: {niter}")
    print(f"total_time: {niter*dt:.8f}")
    print("="*50)

    x = np.arange(nx, dtype=np.float64)*dx
    phi = np.sin(2.0*pi*x)
    phi1 = np.empty_like(phi)
    phi2 = np.empty_like(phi)
    k1 = np.empty_like(phi)
    k2 = np.empty_like(phi)
    k3 = np.empty_like(phi)
    dphi_dx = np.empty_like(phi)

    print("\nwarming up numba compiler...")
    rk3_step(phi, phi1, phi2, k1, k2, k3, dphi_dx, c0, dx, dt)
    phi[:] = np.sin(2.0*pi*x)

    print("\nstarting simulation...")
    start_time = perf_counter()
    for iteration in range(1, niter + 1):
        rk3_step(phi, phi1, phi2, k1, k2, k3, dphi_dx, c0, dx, dt)

        if iteration % 100 == 0:
            if np.isnan(phi).any():
                raise RuntimeError(f"nan detected at iteration {iteration}")

            elapsed = perf_counter() - start_time
            print(
                f"iter {iteration:4d}/{niter} | time: {elapsed:.2f}s | "
                f"phi: [{phi.min():+.6f}, {phi.max():+.6f}]"
            )

    total_time = perf_counter() - start_time
    print("="*50)
    print("simulation_complete!")
    print(f"total_time: {total_time:.8f} seconds")
    print(f"time_per_iteration: {total_time/niter*1000.0:.8f} ms")
    print(f"iterations_per_second: {niter/total_time:.8f}")
    print("="*50)

    t_final = niter*dt
    phi_exact = np.sin(2.0*pi*(x - c0*t_final))
    error = np.abs(phi - phi_exact)
    print("\nnumerical vs analytical:")
    print(f"max_error: {error.max():.6e}")
    print(f"l2_error: {np.sqrt(np.mean(error**2)):.6e}")

    with h5py.File("wave_solution_python_numba.h5", "w") as output_file:
        output_file["x"] = x
        output_file["phi_numerical"] = phi
        output_file["phi_exact"] = phi_exact
        output_file["error"] = error

    print("\nresults written to wave_solution_python_numba.h5")


if __name__ == "__main__":
    main()