#!/usr/bin/env python
"""
wave.py

solves the one-dimensional linear advection equation
    dphi/dt = -c0*dphi/dx
on a periodic domain using a fourth-order central spatial derivative
and third-order runge-kutta time integration.

this serial baseline reports timing and numerical error against the
analytical sine-wave solution. it is intended as a reference for
comparing parallel implementations and optimisation strategies.

output: wave_solution_python.h5
author: dekeract01, 2025
university of southampton
"""

from time import perf_counter

import h5py
import numpy as np

# simulation parameters
c0 = 0.5
dt = 0.0000001
niter = 1000
nx = 2000000
dx = 1.0/nx

print("1d wave equation simulation (python)")
print(f"="*50)
print(f"grid points: {nx}")
print(f"dx: {dx:.8f}")
print(f"dt: {dt:.8f}")
print(f"wave_speed_c0: {c0:.8f}")
print(f"cfl_number: {c0*dt/dx:.8f}")
print(f"number_of_iterations: {niter}")
print(f"total_time: {niter*dt:.8f}")
print(f"="*50)

# Grid setup (periodic domain from 0 to 1)
x0 = np.linspace(0, 1.0 - dx, nx)  # Periodic, so exclude endpoint

# Initial condition: phi = sin(2*pi*x)
phi = np.sin(2.0 * np.pi * x0)

# For 4th order central differences and RK3, we need storage
phi_old = phi.copy()

# 4th order central difference for spatial derivative
def spatial_derivative_4th(f, dx):
    """4th order central difference with periodic boundaries"""
    df = np.zeros_like(f)
    n = len(f)
    
    # Interior points (4th order central)
    for i in range(n):
        df[i] = (-f[(i+2)%n] + 8*f[(i+1)%n] - 8*f[(i-1)%n] + f[(i-2)%n]) / (12.0*dx)
    
    return df

# RHS function: dphi/dt = -c * dphi/dx
def rhs(phi, c0, dx):
    return -c0 * spatial_derivative_4th(phi, dx)

# Runge-Kutta 3 time integration
def rk3_step(phi, c0, dx, dt):
    """Third-order Runge-Kutta (RK3) time step"""
    # Stage 1
    k1 = rhs(phi, c0, dx)
    phi1 = phi + dt * k1
    
    # Stage 2
    k2 = rhs(phi1, c0, dx)
    phi2 = phi + dt * (0.25*k1 + 0.25*k2)
    
    # Stage 3
    k3 = rhs(phi2, c0, dx)
    phi_new = phi + dt * (k1/6.0 + k2/6.0 + 2.0*k3/3.0)
    
    return phi_new

# time stepping loop
print("\nstarting simulation...")
start_time = perf_counter()

for n in range(niter):
    phi = rk3_step(phi, c0, dx, dt)
    
    # nan check and progress report
    if (n + 1) % 100 == 0:
        if np.isnan(phi).any():
            print(f"nan detected at iteration {n + 1}!")
            break
        elapsed = perf_counter() - start_time
        phi_min, phi_max = phi.min(), phi.max()
        print(f"iter {n + 1:4d}/{niter} | time: {elapsed:.2f}s | "
              f"phi: [{phi_min:+.6f}, {phi_max:+.6f}]")

end_time = perf_counter()
total_time = end_time - start_time

print(f"="*50)
print("simulation_complete!")
print(f"total_time: {total_time:.8f} seconds")
print(f"time_per_iteration: {total_time/niter*1000:.8f} ms")
print(f"iterations_per_second: {niter/total_time:.8f}")
print(f"="*50)

# Analytical solution for comparison: phi = sin(2*pi*(x - c*t))
t_final = niter * dt
phi_exact = np.sin(2.0 * np.pi * (x0 - c0*t_final))

# Calculate error
error = np.abs(phi - phi_exact)
print("\nnumerical vs analytical:")
print(f"max_error: {error.max():.6e}")
print(f"l2_error: {np.sqrt(np.mean(error**2)):.6e}")

with h5py.File("wave_solution_python.h5", "w") as output_file:
    output_file["x"] = x0
    output_file["phi_numerical"] = phi
    output_file["phi_exact"] = phi_exact
    output_file["error"] = error

print("\nresults written to wave_solution_python.h5")