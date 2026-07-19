"""
wave.jl

    solves the one-dimensional linear advection equation
        dphi/dt = -c0*dphi/dx
    on a periodic domain using a fourth-order central spatial derivative
    and third-order runge-kutta time integration.

    this serial baseline reports timing and numerical error against the
    analytical wave solution. it is intended as a reference for
    comparing parallel implementations and optimisation strategies.

    output: wave_solution_julia.h5
    author: dekeract01, 2025
    university of southampton
"""

using Printf
using HDF5

# Simulation parameters (matching your OpenSBLI setup)
const c0 = 0.5                    # Wave speed
const dt = 0.0000001              # Time step
const niter = 1000                # Number of iterations
const nx = 2000000                # Grid points
const dx = 1.0/nx                 # Grid spacing

println("1D Wave Equation Simulation (Julia)")
println("-"^50)
println("Grid points: $nx")
@printf("dx: %.6f\n", dx)
@printf("dt: %.6f\n", dt)
println("Wave speed c0: $c0")
@printf("CFL number: %.3f\n", c0*dt/dx)
println("Number of iterations: $niter")
@printf("Total time: %.8f\n", niter*dt)
println("-"^50)

# 4th order central difference for spatial derivative (periodic)
@inline function spatial_derivative_4th!(df, f, dx, n)
    inv_12dx = 1.0 / (12.0 * dx)
    # interior: contiguous stencil, SIMD-friendly
    @inbounds @simd for i in 3:n-2
        df[i] = (-f[i+2] + 8.0*f[i+1] - 8.0*f[i-1] + f[i-2]) * inv_12dx
    end
    # periodic wrap: four edge points
    @inbounds begin
        df[1]   = (-f[3] + 8.0*f[2]   - 8.0*f[n]   + f[n-1]) * inv_12dx
        df[2]   = (-f[4] + 8.0*f[3]   - 8.0*f[1]   + f[n]  ) * inv_12dx
        df[n-1] = (-f[1] + 8.0*f[n]   - 8.0*f[n-2] + f[n-3]) * inv_12dx
        df[n]   = (-f[2] + 8.0*f[1]   - 8.0*f[n-1] + f[n-2]) * inv_12dx
    end
    nothing
end

# RHS function: dphi/dt = -c * dphi/dx
@inline function rhs!(dphi, phi, dphi_dx, c0, dx, n)
    spatial_derivative_4th!(dphi_dx, phi, dx, n)
    @inbounds for i in 1:n
        dphi[i] = -c0*dphi_dx[i]
    end
    nothing
end

# runge-kutta 3 time integration (optimized, in-place)
function rk3_step!(phi, phi1, phi2, k1, k2, k3, dphi_dx, c0, dx, dt, n)
    # Stage 1
    rhs!(k1, phi, dphi_dx, c0, dx, n)
    @inbounds for i in 1:n
        phi1[i] = phi[i] + dt*k1[i]
    end
    
    # Stage 2
    rhs!(k2, phi1, dphi_dx, c0, dx, n)
    @inbounds for i in 1:n
        phi2[i] = phi[i] + dt*(0.25*k1[i] + 0.25*k2[i])
    end
    
    # Stage 3
    rhs!(k3, phi2, dphi_dx, c0, dx, n)
    @inbounds for i in 1:n
        phi[i] = phi[i] + dt*(k1[i]/6.0 + k2[i]/6.0 + 2.0*k3[i]/3.0)
    end
    nothing
end

# Grid setup (periodic domain from 0 to 1)
x0 = collect(range(0.0, step=dx, length=nx))

# Initial condition: phi = sin(2*pi*x)
phi = sin.(2.0 * π .* x0)

# Preallocate arrays for RK3 stages
phi1 = similar(phi)
phi2 = similar(phi)
k1 = similar(phi)
k2 = similar(phi)
k3 = similar(phi)
dphi_dx = similar(phi)

# Warm-up run (for JIT compilation)
println("\nWarming up JIT compiler...")
for n in 1:10
    rk3_step!(phi, phi1, phi2, k1, k2, k3, dphi_dx, c0, dx, dt, nx)
end

# Reset to initial condition
phi .= sin.(2.0 * π .* x0)

# Time stepping loop
println("\nStarting simulation...")
start_time = time()

for n in 1:niter
    rk3_step!(phi, phi1, phi2, k1, k2, k3, dphi_dx, c0, dx, dt, nx)
    
    # NaN check (matching your print_iteration_ops)
    if n % 100 == 0
        if any(isnan, phi)
            println("NaN detected at iteration $n!")
            break
        end
        elapsed = time() - start_time
        phi_min, phi_max = minimum(phi), maximum(phi)
        @printf("Iter %4d/%d | Time: %.2fs | phi: [%+.6f, %+.6f]\n", 
                n, niter, elapsed, phi_min, phi_max)
    end
end

end_time = time()
total_time = end_time - start_time

println("="^50)
println("simulation_complete!")
@printf("total_time: %.8f seconds\n", total_time)
@printf("time_per_iteration: %.8f ms\n", total_time/niter*1000)
@printf("iterations_per_second: %.8f\n", niter/total_time)
println("="^50)

# Analytical solution for comparison: phi = sin(2*pi*(x - c*t))
t_final = niter * dt
phi_exact = sin.(2.0 * π .* (x0 .- c0*t_final))

# Calculate error
error = abs.(phi .- phi_exact)
println("\nNumerical vs Analytical:")
@printf("max_error: %.6e\n", maximum(error))
@printf("l2_error: %.6e\n", sqrt(sum(error.^2) / nx))

# write datasets with the same names used by the other implementations.
h5open("wave_solution_julia.h5", "w") do output_file
    output_file["x"] = x0
    output_file["phi_numerical"] = phi
    output_file["phi_exact"] = phi_exact
    output_file["error"] = error
end
println("\nresults written to wave_solution_julia.h5")