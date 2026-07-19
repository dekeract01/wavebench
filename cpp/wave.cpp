/*
 * wave.cpp
 *
 * solves the one-dimensional linear advection equation
 *     dphi/dt = -c0 * dphi/dx
 * on a periodic domain using a fourth-order central spatial derivative
 * and third-order runge-kutta time integration.
 *
 * this serial baseline reports timing and numerical error against the
 * analytical sine-wave solution. it is intended as a reference for
 * comparing parallel implementations and optimisation strategies.
 *
 * output: wave_solution_cpp.dat
 * author: dekeract01, 2021
 * university of southampton
 */

#include <cstdio>
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>

using vectorfield = std::vector<double>;
using timer_clock = std::chrono::high_resolution_clock;
using seconds = std::chrono::duration<double>;

// simulation parameters
const double c0 = 0.5;            // wave speed
const double dt = 0.0000001;      // time step
const int niter = 1000;           // number of iterations
const int nx = 2000000;           // grid points
const double dx = 1.0 / nx;       // grid spacing
const double pi = 3.14159265358979323846;

// Decimal place control - change these to adjust output precision.
const int decimals_general = 8;   // For dx, dt, c0, etc.
const int decimals_cfl = 8;       // For CFL number
const int decimals_time = 8;      // For total time, time per iter
const int decimals_phi = 6;       // For phi min/max values
const int decimals_error = 6;     // For error output (scientific notation)
const int decimals_file = 8;      // For file output

// 4th order central difference for spatial derivative (periodic)
inline void spatial_derivative_fourth(vectorfield& df, const vectorfield& f, double dx) {
    const int n = f.size();
    const double inv_12dx = 1.0 / (12.0 * dx);

    // interior: pure contiguous stencil, vectorizes cleanly
    for (int i = 2; i < n-2; i++) {
        df[i] = (-f[i+2] + 8.0*f[i+1] - 8.0*f[i-1] + f[i-2])*inv_12dx;
    }

    // periodic wrap: just the 4 edge points
    df[0]   = (-f[2]   + 8.0*f[1]   - 8.0*f[n-1] + f[n-2])*inv_12dx;
    df[1]   = (-f[3]   + 8.0*f[2]   - 8.0*f[0]   + f[n-1])*inv_12dx;
    df[n-2] = (-f[0]   + 8.0*f[n-1] - 8.0*f[n-3] + f[n-4])*inv_12dx;
    df[n-1] = (-f[1]   + 8.0*f[0]   - 8.0*f[n-2] + f[n-3])*inv_12dx;
}

// rhs function: dphi/dt = -c*dphi/dx
inline void rhs(vectorfield& dphi, const vectorfield& phi, vectorfield& dphi_dx,
                double c0, double dx) {
    spatial_derivative_fourth(dphi_dx, phi, dx);

    for (size_t i = 0; i < phi.size(); i++) {
        dphi[i] = -c0*dphi_dx[i];
    }
}

// Runge-Kutta 3 time integration
inline void rk3_step(vectorfield& phi, vectorfield& phi1, vectorfield& phi2,
                     vectorfield& k1, vectorfield& k2, vectorfield& k3,
                     vectorfield& dphi_dx, double c0, double dx, double dt) {
    int n = phi.size();
    
    // stage 1
    rhs(k1, phi, dphi_dx, c0, dx);
    for (int i = 0; i < n; i++) {
        phi1[i] = phi[i] + dt*k1[i];
    }
    
    // stage 2
    rhs(k2, phi1, dphi_dx, c0, dx);
    for (int i = 0; i < n; i++) {
        phi2[i] = phi[i] + dt*(0.25*k1[i] + 0.25*k2[i]);
    }
    
    // stage 3
    rhs(k3, phi2, dphi_dx, c0, dx);
    for (int i = 0; i < n; i++) {
        phi[i] = phi[i] + dt*(k1[i]/6.0 + k2[i]/6.0 + 2.0*k3[i]/3.0);
    }
}

int main() {
    printf("1D wave equation simulation (C++)\n");
    printf("---------------------------------------------------\n");
    printf("grid points: %d\n", nx);
    printf("dx: %.*f\n", decimals_general, dx);
    printf("dt: %.*f\n", decimals_general, dt);
    printf("wave_speed_c0: %.*f\n", decimals_general, c0);
    printf("cfl_number: %.*f\n", decimals_cfl, c0 * dt / dx);
    printf("number_of_iterations: %d\n", niter);
    printf("total_time: %.*f\n", decimals_time, niter * dt);
    printf("---------------------------------------------------\n");
    
    // grid setup (periodic domain from 0 to 1)
    vectorfield x0(nx);
    for (int i = 0; i < nx; i++) {
        x0[i] = i*dx;
    }
    
    // initial condition: phi = sin(2*pi*x)
    vectorfield phi(nx);
    for (int i = 0; i < nx; i++) {
        phi[i] = std::sin(2.0*pi*x0[i]);
    }
    
    // preallocate arrays for RK3 stages
    vectorfield phi1(nx);
    vectorfield phi2(nx);
    vectorfield k1(nx);
    vectorfield k2(nx);
    vectorfield k3(nx);
    vectorfield dphi_dx(nx);
    
    // time stepping loop
    printf("\nstarting simulation...\n");
    auto start_time = timer_clock::now();
    
    
    for (int n = 0; n < niter; n++) {
        rk3_step(phi, phi1, phi2, k1, k2, k3, dphi_dx, c0, dx, dt);
        
        // NaN check and progress report
        if ((n+1) % 100 == 0) {
            bool nan_detected = false;
            for (double val : phi) {
                if (std::isnan(val)) {
                    nan_detected = true;
                    break;
                }
            }
            
            if (nan_detected) {
                printf("NaN detected at iteration %d!\n", n + 1);
                return 1;
            }
            
            auto current_time = timer_clock::now();
            seconds elapsed = current_time - start_time;
            double phi_min = *std::min_element(phi.begin(), phi.end());
            double phi_max = *std::max_element(phi.begin(), phi.end());
            
            printf("iter %4d/%d | time: %.2fs | phi: [%+.*f, %+.*f]\n",
                   n + 1, niter, elapsed.count(), 
                   decimals_phi, phi_min, decimals_phi, phi_max);
        }
    }
    
    auto end_time = timer_clock::now();
    seconds total_time = end_time - start_time;
    
    printf("==================================================\n");
    printf("simulation_complete!\n");
    printf("total_time: %.*f seconds\n", decimals_time, total_time.count());
    printf("time_per_iteration: %.*f ms\n", decimals_time, total_time.count() / niter * 1000.0);
    printf("iterations_per_second: %.1f\n", niter / total_time.count());
    printf("==================================================\n");
    
    // Analytical solution: phi = sin(2*pi*(x - c*t))
    double t_final = niter * dt;
    vectorfield phi_exact(nx);
    vectorfield error(nx);
    
    for (int i = 0; i < nx; i++) {
        phi_exact[i] = std::sin(2.0*pi*(x0[i] - c0*t_final));
        error[i] = std::abs(phi[i] - phi_exact[i]);
    }
    
    // Calculate errors
    double max_error = *std::max_element(error.begin(), error.end());
    double l2_error = 0.0;
    for (double e : error) {
        l2_error += e * e;
    }
    l2_error = std::sqrt(l2_error/nx);
    
    printf("\nnumerical vs Analytical:\n");
    printf("max_error: %.*e\n", decimals_error, max_error);
    printf("l2_error: %.*e\n", decimals_error, l2_error);
    
    // Write results to file
    FILE* outfile = fopen("wave_solution_cpp.dat", "w");
    if (outfile == nullptr) {
        printf("Error: Could not open output file!\n");
        return 1;
    }
    
    fprintf(outfile, "# x, phi_numerical, phi_exact, error\n");
    for (int i = 0; i < nx; i++) {
        fprintf(outfile, "%.*e %.*e %.*e %.*e\n", 
                decimals_file, x0[i], 
                decimals_file, phi[i], 
                decimals_file, phi_exact[i], 
                decimals_file, error[i]);
    }
    fclose(outfile);
    
    printf("\nresults written to wave_solution_cpp.dat\n");
    
    return 0;
}