/*
 * wave_openmp.cpp
 *
 * solves the one-dimensional linear advection equation
 *     dphi/dt = -c0*dphi/dx
 * on a periodic domain using a fourth-order central spatial derivative
 * and third-order runge-kutta time integration.
 *
 * this openmp baseline reports timing and numerical error against the
 * analytical sine-wave solution. it is intended as a reference for
 * comparing parallel implementations and optimisation strategies.
 *
 * output: wave_solution_openmp.h5
 * author: dekeract01, 2021
 * university of southampton
 */

#include <cstdio>
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <hdf5.h>
#include <omp.h>

using vectorfield = std::vector<double>;
using timer_clock = std::chrono::high_resolution_clock;
using seconds = std::chrono::duration<double>;

// simulation parameters
constexpr double c0 = 0.5;
constexpr double dt = 0.0000001;
constexpr int niter = 1000;
constexpr int nx = 2000000;
constexpr double dx = 1.0/nx;
constexpr double pi = 3.14159265358979323846;

// decimal place control
constexpr int decimals_general = 6;
constexpr int decimals_cfl = 3;
constexpr int decimals_time = 8;
constexpr int decimals_phi = 6;
constexpr int decimals_error = 6;

inline void check_hdf5_error(herr_t error, const char* message) {
    if (error < 0) {
        fprintf(stderr, "hdf5 error: %s\n", message);
        exit(EXIT_FAILURE);
    }
}

void write_dataset(hid_t file, const char* name, const vectorfield& values) {
    hsize_t dimensions[] = {values.size()};
    hid_t dataspace = H5Screate_simple(1, dimensions, nullptr);
    if (dataspace < 0) {
        fprintf(stderr, "hdf5 error: could not create dataspace\n");
        exit(EXIT_FAILURE);
    }

    hid_t dataset = H5Dcreate2(file, name, H5T_IEEE_F64LE, dataspace,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0) {
        fprintf(stderr, "hdf5 error: could not create dataset %s\n", name);
        H5Sclose(dataspace);
        exit(EXIT_FAILURE);
    }

    check_hdf5_error(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, values.data()), "could not write dataset");
    check_hdf5_error(H5Dclose(dataset), "could not close dataset");
    check_hdf5_error(H5Sclose(dataspace), "could not close dataspace");
}

// 4th order central difference for spatial derivative (periodic) - OpenMP parallelized
inline void spatial_derivative_fourth(vectorfield& df, const vectorfield& f, double dx) {
    int n = f.size();
    double inv_12dx = 1.0 / (12.0 * dx);
    
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        int im2 = (i - 2 + n) % n;
        int im1 = (i - 1 + n) % n;
        int ip1 = (i + 1) % n;
        int ip2 = (i + 2) % n;
        
        df[i] = (-f[ip2] + 8.0 * f[ip1] - 8.0 * f[im1] + f[im2]) * inv_12dx;
    }
}

// RHS function: dphi/dt = -c * dphi/dx - OpenMP parallelized
inline void rhs(vectorfield& dphi, const vectorfield& phi, vectorfield& dphi_dx,
                double c0, double dx) {
    spatial_derivative_fourth(dphi_dx, phi, dx);
    
    int n = phi.size();
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        dphi[i] = -c0 * dphi_dx[i];
    }
}

// Runge-Kutta 3 time integration - OpenMP parallelized
inline void rk3_step(vectorfield& phi, vectorfield& phi1, vectorfield& phi2,
                     vectorfield& k1, vectorfield& k2, vectorfield& k3,
                     vectorfield& dphi_dx, double c0, double dx, double dt) {
    int n = phi.size();
    
    // Stage 1
    rhs(k1, phi, dphi_dx, c0, dx);
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        phi1[i] = phi[i] + dt * k1[i];
    }
    
    // Stage 2
    rhs(k2, phi1, dphi_dx, c0, dx);
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        phi2[i] = phi[i] + dt * (0.25 * k1[i] + 0.25 * k2[i]);
    }
    
    // Stage 3
    rhs(k3, phi2, dphi_dx, c0, dx);
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        phi[i] = phi[i] + dt * (k1[i] / 6.0 + k2[i] / 6.0 + 2.0 * k3[i] / 3.0);
    }
}

int main() {
    // Get number of threads
    int num_threads = omp_get_max_threads();
    
    printf("1d wave equation simulation (c++ with openmp)\n");
    printf("==================================================\n");
    printf("grid points: %d\n", nx);
    printf("openmp_threads: %d\n", num_threads);
    printf("dx: %.*f\n", decimals_general, dx);
    printf("dt: %.*f\n", decimals_general, dt);
    printf("wave_speed_c0: %.*f\n", decimals_general, c0);
    printf("cfl_number: %.*f\n", decimals_cfl, c0*dt/dx);
    printf("number_of_iterations: %d\n", niter);
    printf("total_time: %.*f\n", decimals_time, niter*dt);
    printf("==================================================\n");
    
    // Grid setup (periodic domain from 0 to 1)
    vectorfield x0(nx);
    #pragma omp parallel for
    for (int i = 0; i < nx; i++) {
        x0[i] = i*dx;
    }
    
    // Initial condition: phi = sin(2*pi*x)
    vectorfield phi(nx);
    #pragma omp parallel for
    for (int i = 0; i < nx; i++) {
        phi[i] = std::sin(2.0*pi*x0[i]);
    }
    
    // Preallocate arrays for RK3 stages
    vectorfield phi1(nx);
    vectorfield phi2(nx);
    vectorfield k1(nx);
    vectorfield k2(nx);
    vectorfield k3(nx);
    vectorfield dphi_dx(nx);
    
    // Time stepping loop
    printf("\nStarting simulation...\n");
    auto start_time = timer_clock::now();
    
    for (int n = 0; n < niter; n++) {
        rk3_step(phi, phi1, phi2, k1, k2, k3, dphi_dx, c0, dx, dt);
        
        // NaN check and progress report
        if ((n + 1) % 100 == 0) {
            bool nan_detected = false;
            #pragma omp parallel for reduction(||:nan_detected)
            for (int i = 0; i < nx; i++) {
                if (std::isnan(phi[i])) {
                    nan_detected = true;
                }
            }
            
            if (nan_detected) {
                printf("nan detected at iteration %d!\n", n+1);
                return 1;
            }
            
            auto current_time = timer_clock::now();
            seconds elapsed = current_time-start_time;
            double phi_min = *std::min_element(phi.begin(), phi.end());
            double phi_max = *std::max_element(phi.begin(), phi.end());
            
                 printf("iter %4d/%d | time: %.2fs | phi: [%+.*f, %+.*f]\n",
                     n+1, niter, elapsed.count(),
                     decimals_phi, phi_min, decimals_phi, phi_max);
        }
    }
    
    auto end_time = timer_clock::now();
    seconds total_time = end_time-start_time;
    
    printf("==================================================\n");
    printf("simulation_complete!\n");
    printf("total_time: %.*f seconds\n", decimals_time, total_time.count());
    printf("time_per_iteration: %.*f ms\n", decimals_time, total_time.count()/niter*1000.0);
    printf("iterations_per_second: %.8f\n", niter/total_time.count());
    printf("==================================================\n");
    
    // Analytical solution: phi = sin(2*pi*(x - c*t))
    double t_final = niter*dt;
    vectorfield phi_exact(nx);
    vectorfield error(nx);
    
    #pragma omp parallel for
    for (int i = 0; i < nx; i++) {
        phi_exact[i] = std::sin(2.0*pi*(x0[i]-c0*t_final));
        error[i] = std::abs(phi[i] - phi_exact[i]);
    }
    
    // Calculate errors
    double max_error = *std::max_element(error.begin(), error.end());
    double l2_error = 0.0;
    
    #pragma omp parallel for reduction(+:l2_error)
    for (int i = 0; i < nx; i++) {
        l2_error += error[i] * error[i];
    }
    l2_error = std::sqrt(l2_error/nx);
    
    printf("\nnumerical vs analytical:\n");
    printf("max_error: %.*e\n", decimals_error, max_error);
    printf("l2_error: %.*e\n", decimals_error, l2_error);
    
    hid_t output_file = H5Fcreate("wave_solution_openmp.h5", H5F_ACC_TRUNC,
                                  H5P_DEFAULT, H5P_DEFAULT);
    if (output_file < 0) {
        fprintf(stderr, "hdf5 error: could not create output file\n");
        return 1;
    }

    write_dataset(output_file, "x", x0);
    write_dataset(output_file, "phi_numerical", phi);
    write_dataset(output_file, "phi_exact", phi_exact);
    write_dataset(output_file, "error", error);
    check_hdf5_error(H5Fclose(output_file), "could not close output file");
    printf("\nresults written to wave_solution_openmp.h5\n");
    
    return 0;
}