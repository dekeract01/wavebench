/*
 * wave_mpi.cpp
 *
 * solves the one-dimensional linear advection equation
 *     dphi/dt = -c0*dphi/dx
 * on a periodic domain using mpi domain decomposition, a fourth-order
 * central spatial derivative, and third-order runge-kutta time integration.
 *
 * this version uses 2-wide ghost cells (matching the 4th-order stencil
 * width) so the derivative is computed with a single branchless loop over
 * the interior. rank 0 gathers and writes the comparison result after the
 * parallel solve completes.
 *
 * output: wave_solution_mpi.h5
 * author: dekeract01, 2023 (ghost-cell rework 2026)
 * university of southampton
 */

#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <mpi.h>
#include <hdf5.h>

using vectorfield = std::vector<double>;

// Simulation parameters
constexpr double C0 = 0.5;           // Wave speed
constexpr double DT = 0.0000001;     // Time step
constexpr int NITER = 1000;          // Number of iterations
constexpr int NX_GLOBAL = 2000000;   // Total grid points
constexpr double DX = 1.0 / NX_GLOBAL;
constexpr double PI = 3.14159265358979323846;

// Ghost cells: 4th-order central stencil needs 2 on each side
constexpr int NG = 2;

// DECIMAL PLACE CONTROL
constexpr int DECIMALS_GENERAL = 6;
constexpr int DECIMALS_CFL = 3;
constexpr int DECIMALS_TIME = 3;
constexpr int DECIMALS_PHI = 6;
constexpr int DECIMALS_ERROR = 6;

inline void check_hdf5_error(herr_t error, const char* message) {
    if (error < 0) {
        fprintf(stderr, "hdf5 error: %s\n", message);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

void write_dataset(hid_t file, const char* name, const vectorfield& values) {
    hsize_t dimensions[] = {values.size()};
    hid_t dataspace = H5Screate_simple(1, dimensions, nullptr);
    if (dataspace < 0) {
        fprintf(stderr, "hdf5 error: could not create dataspace\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    hid_t dataset = H5Dcreate2(file, name, H5T_IEEE_F64LE, dataspace,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0) {
        fprintf(stderr, "hdf5 error: could not create dataset %s\n", name);
        H5Sclose(dataspace);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    check_hdf5_error(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, values.data()), "could not write dataset");
    check_hdf5_error(H5Dclose(dataset), "could not close dataset");
    check_hdf5_error(H5Sclose(dataspace), "could not close dataspace");
}

/*
 * Exchange ghost cells with periodic neighbours.
 *
 * Array layout (size n + 2*NG):
 *   [ 0 .. NG-1 ]              left ghosts   (copies of left neighbour's last NG interior points)
 *   [ NG .. NG+n-1 ]           interior
 *   [ NG+n .. NG+n+NG-1 ]      right ghosts  (copies of right neighbour's first NG interior points)
 */
inline void exchange_boundaries(vectorfield& f, int left_rank, int right_rank, int n) {
    MPI_Status status;

    // send my last NG interior points to the right neighbour,
    // receive my left ghosts from the left neighbour
    MPI_Sendrecv(&f[NG + n - NG], NG, MPI_DOUBLE, right_rank, 0,
                 &f[0],           NG, MPI_DOUBLE, left_rank,  0,
                 MPI_COMM_WORLD, &status);

    // send my first NG interior points to the left neighbour,
    // receive my right ghosts from the right neighbour
    MPI_Sendrecv(&f[NG],     NG, MPI_DOUBLE, left_rank,  1,
                 &f[NG + n], NG, MPI_DOUBLE, right_rank, 1,
                 MPI_COMM_WORLD, &status);
}

// 4th order central difference for spatial derivative.
// Assumes ghost cells are up to date; single branchless loop over interior.
inline void spatial_derivative_4th(vectorfield& df, const vectorfield& f, int n) {
    const double inv_12dx = 1.0 / (12.0 * DX);

    for (int i = NG; i < NG + n; i++) {
        df[i] = (-f[i+2] + 8.0*f[i+1] - 8.0*f[i-1] + f[i-2]) * inv_12dx;
    }
}

// RHS function: dphi/dt = -c0 * dphi/dx (interior only)
inline void rhs(vectorfield& dphi, const vectorfield& f,
                vectorfield& dphi_dx, int n) {
    spatial_derivative_4th(dphi_dx, f, n);
    for (int i = NG; i < NG + n; i++) {
        dphi[i] = -C0 * dphi_dx[i];
    }
}

// RK3 time step with MPI ghost exchange before each stage's derivative
inline void rk3_step(vectorfield& phi, vectorfield& phi1, vectorfield& phi2,
                     vectorfield& k1, vectorfield& k2, vectorfield& k3,
                     vectorfield& dphi_dx, int left_rank, int right_rank, int n) {
    // Stage 1
    exchange_boundaries(phi, left_rank, right_rank, n);
    rhs(k1, phi, dphi_dx, n);
    for (int i = NG; i < NG + n; i++) {
        phi1[i] = phi[i] + DT * k1[i];
    }

    // Stage 2
    exchange_boundaries(phi1, left_rank, right_rank, n);
    rhs(k2, phi1, dphi_dx, n);
    for (int i = NG; i < NG + n; i++) {
        phi2[i] = phi[i] + DT * (0.25 * k1[i] + 0.25 * k2[i]);
    }

    // Stage 3
    exchange_boundaries(phi2, left_rank, right_rank, n);
    rhs(k3, phi2, dphi_dx, n);
    for (int i = NG; i < NG + n; i++) {
        phi[i] = phi[i] + DT * (k1[i] / 6.0 + k2[i] / 6.0 + 2.0 * k3[i] / 3.0);
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Domain decomposition (last rank takes the remainder)
    int nx_per_proc = NX_GLOBAL / nprocs;
    int remainder = NX_GLOBAL % nprocs;

    int nx_local = (rank < nprocs - 1) ? nx_per_proc : nx_per_proc + remainder;
    int i_start = rank * nx_per_proc;

    // Neighbour ranks for periodic boundaries
    int left_rank = (rank - 1 + nprocs) % nprocs;
    int right_rank = (rank + 1) % nprocs;

    if (rank == 0) {
        printf("1D Wave Equation Simulation (C++ with MPI, 2-wide ghost cells)\n");
        printf("==================================================\n");
        printf("Total grid points: %d\n", NX_GLOBAL);
        printf("Number of MPI processes: %d\n", nprocs);
        printf("Points per process: %d\n", nx_per_proc);
        printf("dx: %.*f\n", DECIMALS_GENERAL, DX);
        printf("dt: %.*f\n", DECIMALS_GENERAL, DT);
        printf("Wave speed c0: %.*f\n", DECIMALS_GENERAL, C0);
        printf("CFL number: %.*f\n", DECIMALS_CFL, C0 * DT / DX);
        printf("Number of iterations: %d\n", NITER);
        printf("Total time: %.*f\n", DECIMALS_TIME, NITER * DT);
        printf("==================================================\n");
    }

    // All work arrays share the ghost-padded layout: [NG | interior | NG]
    const int ntot = nx_local + 2 * NG;
    vectorfield phi(ntot, 0.0);
    vectorfield phi1(ntot, 0.0);
    vectorfield phi2(ntot, 0.0);
    vectorfield k1(ntot, 0.0);
    vectorfield k2(ntot, 0.0);
    vectorfield k3(ntot, 0.0);
    vectorfield dphi_dx(ntot, 0.0);

    // Local coordinates (interior only)
    vectorfield x0(nx_local);
    for (int i = 0; i < nx_local; i++) {
        x0[i] = (i_start + i) * DX;
    }

    // Initial condition: phi = sin(2*pi*x) on the interior
    for (int i = 0; i < nx_local; i++) {
        phi[NG + i] = std::sin(2.0 * PI * x0[i]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        printf("\nStarting simulation...\n");
    }

    double start_time = MPI_Wtime();

    for (int n = 0; n < NITER; n++) {
        rk3_step(phi, phi1, phi2, k1, k2, k3, dphi_dx, left_rank, right_rank, nx_local);

        // NaN check and progress report
        if ((n + 1) % 100 == 0) {
            bool nan_detected = false;
            for (int i = NG; i < NG + nx_local; i++) {
                if (std::isnan(phi[i])) {
                    nan_detected = true;
                    break;
                }
            }

            if (nan_detected) {
                printf("NaN detected at iteration %d on rank %d!\n", n + 1, rank);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            if (rank == 0) {
                double elapsed = MPI_Wtime() - start_time;
                double phi_min = *std::min_element(phi.begin() + NG, phi.begin() + NG + nx_local);
                double phi_max = *std::max_element(phi.begin() + NG, phi.begin() + NG + nx_local);
                printf("Iter %4d/%d | Time: %.2fs | phi: [%+.*f, %+.*f]\n",
                       n + 1, NITER, elapsed, DECIMALS_PHI, phi_min, DECIMALS_PHI, phi_max);
            }
        }
    }

    double end_time = MPI_Wtime();
    double total_time = end_time - start_time;

    // Gather interior portions to rank 0
    vectorfield phi_global;
    std::vector<int> recvcounts(nprocs);
    std::vector<int> displs(nprocs);

    if (rank == 0) {
        phi_global.resize(NX_GLOBAL);
        for (int i = 0; i < nprocs; i++) {
            recvcounts[i] = (i < nprocs - 1) ? nx_per_proc : nx_per_proc + remainder;
            displs[i] = (i == 0) ? 0 : displs[i - 1] + recvcounts[i - 1];
        }
    }

    MPI_Gatherv(phi.data() + NG, nx_local, MPI_DOUBLE,
                phi_global.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("==================================================\n");
        printf("Simulation complete!\n");
        printf("Total time: %.*f seconds\n", DECIMALS_TIME, total_time);
        printf("Time per iteration: %.*f ms\n", DECIMALS_TIME, total_time / NITER * 1000.0);
        printf("Iterations per second: %.1f\n", NITER / total_time);
        printf("==================================================\n");

        // Analytical solution: phi = sin(2*pi*(x - c*t))
        double t_final = NITER * DT;
        vectorfield x_global(NX_GLOBAL);
        vectorfield phi_exact(NX_GLOBAL);
        vectorfield error(NX_GLOBAL);

        for (int i = 0; i < NX_GLOBAL; i++) {
            double x = i * DX;
            x_global[i] = x;
            phi_exact[i] = std::sin(2.0 * PI * (x - C0 * t_final));
            error[i] = std::abs(phi_global[i] - phi_exact[i]);
        }

        double max_error = *std::max_element(error.begin(), error.end());
        double l2_error = 0.0;
        for (double e : error) {
            l2_error += e * e;
        }
        l2_error = std::sqrt(l2_error / NX_GLOBAL);

        printf("\nNumerical vs Analytical:\n");
        printf("Max error: %.*e\n", DECIMALS_ERROR, max_error);
        printf("L2 error: %.*e\n", DECIMALS_ERROR, l2_error);

        hid_t output_file = H5Fcreate("wave_solution_mpi.h5", H5F_ACC_TRUNC,
                                      H5P_DEFAULT, H5P_DEFAULT);
        if (output_file < 0) {
            fprintf(stderr, "hdf5 error: could not create output file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        write_dataset(output_file, "x", x_global);
        write_dataset(output_file, "phi_numerical", phi_global);
        write_dataset(output_file, "phi_exact", phi_exact);
        write_dataset(output_file, "error", error);
        check_hdf5_error(H5Fclose(output_file), "could not close output file");
        printf("\nresults written to wave_solution_mpi.h5\n");
    }

    MPI_Finalize();
    return 0;
}