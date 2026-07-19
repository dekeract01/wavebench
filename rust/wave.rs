/*
 * wave.rs
 *
 * solves the one-dimensional linear advection equation
 *     dphi/dt = -c0*dphi/dx
 * on a periodic domain using a fourth-order central spatial derivative
 * and third-order runge-kutta time integration.
 *
 * this serial baseline reports timing and numerical error against the
 * analytical sine-wave solution. it is intended as a reference for
 * comparing parallel implementations and optimisation strategies.
 *
 * output: wave_solution_rust.h5
 * author: dekeract01, 2025
 * university of southampton
 */

use hdf5::File;
use std::f64::consts::PI;
use std::time::Instant;

const C0: f64 = 0.5;           // Wave speed
const DT: f64 = 0.0000001;     // Time step
const NITER: usize = 1000;     // Number of iterations
const NX: usize = 2_000_000;   // Grid points
const DX: f64 = 1.0 / NX as f64; // Grid spacing

fn main() {
    println!("1d wave equation simulation (rust)");
    println!("{}", "=".repeat(50));
    println!("grid points: {}", NX);
    println!("dx: {:.6}", DX);
    println!("dt: {:.6}", DT);
    println!("wave_speed_c0: {:.8}", C0);
    println!("cfl_number: {:.8}", C0 * DT / DX);
    println!("number_of_iterations: {}", NITER);
    println!("total_time: {:.8}", NITER as f64 * DT);
    println!("{}", "=".repeat(50));

    // Grid setup (periodic domain from 0 to 1)
    let x0: Vec<f64> = (0..NX).map(|i| i as f64 * DX).collect();

    // Initial condition: phi = sin(2*pi*x)
    let mut phi: Vec<f64> = x0.iter().map(|&x| (2.0 * PI * x).sin()).collect();

    // Preallocate arrays for RK3 stages
    let mut phi1 = vec![0.0; NX];
    let mut phi2 = vec![0.0; NX];
    let mut k1 = vec![0.0; NX];
    let mut k2 = vec![0.0; NX];
    let mut k3 = vec![0.0; NX];
    let mut dphi_dx = vec![0.0; NX];

    // Time stepping loop
    println!("\nstarting simulation...");
    let start_time = Instant::now();

    for n in 0..NITER {
        rk3_step(
            &mut phi,
            &mut phi1,
            &mut phi2,
            &mut k1,
            &mut k2,
            &mut k3,
            &mut dphi_dx,
        );

        // NaN check and progress report
        if (n + 1) % 100 == 0 {
            if phi.iter().any(|&x| x.is_nan()) {
                println!("nan detected at iteration {}!", n + 1);
                std::process::exit(1);
            }

            let elapsed = start_time.elapsed().as_secs_f64();
            let phi_min = phi.iter().cloned().fold(f64::INFINITY, f64::min);
            let phi_max = phi.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            println!(
                "iter {:4}/{} | time: {:.2}s | phi: [{:+.6}, {:+.6}]",
                n + 1,
                NITER,
                elapsed,
                phi_min,
                phi_max
            );
        }
    }

    let total_time = start_time.elapsed().as_secs_f64();

    println!("{}", "=".repeat(50));
    println!("simulation_complete!");
    println!("total_time: {:.8} seconds", total_time);
    println!("time_per_iteration: {:.8} ms", total_time / NITER as f64 * 1000.0);
    println!("iterations_per_second: {:.8}", NITER as f64 / total_time);
    println!("{}", "=".repeat(50));

    // Analytical solution: phi = sin(2*pi*(x - c*t))
    let t_final = NITER as f64 * DT;
    let phi_exact: Vec<f64> = x0.iter()
        .map(|&x| (2.0 * PI * (x - C0 * t_final)).sin())
        .collect();

    // Calculate errors
    let error: Vec<f64> = phi.iter()
        .zip(phi_exact.iter())
        .map(|(&p, &pe)| (p - pe).abs())
        .collect();

    let max_error = error.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let l2_error = (error.iter().map(|&e| e * e).sum::<f64>() / NX as f64).sqrt();

    println!("\nnumerical vs analytical:");
    println!("max_error: {:.6e}", max_error);
    println!("l2_error: {:.6e}", l2_error);

    write_hdf5_solution(&x0, &phi, &phi_exact, &error)
        .expect("could not write wave_solution_rust.h5");
    println!("\nresults written to wave_solution_rust.h5");
}

fn write_hdf5_solution(
    x: &[f64],
    phi_numerical: &[f64],
    phi_exact: &[f64],
    error: &[f64],
) -> hdf5::Result<()> {
    let output_file = File::create("wave_solution_rust.h5")?;
    output_file.new_dataset_builder().with_data(x).create("x")?;
    output_file
        .new_dataset_builder()
        .with_data(phi_numerical)
        .create("phi_numerical")?;
    output_file
        .new_dataset_builder()
        .with_data(phi_exact)
        .create("phi_exact")?;
    output_file
        .new_dataset_builder()
        .with_data(error)
        .create("error")?;
    Ok(())
}

// 4th order central difference for spatial derivative (periodic)
#[inline]
fn spatial_derivative_4th(df: &mut [f64], f: &[f64]) {
    let n = f.len();
    let inv_12dx = 1.0 / (12.0 * DX);

    // interior: windows give the optimizer provably in-bounds, contiguous access
    for (df_i, w) in df[2..n-2].iter_mut().zip(f.windows(5)) {
        *df_i = (-w[4] + 8.0 * w[3] - 8.0 * w[1] + w[0]) * inv_12dx;
    }

    // periodic wrap: the four edge points
    df[0]   = (-f[2] + 8.0 * f[1]     - 8.0 * f[n-1] + f[n-2]) * inv_12dx;
    df[1]   = (-f[3] + 8.0 * f[2]     - 8.0 * f[0]   + f[n-1]) * inv_12dx;
    df[n-2] = (-f[0] + 8.0 * f[n-1]   - 8.0 * f[n-3] + f[n-4]) * inv_12dx;
    df[n-1] = (-f[1] + 8.0 * f[0]     - 8.0 * f[n-2] + f[n-3]) * inv_12dx;
}

// RHS function: dphi/dt = -c * dphi/dx
#[inline]
fn rhs(dphi: &mut [f64], phi: &[f64], dphi_dx: &mut [f64]) {
    spatial_derivative_4th(dphi_dx, phi);
    for i in 0..phi.len() {
        dphi[i] = -C0 * dphi_dx[i];
    }
}

// Runge-Kutta 3 time integration
#[inline]
fn rk3_step(
    phi: &mut [f64],
    phi1: &mut [f64],
    phi2: &mut [f64],
    k1: &mut [f64],
    k2: &mut [f64],
    k3: &mut [f64],
    dphi_dx: &mut [f64],
) {
    let n = phi.len();

    // Stage 1
    rhs(k1, phi, dphi_dx);
    for i in 0..n {
        phi1[i] = phi[i] + DT * k1[i];
    }

    // Stage 2
    rhs(k2, phi1, dphi_dx);
    for i in 0..n {
        phi2[i] = phi[i] + DT * (0.25 * k1[i] + 0.25 * k2[i]);
    }

    // Stage 3
    rhs(k3, phi2, dphi_dx);
    for i in 0..n {
        phi[i] = phi[i] + DT * (k1[i] / 6.0 + k2[i] / 6.0 + 2.0 * k3[i] / 3.0);
    }
}