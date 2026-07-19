/*
 * wave_opencl.cpp
 *
 * solves the one-dimensional linear advection equation
 *     dphi/dt = -c0*dphi/dx
 * on an opencl device using a fourth-order central spatial derivative
 * and third-order runge-kutta time integration.
 *
 * this single-precision gpu baseline reports timing and numerical error
 * against the analytical sine-wave solution. it is intended as a reference
 * for comparing parallel implementations and optimisation strategies.
 *
 * output: wave_solution_opencl.h5
 * author: dekeract01, 2023
 * university of southampton
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <string>
#include <hdf5.h>

#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

using vectorfield = std::vector<float>;
using timer_clock = std::chrono::high_resolution_clock;
using seconds = std::chrono::duration<double>;

// simulation parameters - float is used for apple gpu compatibility.
constexpr float c0 = 0.5f;
constexpr float dt = 0.0000001f;
constexpr int niter = 1000;
constexpr int nx = 2000000;
constexpr float dx = 1.0f/nx;
constexpr float pi = 3.14159265358979323846f;

// decimal place control
constexpr int decimals_general = 6;
constexpr int decimals_cfl = 3;
constexpr int decimals_time = 8;
constexpr int decimals_phi = 6;
constexpr int decimals_error = 6;

inline void check_cl_error(cl_int error, const char* message){
    if (error != CL_SUCCESS){
        fprintf(stderr, "opencl error (%d): %s\n", error, message);
        exit(EXIT_FAILURE);
        }
    }

std::string read_kernel_source(const char* filename){
    std::ifstream input(filename, std::ios::binary);
    if (!input){
        fprintf(stderr, "error: could not open opencl kernel file %s\n", filename);
        exit(EXIT_FAILURE);
        }

        return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
    }

inline void check_hdf5_error(herr_t error, const char* message){
    if (error < 0){
        fprintf(stderr, "hdf5 error: %s\n", message);
        exit(EXIT_FAILURE);
        }
    }

void write_dataset(hid_t file, const char* name, const vectorfield& values){
    hsize_t dimensions[] = {values.size()};
    hid_t dataspace = H5Screate_simple(1, dimensions, nullptr);
    if (dataspace < 0){
        fprintf(stderr, "hdf5 error: could not create dataspace\n");
        exit(EXIT_FAILURE);
        }

    hid_t dataset = H5Dcreate2(file, name, H5T_IEEE_F32LE, dataspace,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0){
        fprintf(stderr, "hdf5 error: could not create dataset %s\n", name);
        H5Sclose(dataspace);
        exit(EXIT_FAILURE);
        }

    check_hdf5_error(H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, values.data()), "could not write dataset");
    check_hdf5_error(H5Dclose(dataset), "could not close dataset");
    check_hdf5_error(H5Sclose(dataspace), "could not close dataspace");
    }

int main(){
    cl_int err;
    
    printf("1d wave equation simulation (opencl - single precision)\n");
    printf("-------------------------------------------------------\n");
    
    // get platform
    cl_platform_id platform;
    err = clGetPlatformIDs(1, &platform, NULL);
    check_cl_error(err, "failed to get platform");
    
    // Get device (try GPU first, fallback to CPU)
    cl_device_id device;
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
    if (err != CL_SUCCESS){
        err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, NULL);
        check_cl_error(err, "failed to get device");
    }
    
    // Get device name
    char device_name[128];
    clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_name), device_name, NULL);
    printf("using device: %s\n", device_name);
    
    // Create context
    cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    check_cl_error(err, "failed to create context");
    
    // Create command queue
    cl_command_queue queue = clCreateCommandQueue(context, device, 0, &err);
    check_cl_error(err, "failed to create command queue");
    
    // create and build program
    std::string kernel_source = read_kernel_source("wave_opencl.cl");
    const char* kernel_source_pointer = kernel_source.c_str();
    size_t kernel_source_size = kernel_source.size();
    cl_program program = clCreateProgramWithSource(
        context, 1, &kernel_source_pointer, &kernel_source_size, &err
    );
    check_cl_error(err, "failed to create program");
    
    err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    if (err != CL_SUCCESS){
        size_t log_size;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
        char* log = (char*)malloc(log_size);
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
        fprintf(stderr, "build error:\n%s\n", log);
        free(log);
        exit(EXIT_FAILURE);
    }
    
    // Create kernels
    cl_kernel kernel_deriv = clCreateKernel(program, "spatial_derivative_4th", &err);
    check_cl_error(err, "failed to create derivative kernel");
    
    cl_kernel kernel_rhs = clCreateKernel(program, "compute_rhs", &err);
    check_cl_error(err, "failed to create rhs kernel");
    
    cl_kernel kernel_stage1 = clCreateKernel(program, "rk3_stage1", &err);
    check_cl_error(err, "failed to create stage1 kernel");
    
    cl_kernel kernel_stage2 = clCreateKernel(program, "rk3_stage2", &err);
    check_cl_error(err, "failed to create stage2 kernel");
    
    cl_kernel kernel_stage3 = clCreateKernel(program, "rk3_stage3", &err);
    check_cl_error(err, "failed to create stage3 kernel");
    
    printf("grid points: %d\n", nx);
    printf("dx: %.*f\n", decimals_general, dx);
    printf("dt: %.*f\n", decimals_general, dt);
    printf("wave_speed_c0: %.*f\n", decimals_general, c0);
    printf("cfl_number: %.*f\n", decimals_cfl, c0*dt/dx);
    printf("number_of_iterations: %d\n", niter);
    printf("total_time: %.*f\n", decimals_time, niter*dt);
    printf("-------------------------------------------------------\n");
    
    // Initialize host arrays
    vectorfield x0(nx);
    vectorfield phi_host(nx);
    
    for (int i = 0; i < nx; i++){
        x0[i] = i*dx;
        phi_host[i] = std::sin(2.0f*pi*x0[i]);
    }
    
    // Create device buffers
    cl_mem d_phi = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                  nx*sizeof(float), phi_host.data(), &err);
    check_cl_error(err, "failed to create phi buffer");
    
    cl_mem d_phi1 = clCreateBuffer(context, CL_MEM_READ_WRITE, nx*sizeof(float), NULL, &err);
    check_cl_error(err, "failed to create phi1 buffer");
    
    cl_mem d_phi2 = clCreateBuffer(context, CL_MEM_READ_WRITE, nx*sizeof(float), NULL, &err);
    check_cl_error(err, "failed to create phi2 buffer");
    
    cl_mem d_k1 = clCreateBuffer(context, CL_MEM_READ_WRITE, nx*sizeof(float), NULL, &err);
    check_cl_error(err, "failed to create k1 buffer");
    
    cl_mem d_k2 = clCreateBuffer(context, CL_MEM_READ_WRITE, nx*sizeof(float), NULL, &err);
    check_cl_error(err, "failed to create k2 buffer");
    
    cl_mem d_k3 = clCreateBuffer(context, CL_MEM_READ_WRITE, nx*sizeof(float), NULL, &err);
    check_cl_error(err, "failed to create k3 buffer");
    
    cl_mem d_dphi_dx = clCreateBuffer(context, CL_MEM_READ_WRITE, nx*sizeof(float), NULL, &err);
    check_cl_error(err, "failed to create dphi_dx buffer");
    
    // Set constant kernel arguments
    float inv_12dx = 1.0f/(12.0f*dx);
    
    // Work group size
    size_t global_size = nx;
    size_t local_size = 64;
    if (global_size % local_size != 0){
        global_size = ((nx+local_size-1)/local_size)*local_size;
    }
    
    printf("\nstarting simulation...\n");
    auto start_time = timer_clock::now();
    
    // time stepping loop
    for (int n = 0; n < niter; n++){
        // Stage 1
        clSetKernelArg(kernel_deriv, 0, sizeof(cl_mem), &d_phi);
        clSetKernelArg(kernel_deriv, 1, sizeof(cl_mem), &d_dphi_dx);
        clSetKernelArg(kernel_deriv, 2, sizeof(float), &inv_12dx);
        clSetKernelArg(kernel_deriv, 3, sizeof(int), &nx);
        clEnqueueNDRangeKernel(queue, kernel_deriv, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
        
        clSetKernelArg(kernel_rhs, 0, sizeof(cl_mem), &d_dphi_dx);
        clSetKernelArg(kernel_rhs, 1, sizeof(cl_mem), &d_k1);
        clSetKernelArg(kernel_rhs, 2, sizeof(float), &c0);
        clSetKernelArg(kernel_rhs, 3, sizeof(int), &nx);
        clEnqueueNDRangeKernel(queue, kernel_rhs, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
        
        clSetKernelArg(kernel_stage1, 0, sizeof(cl_mem), &d_phi);
        clSetKernelArg(kernel_stage1, 1, sizeof(cl_mem), &d_k1);
        clSetKernelArg(kernel_stage1, 2, sizeof(cl_mem), &d_phi1);
        clSetKernelArg(kernel_stage1, 3, sizeof(float), &dt);
        clSetKernelArg(kernel_stage1, 4, sizeof(int), &nx);
        clEnqueueNDRangeKernel(queue, kernel_stage1, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
        
        // Stage 2
        clSetKernelArg(kernel_deriv, 0, sizeof(cl_mem), &d_phi1);
        clEnqueueNDRangeKernel(queue, kernel_deriv, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
        
        clSetKernelArg(kernel_rhs, 1, sizeof(cl_mem), &d_k2);
        clEnqueueNDRangeKernel(queue, kernel_rhs, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
        
        clSetKernelArg(kernel_stage2, 0, sizeof(cl_mem), &d_phi);
        clSetKernelArg(kernel_stage2, 1, sizeof(cl_mem), &d_k1);
        clSetKernelArg(kernel_stage2, 2, sizeof(cl_mem), &d_k2);
        clSetKernelArg(kernel_stage2, 3, sizeof(cl_mem), &d_phi2);
        clSetKernelArg(kernel_stage2, 4, sizeof(float), &dt);
        clSetKernelArg(kernel_stage2, 5, sizeof(int), &nx);
        clEnqueueNDRangeKernel(queue, kernel_stage2, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
        
        // Stage 3
        clSetKernelArg(kernel_deriv, 0, sizeof(cl_mem), &d_phi2);
        clEnqueueNDRangeKernel(queue, kernel_deriv, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
        
        clSetKernelArg(kernel_rhs, 1, sizeof(cl_mem), &d_k3);
        clEnqueueNDRangeKernel(queue, kernel_rhs, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
        
        clSetKernelArg(kernel_stage3, 0, sizeof(cl_mem), &d_phi);
        clSetKernelArg(kernel_stage3, 1, sizeof(cl_mem), &d_k1);
        clSetKernelArg(kernel_stage3, 2, sizeof(cl_mem), &d_k2);
        clSetKernelArg(kernel_stage3, 3, sizeof(cl_mem), &d_k3);
        clSetKernelArg(kernel_stage3, 4, sizeof(float), &dt);
        clSetKernelArg(kernel_stage3, 5, sizeof(int), &nx);
        clEnqueueNDRangeKernel(queue, kernel_stage3, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
        
        // Progress report
        if ((n + 1) % 100 == 0){
            clEnqueueReadBuffer(queue, d_phi, CL_TRUE, 0, nx*sizeof(float),
                              phi_host.data(), 0, NULL, NULL);
            
            bool nan_detected = false;
            for (float val : phi_host){
                if (std::isnan(val)){
                    nan_detected = true;
                    break;
                }
            }
            
            if (nan_detected){
                printf("nan detected at iteration %d!\n", n+1);
                return 1;
            }
            
            auto current_time = timer_clock::now();
            seconds elapsed = current_time-start_time;
            float phi_min = *std::min_element(phi_host.begin(), phi_host.end());
            float phi_max = *std::max_element(phi_host.begin(), phi_host.end());
            
                 printf("iter %4d/%d | time: %.2fs | phi: [%+.*f, %+.*f]\n",
                     n+1, niter, elapsed.count(),
                     decimals_phi, phi_min, decimals_phi, phi_max);
        }
    }
    
    clEnqueueReadBuffer(queue, d_phi, CL_TRUE, 0, nx*sizeof(float),
                       phi_host.data(), 0, NULL, NULL);
    
    auto end_time = timer_clock::now();
    seconds total_time = end_time-start_time;
    
    printf("-------------------------------------------------------\n");
    printf("simulation_complete!\n");
    printf("total_time: %.*f seconds\n", decimals_time, total_time.count());
    printf("time_per_iteration: %.*f ms\n", decimals_time, total_time.count()/niter*1000.0);
    printf("iterations_per_second: %.1f\n", niter/total_time.count());
    printf("-------------------------------------------------------\n");
    
    // Analytical solution
    float t_final = niter*dt;
    vectorfield phi_exact(nx);
    vectorfield error(nx);
    
    for (int i = 0; i < nx; i++){
        phi_exact[i] = std::sin(2.0f*pi*(x0[i]-c0*t_final));
        error[i] = std::abs(phi_host[i] - phi_exact[i]);
    }
    
    float max_error = *std::max_element(error.begin(), error.end());
    float l2_error = 0.0f;
    for (float e : error){
        l2_error += e * e;
    }
    l2_error = std::sqrt(l2_error/nx);
    
    printf("\nnumerical vs analytical:\n");
    printf("max_error: %.*e\n", decimals_error, max_error);
    printf("l2_error: %.*e\n", decimals_error, l2_error);
    
    hid_t output_file = H5Fcreate("wave_solution_opencl.h5", H5F_ACC_TRUNC,
                                  H5P_DEFAULT, H5P_DEFAULT);
    if (output_file < 0){
        fprintf(stderr, "hdf5 error: could not create output file\n");
        return 1;
    }

    write_dataset(output_file, "x", x0);
    write_dataset(output_file, "phi_numerical", phi_host);
    write_dataset(output_file, "phi_exact", phi_exact);
    write_dataset(output_file, "error", error);
    check_hdf5_error(H5Fclose(output_file), "could not close output file");
    printf("\nresults written to wave_solution_opencl.h5\n");
    
    // Cleanup
    clReleaseMemObject(d_phi);
    clReleaseMemObject(d_phi1);
    clReleaseMemObject(d_phi2);
    clReleaseMemObject(d_k1);
    clReleaseMemObject(d_k2);
    clReleaseMemObject(d_k3);
    clReleaseMemObject(d_dphi_dx);
    clReleaseKernel(kernel_deriv);
    clReleaseKernel(kernel_rhs);
    clReleaseKernel(kernel_stage1);
    clReleaseKernel(kernel_stage2);
    clReleaseKernel(kernel_stage3);
    clReleaseProgram(program);
    clReleaseCommandQueue(queue);
    clReleaseContext(context);
    
    return 0;
}