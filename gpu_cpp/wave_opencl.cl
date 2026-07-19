__kernel void spatial_derivative_4th(
    __global const float* f,
    __global float* df,
    const float inv_12dx,
    const int n){
    int i = get_global_id(0);

    if (i < n){
        int im2 = (i-2+n) % n;
        int im1 = (i-1+n) % n;
        int ip1 = (i+1) % n;
        int ip2 = (i+2) % n;

        df[i] = (-f[ip2]+8.0f*f[ip1]-8.0f*f[im1]+f[im2])*inv_12dx;
        }
    }

__kernel void compute_rhs(
    __global const float* dphi_dx,
    __global float* dphi,
    const float c0,
    const int n){
    int i = get_global_id(0);

    if (i < n) {
        dphi[i] = -c0*dphi_dx[i];
        }
    }

__kernel void rk3_stage1(
    __global const float* phi,
    __global const float* k1,
    __global float* phi1,
    const float dt,
    const int n){
    int i = get_global_id(0);

    if (i < n) {
        phi1[i] = phi[i]+dt*k1[i];
        }
    }

__kernel void rk3_stage2(
    __global const float* phi,
    __global const float* k1,
    __global const float* k2,
    __global float* phi2,
    const float dt,
    const int n){
    int i = get_global_id(0);

    if (i < n) {
        phi2[i] = phi[i]+dt*(0.25f*k1[i]+0.25f*k2[i]);
        }
    }

__kernel void rk3_stage3(
    __global float* phi,
    __global const float* k1,
    __global const float* k2,
    __global const float* k3,
    const float dt,
    const int n){
    int i = get_global_id(0);

    if (i < n) {
        phi[i] = phi[i]+dt*(k1[i]/6.0f+k2[i]/6.0f+2.0f*k3[i]/3.0f);
        }
    }