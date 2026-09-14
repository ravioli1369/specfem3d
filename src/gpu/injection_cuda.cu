/*
 !=====================================================================
 !
 !                         S p e c f e m 3 D
 !                         -----------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                              CNRS, France
 !                       and Princeton University, USA
 !                 (there are currently many more authors!)
 !                           (c) October 2017
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 3 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================
*/

// Type-5 (Rayleigh/Love runtime-eval) surface-wave injection, evaluated on the device.
//
// Stage 1 (couple_with_injection.f90) precomputes, once, the time-invariant coefficients
// ray_Vs/ray_Vq/ray_Ts/ray_Tq and the per-point time offset ray_t0p, and lists only the
// ray_nact boundary points above the table's deepest sample in ray_act -- everything else
// stays at its zero-initialised value for the whole run. This file uploads those arrays once
// (prepare_rayleigh_injection_device) and evaluates the wavelet directly into d_veloc_inj /
// d_tract_inj every step (compute_rayleigh_injection_cuda), removing the host round trip
// stage 1 still needed: no host buffer, no per-step allocation, no H2D copy of the field.
//
// ray_t0p is uploaded and summed in DOUBLE: tau = tnow + t0p cancels against terms far
// larger than itself (see couple_with_injection.f90), so forming it in single would cost
// about 1.5 decimal digits of wavelet phase, same as on the host.

#include "mesh_constants_gpu.h"


/* ----------------------------------------------------------------------------------------------- */

__global__ void rayleigh_injection_kernel(const int* __restrict__ d_ray_act,
                                          const double* __restrict__ d_ray_t0p,
                                          const realw* __restrict__ d_ray_Vs,
                                          const realw* __restrict__ d_ray_Vq,
                                          const realw* __restrict__ d_ray_Ts,
                                          const realw* __restrict__ d_ray_Tq,
                                          realw* d_veloc_inj,
                                          realw* d_tract_inj,
                                          const int nact,
                                          const double tnow_d,
                                          const double om_d,
                                          const realw om,
                                          const realw ig2) {

  int bx = blockIdx.y*gridDim.x + blockIdx.x;
  int k = threadIdx.x + bx*blockDim.x;

  if (k >= nact) return;

  // tau formed in double, arg rounded to single once -- mirrors compute_rayleigh_field()
  double tau_d = tnow_d + d_ray_t0p[k];
  realw arg = (realw)(om_d * tau_d);

  realw env = exp(-(realw)0.5*(arg*arg)*ig2);
  realw ca = cos(arg);
  realw sa = sin(arg);

  realw s  = env*ca;
  realw q  = env*sa;
  realw w1 = -om*arg*ig2;
  realw sp = env*(w1*ca - om*sa);
  realw qp = env*(w1*sa + om*ca);

  // d_ray_act stores Fortran 1-based boundary-point indices (ray_act is also used
  // by the CPU fallback compute_rayleigh_field, which needs it 1-based); the flat
  // (NDIM,npt) layout the Stacey kernel reads is 0-based, hence the -1 here.
  int ipt = d_ray_act[k] - 1;
  int i0 = 3*ipt;
  int k0 = 3*k;

  d_veloc_inj[i0  ] = d_ray_Vs[k0  ]*sp + d_ray_Vq[k0  ]*qp;
  d_veloc_inj[i0+1] = d_ray_Vs[k0+1]*sp + d_ray_Vq[k0+1]*qp;
  d_veloc_inj[i0+2] = d_ray_Vs[k0+2]*sp + d_ray_Vq[k0+2]*qp;

  d_tract_inj[i0  ] = d_ray_Ts[k0  ]*s + d_ray_Tq[k0  ]*q;
  d_tract_inj[i0+1] = d_ray_Ts[k0+1]*s + d_ray_Tq[k0+1]*q;
  d_tract_inj[i0+2] = d_ray_Ts[k0+2]*s + d_ray_Tq[k0+2]*q;
}

/* ----------------------------------------------------------------------------------------------- */

// uploads the precomputed coefficients once and zeroes d_veloc_inj/d_tract_inj: points below
// the table (not listed in ray_act) must stay at exactly zero for the whole run, matching the
// CPU path where Veloc_specfem/Tract_specfem are zero-initialised and never touched for them.

extern EXTERN_LANG
void FC_FUNC_(prepare_rayleigh_injection_device,
              PREPARE_RAYLEIGH_INJECTION_DEVICE)(long* Mesh_pointer,
                                                 int* nact,
                                                 int* ray_act,
                                                 double* ray_t0p,
                                                 realw* ray_Vs, realw* ray_Vq,
                                                 realw* ray_Ts, realw* ray_Tq,
                                                 int* size_face) {

  TRACE("prepare_rayleigh_injection_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int n = *nact;
  mp->ray_nact = n;

  gpuMalloc_int((void**)&mp->d_ray_act,(size_t)n);
  gpuMemcpy_todevice_int(mp->d_ray_act,ray_act,(size_t)n);

#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_ray_t0p,(size_t)n*sizeof(double)),9100);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ray_t0p,ray_t0p,(size_t)n*sizeof(double),cudaMemcpyHostToDevice),9101);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMalloc((void**)&mp->d_ray_t0p,(size_t)n*sizeof(double)),9100);
    print_HIP_error_if_any(hipMemcpy(mp->d_ray_t0p,ray_t0p,(size_t)n*sizeof(double),hipMemcpyHostToDevice),9101);
  }
#endif

  gpuMalloc_realw((void**)&mp->d_ray_Vs,(size_t)NDIM*n);
  gpuMalloc_realw((void**)&mp->d_ray_Vq,(size_t)NDIM*n);
  gpuMalloc_realw((void**)&mp->d_ray_Ts,(size_t)NDIM*n);
  gpuMalloc_realw((void**)&mp->d_ray_Tq,(size_t)NDIM*n);
  gpuMemcpy_todevice_realw(mp->d_ray_Vs,ray_Vs,(size_t)NDIM*n);
  gpuMemcpy_todevice_realw(mp->d_ray_Vq,ray_Vq,(size_t)NDIM*n);
  gpuMemcpy_todevice_realw(mp->d_ray_Ts,ray_Ts,(size_t)NDIM*n);
  gpuMemcpy_todevice_realw(mp->d_ray_Tq,ray_Tq,(size_t)NDIM*n);

  // d_veloc_inj/d_tract_inj already exist (allocated in prepare_constants_device); zero them
  // once here so inactive points stay at zero for the whole run.
#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMemset(mp->d_veloc_inj,0,(size_t)(*size_face)*sizeof(realw)),9102);
    print_CUDA_error_if_any(cudaMemset(mp->d_tract_inj,0,(size_t)(*size_face)*sizeof(realw)),9103);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMemset(mp->d_veloc_inj,0,(size_t)(*size_face)*sizeof(realw)),9102);
    print_HIP_error_if_any(hipMemset(mp->d_tract_inj,0,(size_t)(*size_face)*sizeof(realw)),9103);
  }
#endif

  GPU_ERROR_CHECKING("prepare_rayleigh_injection_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_rayleigh_injection_cuda,
              COMPUTE_RAYLEIGH_INJECTION_CUDA)(long* Mesh_pointer,
                                               double* tnow_d,
                                               double* om_d,
                                               realw* om,
                                               realw* ig2) {

  TRACE("compute_rayleigh_injection_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int nact = mp->ray_nact;
  if (nact == 0) return;

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)nact)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

#ifdef USE_CUDA
  if (run_cuda){
    rayleigh_injection_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ray_act,
                                                                     mp->d_ray_t0p,
                                                                     mp->d_ray_Vs,
                                                                     mp->d_ray_Vq,
                                                                     mp->d_ray_Ts,
                                                                     mp->d_ray_Tq,
                                                                     mp->d_veloc_inj,
                                                                     mp->d_tract_inj,
                                                                     nact,
                                                                     *tnow_d,
                                                                     *om_d,
                                                                     *om,
                                                                     *ig2);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(rayleigh_injection_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                       mp->d_ray_act,
                       mp->d_ray_t0p,
                       mp->d_ray_Vs,
                       mp->d_ray_Vq,
                       mp->d_ray_Ts,
                       mp->d_ray_Tq,
                       mp->d_veloc_inj,
                       mp->d_tract_inj,
                       nact,
                       *tnow_d,
                       *om_d,
                       *om,
                       *ig2);
  }
#endif

  GPU_ERROR_CHECKING("rayleigh_injection_kernel");
}
