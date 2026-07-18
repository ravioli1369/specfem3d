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

// Newtonian-noise gravity-perturbation integral on the device.
//
// This replaces the per-output-step transfer of the full displacement field back to the
// host (NDIM*NGLOB_AB reals, ~761 MB on production runs) followed by a CPU reduction. The
// time-invariant per-node per-station weights w3/w5 (built once in gravity_init, folding in
// G, rho0_wm and the 1/Rg**3, 3/Rg**5 geometric factors) and the node coordinates are
// uploaded once by prepare_gravity_device(); compute_gravity_cuda() then reduces over
// d_displ directly on the GPU each output step and returns only the per-station E/N/Z sums.

#include "mesh_constants_gpu.h"


/* ----------------------------------------------------------------------------------------------- */

// block-strided reduction over NGLOB for a single gravity station.
// each thread accumulates its node contribution in DOUBLE precision (the CPU sum is
// single-precision and sequential; double accumulation here is strictly more accurate and
// is the one intended numerical change), then a shared-memory tree reduction produces three
// double partials per thread-block (E,N,Z), matching the CPU per-node formula term-for-term:
//   dotP = dx*u1 + dy*u2 + dz*u3
//   E   += w3*u1 - w5*dx*dotP   (and N with u2/dy, Z with u3/dz)
// the per-node term is formed in realw (as on the CPU) and only the accumulation is double.

__global__ void gravity_reduction_kernel(realw* d_displ,
                                         realw* d_grav_w3_col,
                                         realw* d_grav_w5_col,
                                         realw* d_grav_xstore,
                                         realw* d_grav_ystore,
                                         realw* d_grav_zstore,
                                         realw xstat, realw ystat, realw zstat,
                                         int nglob,
                                         double* d_partE, double* d_partN, double* d_partZ) {

  // shared-memory scratch for the per-block tree reduction (one slot per thread, per component)
  __shared__ double sE[BLOCKSIZE_TRANSFER];
  __shared__ double sN[BLOCKSIZE_TRANSFER];
  __shared__ double sZ[BLOCKSIZE_TRANSFER];

  unsigned int tid = threadIdx.x;
  unsigned int bx = blockIdx.y*gridDim.x + blockIdx.x;
  unsigned int i = tid + bx*blockDim.x;

  // per-thread double accumulators (zero for padding threads with i >= nglob)
  double eE = 0.0;
  double eN = 0.0;
  double eZ = 0.0;

  if (i < nglob){
    // node-to-station offset (kept inline, as on the CPU) and the displacement at this node.
    // displ is stored as displ(NDIM,iglob) column-major => x,y,z at d_displ[3*i..3*i+2].
    realw dx = d_grav_xstore[i] - xstat;
    realw dy = d_grav_ystore[i] - ystat;
    realw dz = d_grav_zstore[i] - zstat;

    realw u1 = d_displ[i*3];
    realw u2 = d_displ[i*3+1];
    realw u3 = d_displ[i*3+2];

    realw dotP = dx*u1 + dy*u2 + dz*u3;

    realw w3 = d_grav_w3_col[i];
    realw w5 = d_grav_w5_col[i];

    eE = (double)(w3*u1 - w5*dx*dotP);
    eN = (double)(w3*u2 - w5*dy*dotP);
    eZ = (double)(w3*u3 - w5*dz*dotP);
  }

  sE[tid] = eE;
  sN[tid] = eN;
  sZ[tid] = eZ;

  __syncthreads();

  // tree reduction in shared memory (blockDim.x is a power of two: BLOCKSIZE_TRANSFER)
  for (unsigned int s = blockDim.x/2; s > 0; s >>= 1){
    if (tid < s){
      sE[tid] += sE[tid + s];
      sN[tid] += sN[tid + s];
      sZ[tid] += sZ[tid + s];
    }
    __syncthreads();
  }

  // one partial per block, per component
  if (tid == 0){
    d_partE[bx] = sE[0];
    d_partN[bx] = sN[0];
    d_partZ[bx] = sZ[0];
  }
}


/* ----------------------------------------------------------------------------------------------- */

// uploads the time-invariant gravity weights + node coordinates to the device once,
// and stores the (small) station coordinate arrays on the mesh struct.

extern EXTERN_LANG
void FC_FUNC_(prepare_gravity_device,
              PREPARE_GRAVITY_DEVICE)(long* Mesh_pointer,
                                      realw* w3, realw* w5,
                                      realw* xstore, realw* ystore, realw* zstore,
                                      int* NGLOB_AB, int* nstat,
                                      realw* xstat, realw* ystat, realw* zstat) {

  TRACE("prepare_gravity_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int nglob = *NGLOB_AB;
  int ns = *nstat;

  mp->gravity_nstat = ns;

  // node coordinates (constant in time)
  gpuMalloc_realw((void**)&mp->d_grav_xstore,(size_t)nglob);
  gpuMalloc_realw((void**)&mp->d_grav_ystore,(size_t)nglob);
  gpuMalloc_realw((void**)&mp->d_grav_zstore,(size_t)nglob);
  gpuMemcpy_todevice_realw(mp->d_grav_xstore,xstore,(size_t)nglob);
  gpuMemcpy_todevice_realw(mp->d_grav_ystore,ystore,(size_t)nglob);
  gpuMemcpy_todevice_realw(mp->d_grav_zstore,zstore,(size_t)nglob);

  // per-node per-station weights, stored NGLOB_AB x nstat (fortran column-major):
  // element (iglob,istat) is at linear index istat*nglob + iglob.
  gpuMalloc_realw((void**)&mp->d_grav_w3,(size_t)nglob*ns);
  gpuMalloc_realw((void**)&mp->d_grav_w5,(size_t)nglob*ns);
  gpuMemcpy_todevice_realw(mp->d_grav_w3,w3,(size_t)nglob*ns);
  gpuMemcpy_todevice_realw(mp->d_grav_w5,w5,(size_t)nglob*ns);

  // station coordinates are tiny; keep host copies on the struct and pass them as scalar
  // kernel arguments per launch (no device array needed for these).
  mp->h_grav_xstat = (realw*) malloc(ns*sizeof(realw));
  mp->h_grav_ystat = (realw*) malloc(ns*sizeof(realw));
  mp->h_grav_zstat = (realw*) malloc(ns*sizeof(realw));
  if (mp->h_grav_xstat == NULL || mp->h_grav_ystat == NULL || mp->h_grav_zstat == NULL){
    exit_on_error("Error allocating host station coordinate arrays in prepare_gravity_device");
  }
  memcpy(mp->h_grav_xstat,xstat,ns*sizeof(realw));
  memcpy(mp->h_grav_ystat,ystat,ns*sizeof(realw));
  memcpy(mp->h_grav_zstat,zstat,ns*sizeof(realw));

  GPU_ERROR_CHECKING("prepare_gravity_device");
}


/* ----------------------------------------------------------------------------------------------- */

// computes, for every gravity station, this rank's partial E/N/Z gravity-perturbation sums
// by reducing over d_displ on the device. returns the three per-station scalars (double
// accumulated on the device, cast back to CUSTOM_REAL) so the caller can MPI-reduce and store
// them; the full displacement field never leaves the GPU.

extern EXTERN_LANG
void FC_FUNC_(compute_gravity_cuda,
              COMPUTE_GRAVITY_CUDA)(long* Mesh_pointer,
                                    realw* h_grav_E,
                                    realw* h_grav_N,
                                    realw* h_grav_Z) {

  TRACE("compute_gravity_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int nglob = mp->NGLOB_AB;
  int nstat = mp->gravity_nstat;

  // reduction geometry: one thread per node, grid padded up to a multiple of the block size
  // (mirrors get_norm_elastic_from_device / get_maximum_vector_kernel).
  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)nglob)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  int num_blocks = num_blocks_x*num_blocks_y;

  // one contiguous double partials buffer holding [E block | N block | Z block]
  double* d_part = NULL;
  double* h_part = (double*) calloc(3*num_blocks,sizeof(double));
  if (!h_part){ exit_on_error("Error allocating temporary host gravity partials array"); }

#ifdef USE_CUDA
  if (run_cuda){
    print_CUDA_error_if_any(cudaMalloc((void**)&d_part,3*num_blocks*sizeof(double)),9000);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    print_HIP_error_if_any(hipMalloc((void**)&d_part,3*num_blocks*sizeof(double)),9000);
  }
#endif

  double* d_partE = d_part;
  double* d_partN = d_part + num_blocks;
  double* d_partZ = d_part + 2*num_blocks;

  // one launch per station over the identical NGLOB range as the CPU (no free-surface special-casing)
  for (int istat = 0; istat < nstat; istat++){

    realw xstat = mp->h_grav_xstat[istat];
    realw ystat = mp->h_grav_ystat[istat];
    realw zstat = mp->h_grav_zstat[istat];

    // column istat of the NGLOB x nstat weight arrays
    realw* w3_col = mp->d_grav_w3 + (size_t)istat*nglob;
    realw* w5_col = mp->d_grav_w5 + (size_t)istat*nglob;

#ifdef USE_CUDA
    if (run_cuda){
      gravity_reduction_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_displ,
                                                                      w3_col,w5_col,
                                                                      mp->d_grav_xstore,
                                                                      mp->d_grav_ystore,
                                                                      mp->d_grav_zstore,
                                                                      xstat,ystat,zstat,
                                                                      nglob,
                                                                      d_partE,d_partN,d_partZ);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(gravity_reduction_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                         mp->d_displ,
                         w3_col,w5_col,
                         mp->d_grav_xstore,
                         mp->d_grav_ystore,
                         mp->d_grav_zstore,
                         xstat,ystat,zstat,
                         nglob,
                         d_partE,d_partN,d_partZ);
    }
#endif

    GPU_ERROR_CHECKING("gravity_reduction_kernel");

    // explicitly waits for the stream to finish before copying the partials back
    gpuStreamSynchronize(mp->compute_stream);

#ifdef USE_CUDA
    if (run_cuda){
      print_CUDA_error_if_any(cudaMemcpy(h_part,d_part,3*num_blocks*sizeof(double),cudaMemcpyDeviceToHost),9010);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      print_HIP_error_if_any(hipMemcpy(h_part,d_part,3*num_blocks*sizeof(double),hipMemcpyDeviceToHost),9010);
    }
#endif

    // final sum of block partials in double precision, then cast back to CUSTOM_REAL
    double sumE = 0.0;
    double sumN = 0.0;
    double sumZ = 0.0;
    for (int i = 0; i < num_blocks; i++){
      sumE += h_part[i];
      sumN += h_part[num_blocks + i];
      sumZ += h_part[2*num_blocks + i];
    }

    h_grav_E[istat] = (realw) sumE;
    h_grav_N[istat] = (realw) sumN;
    h_grav_Z[istat] = (realw) sumZ;
  }

  free(h_part);
  gpuFree(d_part);

  GPU_ERROR_CHECKING("compute_gravity_cuda");
}
