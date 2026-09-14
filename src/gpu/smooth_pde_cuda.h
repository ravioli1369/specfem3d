#ifndef SMOOTH_PDE_CUDA_H
#define SMOOTH_PDE_CUDA_H

#include "mesh_constants_gpu.h"

typedef struct Smooth_pde_data_ {
  field * d_dat_smooth_glob; // data to be smoothed, on device
  field * d_ddat_smooth_glob; // data update, on device
  field * d_send_buffer; // send buffer, on device
  int size_mpi_buffer_smooth;
  realw * d_rvol;
  int * d_phase_ispec_inner;
  int num_phase_ispec;
  int * d_CPML_to_spec;
  int NSPEC_CPML;
  realw cv;
  realw ch;
} Smooth_pde_data;

#endif
