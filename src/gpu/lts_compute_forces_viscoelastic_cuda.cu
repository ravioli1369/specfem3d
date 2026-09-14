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

#include "mesh_constants_gpu.h"

// generalized LTS viscoelastic forces kernel launcher
void lts_compute_forces_viscoelastic_kernel(long* Mesh_pointer,
                                            int num_elements,
                                            int* d_element_list,
                                            realw* d_displ,
                                            realw* d_accel) {

  TRACE("lts_compute_forces_viscoelastic_kernel");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_elements,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // kernel timing
  gpu_event start,stop;
  if (GPU_TIMING){ start_timing_gpu(&start,&stop); }

#ifdef USE_CUDA
  if (run_cuda){
    compute_forces_viscoelastic_cuda_lts_kernel<<<grid,threads,0,mp->compute_stream>>>(num_elements,
                                                                                       mp->d_ibool,
                                                                                       mp->d_irregular_element_number,
                                                                                       d_displ,
                                                                                       d_accel,
                                                                                       mp->d_xix, mp->d_xiy, mp->d_xiz,
                                                                                       mp->d_etax, mp->d_etay, mp->d_etaz,
                                                                                       mp->d_gammax, mp->d_gammay, mp->d_gammaz,
                                                                                       mp->xix_regular, mp->jacobian_regular,
                                                                                       mp->d_hprime_xx,
                                                                                       mp->d_hprimewgll_xx,
                                                                                       mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                       mp->d_kappav, mp->d_muv,
                                                                                       d_element_list);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(compute_forces_viscoelastic_cuda_lts_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                                       num_elements,
                                                                                       mp->d_ibool,
                                                                                       mp->d_irregular_element_number,
                                                                                       d_displ,
                                                                                       d_accel,
                                                                                       mp->d_xix, mp->d_xiy, mp->d_xiz,
                                                                                       mp->d_etax, mp->d_etay, mp->d_etaz,
                                                                                       mp->d_gammax, mp->d_gammay, mp->d_gammaz,
                                                                                       mp->xix_regular, mp->jacobian_regular,
                                                                                       mp->d_hprime_xx,
                                                                                       mp->d_hprimewgll_xx,
                                                                                       mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                                                       mp->d_kappav, mp->d_muv,
                                                                                       d_element_list);
  }
#endif

  // kernel timing
  if (GPU_TIMING){ stop_timing_gpu(&start,&stop,"compute_forces_viscoelastic_cuda_lts_kernel"); }

  GPU_ERROR_CHECKING("compute_forces_viscoelastic_cuda_lts_kernel");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_forces_viscoelastic_cuda_lts,
              COMPUTE_FORCES_VISCOELASTIC_CUDA_LTS)(long* Mesh_pointer,
                                                    int* ilevel_f,
                                                    int* iphase_f) {
  TRACE("compute_forces_viscoelastic_cuda_lts");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int ilevel = *ilevel_f;
  int iphase = *iphase_f;

  int num_elements = mp->lts_num_element_list[(iphase-1) + 2*(ilevel-1)];

  // checks if anything to do
  if (num_elements == 0) return;

  // pointer to current section of lts displacement
  realw* d_displ = &(mp->d_lts_displ_p[NDIM*mp->NGLOB_AB*(ilevel-1)]);

  // "compute" elements only w/out boundary elements
  // pointer to current section of element_list(NSPEC_AB,2,num_p_level)
  int* d_element_list = &(mp->d_lts_element_list[mp->NSPEC_AB*(iphase-1) + mp->NSPEC_AB*2*(ilevel-1)]);

  // launch kernel with list of boundary elements and current tmp_displ
  lts_compute_forces_viscoelastic_kernel(Mesh_pointer,
                                         num_elements,
                                         d_element_list,
                                         d_displ,
                                         mp->d_accel);
}

/* ----------------------------------------------------------------------------------------------- */

// transfer functions

/* ----------------------------------------------------------------------------------------------- */

// not used yet..

//extern EXTERN_LANG
//void FC_FUNC_(copy_field_to_device,
//              COPY_FIELD_TO_DEVICE)(realw* field,
//                                    long* Mesh_pointer) {
//  TRACE("copy_field_to_device");
//  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
//
//  gpuMemcpy_todevice_realw(mp->d_displ, field, NDIM * mp->NGLOB_AB);
//}

/* ----------------------------------------------------------------------------------------------- */

// not used yet..

//extern EXTERN_LANG
//void FC_FUNC_(copy_p_fields_to_device,
//              COPY_P_FIELDS_TO_DEVICE)(realw* field_p,
//                                       realw* field_p2,
//                                       long* Mesh_pointer) {
//  TRACE("copy_p_fields_to_device");
//  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
//
//  gpuMemcpy_todevice_realw(mp->d_lts_displ_p, field_p, NDIM * mp->NGLOB_AB * mp->lts_num_p_level);
//  gpuMemcpy_todevice_realw(mp->d_lts_veloc_p, field_p2, NDIM * mp->NGLOB_AB * mp->lts_num_p_level);
//}

/* ----------------------------------------------------------------------------------------------- */

// not used yet..

//extern EXTERN_LANG
//void FC_FUNC_(copy_accelfield_to_device,
//              COPY_ACCELFIELD_TO_DEVICE)(realw* field,
//                                         long* Mesh_pointer) {
//  TRACE("copy_accelfield_to_device");
//  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
//
//  gpuMemcpy_todevice_realw(mp->d_accel, field, NDIM * mp->NGLOB_AB);
//}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(copy_accelfield_from_device,
              COPY_ACCELFIELD_FROM_DEVICE)(realw* field,
                                           long* Mesh_pointer) {
  TRACE("copy_accelfield_from_device");
  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  gpuMemcpy_tohost_realw(field,mp->d_accel,NDIM * mp->NGLOB_AB);
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(copy_lts_fields_from_device,
              COPY_LTS_FIELDS_FROM_DEVICE)(realw* displ_p,
                                           realw* veloc_p,
                                           realw* accel,
                                           long* Mesh_pointer) {
  TRACE("copy_lts_fields_from_device");
  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  gpuMemcpy_tohost_realw(displ_p, mp->d_lts_displ_p, NDIM * mp->NGLOB_AB * mp->lts_num_p_level);
  gpuMemcpy_tohost_realw(veloc_p, mp->d_lts_veloc_p, NDIM * mp->NGLOB_AB * mp->lts_num_p_level);
  gpuMemcpy_tohost_realw(accel, mp->d_accel, NDIM * mp->NGLOB_AB);
}


/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(copy_lts_fields_to_device,
              COPY_lts_FIELDS_TO_DEVICE)(realw* displ_p,
                                         realw* veloc_p,
                                         realw* accel,
                                         long* Mesh_pointer) {
  TRACE("copy_lts_fields_to_device");
  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  gpuMemcpy_todevice_realw(mp->d_lts_displ_p, displ_p, NDIM * mp->NGLOB_AB * mp->lts_num_p_level);
  gpuMemcpy_todevice_realw(mp->d_lts_veloc_p, veloc_p, NDIM * mp->NGLOB_AB * mp->lts_num_p_level);
  gpuMemcpy_todevice_realw(mp->d_accel, accel, NDIM * mp->NGLOB_AB);
}

/* ----------------------------------------------------------------------------------------------- */

// not used yet..

//extern EXTERN_LANG
//void FC_FUNC_(copy_pfields_global_to_device,
//              COPY_PFIELDS_GLOBAL_TO_DEVICE)(realw* displ_p, realw* veloc_p, realw* accel, long* Mesh_pointer) {
//  TRACE("copy_pfields_global_to_device");
//  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
//
//  gpuMemcpy_todevice_realw(mp->d_displ, displ_p, NDIM * mp->NGLOB_AB);
//  gpuMemcpy_todevice_realw(mp->d_veloc, veloc_p, NDIM * mp->NGLOB_AB);
//  gpuMemcpy_todevice_realw(mp->d_accel, accel, NDIM * mp->NGLOB_AB);
//}

/* ----------------------------------------------------------------------------------------------- */

// not used yet..

//extern EXTERN_LANG
//void FC_FUNC_(copy_global_fields_to_device,
//              COPY_GLOBAL_FIELDS_TO_DEVICE)(realw* displ, realw* veloc, long* Mesh_pointer) {
//  TRACE("copy_global_fields_to_device");
//  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
//
//  gpuMemcpy_todevice_realw(mp->d_displ, displ, NDIM * mp->NGLOB_AB);
//  gpuMemcpy_todevice_realw(mp->d_veloc, veloc, NDIM * mp->NGLOB_AB);
//}

/* ----------------------------------------------------------------------------------------------- */

// boundary updates

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(lts_boundary_contribution_cuda,
              LTS_BOUNDARY_CONTRIBUTION_CUDA)(long* Mesh_pointer,
                                              int* num_points_f,
                                              int* num_elements_f,
                                              int* iphase_f,
                                              int* ilevel_f) {
  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  // checks if nothing to do, skip
  int num_points = *num_points_f;
  if (num_points == 0) return;

  int num_elements = *num_elements_f;
  if (num_elements == 0) return;

  int iphase = *iphase_f;
  int ilevel = *ilevel_f;

  // setup threads/blocks
  int blocksize = BLOCKSIZE_KERNEL1;

  int size_padded = ((int)ceil(((double)num_points)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // TODO: LTS w/ attenuation?
  // velocity used for attenuation and dynamic rupture...
  if (mp->attenuation) {printf("Error: lts_boundary_contribution_cuda() not supported yet for attenuation\n"); exit(1);}

  // see note in routine lts_boundary_contribution().
  // the boundary field d_lts_displ_tmp can be initialized once at beginning during the preparation.
  //gpuMemset_realw(mp->d_lts_displ_tmp,NDIM * mp->NGLOB_AB,0);

  // pointer to boundary nodes
  int* ibool_from = &(mp->d_lts_boundary_node[mp->NGLOB_AB*(iphase-1) + 2*mp->NGLOB_AB*(ilevel-1)]);
  int* ilevel_from = &(mp->d_lts_boundary_ilevel_from[mp->NGLOB_AB*(iphase-1) + 2*mp->NGLOB_AB*(ilevel-1)]);

  // launch kernel to setup temporary boundary array
#ifdef USE_CUDA
  if (run_cuda){
    setup_lts_boundary_array<<<grid,threads,0,mp->compute_stream>>>(ibool_from,
                                                                    ilevel_from,
                                                                    num_points,
                                                                    mp->d_lts_displ_p,
                                                                    mp->d_lts_displ_tmp,
                                                                    mp->d_lts_accel_tmp,
                                                                    mp->NGLOB_AB);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(setup_lts_boundary_array, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                    ibool_from,
                                                                    ilevel_from,
                                                                    num_points,
                                                                    mp->d_lts_displ_p,
                                                                    mp->d_lts_displ_tmp,
                                                                    mp->d_lts_accel_tmp,
                                                                    mp->NGLOB_AB);
  }
#endif
  GPU_ERROR_CHECKING("setup_lts_boundary_array");

  // pointer to current element list
  int* d_element_list = &(mp->d_lts_boundary_ispec[mp->NSPEC_AB*(iphase-1) + mp->NSPEC_AB*2*(ilevel-1)]);

  // launch kernel with list of boundary elements and current tmp_displ
  lts_compute_forces_viscoelastic_kernel(Mesh_pointer,
                                         num_elements,
                                         d_element_list,
                                         mp->d_lts_displ_tmp,
                                         mp->d_lts_accel_tmp);

  // fills element values back into "global" array
  // TODO: check if d_lts_accel_tmp is really needed or if the contributions could be directly added to d_accel
  //       in routine lts_compute_forces_viscoelastic_kernel() above by passing mp->d_accel instead of mp->d_lts_accel_tmp.
  //       at the moment, this mimicks the same handling as done on the CPU routine.
#ifdef USE_CUDA
  if (run_cuda){
    add_lts_boundary_contribution<<<grid,threads,0,mp->compute_stream>>>(ibool_from,
                                                                         num_points,
                                                                         mp->d_lts_accel_tmp,
                                                                         mp->d_accel);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(add_lts_boundary_contribution, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                         ibool_from,
                                                                         num_points,
                                                                         mp->d_lts_accel_tmp,
                                                                         mp->d_accel);  }
#endif
  GPU_ERROR_CHECKING("add_lts_boundary_contribution");
}


/* ----------------------------------------------------------------------------------------------- */

// LTS Newmark updates

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(lts_newmark_update_cuda,
              LTS_NEWMARK_UPDATE_CUDA)(long* Mesh_pointer,
                                       realw* deltat_lts_f,
                                       int* ilevel_f,
                                       int* step_m_f,
                                       int* lts_current_m,
                                       int* p_level_iglob_start,
                                       int* p_level_iglob_end,
                                       int* num_p_level_coarser_to_update) {

  TRACE("lts_newmark_update_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  realw deltat_lts = *deltat_lts_f;
  int ilevel = *ilevel_f;
  int step_m = *step_m_f;

  int blocksize = BLOCKSIZE_KERNEL1;

  int num_blocks_x,num_blocks_y;
  int num_points,size_padded;
  dim3 grid,threads;

  // safety check
  if (p_level_iglob_start[0] != 1) {printf("ASSERT ERROR: p_level_iglob[0] must start at 1\n"); exit(1);}

  // collect accelerations for seismogram & shakemovie outputs
  int collect_accel = 0; // as boolean 0==false/1==true
  if (mp->lts_use_accel_collected){
    // collects acceleration B u_n+1 from this current level ilevel
    // this is computed in the first local iteration (m==1) where the initial condition sets u_0 = u_n+1
    if (ilevel < mp->lts_num_p_level) {
      int next_level_m = lts_current_m[ilevel+1-1]; // fortran uses lts_current_m(ilevel+1), in C indexing starts at 0 -> lts_current_m[ilevel+1-1]
      // only store if accel was computed the very first time this level was called
      if (step_m == 1 && next_level_m == 1) collect_accel = 1;
    } else {
      // coarsest p-level (ilevel == num_p_level)
      collect_accel = 1;
    }
  }

  // ---- updates to P nodes -----------------------------------
  // contiguous range of p-level nodes
  // all up to current level

  // compute updates for DOFs under this ilevel (skip for ilevel==1)
  if (ilevel > 1) {
    num_points = p_level_iglob_end[ilevel-2];
    if (num_points > 0) {
      size_padded = ((int)ceil(((double)num_points)/((double)blocksize)))*blocksize;
      get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

      grid = dim3(num_blocks_x,num_blocks_y);
      threads = dim3(blocksize,1,1);

#ifdef USE_CUDA
      if (run_cuda){
        lts_newmark_update_kernel_P<<<grid,threads,0,mp->compute_stream>>>(mp->d_lts_displ_p,
                                                                           mp->d_lts_veloc_p,
                                                                           mp->d_accel,
                                                                           mp->d_lts_rmassxyz,
                                                                           mp->d_lts_rmassxyz_mod,
                                                                           mp->d_lts_cmassxyz,
                                                                           ilevel,
                                                                           step_m,
                                                                           mp->NGLOB_AB,
                                                                           num_points,
                                                                           mp->lts_num_p_level,deltat_lts,
                                                                           collect_accel,mp->d_lts_accel_collected,mp->d_lts_mask_ibool_collected);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(lts_newmark_update_kernel_P, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                           mp->d_lts_displ_p,
                                                                           mp->d_lts_veloc_p,
                                                                           mp->d_accel,
                                                                           mp->d_lts_rmassxyz,
                                                                           mp->d_lts_rmassxyz_mod,
                                                                           mp->d_lts_cmassxyz,
                                                                           ilevel,
                                                                           step_m,
                                                                           mp->NGLOB_AB,
                                                                           num_points,
                                                                           mp->lts_num_p_level,deltat_lts,
                                                                           collect_accel,mp->d_lts_accel_collected,mp->d_lts_mask_ibool_collected);
      }
#endif
      GPU_ERROR_CHECKING("lts_newmark_update_kernel_P");
    }
  } // ilevel > 1

  // start index of current p-level
  int is = p_level_iglob_start[ilevel-1];
  // end index of current p-level
  int ie = p_level_iglob_end[ilevel-1];

  if (ie - is > 0) {
    // size has (+1)
    // in case ie == is -> single node
    // in case ie > is  -> fortran range (is:ie) has ie-is+1 elements, for example: (1:5) in C becomes [0,1,2,3,4] with size == 5 - 1 + 1
    num_points = ie - is + 1;

    size_padded = ((int)ceil(((double)num_points)/((double)blocksize)))*blocksize;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

#ifdef USE_CUDA
    if (run_cuda){
      lts_newmark_update_kernel_P2<<<grid,threads,0,mp->compute_stream>>>(mp->d_lts_displ_p,
                                                                          mp->d_lts_veloc_p,
                                                                          mp->d_accel,
                                                                          mp->d_lts_rmassxyz,
                                                                          mp->d_lts_rmassxyz_mod,
                                                                          ilevel,
                                                                          step_m,
                                                                          mp->NGLOB_AB,
                                                                          is,
                                                                          num_points,
                                                                          mp->lts_num_p_level,deltat_lts,
                                                                          collect_accel,mp->d_lts_accel_collected,mp->d_lts_mask_ibool_collected);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(lts_newmark_update_kernel_P2, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                          mp->d_lts_displ_p,
                                                                          mp->d_lts_veloc_p,
                                                                          mp->d_accel,
                                                                          mp->d_lts_rmassxyz,
                                                                          mp->d_lts_rmassxyz_mod,
                                                                          ilevel,
                                                                          step_m,
                                                                          mp->NGLOB_AB,
                                                                          is,
                                                                          num_points,
                                                                          mp->lts_num_p_level,deltat_lts,
                                                                          collect_accel,mp->d_lts_accel_collected,mp->d_lts_mask_ibool_collected);
    }
#endif
    GPU_ERROR_CHECKING("lts_newmark_update_kernel_P2");
  }

  // ---------------- updates to R ----------------------------------------
  // p-level boundary updates
  // coarsest level is already finished from above
  if (ilevel < mp->lts_num_p_level) {
    num_points = num_p_level_coarser_to_update[ilevel-1];
    if (num_points > 0) {
      size_padded = ((int)ceil(((double)num_points)/((double)blocksize)))*blocksize;
      get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

      grid = dim3(num_blocks_x,num_blocks_y);
      threads = dim3(blocksize,1,1);

#ifdef USE_CUDA
      if (run_cuda){
        lts_newmark_update_kernel_R<<<grid,threads,0,mp->compute_stream>>>(mp->d_lts_displ_p,
                                                                           mp->d_lts_veloc_p,
                                                                           mp->d_accel,
                                                                           mp->d_lts_rmassxyz,
                                                                           mp->d_lts_rmassxyz_mod,
                                                                           ilevel,
                                                                           step_m,
                                                                           mp->NGLOB_AB,mp->lts_num_p_level,deltat_lts,
                                                                           num_p_level_coarser_to_update[ilevel-1],
                                                                           &(mp->d_lts_p_level_coarser_to_update[mp->NGLOB_AB*(ilevel-1)]),
                                                                           mp->d_lts_iglob_p_refine,
                                                                           collect_accel,mp->d_lts_accel_collected,mp->d_lts_mask_ibool_collected);
      }
#endif
#ifdef USE_HIP
      if (run_hip){
        hipLaunchKernelGGL(lts_newmark_update_kernel_R, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                           mp->d_lts_displ_p,
                                                                           mp->d_lts_veloc_p,
                                                                           mp->d_accel,
                                                                           mp->d_lts_rmassxyz,
                                                                           mp->d_lts_rmassxyz_mod,
                                                                           ilevel,
                                                                           step_m,
                                                                           mp->NGLOB_AB,mp->lts_num_p_level,deltat_lts,
                                                                           num_p_level_coarser_to_update[ilevel-1],
                                                                           &(mp->d_lts_p_level_coarser_to_update[mp->NGLOB_AB*(ilevel-1)]),
                                                                           mp->d_lts_iglob_p_refine,
                                                                           collect_accel,mp->d_lts_accel_collected,mp->d_lts_mask_ibool_collected);
      }
#endif
      GPU_ERROR_CHECKING("lts_newmark_update_kernel_R");
    }
  }

  // updates global main wavefields
  if (ilevel == mp->lts_num_p_level) {
    // each thread is one global point and updates it 3 components, i.e., uses NGLOB as size - not NDIM * NGLOB
    size_padded = ((int)ceil(((double)mp->NGLOB_AB)/((double)blocksize)))*blocksize;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

#ifdef USE_CUDA
    if (run_cuda){
      lts_newmark_update_veloc_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_lts_veloc_p,
                                                                             mp->d_veloc,
                                                                             mp->d_accel,
                                                                             mp->NGLOB_AB,
                                                                             mp->lts_num_p_level,
                                                                             mp->lts_use_accel_collected,mp->d_lts_accel_collected,mp->d_lts_mask_ibool_collected);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(lts_newmark_update_veloc_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                             mp->d_lts_veloc_p,
                                                                             mp->d_veloc,
                                                                             mp->d_accel,
                                                                             mp->NGLOB_AB,
                                                                             mp->lts_num_p_level,
                                                                             mp->lts_use_accel_collected,mp->d_lts_accel_collected,mp->d_lts_mask_ibool_collected);
    }
#endif
    GPU_ERROR_CHECKING("lts_newmark_update_global_wavefields_kernel");
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(lts_newmark_update_displ_cuda,
              LTS_NEWMARK_UPDATE_displ_CUDA)(long* Mesh_pointer) {

  TRACE("lts_newmark_update_displ_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int blocksize = BLOCKSIZE_KERNEL1;

  // each thread is one global point and updates it 3 components, i.e., uses NGLOB as size - not NDIM * NGLOB
  int size_padded = ((int)ceil(((double)mp->NGLOB_AB)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

#ifdef USE_CUDA
  if (run_cuda){
    lts_newmark_update_displ_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_lts_displ_p,
                                                                           mp->d_displ,
                                                                           mp->NGLOB_AB,
                                                                           mp->lts_num_p_level);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    hipLaunchKernelGGL(lts_newmark_update_displ_kernel, dim3(grid), dim3(threads), 0, mp->compute_stream,
                                                                           mp->d_lts_displ_p,
                                                                           mp->d_displ,
                                                                           mp->NGLOB_AB,
                                                                           mp->lts_num_p_level);
  }
#endif
  GPU_ERROR_CHECKING("lts_newmark_update_displ_kernel");
}



/* ----------------------------------------------------------------------------------------------- */

// helper functions

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(lts_set_accel_zero_cuda,
              LTS_SET_ACCEL_ZERO_CUDA)(long* Mesh_pointer) {
  TRACE("lts_set_accel_zero_cuda");
  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  gpuMemset_realw(mp->d_accel,NDIM * mp->NGLOB_AB,0);
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(lts_set_accel_zero_level_cuda,
              LTS_SET_ACCEL_ZERO_LEVEL_CUDA)(long* Mesh_pointer,
                                             int* ilevel_f,
                                             int* p_level_iglob_end,
                                             int* num_p_level_coarser_to_update) {
  TRACE("lts_set_accel_zero_level_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int ilevel = *ilevel_f;

  int num_blocks_x, num_blocks_y;
  int size,size_padded,num_points;

  int blocksize = BLOCKSIZE_KERNEL1;

  // gpuMemset_realw(&mp->d_lts_displ_p[3*mp->NGLOB_AB*(ilevel-1)],NDIM*mp->NGLOB_AB,0);

  size = p_level_iglob_end[ilevel-1];
  if (size > 0) {
    size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 gridP(num_blocks_x,num_blocks_y);
    dim3 threadsP(blocksize,1,1);

#ifdef USE_CUDA
    if (run_cuda){
      zero_accel_kernel_P<<<gridP,threadsP,0,mp->compute_stream>>>(mp->d_accel,
                                                                   size);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(zero_accel_kernel_P, dim3(gridP), dim3(threadsP), 0, mp->compute_stream,
                                                                   mp->d_accel,
                                                                   size);
    }
#endif
    GPU_ERROR_CHECKING("zero_accel_kernel_P");
  }

  num_points = num_p_level_coarser_to_update[ilevel-1];
  if (num_points > 0) {
    size_padded = ((int)ceil(((double)num_points)/((double)blocksize)))*blocksize;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 gridR(num_blocks_x,num_blocks_y);
    dim3 threadsR(blocksize,1,1);

#ifdef USE_CUDA
    if (run_cuda){
      zero_accel_kernel_R<<<gridR,threadsR,0,mp->compute_stream>>>(mp->d_accel,
                                                                   num_points,
                                                                   &(mp->d_lts_p_level_coarser_to_update[mp->NGLOB_AB*(ilevel-1)]) );
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(zero_accel_kernel_R, dim3(gridR), dim3(threadsR), 0, mp->compute_stream,
                                                                   mp->d_accel,
                                                                   num_points,
                                                                   &(mp->d_lts_p_level_coarser_to_update[mp->NGLOB_AB*(ilevel-1)]) );
    }
#endif
    GPU_ERROR_CHECKING("zero_accel_kernel_R");
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(lts_set_finer_initial_condition_cuda,
              LTS_SET_FINER_INITIAL_CONDITION_CUDA)(long* Mesh_pointer,
                                                    int* ilevel_f,
                                                    int* p_level_iglob_end,
                                                    int* num_p_level_coarser_to_update) {
  TRACE("lts_set_finer_initial_condition_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int ilevel = *ilevel_f;

  int num_blocks_x, num_blocks_y;
  int size,size_padded,num_points;

  // P matrix
  int blocksize = BLOCKSIZE_KERNEL1;

  size = p_level_iglob_end[ilevel-1];
  if (size > 0) {
    size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 gridP(num_blocks_x,num_blocks_y);
    dim3 threadsP(blocksize,1,1);

#ifdef USE_CUDA
    if (run_cuda){
      set_finer_initial_condition_fast_kernel_P<<<gridP,threadsP,0,mp->compute_stream>>>(mp->d_lts_displ_p,
                                                                                         ilevel,
                                                                                         mp->NGLOB_AB,
                                                                                         size);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(set_finer_initial_condition_fast_kernel_P, dim3(gridP), dim3(threadsP), 0, mp->compute_stream,
                                                                                         mp->d_lts_displ_p,
                                                                                         ilevel,
                                                                                         mp->NGLOB_AB,
                                                                                         size);
    }
#endif
    GPU_ERROR_CHECKING("set_finer_initial_condition_fast_kernel_P");
  }

  // R matrix
  num_points = num_p_level_coarser_to_update[ilevel-1];
  if (num_points > 0) {
    size_padded = ((int)ceil(((double)num_points)/((double)blocksize)))*blocksize;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 gridR(num_blocks_x,num_blocks_y);
    dim3 threadsR(blocksize,1,1);

#ifdef USE_CUDA
    if (run_cuda){
      set_finer_initial_condition_fast_kernel_R<<<gridR,threadsR,0,mp->compute_stream>>>(mp->d_lts_displ_p,
                                                                                         ilevel,
                                                                                         mp->NGLOB_AB,
                                                                                         num_points,
                                                                                         &(mp->d_lts_p_level_coarser_to_update[mp->NGLOB_AB*(ilevel-1)]) );
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      hipLaunchKernelGGL(set_finer_initial_condition_fast_kernel_R, dim3(gridR), dim3(threadsR), 0, mp->compute_stream,
                                                                                         mp->d_lts_displ_p,
                                                                                         ilevel,
                                                                                         mp->NGLOB_AB,
                                                                                         num_points,
                                                                                         &(mp->d_lts_p_level_coarser_to_update[mp->NGLOB_AB*(ilevel-1)]) );
    }
#endif
    GPU_ERROR_CHECKING("set_finer_initial_condition_fast_kernel_R");
  }
}


