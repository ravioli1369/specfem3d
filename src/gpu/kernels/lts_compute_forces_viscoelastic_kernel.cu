/*
!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!    Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                             CNRS, France
!                      and Princeton University, USA
!                (there are currently many more authors!)
!                          (c) October 2017
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


#ifdef USE_TEXTURES_CONSTANTS
realw_texture d_hprime_xx_tex;
#endif


/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_forces_viscoelastic_cuda_lts_kernel (const int nb_blocks_to_compute,
                                                             const int* d_ibool,
                                                             const int* d_irregular_element_number,
                                                             realw_const_p d_displ,
                                                             realw* d_accel,
                                                             realw_const_p d_xix, realw_const_p d_xiy, realw_const_p d_xiz,
                                                             realw_const_p d_etax, realw_const_p d_etay, realw_const_p d_etaz,
                                                             realw_const_p d_gammax, realw_const_p d_gammay, realw_const_p d_gammaz,
                                                             const realw xix_regular, const realw jacobian_regular,
                                                             realw_const_p d_hprime_xx,
                                                             realw_const_p d_hprimewgll_xx,
                                                             realw_const_p d_wgllwgll_xy,
                                                             realw_const_p d_wgllwgll_xz,
                                                             realw_const_p d_wgllwgll_yz,
                                                             realw_const_p d_kappav, realw_const_p d_muv,
                                                             const int* element_list) {

  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int active,offset;
  int iglob = 0;
  int working_element, ispec_irreg;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  __shared__ realw s_tempx1[NGLL3];
  __shared__ realw s_tempx2[NGLL3];
  __shared__ realw s_tempx3[NGLL3];

  __shared__ realw s_tempy1[NGLL3];
  __shared__ realw s_tempy2[NGLL3];
  __shared__ realw s_tempy3[NGLL3];

  __shared__ realw s_tempz1[NGLL3];
  __shared__ realw s_tempz2[NGLL3];
  __shared__ realw s_tempz3[NGLL3];

  __shared__ realw sh_hprime[NGLL2];
  __shared__ realw sh_hprimewgll[NGLL2];

  // use only NGLL^3 = 125 active threads, plus 3 inactive/ghost
  // threads, because we used memory padding from NGLL^3 = 125 to 128
  // to get coalescent memory accesses
  active = (tx < NGLL3 && bx < nb_blocks_to_compute) ? 1:0;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (active) {
    // spectral-element id
    working_element = element_list[bx] - 1; // note: fortran 1-indexing

    ispec_irreg = d_irregular_element_number[working_element] - 1;

    // local padded index
    offset = working_element*NGLL3_PADDED + tx;

    // global index
    iglob = d_ibool[offset] - 1;

    // changing iglob indexing to match fortran row changes fast style
    s_dummyx_loc[tx] = d_displ[iglob*3];
    s_dummyy_loc[tx] = d_displ[iglob*3 + 1];
    s_dummyz_loc[tx] = d_displ[iglob*3 + 2];

    if (tx < NGLL2) {
#ifdef USE_TEXTURES_CONSTANTS
      sh_hprime[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
      sh_hprimewgll[tx] = tex1Dfetch(d_hprimewgll_xx_tex,tx);
#else
      sh_hprime[tx] = d_hprime_xx[tx];
      sh_hprimewgll[tx] = d_hprimewgll_xx[tx];
#endif
    }
  }

  // synchronize all the threads (one thread for each of the NGLL
  // grid points of the current spectral element) because we need
  // the whole element to be ready in order to be able to compute
  // the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active) {
    tempx1l = 0.f;
    tempx2l = 0.f;
    tempx3l = 0.f;

    tempy1l = 0.f;
    tempy2l = 0.f;
    tempy3l = 0.f;

    tempz1l = 0.f;
    tempz2l = 0.f;
    tempz3l = 0.f;

    for (int l=0;l<NGLLX;l++) {
      realw hp1 = sh_hprime[l*NGLLX+I];
      int index1 = K*NGLL2+J*NGLLX+l;
      tempx1l += s_dummyx_loc[index1]*hp1;
      tempy1l += s_dummyy_loc[index1]*hp1;
      tempz1l += s_dummyz_loc[index1]*hp1;

      //assumes that hprime_xx = hprime_yy = hprime_zz
      realw hp2 = sh_hprime[l*NGLLX+J];
      int index2 = K*NGLL2+l*NGLLX+I;
      tempx2l += s_dummyx_loc[index2]*hp2;
      tempy2l += s_dummyy_loc[index2]*hp2;
      tempz2l += s_dummyz_loc[index2]*hp2;

      realw hp3 = sh_hprime[l*NGLLX+K];
      int index3 = l*NGLL2+J*NGLLX+I;
      tempx3l += s_dummyx_loc[index3]*hp3;
      tempy3l += s_dummyy_loc[index3]*hp3;
      tempz3l += s_dummyz_loc[index3]*hp3;
    }

    // compute derivatives of ux, uy and uz with respect to x, y and z
    if (ispec_irreg >= 0){
      xixl = d_xix[offset];
      xiyl = d_xiy[offset];
      xizl = d_xiz[offset];
      etaxl = d_etax[offset];
      etayl = d_etay[offset];
      etazl = d_etaz[offset];
      gammaxl = d_gammax[offset];
      gammayl = d_gammay[offset];
      gammazl = d_gammaz[offset];

      jacobianl = 1.0f / (xixl*(etayl*gammazl-etazl*gammayl)-
                          xiyl*(etaxl*gammazl-etazl*gammaxl)+
                          xizl*(etaxl*gammayl-etayl*gammaxl));

      duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
      duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
      duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

      duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
      duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
      duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

      duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
      duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
      duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

    } else {
      // regular element
      jacobianl = jacobian_regular;

      duxdxl = xix_regular*tempx1l;
      duxdyl = xix_regular*tempx2l;
      duxdzl = xix_regular*tempx3l;

      duydxl = xix_regular*tempy1l;
      duydyl = xix_regular*tempy2l;
      duydzl = xix_regular*tempy3l;

      duzdxl = xix_regular*tempz1l;
      duzdyl = xix_regular*tempz2l;
      duzdzl = xix_regular*tempz3l;
    }

    // precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl;
    duxdxl_plus_duzdzl = duxdxl + duzdzl;
    duydyl_plus_duzdzl = duydyl + duzdzl;
    duxdyl_plus_duydxl = duxdyl + duydxl;
    duzdxl_plus_duxdzl = duzdxl + duxdzl;
    duzdyl_plus_duydzl = duzdyl + duydzl;

    // compute elements with an elastic isotropic rheology
    kappal = d_kappav[offset];
    mul = d_muv[offset];

    // isotropic case
    lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
    lambdal = lambdalplus2mul - 2.0f * mul;

    // compute the six components of the stress tensor sigma
    sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
    sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
    sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

    sigma_xy = mul*duxdyl_plus_duydxl;
    sigma_xz = mul*duzdxl_plus_duxdzl;
    sigma_yz = mul*duzdyl_plus_duydzl;

    // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
    sigma_yx = sigma_xy;
    sigma_zx = sigma_xz;
    sigma_zy = sigma_yz;

    // form dot product with test vector, non-symmetric form
    if (ispec_irreg >= 0){
      s_tempx1[tx] = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl);
      s_tempy1[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl);
      s_tempz1[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl);

      s_tempx2[tx] = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl);
      s_tempy2[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl);
      s_tempz2[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl);

      s_tempx3[tx] = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl);
      s_tempy3[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl);
      s_tempz3[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl);
    } else {
      // regular element
      s_tempx1[tx] = jacobian_regular * (sigma_xx*xix_regular);
      s_tempy1[tx] = jacobian_regular * (sigma_xy*xix_regular);
      s_tempz1[tx] = jacobian_regular * (sigma_xz*xix_regular);

      s_tempx2[tx] = jacobian_regular * (sigma_yx*xix_regular);
      s_tempy2[tx] = jacobian_regular * (sigma_yy*xix_regular);
      s_tempz2[tx] = jacobian_regular * (sigma_yz*xix_regular);

      s_tempx3[tx] = jacobian_regular * (sigma_zx*xix_regular);
      s_tempy3[tx] = jacobian_regular * (sigma_zy*xix_regular);
      s_tempz3[tx] = jacobian_regular * (sigma_zz*xix_regular);
    }
  }

  // synchronize all the threads (one thread for each of the NGLL
  // grid points of the current spectral element) because we need
  // the whole element to be ready in order to be able to compute
  // the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active) {
    tempx1l = 0.f;
    tempy1l = 0.f;
    tempz1l = 0.f;

    tempx2l = 0.f;
    tempy2l = 0.f;
    tempz2l = 0.f;

    tempx3l = 0.f;
    tempy3l = 0.f;
    tempz3l = 0.f;

    for (int l=0;l<NGLLX;l++) {
      realw hp1 = sh_hprimewgll[I*NGLLX+l];
      int index1 = K*NGLL2+J*NGLLX+l;
      tempx1l += s_tempx1[index1]*hp1;
      tempy1l += s_tempy1[index1]*hp1;
      tempz1l += s_tempz1[index1]*hp1;

      // assumes hprimewgll_xx == hprimewgll_yy == hprimewgll_zz
      realw hp2 = sh_hprimewgll[J*NGLLX+l];
      int index2 = K*NGLL2+l*NGLLX+I;
      tempx2l += s_tempx2[index2]*hp2;
      tempy2l += s_tempy2[index2]*hp2;
      tempz2l += s_tempz2[index2]*hp2;

      realw hp3 = sh_hprimewgll[K*NGLLX+l];
      int index3 = l*NGLL2+J*NGLLX+I;
      tempx3l += s_tempx3[index3]*hp3;
      tempy3l += s_tempy3[index3]*hp3;
      tempz3l += s_tempz3[index3]*hp3;
    }

    realw fac1 = d_wgllwgll_yz[K*NGLLX+J];
    realw fac2 = d_wgllwgll_xz[K*NGLLX+I];
    realw fac3 = d_wgllwgll_xy[J*NGLLX+I];

    realw sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
    realw sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
    realw sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

    atomicAdd(&d_accel[iglob*3], sum_terms1);
    atomicAdd(&d_accel[iglob*3+1], sum_terms2);
    atomicAdd(&d_accel[iglob*3+2], sum_terms3);
  } // if (active)
}


/* ----------------------------------------------------------------------------------------------- */
// boundary kernel
/* ----------------------------------------------------------------------------------------------- */

__global__ void setup_lts_boundary_array(const int* ibool_from,
                                         const int* ilevel_from,
                                         int num_points,
                                         const realw* displ_p,
                                         realw* displ_tmp,
                                         realw* accel_tmp,
                                         int nglob_ab) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  if (id < num_points) {
    // put appropriate values into displ_tmp
    // displ_tmp is full sized, so we just copy them into the same iglob locations
    // note (-1) due to fortran 1-indexing
    int iglob = ibool_from[id] - 1;
    int iilevel = ilevel_from[id] - 1;

    int offset = NDIM * nglob_ab * iilevel;

    displ_tmp[iglob*3]   = displ_p[iglob*3 + offset];
    displ_tmp[iglob*3+1] = displ_p[iglob*3+1 + offset];
    displ_tmp[iglob*3+2] = displ_p[iglob*3+2 + offset];

    // zeros out accel
    accel_tmp[iglob*3]   = 0.f;
    accel_tmp[iglob*3+1] = 0.f;
    accel_tmp[iglob*3+2] = 0.f;

    // if (threadIdx.x == 63 && blockIdx.x == 0) {
    //   printf("id=%d,ibool=%d,ilevel=%d\n",id,ibool_from[id]-1,ilevel_from[id]-1);
    // }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void add_lts_boundary_contribution(const int* ibool_from,
                                              int num_points,
                                              const realw* accel_tmp,
                                              realw* accel) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  if (id < num_points) {
    // add appropriate boundary values into accel
    int iglob = ibool_from[id] - 1;  // (-1) due to fortran 1-indexing

    accel[iglob*3]   += accel_tmp[iglob*3];
    accel[iglob*3+1] += accel_tmp[iglob*3+1];
    accel[iglob*3+2] += accel_tmp[iglob*3+2];
  }
}


/* ----------------------------------------------------------------------------------------------- */
// LTS Newmark update kernels
/* ----------------------------------------------------------------------------------------------- */

/* not used...

__global__ void lts_newmark_update_kernel_slow(realw* displ,
                                               realw* veloc,
                                               realw* accel,
                                               realw* rmassx,
                                               int ilevel,
                                               int step_m,
                                               int NGLOB_AB,
                                               int num_p_level,
                                               realw deltat_lts,
                                               realw* displ_global,
                                               realw* veloc_global) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < NDIM * NGLOB_AB) {
    accel[id] = rmassx[id/3] * accel[id];

    int index = id + NDIM * NGLOB_AB*(ilevel-1);
    int indexm1 = id + NDIM * NGLOB_AB*(ilevel-2);

    if (step_m == 1 && ilevel < num_p_level) {
      veloc[index] = deltat_lts/2.0 * accel[id];
      if (ilevel > 1) veloc[index] += 1.0/deltat_lts * (displ[indexm1] - displ[index]);
    } else {
      // after first step
      veloc[index] += deltat_lts * accel[id];
      if (ilevel > 1) veloc[index] += 2.0/deltat_lts * (displ[indexm1] - displ[index]);
    }

    // finalize step
    displ[index] = displ[index] + deltat_lts * veloc[index];

    // finalize global step and setup next round of fine steps
    if (ilevel == num_p_level) {
      displ_global[id] = displ[index];
      veloc_global[id] = veloc[index];
    }
  }
}
*/

/* ----------------------------------------------------------------------------------------------- */

__global__ void lts_newmark_update_kernel_P(realw* displ,
                                            realw* veloc,
                                            realw* accel,
                                            realw* rmassxyz,
                                            realw* rmassxyz_mod,
                                            realw* cmassxyz,
                                            int ilevel,
                                            int step_m,
                                            int NGLOB_AB,
                                            int NGLOB_PLEVELS_END,
                                            int num_p_level,
                                            realw deltat_lts,
                                            int collect_accel,
                                            realw* accel_collected,
                                            int* mask_ibool_collected) {

  int iglob = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // P node updates (this level and finer)
  if (iglob < NGLOB_PLEVELS_END) {
    // acceleration update
    if (ilevel < num_p_level) {
      accel[iglob*3]   = rmassxyz[iglob*3]   * accel[iglob*3];
      accel[iglob*3+1] = rmassxyz[iglob*3+1] * accel[iglob*3+1];
      accel[iglob*3+2] = rmassxyz[iglob*3+2] * accel[iglob*3+2];
    } else {
      accel[iglob*3]   = rmassxyz_mod[iglob*3]   * accel[iglob*3];
      accel[iglob*3+1] = rmassxyz_mod[iglob*3+1] * accel[iglob*3+1];
      accel[iglob*3+2] = rmassxyz_mod[iglob*3+2] * accel[iglob*3+2];
    }

    // velocity update
    int index = iglob*3 + NDIM * NGLOB_AB * (ilevel-1);
    int indexm1 = iglob*3 + NDIM * NGLOB_AB * (ilevel-2);

    if (step_m == 1 && ilevel < num_p_level) {
      // 1st local time step for finer p-levels
      realw fac = 0.5f * deltat_lts;
      veloc[index]   = fac * accel[iglob*3];
      veloc[index+1] = fac * accel[iglob*3+1];
      veloc[index+2] = fac * accel[iglob*3+2];
      if (ilevel > 1) {
        realw deltat_inv = 1.0f/deltat_lts;
        veloc[index]   += deltat_inv * (displ[indexm1]   - displ[index]);
        veloc[index+1] += deltat_inv * (displ[indexm1+1] - displ[index+1]);
        veloc[index+2] += deltat_inv * (displ[indexm1+2] - displ[index+2]);
      }
    } else {
      // after first step
      veloc[index]   += deltat_lts * accel[iglob*3];
      veloc[index+1] += deltat_lts * accel[iglob*3+1];
      veloc[index+2] += deltat_lts * accel[iglob*3+2];
      if (ilevel < num_p_level) {
        // finer p-levels
        if (ilevel > 1) {
          realw fac2 = 2.0f/deltat_lts;
          veloc[index]   += fac2 * (displ[indexm1]   - displ[index]);
          veloc[index+1] += fac2 * (displ[indexm1+1] - displ[index+1]);
          veloc[index+2] += fac2 * (displ[indexm1+2] - displ[index+2]);
        }
      } else {
        // coarsest p-level
        // includes absorbing boundary term
        realw fac2 = 2.0f/deltat_lts;
        veloc[index]   += 1.0/(1.0 + rmassxyz[iglob*3]   * cmassxyz[iglob*3])   * fac2 * (displ[indexm1]   - displ[index]);
        veloc[index+1] += 1.0/(1.0 + rmassxyz[iglob*3+1] * cmassxyz[iglob*3+1]) * fac2 * (displ[indexm1+1] - displ[index+1]);
        veloc[index+2] += 1.0/(1.0 + rmassxyz[iglob*3+2] * cmassxyz[iglob*3+2]) * fac2 * (displ[indexm1+2] - displ[index+2]);
      }
    }

    // displacement update
    displ[index]   = displ[index]   + deltat_lts * veloc[index];
    displ[index+1] = displ[index+1] + deltat_lts * veloc[index+1];
    displ[index+2] = displ[index+2] + deltat_lts * veloc[index+2];

    // collects accel wavefield for (possible) seismogram outputs or shakemaps
    if (collect_accel){
      if (! mask_ibool_collected[iglob]){
        if (ilevel < num_p_level){
          // uses same mass matrix as on coarsest level
          accel_collected[iglob*3]   = rmassxyz_mod[iglob*3]   / rmassxyz[iglob*3]   * accel[iglob*3];
          accel_collected[iglob*3+1] = rmassxyz_mod[iglob*3+1] / rmassxyz[iglob*3+1] * accel[iglob*3+1];
          accel_collected[iglob*3+2] = rmassxyz_mod[iglob*3+2] / rmassxyz[iglob*3+2] * accel[iglob*3+2];
        } else {
          // coarsest level (ilevel == num_p_level)
          accel_collected[iglob*3]   = accel[iglob*3];
          accel_collected[iglob*3+1] = accel[iglob*3+1];
          accel_collected[iglob*3+2] = accel[iglob*3+2];
        }
        mask_ibool_collected[iglob] = 1; // mask
      }
    }  // collect_accel
  }
}

/* ----------------------------------------------------------------------------------------------- */

// updates to this level only
__global__ void lts_newmark_update_kernel_P2(realw* displ,
                                             realw* veloc,
                                             realw* accel,
                                             realw* rmassxyz,
                                             realw* rmassxyz_mod,
                                             int ilevel,
                                             int step_m,
                                             int NGLOB_AB,
                                             int NGLOB_PLEVELS_START,
                                             int NGLOB_PLEVELS,
                                             int num_p_level,
                                             realw deltat_lts,
                                             int collect_accel,
                                             realw* accel_collected,
                                             int* mask_ibool_collected) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // P node updates (this level only)
  if (id < NGLOB_PLEVELS) {
    int iglob = id + NGLOB_PLEVELS_START - 1;

    // acceleration update
    if (ilevel < num_p_level) {
      accel[iglob*3]   = rmassxyz[iglob*3]   * accel[iglob*3];
      accel[iglob*3+1] = rmassxyz[iglob*3+1] * accel[iglob*3+1];
      accel[iglob*3+2] = rmassxyz[iglob*3+2] * accel[iglob*3+2];
    } else {
      // includes absorbing boundary term in mass matrix
      accel[iglob*3]   = rmassxyz_mod[iglob*3]   * accel[iglob*3];
      accel[iglob*3+1] = rmassxyz_mod[iglob*3+1] * accel[iglob*3+1];
      accel[iglob*3+2] = rmassxyz_mod[iglob*3+2] * accel[iglob*3+2];
    }

    // velocity update
    int index = iglob*3 + NDIM * NGLOB_AB * (ilevel-1);
    int indexm1 = iglob*3 + NDIM * NGLOB_AB * (ilevel-2);

    if (step_m == 1 && ilevel < num_p_level) {
      // 1st local time step for finer p-levels
      realw fac = 0.5f * deltat_lts;
      veloc[index]   = fac * accel[iglob*3];
      veloc[index+1] = fac * accel[iglob*3+1];
      veloc[index+2] = fac * accel[iglob*3+2];
      if (ilevel > 1) {
        realw deltat_inv = 1.0f/deltat_lts;
        veloc[index]   += deltat_inv * displ[indexm1];
        veloc[index+1] += deltat_inv * displ[indexm1+1];
        veloc[index+2] += deltat_inv * displ[indexm1+2];
      }
    } else {
      // after first step
      veloc[index]   += deltat_lts * accel[iglob*3];
      veloc[index+1] += deltat_lts * accel[iglob*3+1];
      veloc[index+2] += deltat_lts * accel[iglob*3+2];
      if (ilevel > 1) {
        realw fac2 = 2.0f/deltat_lts;
        veloc[index]   += fac2 * displ[indexm1];
        veloc[index+1] += fac2 * displ[indexm1+1];
        veloc[index+2] += fac2 * displ[indexm1+2];
      }
    }

    // displacement update
    displ[index]   = displ[index]   + deltat_lts * veloc[index];
    displ[index+1] = displ[index+1] + deltat_lts * veloc[index+1];
    displ[index+2] = displ[index+2] + deltat_lts * veloc[index+2];

    // collects accel wavefield for (possible) seismogram outputs or shakemaps
    if (collect_accel){
      if (! mask_ibool_collected[iglob]){
        if (ilevel < num_p_level){
          // uses same mass matrix as on coarsest level
          accel_collected[iglob*3]   = rmassxyz_mod[iglob*3]   / rmassxyz[iglob*3]   * accel[iglob*3];
          accel_collected[iglob*3+1] = rmassxyz_mod[iglob*3+1] / rmassxyz[iglob*3+1] * accel[iglob*3+1];
          accel_collected[iglob*3+2] = rmassxyz_mod[iglob*3+2] / rmassxyz[iglob*3+2] * accel[iglob*3+2];
        } else {
          // coarsest level (ilevel == num_p_level)
          accel_collected[iglob*3]   = accel[iglob*3];
          accel_collected[iglob*3+1] = accel[iglob*3+1];
          accel_collected[iglob*3+2] = accel[iglob*3+2];
        }
        mask_ibool_collected[iglob] = 1; // mask
      }
    } // collect_accel
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void lts_newmark_update_kernel_R(realw* displ,
                                            realw* veloc,
                                            realw* accel,
                                            realw* rmassxyz,
                                            realw* rmassxyz_mod,
                                            int ilevel,
                                            int step_m,
                                            int NGLOB_AB,
                                            int num_p_level,
                                            realw deltat_lts,
                                            int num_p_level_coarser_to_update,
                                            int* p_level_coarser_to_update,
                                            int* iglob_p_refine,
                                            int collect_accel,
                                            realw* accel_collected,
                                            int* mask_ibool_collected) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // R node updates (this level and finer)
  //
  // this kernel lts_newmark_update_kernel_R() only gets called for ilevel < num_p_level
  if (id < num_p_level_coarser_to_update) {
    int iglob = p_level_coarser_to_update[id] - 1;

    // acceleration update
    //
    // if coarser node is p==1, then need to use modified mass matrix
    // (for absorbing boundaries)
    // if (iglob*3 + 2 == 82122*3+2) printf("accelR[%d]=%g,rmassxyz[idx]=%g\n",3*128881+1,accel[iglob*3 + 2],rmassxyz[iglob*3 + 2]);
    if (iglob_p_refine[iglob] > 1) {
      accel[iglob*3]   = rmassxyz[iglob*3]   * accel[iglob*3];
      accel[iglob*3+1] = rmassxyz[iglob*3+1] * accel[iglob*3+1];
      accel[iglob*3+2] = rmassxyz[iglob*3+2] * accel[iglob*3+2];
    } else {
      accel[iglob*3]   = rmassxyz_mod[iglob*3]   * accel[iglob*3];
      accel[iglob*3+1] = rmassxyz_mod[iglob*3+1] * accel[iglob*3+1];
      accel[iglob*3+2] = rmassxyz_mod[iglob*3+2] * accel[iglob*3+2];
    }

    int index = iglob*3 + NDIM * NGLOB_AB * (ilevel-1);
    int indexm1 = iglob*3 + NDIM * NGLOB_AB * (ilevel-2);

    // velocity update
    if (step_m == 1 && ilevel < num_p_level) {
      realw fac = 0.5f * deltat_lts;
      veloc[index]   = fac * accel[iglob*3];
      veloc[index+1] = fac * accel[iglob*3+1];
      veloc[index+2] = fac * accel[iglob*3+2];
      if (ilevel > 1) {
        realw deltat_inv = 1.0f/deltat_lts;
        veloc[index]   += deltat_inv * (displ[indexm1]);
        veloc[index+1] += deltat_inv * (displ[indexm1+1]);
        veloc[index+2] += deltat_inv * (displ[indexm1+2]);
      }
    } else {
      // after first step
      veloc[index]   += deltat_lts * accel[iglob*3];
      veloc[index+1] += deltat_lts * accel[iglob*3+1];
      veloc[index+2] += deltat_lts * accel[iglob*3+2];
      if (ilevel > 1) {
        realw fac2 = 2.0f/deltat_lts;
        veloc[index]   += fac2 * (displ[indexm1]);
        veloc[index+1] += fac2 * (displ[indexm1+1]);
        veloc[index+2] += fac2 * (displ[indexm1+2]);
      }
    }

    // displacement update
    displ[index]   = displ[index]   + deltat_lts * veloc[index];
    displ[index+1] = displ[index+1] + deltat_lts * veloc[index+1];
    displ[index+2] = displ[index+2] + deltat_lts * veloc[index+2];

    // collects accel wavefield for (possible) seismogram outputs or shakemaps
    if (collect_accel){
      if (! mask_ibool_collected[iglob]){
        if (iglob_p_refine[iglob] > 1) {
          // node belongs to finer p-levels
          // uses same mass matrix as on coarsest level
          accel_collected[iglob*3]   = rmassxyz_mod[iglob*3]   / rmassxyz[iglob*3]   * accel[iglob*3];
          accel_collected[iglob*3+1] = rmassxyz_mod[iglob*3+1] / rmassxyz[iglob*3+1] * accel[iglob*3+1];
          accel_collected[iglob*3+2] = rmassxyz_mod[iglob*3+2] / rmassxyz[iglob*3+2] * accel[iglob*3+2];
        } else {
          // node belongs to coarsest level
          accel_collected[iglob*3]   = accel[iglob*3];
          accel_collected[iglob*3+1] = accel[iglob*3+1];
          accel_collected[iglob*3+2] = accel[iglob*3+2];
        }
        mask_ibool_collected[iglob] = 1; // mask
      }
    }  // collect_accel
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void lts_newmark_update_veloc_kernel(realw* veloc,
                                                realw* veloc_global,
                                                realw* accel_global,
                                                int nglob,
                                                int num_p_level,
                                                int use_accel_collected,
                                                realw* accel_collected,
                                                int* mask_ibool_collected) {

  int iglob = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  // size == NGLOB_AB
  if (iglob < nglob) {
    // veloc_p(NDIM,NGLOB,num_p_level) indexing
    int index = iglob*3 + NDIM * nglob * (num_p_level-1);

    // update velocity with last p-level (ilevel == num_p_level) velocity values
    veloc_global[iglob*3]   = veloc[index];
    veloc_global[iglob*3+1] = veloc[index+1];
    veloc_global[iglob*3+2] = veloc[index+2];

    // sets accel to collected accel wavefield for (possible) seismogram outputs or shakemaps
    if (use_accel_collected){
      accel_global[iglob*3]   = accel_collected[iglob*3];
      accel_global[iglob*3+1] = accel_collected[iglob*3+1];
      accel_global[iglob*3+2] = accel_collected[iglob*3+2];
      // re-sets flags on global nodes (for next lts stepping)
      mask_ibool_collected[iglob] = 0;
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void lts_newmark_update_displ_kernel(realw* displ,
                                                realw* displ_global,
                                                int nglob,
                                                int num_p_level) {

  int iglob = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  // size == NDIM * NGLOB_AB
  if (iglob < nglob) {
    // displ_p(NDIM,NGLOB,num_p_level) indexing
    int index = iglob*3 + NDIM * nglob * (num_p_level-1);

    // update displacement: time step displ_global(n) -> displ_global(n+1)
    //                      see comment in lts_newmark_update_displ() routine about updating displ/veloc/accel
    displ_global[iglob*3]   = displ[index];
    displ_global[iglob*3+1] = displ[index+1];
    displ_global[iglob*3+2] = displ[index+2];
  }
}


/* ----------------------------------------------------------------------------------------------- */
// zeroing kernels
/* ----------------------------------------------------------------------------------------------- */

__global__ void zero_accel_kernel_P(realw* accel,
                                    int NGLOB_PLEVELS_END) {

  int iglob = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // get nodes from the coarser level
  if (iglob < NGLOB_PLEVELS_END) {
    accel[iglob*3] = 0.f;
    accel[iglob*3+1] = 0.f;
    accel[iglob*3+2] = 0.f;
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void zero_accel_kernel_R(realw* accel,
                                    int num_p_level_coarser_to_update,
                                    int* p_level_coarser_to_update) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // P node updates (this level and finer)
  if (id < num_p_level_coarser_to_update) {
    int iglob = p_level_coarser_to_update[id] - 1;

    accel[iglob*3] = 0.f;
    accel[iglob*3+1] = 0.f;
    accel[iglob*3+2] = 0.f;
  }
}


/* ----------------------------------------------------------------------------------------------- */
// initial conditions kernels
/* ----------------------------------------------------------------------------------------------- */

__global__ void set_finer_initial_condition_fast_kernel_P(realw* displ_p,
                                                          int ilevel,
                                                          int NGLOB_AB,
                                                          int NGLOB_PLEVELS_END) {

  int iglob = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // get nodes from the coarser level
  if (iglob < NGLOB_PLEVELS_END) {
    int index = iglob*3 + NDIM * NGLOB_AB * (ilevel-1);
    int indexp1 = iglob*3 + NDIM * NGLOB_AB * (ilevel);

    displ_p[index]   = displ_p[indexp1];
    displ_p[index+1] = displ_p[indexp1+1];
    displ_p[index+2] = displ_p[indexp1+2];
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void set_finer_initial_condition_fast_kernel_R(realw* displ_p,
                                                          int ilevel,
                                                          int NGLOB_AB,
                                                          int num_p_level_coarser_to_update,
                                                          int* p_level_coarser_to_update) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  // P node updates (this level and finer)
  if (id < num_p_level_coarser_to_update) {
    int iglob = p_level_coarser_to_update[id]-1;
    int index = NDIM * iglob + NDIM * NGLOB_AB * (ilevel-1);

    displ_p[index] = 0.f;
    displ_p[index+1] = 0.f;
    displ_p[index+2] = 0.f;
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void set_finer_initial_condition_kernel(realw* displ_p,
                                                   int ilevel,
                                                   int NGLOB_AB) {

  int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if (id < NDIM * NGLOB_AB) {
    int index = id + NDIM * NGLOB_AB * (ilevel-1);
    int indexp1 = id + NDIM * NGLOB_AB * (ilevel);

    displ_p[index] = displ_p[indexp1];
  }
}
