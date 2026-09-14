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


__global__ void compute_add_sources_kernel(realw* accel,
                                           int* d_ibool,
                                           realw* sourcearrays,
                                           field* stf_pre_compute,
                                           int myrank,
                                           int* islice_selected_source,
                                           int* ispec_selected_source,
                                           int* ispec_is_elastic,
                                           int NSOURCES) {
  int i = threadIdx.x;
  int j = threadIdx.y;
  int k = threadIdx.z;

  int isource  = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int ispec,iglob;

  if (isource < NSOURCES) { // when NSOURCES > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.

    if (myrank == islice_selected_source[isource]) {

      ispec = ispec_selected_source[isource]-1;

      if (ispec_is_elastic[ispec]) {

        iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

        realw stf = (realw) stf_pre_compute[isource];

        realw stf_x = sourcearrays[INDEX5(NSOURCES,NDIM,NGLLX,NGLLX,isource,0,i,j,k)] * stf;
        realw stf_y = sourcearrays[INDEX5(NSOURCES,NDIM,NGLLX,NGLLX,isource,1,i,j,k)] * stf;
        realw stf_z = sourcearrays[INDEX5(NSOURCES,NDIM,NGLLX,NGLLX,isource,2,i,j,k)] * stf;

        atomicAdd(&accel[iglob*3],stf_x);
        atomicAdd(&accel[iglob*3+1],stf_y);
        atomicAdd(&accel[iglob*3+2],stf_z);
      }
    }
  }

}


