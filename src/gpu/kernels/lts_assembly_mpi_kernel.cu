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



/* ----------------------------------------------------------------------------------------------- */
// prepare/assemble kernels
/* ----------------------------------------------------------------------------------------------- */

__global__ void prepare_reduced_boundary_lts_accel_on_device(realw* d_accel,
                                                             realw* d_send_accel_buffer,
                                                             int num_interface_p_refine_boundary,
                                                             int* interface_p_refine_boundary,
                                                             int max_nibool_interfaces_boundary,
                                                             int ilevel) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // go through list of nodes on p-level interface
  if(id < num_interface_p_refine_boundary) {
    int iglob = interface_p_refine_boundary[id + (ilevel-1)*max_nibool_interfaces_boundary]-1;

    d_send_accel_buffer[3*id] = d_accel[3*iglob];
    d_send_accel_buffer[3*id+1] = d_accel[3*iglob+1];
    d_send_accel_buffer[3*id+2] = d_accel[3*iglob+2];
  }

}

/* ----------------------------------------------------------------------------------------------- */

__global__ void assemble_reduced_boundary_lts_accel_on_device(realw* d_accel,
                                                              realw* d_send_accel_buffer,
                                                              int num_interface_p_refine_boundary,
                                                              int* interface_p_refine_boundary,
                                                              int max_nibool_interfaces_boundary,
                                                              int ilevel) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // go through list of nodes on p-level interface
  if(id < num_interface_p_refine_boundary) {
    int iglob = interface_p_refine_boundary[id + (ilevel-1)*max_nibool_interfaces_boundary]-1;

    atomicAdd(&d_accel[3*iglob],d_send_accel_buffer[3*id]);
    atomicAdd(&d_accel[3*iglob+1],d_send_accel_buffer[3*id+1]);
    atomicAdd(&d_accel[3*iglob+2],d_send_accel_buffer[3*id+2]);
  }
}

/* ----------------------------------------------------------------------------------------------- */


__global__ void assemble_boundary_lts_accel_on_device(realw* d_accel,
                                                      realw* d_send_accel_buffer,
                                                      int num_interfaces_ext_mesh,
                                                      int* num_interface_p_refine_ibool,
                                                      int* interface_p_refine_ibool,
                                                      int max_nibool_interfaces_ext_mesh,
                                                      int ilevel) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  for( int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    int num_interface_iglob = num_interface_p_refine_ibool[iinterface + num_interfaces_ext_mesh*(ilevel-1)];
    if(id < num_interface_iglob) {
      int iglob = interface_p_refine_ibool[id + max_nibool_interfaces_ext_mesh*(iinterface + num_interfaces_ext_mesh*(ilevel-1))]-1;
      int iglob_send = id + iinterface*max_nibool_interfaces_ext_mesh;

      atomicAdd(&d_accel[3*iglob],d_send_accel_buffer[3*iglob_send]);
      atomicAdd(&d_accel[3*iglob+1],d_send_accel_buffer[3*iglob_send+1]);
      atomicAdd(&d_accel[3*iglob+2],d_send_accel_buffer[3*iglob_send+2]);
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */


__global__ void prepare_boundary_lts_accel_on_device(realw* d_accel,
                                                     realw* d_send_accel_buffer,
                                                     int num_interfaces_ext_mesh,
                                                     int* num_interface_p_refine_ibool,
                                                     int* interface_p_refine_ibool,
                                                     int max_nibool_interfaces_ext_mesh,
                                                     int ilevel) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  for( int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    int num_interface_iglob = num_interface_p_refine_ibool[iinterface + num_interfaces_ext_mesh*(ilevel-1)];
    if(id < num_interface_iglob) {
      int iglob = interface_p_refine_ibool[id + max_nibool_interfaces_ext_mesh*(iinterface + num_interfaces_ext_mesh*(ilevel-1))]-1;
      int iglob_send = id + iinterface*max_nibool_interfaces_ext_mesh;

      d_send_accel_buffer[3*iglob_send] = d_accel[3*iglob];
      d_send_accel_buffer[3*iglob_send+1] = d_accel[3*iglob+1];
      d_send_accel_buffer[3*iglob_send+2] = d_accel[3*iglob+2];
    }
  }
}
