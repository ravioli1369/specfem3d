!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

!! Solving the wavefield discontinuity problem with a non-split-node
!! scheme
!! Tianshi Liu, 2023.5
module wavefield_discontinuity_solver

  use constants, only: CUSTOM_REAL

  implicit none

  !! ispec_to_elem_wd(NSPEC_AB)
  !! ispec_to_elem_wd(ispec) = ispec_wd (0 if element not belong to boundary)
  !! read from solver database
  integer, dimension(:), allocatable :: ispec_to_elem_wd

  !! number of distinct GLL points on the boundary
  !! read from solver database
  integer :: nglob_wd

  !! number of elements on the inner side of the boundary
  !! read from solver database
  integer :: nspec_wd

  !! ibool_wd(NGLLX, NGLLY, NGLLZ, nspec_wd)
  !! ibool_wd(i,j,k,ispec_wd) = iglob_wd (0 if point not on boundary)
  !! read from solver database
  integer, dimension(:,:,:,:), allocatable :: ibool_wd

  !! boundary_to_iglob_wd(nglob_wd)
  !! boundary_to_iglob_wd(iglob_wd) = iglob
  !! read from solver database
  integer, dimension(:), allocatable :: boundary_to_iglob_wd

  !! mass_in_wd(nglob_wd)
  !! mass matrix on the inner side of the boundary
  !! note that it is not assembled over processors
  !! read from solver database
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: mass_in_wd

  !! number of faces on the boundary
  !! read from solver database
  integer :: nfaces_wd

  !! face_ijk_wd(NDIM, NGLLSQUARE, nfaces_wd)
  !! read from solver database
  integer, dimension(:,:,:), allocatable :: face_ijk_wd

  !! face_ispec_wd(nfaces_wd)
  !! read from solver database
  integer, dimension(:), allocatable :: face_ispec_wd

  !! face_normal_wd(NDIM, NGLLSQUARE, nfaces_wd)
  !! read from solver database
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: face_normal_wd

  !! face_jacobian2Dw_wd(NGLLSQUARE, nfaces_wd)
  !! read from solver database
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: face_jacobian2dw_wd

  !! displ_wd(NDIM, nglob_wd)
  !! displacement discontinuity condition at current time step
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: displ_wd

  !! accel_wd(NDIM, nglob_wd)
  !! acceleration discontinuity condition at current time step
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_wd

  !! traction_wd(NDIM, NGLLSQUARE, nfaces_wd)
  !! traction discontinuity condition at current time step
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: traction_wd

contains

  subroutine prepare_timerun_wavefield_discontinuity()
  use specfem_par, only: IS_WAVEFIELD_DISCONTINUITY
  implicit none
  if (IS_WAVEFIELD_DISCONTINUITY) then
    call read_mesh_databases_wavefield_discontinuity()
    call open_wavefield_discontinuity_file()
  endif
  end subroutine prepare_timerun_wavefield_discontinuity

  subroutine read_mesh_databases_wavefield_discontinuity()
  use constants, only: IFILE_WAVEFIELD_DISCONTINUITY, &
                       FNAME_WAVEFIELD_DISCONTINUITY_DATABASE
  use specfem_par, only: prname, &
                         NSPEC_AB, NGLLX, NGLLY, NGLLZ, NDIM, NGLLSQUARE
  implicit none
  open(unit=IFILE_WAVEFIELD_DISCONTINUITY, &
       file = trim(prname)//trim(FNAME_WAVEFIELD_DISCONTINUITY_DATABASE), &
       action='read', form='unformatted')

  allocate(ispec_to_elem_wd(NSPEC_AB))
  read(IFILE_WAVEFIELD_DISCONTINUITY) ispec_to_elem_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) nglob_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) nspec_wd
  allocate(ibool_wd(NGLLX, NGLLY, NGLLZ, nspec_wd), &
           boundary_to_iglob_wd(nglob_wd), &
           mass_in_wd(nglob_wd))
  read(IFILE_WAVEFIELD_DISCONTINUITY) ibool_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) boundary_to_iglob_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) mass_in_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) nfaces_wd
  allocate(face_ijk_wd(NDIM, NGLLSQUARE, nfaces_wd), &
           face_ispec_wd(nfaces_wd), &
           face_normal_wd(NDIM, NGLLSQUARE, nfaces_wd), &
           face_jacobian2Dw_wd(NGLLSQUARE, nfaces_wd))
  read(IFILE_WAVEFIELD_DISCONTINUITY) face_ijk_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) face_ispec_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) face_normal_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) face_jacobian2Dw_wd
  allocate(displ_wd(NDIM, nglob_wd), &
           accel_wd(NDIM, nglob_wd), &
           traction_wd(NDIM, NGLLSQUARE, nfaces_wd))
  end subroutine read_mesh_databases_wavefield_discontinuity

  subroutine open_wavefield_discontinuity_file()
  use specfem_par, only: prname
  use constants, only: IFILE_WAVEFIELD_DISCONTINUITY
  implicit none
  open(unit=IFILE_WAVEFIELD_DISCONTINUITY, &
       file=trim(prname)//'wavefield_discontinuity.bin', &
       status='old',action='read',form='unformatted')
  end subroutine open_wavefield_discontinuity_file

  subroutine read_wavefield_discontinuity_file()
  use constants, only: IFILE_WAVEFIELD_DISCONTINUITY
  implicit none
  read(IFILE_WAVEFIELD_DISCONTINUITY) displ_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) accel_wd
  read(IFILE_WAVEFIELD_DISCONTINUITY) traction_wd
  end subroutine read_wavefield_discontinuity_file

  subroutine transfer_wavefield_discontinuity_to_GPU()
  use specfem_par, only: Mesh_pointer, NDIM, NGLLSQUARE
  implicit none
  call transfer_wavefield_discontinuity_to_device(nglob_wd*NDIM, &
                                                  nfaces_wd*NDIM*NGLLSQUARE, &
                                                  displ_wd, accel_wd, &
                                                  traction_wd, Mesh_pointer)
  end subroutine transfer_wavefield_discontinuity_to_GPU

  subroutine prepare_wavefield_discontinuity_GPU()
  use specfem_par, only: Mesh_pointer
  implicit none
  call prepare_wavefield_discontinuity_device(Mesh_pointer, ispec_to_elem_wd, &
                                                nglob_wd, nspec_wd, ibool_wd, &
                                                boundary_to_iglob_wd, &
                                                mass_in_wd, &
                                                nfaces_wd, face_ijk_wd, &
                                                face_ispec_wd, face_normal_wd, &
                                                face_jacobian2dw_wd)
  end subroutine prepare_wavefield_discontinuity_GPU

  subroutine finalize_wavefield_discontinuity()
  use constants, only: IFILE_WAVEFIELD_DISCONTINUITY
  implicit none
  close(IFILE_WAVEFIELD_DISCONTINUITY)
  deallocate(ispec_to_elem_wd, ibool_wd, boundary_to_iglob_wd, mass_in_wd, &
             face_ijk_wd, face_ispec_wd, face_normal_wd, face_jacobian2Dw_wd, &
             displ_wd, accel_wd, traction_wd)
  end subroutine finalize_wavefield_discontinuity

  subroutine add_displacement_discontinuity_element(ispec, dummyx_loc, &
                                                  dummyy_loc, dummyz_loc)
  use specfem_par, only: CUSTOM_REAL,NGLLX, NGLLY, NGLLZ
  implicit none
  integer, intent(in) :: ispec
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ) :: dummyx_loc, &
                        dummyy_loc, dummyz_loc
  integer :: ispec_wd, i, j, k, iglob_wd
  ispec_wd = ispec_to_elem_wd(ispec)
  if (ispec_wd /= 0) then
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob_wd = ibool_wd(i,j,k,ispec_wd)
          if (iglob_wd /= 0) then
            dummyx_loc(i,j,k) = dummyx_loc(i,j,k) + displ_wd(1, iglob_wd)
            dummyy_loc(i,j,k) = dummyy_loc(i,j,k) + displ_wd(2, iglob_wd)
            dummyz_loc(i,j,k) = dummyz_loc(i,j,k) + displ_wd(3, iglob_wd)
          endif
        enddo
      enddo
    enddo
  endif
  end subroutine add_displacement_discontinuity_element

  subroutine add_traction_discontinuity(accel, nglob)
  use specfem_par, only: CUSTOM_REAL, NGLLSQUARE, &
                         ibool, NDIM
  !use specfem_par_elastic, only: accel
  implicit none
  integer :: iglob_wd, iglob, ispec, i, j, k, iface_wd, igll, nglob
  real(kind=CUSTOM_REAL) :: jacobianw
  real(kind=CUSTOM_REAL) :: accel(NDIM, nglob)
  do iglob_wd = 1, nglob_wd
    iglob = boundary_to_iglob_wd(iglob_wd)
    accel(:,iglob) = accel(:,iglob) - &
                     accel_wd(:,iglob_wd) * mass_in_wd(iglob_wd)
  enddo
  do iface_wd = 1, nfaces_wd
    do igll = 1, NGLLSQUARE
      i = face_ijk_wd(1, igll, iface_wd)
      j = face_ijk_wd(2, igll, iface_wd)
      k = face_ijk_wd(3, igll, iface_wd)
      ispec = face_ispec_wd(iface_wd)
      iglob = ibool(i,j,k,ispec)
      jacobianw = face_jacobian2Dw_wd(igll, iface_wd)
      accel(:,iglob) = accel(:,iglob) + &
                     traction_wd(:,igll,iface_wd) * jacobianw
    enddo
  enddo
  end subroutine add_traction_discontinuity

  subroutine add_traction_discontinuity_GPU()
  use specfem_par, only: Mesh_pointer
  implicit none
  call wavefield_discontinuity_add_traction_cuda(nglob_wd, nfaces_wd, &
                                                 Mesh_pointer)
  end subroutine add_traction_discontinuity_GPU
end module wavefield_discontinuity_solver
