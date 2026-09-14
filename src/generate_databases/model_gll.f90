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

!--------------------------------------------------------------------------------------------------
!
! GLL
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!
!--------------------------------------------------------------------------------------------------

  subroutine model_gll(myrank,nspec,LOCAL_PATH)

  use constants, only: NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,IMAIN,MAX_STRING_LEN,IIN

  use generate_databases_par, only: ANISOTROPY, ATTENUATION

  use create_regions_mesh_ext_par, only: rhostore,kappastore,mustore,rho_vp,rho_vs,qkappa_attenuation_store,qmu_attenuation_store, &
                                         c11store,c12store,c13store,c14store,c15store,c16store, &
                                         c22store,c23store,c24store,c25store,c26store, &
                                         c33store,c34store,c35store,c36store, &
                                         c44store,c45store,c46store, &
                                         c55store,c56store, &
                                         c66store

  use shared_parameters, only: ADIOS_FOR_MESH,HDF5_ENABLED

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read
  integer :: ier
  character(len=MAX_STRING_LEN) :: prname_lp,filename

  ! select routine for file i/o format
  if (ADIOS_FOR_MESH) then
    ! ADIOS
    call model_gll_adios(myrank,nspec,LOCAL_PATH)
    ! all done
    return
  else if (HDF5_ENABLED) then
    ! not implemented yet
    stop 'HDF5_ENABLED not supported yet for model_gll() routine, please return without flag...'
  else
    ! default binary
    ! implemented here below, continue
    continue
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     using GLL model from: ',trim(LOCAL_PATH)
  endif

  ! processors name
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)// '/' //'proc',myrank,'_'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! if only vp structure is available (as is often the case in exploration seismology),
  !!! use lines for vp only

  ! density
  allocate(rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 647')
  if (ier /= 0) stop 'error allocating array rho_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: rho.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'rho.bin'
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading rho.bin file'
  endif

  read(IIN) rho_read
  close(IIN)

  ! vp
  allocate(vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 648')
  if (ier /= 0) stop 'error allocating array vp_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: vp.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'vp.bin'
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vp.bin file'
  endif

  read(IIN) vp_read
  close(IIN)

  ! vs
  allocate(vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 649')
  if (ier /= 0) stop 'error allocating array vs_read'

  ! user output
  if (myrank == 0) write(IMAIN,*) '     reading in: vs.bin'

  filename = prname_lp(1:len_trim(prname_lp))//'vs.bin'
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(filename)
    stop 'error reading vs.bin file'
  endif

  read(IIN) vs_read
  close(IIN)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! in cases where density structure is not given
  !!! modify according to your desire

  !  rho_read = 1000.0
  !  where ( mustore > 100.0 )  &
  !           rho_read = (1.6612 * (vp_read / 1000.0)     &
  !                      -0.4720 * (vp_read / 1000.0)**2  &
  !                      +0.0671 * (vp_read / 1000.0)**3  &
  !                      -0.0043 * (vp_read / 1000.0)**4  &
  !                      +0.000106*(vp_read / 1000.0)**5)*1000.0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! in cases where shear wavespeed structure is not given
  !!! modify according to your desire

  !   vs_read = 0.0
  !   where ( mustore > 100.0 )       vs_read = vp_read / sqrt(3.0)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! update arrays that will be saved and used in the solver xspecfem3D
  !!! the following part is neccessary if you uncommented something above

  ! density
  rhostore(:,:,:,:) = rho_read(:,:,:,:)

  ! bulk moduli: kappa = rho * (vp**2 - 4/3 vs**2)
  kappastore(:,:,:,:) = rhostore(:,:,:,:) * ( vp_read(:,:,:,:) * vp_read(:,:,:,:) &
                                              - FOUR_THIRDS * vs_read(:,:,:,:) * vs_read(:,:,:,:) )

  ! shear moduli: mu = rho * vs**2
  mustore(:,:,:,:) = rhostore(:,:,:,:) * vs_read(:,:,:,:) * vs_read(:,:,:,:)

  ! products rho*vp and rho*vs (used to speed up absorbing boundaries contributions)
  rho_vp(:,:,:,:) = rhostore(:,:,:,:) * vp_read(:,:,:,:)
  rho_vs(:,:,:,:) = rhostore(:,:,:,:) * vs_read(:,:,:,:)

  if (ANISOTROPY) then

    ! c11
    if (myrank == 0) write(IMAIN,*) '     reading in: c11.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c11.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c11.bin file'
    endif

    read(IIN) c11store
    close(IIN)

    ! c12
    if (myrank == 0) write(IMAIN,*) '     reading in: c12.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c12.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c12.bin file'
    endif

    read(IIN) c12store
    close(IIN)

    ! c13
    if (myrank == 0) write(IMAIN,*) '     reading in: c13.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c13.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c13.bin file'
    endif

    read(IIN) c13store
    close(IIN)

    ! c14
    if (myrank == 0) write(IMAIN,*) '     reading in: c14.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c14.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c14.bin file'
    endif

    read(IIN) c14store
    close(IIN)

    ! c15
    if (myrank == 0) write(IMAIN,*) '     reading in: c15.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c15.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c15.bin file'
    endif

    read(IIN) c15store
    close(IIN)

    ! c16
    if (myrank == 0) write(IMAIN,*) '     reading in: c16.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c16.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c16.bin file'
    endif

    read(IIN) c16store
    close(IIN)

    ! c22
    if (myrank == 0) write(IMAIN,*) '     reading in: c22.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c22.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c22.bin file'
    endif

    read(IIN) c22store
    close(IIN)

    ! c23
    if (myrank == 0) write(IMAIN,*) '     reading in: c23.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c23.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c23.bin file'
    endif

    read(IIN) c23store
    close(IIN)

    ! c24
    if (myrank == 0) write(IMAIN,*) '     reading in: c24.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c24.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c24.bin file'
    endif

    read(IIN) c24store
    close(IIN)

    ! c25
    if (myrank == 0) write(IMAIN,*) '     reading in: c25.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c25.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c25.bin file'
    endif

    read(IIN) c25store
    close(IIN)

    ! c26
    if (myrank == 0) write(IMAIN,*) '     reading in: c26.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c26.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c26.bin file'
    endif

    read(IIN) c26store
    close(IIN)

    ! c33
    if (myrank == 0) write(IMAIN,*) '     reading in: c33.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c33.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c33.bin file'
    endif

    read(IIN) c33store
    close(IIN)

    ! c34
    if (myrank == 0) write(IMAIN,*) '     reading in: c34.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c34.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c34.bin file'
    endif

    read(IIN) c34store
    close(IIN)

    ! c35
    if (myrank == 0) write(IMAIN,*) '     reading in: c35.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c35.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c35.bin file'
    endif

    read(IIN) c35store
    close(IIN)

    ! c36
    if (myrank == 0) write(IMAIN,*) '     reading in: c36.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c36.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c36.bin file'
    endif

    read(IIN) c36store
    close(IIN)

    ! c44
    if (myrank == 0) write(IMAIN,*) '     reading in: c44.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c44.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c44.bin file'
    endif

    read(IIN) c44store
    close(IIN)

    ! c45
    if (myrank == 0) write(IMAIN,*) '     reading in: c45.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c45.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c45.bin file'
    endif

    read(IIN) c45store
    close(IIN)

    ! c46
    if (myrank == 0) write(IMAIN,*) '     reading in: c46.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c46.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c46.bin file'
    endif

    read(IIN) c46store
    close(IIN)

    ! c55
    if (myrank == 0) write(IMAIN,*) '     reading in: c55.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c55.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c55.bin file'
    endif

    read(IIN) c55store
    close(IIN)

    ! c56
    if (myrank == 0) write(IMAIN,*) '     reading in: c56.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c56.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c56.bin file'
    endif

    read(IIN) c56store
    close(IIN)

    ! c66
    if (myrank == 0) write(IMAIN,*) '     reading in: c66.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'c66.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading c66.bin file'
    endif

    read(IIN) c66store
    close(IIN)

  endif

  ! gets attenuation arrays from files
  if (ATTENUATION) then
    ! shear attenuation
    ! user output
    if (myrank == 0) write(IMAIN,*) '     reading in: qmu.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'qmu.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      stop 'Error reading qmu.bin file'
    endif

    read(IIN) qmu_attenuation_store
    close(IIN)

    ! bulk attenuation
    ! user output
    if (myrank == 0) write(IMAIN,*) '     reading in: qkappa.bin'

    filename = prname_lp(1:len_trim(prname_lp))//'qkappa.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'error opening file: ',trim(filename)
      stop 'error reading qkappa.bin file'
    endif

    read(IIN) qkappa_attenuation_store
    close(IIN)
  endif

  ! free memory
  deallocate(rho_read,vp_read,vs_read)

  end subroutine model_gll

