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

! This module outputs gravity field at required locations
!
! Authors:
! Surendra Somala (Caltech) surendra@caltech.edu - 2013
! with advice from Jan Harms and Pablo Ampuero (ampuero@gps.caltech.edu)

module gravity_perturbation

  use constants
  ! note: instead of using
  !         use constants, only: CUSTOM_REAL
  !       leads to an internal compiler error for gfortran:
  !         internal compiler error: in gfc_typenode_for_spec, at Fortran/trans-types.c:1086
  !       just use without the only-specifier seems to work.

  implicit none

  private

  integer :: nstat,ntimgap,nstat_local
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstat,ystat,zstat
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accE,accN,accZ
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rho0_wm
  ! time-invariant per-node per-station gravity integral weights (precomputed in gravity_init)
  ! w3 = G_const*rho0_wm/Rg**3 ; w5 = 3*G_const*rho0_wm/Rg**5
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: w3,w5

  logical, save :: GRAVITY_SIMULATION = .false.

  public :: gravity_init, gravity_init_device, gravity_timeseries, gravity_output, GRAVITY_SIMULATION

contains

!=====================================================================

  subroutine gravity_init()

  use constants, only: IMAIN,IIN_G,myrank,NGLLX,NGLLY,NGLLZ,GRAV

  use specfem_par, only: NGLOB_AB, NSTEP, NSPEC_AB, &
       rhostore, &
       xstore, ystore, zstore, &
       xigll, yigll, zigll, &
       wxgll, wygll, wzgll, &
       NGNOD, ibool, GPU_MODE, &
       jacobianstore, irregular_element_number, jacobian_regular

  implicit none

  ! local parameters
  ! same gravitational constant (cast to CUSTOM_REAL) as used in gravity_timeseries
  real(kind=CUSTOM_REAL), parameter :: G_const = GRAV
  ! per-node station distance used only while precomputing the time-invariant weights
  real(kind=CUSTOM_REAL) :: Rg
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ):: rho_elem
  ! jacobian at the current GLL point, taken from the mesher's stored values (see loop below)
  real(kind=CUSTOM_REAL) :: jacobianl
  integer :: ispec_irreg
  double precision :: Jac3D
  ! coordinates of the control points
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)
  integer :: ia
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz,iax,iay,iaz
  integer :: nstep_grav
  integer :: i,j,k,iglob,ispec,istat,ier

  ! opens gravity parameter file
  open(unit=IIN_G,file='./DATA/gravity_stations',status='old',iostat=ier)

  ! checks if file exists
  if (ier /= 0) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  no gravity stations'
      call flush_IMAIN()
    endif
    ! nothing to do
    return
  endif

  ! sets gravity flag
  GRAVITY_SIMULATION = .true.

  ! reads number of stations
  read(IIN_G,*) nstat,ntimgap

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  incorporating gravity field simulation'
    write(IMAIN,*) '    gravity stations: ',nstat
    call flush_IMAIN()
  endif

  allocate(xstat(nstat),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2234')
  allocate(ystat(nstat),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2235')
  allocate(zstat(nstat),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2236')
  do istat = 1,nstat
     read(IIN_G,*) xstat(istat),ystat(istat),zstat(istat)
  enddo
  close(IIN_G)

  nstep_grav = floor(dble(NSTEP)/dble(ntimgap))
  allocate(accE(nstep_grav,nstat),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2237')
  allocate(accN(nstep_grav,nstat),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2238')
  allocate(accZ(nstep_grav,nstat),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2239')

  accE = 0._CUSTOM_REAL
  accN = 0._CUSTOM_REAL
  accZ = 0._CUSTOM_REAL

  allocate(rho0_wm(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2240')
  rho0_wm = 0._CUSTOM_REAL

  call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)

  ! define coordinates of the control points of the element
  do ia = 1,NGNOD

    if (iaddx(ia) == 0) then
      iax(ia) = 1
    else if (iaddx(ia) == 1) then
      iax(ia) = (NGLLX+1)/2
    else if (iaddx(ia) == 2) then
      iax(ia) = NGLLX
    else
      call exit_MPI(myrank,'incorrect value of iaddx')
    endif

    if (iaddy(ia) == 0) then
      iay(ia) = 1
    else if (iaddy(ia) == 1) then
      iay(ia) = (NGLLY+1)/2
    else if (iaddy(ia) == 2) then
      iay(ia) = NGLLY
    else
      call exit_MPI(myrank,'incorrect value of iaddy')
    endif

    if (iaddz(ia) == 0) then
      iaz(ia) = 1
    else if (iaddz(ia) == 1) then
      iaz(ia) = (NGLLZ+1)/2
    else if (iaddz(ia) == 2) then
      iaz(ia) = NGLLZ
    else
      call exit_MPI(myrank,'incorrect value of iaddz')
    endif

  enddo

  do ispec = 1,NSPEC_AB

    rho_elem = rhostore(:,:,:,ispec)
    ! rho0_wm(iglob) = sum over elements of rho * jacobian * GLL weight at each node -- exactly
    ! the (un-inverted, pre-Stacey) elastic mass matrix. Reuse the jacobian the mesher already
    ! stored (jacobianstore for deformed/irregular elements, the constant jacobian_regular
    ! otherwise) instead of recomputing it per GLL point via recompute_jacobian_gravity. That
    ! per-node double-precision recompute dominated the ~20 min "preparing gravity" setup on the
    ! 3.4M-element fine box; this mirrors define_mass_matrices_elastic, so rho0_wm is unchanged.
    ispec_irreg = irregular_element_number(ispec)
    if (ispec_irreg == 0) jacobianl = jacobian_regular

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)
          rho0_wm(iglob) = rho0_wm(iglob) &
               + real(dble(jacobianl) * wxgll(i)*wygll(j)*wzgll(k) * dble(rho_elem(i,j,k)), kind=CUSTOM_REAL)
        enddo
      enddo
    enddo

  enddo

  ! Time-invariant per-node per-station weights (w3 = G*rho0_wm/Rg^3, w5 = 3*G*rho0_wm/Rg^5)
  ! for the CPU hot loop in gravity_timeseries. On the GPU path they are NOT built here: the
  ! reduction kernel rebuilds them per output step from rho0_wm + node/station coordinates
  ! (rho0_wm is uploaded in gravity_init_device), which avoids the NGLOB_AB x nstat host and
  ! device storage (~6 GB on the production fine box) at the cost of a few flops on a kernel
  ! that only runs every ntimgap steps.
  if (.not. GPU_MODE) then
    allocate(w3(NGLOB_AB,nstat),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2243')
    allocate(w5(NGLOB_AB,nstat),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2244')
    w3 = 0._CUSTOM_REAL
    w5 = 0._CUSTOM_REAL

    do istat = 1,nstat
      do iglob = 1,NGLOB_AB
        Rg = sqrt((xstore(iglob)-xstat(istat))**2 &
                + (ystore(iglob)-ystat(istat))**2 &
                + (zstore(iglob)-zstat(istat))**2)
        w3(iglob,istat) = G_const*rho0_wm(iglob)/Rg**3
        w5(iglob,istat) = 3._CUSTOM_REAL*G_const*rho0_wm(iglob)/Rg**5
      enddo
    enddo
  endif

  end subroutine gravity_init

!=====================================================================

  subroutine gravity_init_device()

! uploads the time-invariant gravity integral data to the device once, so that
! gravity_timeseries can run the whole per-output-step reduction on the GPU.
! must be called after gravity_init has filled rho0_wm/xstat/... AND after the GPU
! mesh has been set up (Mesh_pointer valid), i.e. from prepare_GPU.

  use specfem_par, only: NGLOB_AB, xstore, ystore, zstore, Mesh_pointer

  implicit none

  real(kind=CUSTOM_REAL) :: G_const

  ! nothing to do if there are no gravity stations
  if (.not. GRAVITY_SIMULATION) return

  ! same gravitational constant (cast to CUSTOM_REAL) that the device kernel folds with rho0_wm
  G_const = real(GRAV, kind=CUSTOM_REAL)

  ! rho0_wm : NGLOB_AB station-independent per-node mass weight; the device reduction kernel
  !           rebuilds w3 = G*rho0_wm/Rg**3 and w5 = 3*G*rho0_wm/Rg**5 per output step.
  ! xstore/ystore/zstore : node coordinates (constant in time)
  ! xstat/ystat/zstat : station coordinates (kept on host inside the device struct)
  call prepare_gravity_device(Mesh_pointer, rho0_wm, G_const, xstore, ystore, zstore, &
                              NGLOB_AB, nstat, xstat, ystat, zstat)

  end subroutine gravity_init_device

!=====================================================================

! recompute 3D jacobian at a given point for a 8-node element : modified from recompute_jacobian

  subroutine recompute_jacobian_gravity(xelm,yelm,zelm,xi,eta,gamma,jacobian)

  use constants, only: NDIM,ONE,ZERO
  use specfem_par, only: NGNOD

  implicit none

  double precision :: x,y,z
  double precision :: xi,eta,gamma,jacobian

! coordinates of the control points
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

! 3D shape functions and their derivatives at receiver
  double precision :: shape3D(NGNOD)
  double precision :: dershape3D(NDIM,NGNOD)

  double precision :: xxi,yxi,zxi
  double precision :: xeta,yeta,zeta
  double precision :: xgamma,ygamma,zgamma
  double precision :: ra1,ra2,rb1,rb2,rc1,rc2

  integer :: ia

! for 8-node element
  double precision, parameter :: ONE_EIGHTH = 0.125d0

! recompute jacobian for any (xi,eta,gamma) point, not necessarily a GLL point

! check that the parameter file is correct
  if (NGNOD /= 8) &
       stop 'elements should have 8  control nodes'

! ***
! *** create the 3D shape functions and the Jacobian for an 8-node element
! ***

!--- case of an 8-node 3D element (Dhatt-Touzot p. 115)

  ra1 = ONE + xi
  ra2 = ONE - xi

  rb1 = ONE + eta
  rb2 = ONE - eta

  rc1 = ONE + gamma
  rc2 = ONE - gamma

  shape3D(1) = ONE_EIGHTH*ra2*rb2*rc2
  shape3D(2) = ONE_EIGHTH*ra1*rb2*rc2
  shape3D(3) = ONE_EIGHTH*ra1*rb1*rc2
  shape3D(4) = ONE_EIGHTH*ra2*rb1*rc2
  shape3D(5) = ONE_EIGHTH*ra2*rb2*rc1
  shape3D(6) = ONE_EIGHTH*ra1*rb2*rc1
  shape3D(7) = ONE_EIGHTH*ra1*rb1*rc1
  shape3D(8) = ONE_EIGHTH*ra2*rb1*rc1

  dershape3D(1,1) = - ONE_EIGHTH*rb2*rc2
  dershape3D(1,2) = ONE_EIGHTH*rb2*rc2
  dershape3D(1,3) = ONE_EIGHTH*rb1*rc2
  dershape3D(1,4) = - ONE_EIGHTH*rb1*rc2
  dershape3D(1,5) = - ONE_EIGHTH*rb2*rc1
  dershape3D(1,6) = ONE_EIGHTH*rb2*rc1
  dershape3D(1,7) = ONE_EIGHTH*rb1*rc1
  dershape3D(1,8) = - ONE_EIGHTH*rb1*rc1

  dershape3D(2,1) = - ONE_EIGHTH*ra2*rc2
  dershape3D(2,2) = - ONE_EIGHTH*ra1*rc2
  dershape3D(2,3) = ONE_EIGHTH*ra1*rc2
  dershape3D(2,4) = ONE_EIGHTH*ra2*rc2
  dershape3D(2,5) = - ONE_EIGHTH*ra2*rc1
  dershape3D(2,6) = - ONE_EIGHTH*ra1*rc1
  dershape3D(2,7) = ONE_EIGHTH*ra1*rc1
  dershape3D(2,8) = ONE_EIGHTH*ra2*rc1

  dershape3D(3,1) = - ONE_EIGHTH*ra2*rb2
  dershape3D(3,2) = - ONE_EIGHTH*ra1*rb2
  dershape3D(3,3) = - ONE_EIGHTH*ra1*rb1
  dershape3D(3,4) = - ONE_EIGHTH*ra2*rb1
  dershape3D(3,5) = ONE_EIGHTH*ra2*rb2
  dershape3D(3,6) = ONE_EIGHTH*ra1*rb2
  dershape3D(3,7) = ONE_EIGHTH*ra1*rb1
  dershape3D(3,8) = ONE_EIGHTH*ra2*rb1


! compute coordinates and jacobian matrix
  x = ZERO
  y = ZERO
  z = ZERO
  xxi = ZERO
  xeta = ZERO
  xgamma = ZERO
  yxi = ZERO
  yeta = ZERO
  ygamma = ZERO
  zxi = ZERO
  zeta = ZERO
  zgamma = ZERO

  do ia = 1,NGNOD
    x = x+shape3D(ia)*xelm(ia)
    y = y+shape3D(ia)*yelm(ia)
    z = z+shape3D(ia)*zelm(ia)

    xxi = xxi+dershape3D(1,ia)*xelm(ia)
    xeta = xeta+dershape3D(2,ia)*xelm(ia)
    xgamma = xgamma+dershape3D(3,ia)*xelm(ia)
    yxi = yxi+dershape3D(1,ia)*yelm(ia)
    yeta = yeta+dershape3D(2,ia)*yelm(ia)
    ygamma = ygamma+dershape3D(3,ia)*yelm(ia)
    zxi = zxi+dershape3D(1,ia)*zelm(ia)
    zeta = zeta+dershape3D(2,ia)*zelm(ia)
    zgamma = zgamma+dershape3D(3,ia)*zelm(ia)
  enddo

  jacobian = xxi*(yeta*zgamma-ygamma*zeta) - xeta*(yxi*zgamma-ygamma*zxi) + xgamma*(yxi*zeta-yeta*zxi)

  if (jacobian <= ZERO) stop '3D Jacobian undefined'


  end subroutine recompute_jacobian_gravity

!=====================================================================

  subroutine gravity_timeseries()

  use specfem_par, only: xstore, ystore, zstore, it, NGLOB_AB, GPU_MODE, Mesh_pointer
  use specfem_par_elastic, only: displ

  implicit none

  ! local parameters
  ! CPU-path node-sized scratch (only allocated on the non-GPU path)
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: accEdV,accNdV,accZdV
  real(kind=CUSTOM_REAL) :: E_local,N_local,Z_local,E_all,N_all,Z_all
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dotP
  ! per-station partial (this-rank) gravity-perturbation sums returned by the device reduction
  real(kind=CUSTOM_REAL), dimension(nstat) :: grav_E,grav_N,grav_Z
  integer :: istat, it_grav, ier

  if (mod(it,ntimgap) == 0) then
    it_grav = nint(dble(it)/dble(ntimgap))

    if (GPU_MODE) then
      ! GPU path: reduce over d_displ on the device using the precomputed w3/w5 weights.
      ! this avoids transferring the full displacement field (NDIM*NGLOB_AB reals) back to
      ! the host each output step; only the nstat per-station scalars come back.
      ! grav_E/grav_N/grav_Z hold this rank's partial sums (double-accumulated on the device,
      ! returned as CUSTOM_REAL) before the MPI reduction below.
      call compute_gravity_cuda(Mesh_pointer, grav_E, grav_N, grav_Z)
      do istat = 1,nstat
        call sum_all_all_cr(grav_E(istat),E_all)
        accE(it_grav,istat) = E_all
        call sum_all_all_cr(grav_N(istat),N_all)
        accN(it_grav,istat) = N_all
        call sum_all_all_cr(grav_Z(istat),Z_all)
        accZ(it_grav,istat) = Z_all
      enddo
    else
      ! CPU path (unchanged from Task 2)
      allocate(dotP(NGLOB_AB),accEdV(NGLOB_AB),accNdV(NGLOB_AB),accZdV(NGLOB_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2242')

      ! hot loop: multiply-add only (no sqrt/pow/divide per time step).
      ! the geometric factors G*rho0_wm/Rg**3 and 3*G*rho0_wm/Rg**5 are precomputed
      ! once in gravity_init as w3(:,istat)/w5(:,istat); only displ changes in time.
      ! the (xstore-xstat) etc. subtractions are kept inline (cheap, and needed for dotP).
      do istat = 1,nstat
        dotP = (xstore-xstat(istat))*displ(1,:)+(ystore-ystat(istat))*displ(2,:)+(zstore-zstat(istat))*displ(3,:)

        accEdV = w3(:,istat)*displ(1,:)-w5(:,istat)*(xstore-xstat(istat))*dotP
        E_local = sum(accEdV(:))
        call sum_all_all_cr(E_local,E_all)
        accE(it_grav,istat) = E_all

        accNdV = w3(:,istat)*displ(2,:)-w5(:,istat)*(ystore-ystat(istat))*dotP
        N_local = sum(accNdV(:))
        call sum_all_all_cr(N_local,N_all)
        accN(it_grav,istat) = N_all

        accZdV = w3(:,istat)*displ(3,:)-w5(:,istat)*(zstore-zstat(istat))*dotP
        Z_local = sum(accZdV(:))
        call sum_all_all_cr(Z_local,Z_all)
        accZ(it_grav,istat) = Z_all
      enddo
    endif
  endif

  end subroutine gravity_timeseries

!=====================================================================

  subroutine gravity_output()

  use constants, only: IOUT,OUTPUT_FILES,myrank
  use specfem_par, only: NPROC,NSTEP,DT

  implicit none

  integer :: isample,istat,nstep_grav
  character(len=MAX_STRING_LEN) :: sisname

  nstep_grav = floor(dble(NSTEP)/dble(ntimgap))
  nstat_local = nint(dble(nstat)/dble(NPROC))

  do istat = 1,nstat
    if (istat < myrank*nstat_local+1 .or. istat > (myrank+1)*nstat_local) cycle
    write(sisname,"(a,I0,a)") trim(OUTPUT_FILES)//'/stat', istat, '.grav'
    open(unit=IOUT,file=sisname,status='replace')
    do isample = 1,nstep_grav
      write(IOUT,*) isample*DT*ntimgap, accE(isample,istat),accN(isample,istat),accZ(isample,istat)
    enddo
    close(IOUT)
  enddo

  if (myrank == 0) then ! left-over stations
    do istat = NPROC*nstat_local,nstat
      write(sisname,"(a,I0,a)") trim(OUTPUT_FILES)//'/stat', istat, '.grav'
      open(unit=IOUT,file=sisname,status='replace')
      do isample = 1,nstep_grav
        write(IOUT,*) isample*DT*ntimgap, accE(isample,istat),accN(isample,istat),accZ(isample,istat)
      enddo
      close(IOUT)
    enddo
  endif

  end subroutine gravity_output

!=====================================================================


end module gravity_perturbation
