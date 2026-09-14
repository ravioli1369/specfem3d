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
! VS30 model
!
! USGS VS30 layer velocities
! https://earthquake.usgs.gov/data/vs30/
!
! assigns Vs30 layer velocities to top 30-m layer if the data is found in DATA/interface_vs30.dat (optional)
! created by script ./utils/scripts/run_get_USGS_Vs30.py
!
! derives Vp and density by scaling relationship from Brocher (2005)
!
!--------------------------------------------------------------------------------------------------


  subroutine model_vs30_setup_check()

! checks if file exists and reads Vs30 model interface file

  use constants, only: myrank,IIN,IMAIN,HUGEVAL,CUSTOM_REAL,MAX_STRING_LEN,IN_DATA_FILES

  use create_regions_mesh_ext_par, only: USE_MODEL_LAYER_VS30, &
    interface_model_vs30,npx_interface_vs30,npy_interface_vs30, &
    orig_x_interface_vs30,orig_y_interface_vs30, &
    spacing_x_interface_vs30,spacing_y_interface_vs30, &
    SUPPRESS_UTM_PROJECTION_VS30

  implicit none

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename
  logical :: file_exists

  ! file reading
  logical :: SUPPRESS_UTM_PROJECTION
  integer :: ix,iy,ier
  integer :: npx_interface,npy_interface
  double precision :: orig_x_interface,orig_y_interface
  double precision :: spacing_x_interface,spacing_y_interface
  character(len=MAX_STRING_LEN) :: string_read
  character(len=12) :: str_unit
  character(len=12),parameter :: str_unit_m = "(m)",str_unit_deg = "(deg)"
  ! conversion factor:
  ! at equator earth circumference 40,075,161.2 divided by 360.0 degree
  double precision, parameter :: DEGREE_TO_METERS = 111319.8922222222d0

  ! model values
  real(kind=CUSTOM_REAL) :: vs,min_vs,max_vs

  ! initializes
  USE_MODEL_LAYER_VS30 = .false.

  ! checks if USGS Vs30 model interface is provided in folder USGS_VS30/
  !
  ! the script ./utils/scripts/run_get_USGS_Vs30.py can create an interface file for a given region
  ! and stores it in folder USGS_VS30/ as 'interface_vs30.dat'
  !
  filename = trim(IN_DATA_FILES) // 'interface_vs30.dat'

  inquire(file=trim(filename),exist=file_exists)

  ! checks if anything to do
  if (.not. file_exists) return

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) 'Found Vs30 interface file ',trim(filename)
    write(IMAIN,*)
    write(IMAIN,*) '     Reading model parameters from Vs30 interface file ...'
    call flush_IMAIN()
  endif

  ! open file
  open(IIN,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    print *
    print *,'For some reason, inquiring this file was positive, however the file cannot be opened.'
    print *,'Please check your file path and if this file is readable.'
    print *
    print *,'will continue without the Vs30 layer...'
    ! nothing to do
    return
  endif

  ! read in file
  ! this is similar as for (topographic) interface files in src/meshfem3D/create_interfaces_mesh.f90

  ! reads first header line
  call read_next_line(IIN,string_read)

  ! interface definition
  ! # SUPPRESS_UTM_PROJECTION  NXI  NETA  LONG_MIN  LAT_MIN  SPACING_XI  SPACING_ETA
  read(string_read,*,iostat=ier) SUPPRESS_UTM_PROJECTION,npx_interface,npy_interface, &
                                 orig_x_interface,orig_y_interface, &
                                 spacing_x_interface,spacing_y_interface

  ! user output
  if (myrank == 0) then
    if (SUPPRESS_UTM_PROJECTION) then
      str_unit = str_unit_m
    else
      str_unit = str_unit_deg
    endif
    write(IMAIN,*) '       Vs30 interface file   : ',trim(filename)
    write(IMAIN,*)
    write(IMAIN,*) '       number of points x/y = ',npx_interface,npy_interface
    write(IMAIN,*) '       origin x/y     ',trim(str_unit),' = ', sngl(orig_x_interface),sngl(orig_y_interface)
    write(IMAIN,*) '       spacing x/y    ',trim(str_unit),' = ',sngl(spacing_x_interface),sngl(spacing_y_interface)
    if (.not. SUPPRESS_UTM_PROJECTION) then
      write(IMAIN,*) '                       ',trim(str_unit_m),' = ', &
           sngl(spacing_x_interface*DEGREE_TO_METERS),sngl(spacing_y_interface*DEGREE_TO_METERS)
    endif
    write(IMAIN,*)
    write(IMAIN,*) '       dimension x-direction ',trim(str_unit),' = ',sngl(orig_x_interface),'/', &
                          sngl(orig_x_interface + (npx_interface-1)*spacing_x_interface)
    write(IMAIN,*) '       dimension y-direction ',trim(str_unit),' = ',sngl(orig_y_interface),'/', &
                          sngl(orig_y_interface + (npy_interface-1)*spacing_y_interface)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! stores interface dimension
  SUPPRESS_UTM_PROJECTION_VS30 = SUPPRESS_UTM_PROJECTION
  npx_interface_vs30 = npx_interface
  npy_interface_vs30 = npy_interface
  orig_x_interface_vs30 = orig_x_interface
  orig_y_interface_vs30 = orig_y_interface
  spacing_x_interface_vs30 = spacing_x_interface
  spacing_y_interface_vs30 = spacing_y_interface

  ! allocate interface
  allocate(interface_model_vs30(npx_interface,npy_interface),stat=ier)
  if (ier /= 0) stop 'Error allocating Vs30 interface'
  interface_model_vs30(:,:) = 0.0_CUSTOM_REAL

  ! read Vs30 data records
  ! statistics
  max_vs = - HUGEVAL
  min_vs = HUGEVAL

  ! reads in interface points
  do iy = 1,npy_interface
    do ix = 1,npx_interface
      ! reading in data records
      call read_next_line(IIN,string_read)

      read(string_read,*,iostat=ier) vs
      if (ier /= 0) stop 'Error reading interface vs value'

      ! stores values for interpolation
      interface_model_vs30(ix,iy) = vs

      ! statstics
      if (vs < min_vs) min_vs = vs
      if (vs > max_vs) max_vs = vs
    enddo
  enddo
  close(IIN)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '       original interface Vs30 model: vs min/max = ',min_vs,max_vs,'(m/s)'
    write(IMAIN,*) '       done interface reading'
    write(IMAIN,*)
    write(IMAIN,*) '       using Vs30 model values for GLL points in top layer of 30 m depth'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! set model flag
  USE_MODEL_LAYER_VS30 = .true.

contains

    subroutine read_next_line(unit_in,string_read)

    use constants, only: MAX_STRING_LEN
    implicit none

    integer :: unit_in
    character(len=MAX_STRING_LEN) :: string_read

    integer :: ier

    do
      read(unit=unit_in,fmt="(a)",iostat=ier) string_read
      if (ier /= 0) stop 'Error while reading Vs30 interface file'

      ! suppress leading white spaces, if any
      string_read = trim(adjustl(string_read))

      ! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
      if (index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)

      ! reads next line if empty
      if (len_trim(string_read) == 0) cycle

      ! exit loop when we find the first line that is not a comment or a white line
      if (string_read(1:1) /= '#') exit
    enddo

    ! suppress trailing white spaces, if any
    string_read = string_read(1:len_trim(string_read))

    ! suppress trailing comments, if any
    if (index(string_read,'#') > 0) string_read = string_read(1:index(string_read,'#')-1)

    ! suppress leading and trailing white spaces again, if any, after having suppressed the leading junk
    string_read = adjustl(string_read)
    string_read = string_read(1:len_trim(string_read))

    end subroutine read_next_line

  end subroutine model_vs30_setup_check

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_vs30(xmesh,ymesh,zmesh,rho,vp,vs,idomain_id)

! assigns Vs30 velocities (vs and derived vp and rho) to top 30-m layer

  use constants, only: CUSTOM_REAL,HUGEVAL,IDOMAIN_ELASTIC,IUTM2LONGLAT

  use shared_parameters, only: SUPPRESS_UTM_PROJECTION

  use generate_databases_par, only: nspec => NSPEC_AB,ibool

  use create_regions_mesh_ext_par

  implicit none

  ! GLL point
  double precision, intent(in) :: xmesh,ymesh,zmesh
  ! density, Vp and Vs
  real(kind=CUSTOM_REAL),intent(inout) :: rho,vp,vs
  ! acoustic/elastic/.. domain flag ( 1 = acoustic / 2 = elastic / ... )
  integer, intent(in) :: idomain_id

  ! local parameters
  double precision :: x,y,z
  real(kind=CUSTOM_REAL) :: x_target,y_target
  real(kind=CUSTOM_REAL) :: depth
  real(kind=CUSTOM_REAL) :: elevation,distmin

  ! interface
  logical :: SUPPRESS_UTM_PROJECTION_COPY
  integer :: icornerlat,icornerlong
  double precision :: lat,long
  double precision :: long_corner,lat_corner,ratio_xi,ratio_eta

  real(kind=CUSTOM_REAL) :: vs_interface,vp_interface,rho_interface

  ! only applies to elastic domains
  if (idomain_id /= IDOMAIN_ELASTIC) return

  ! GLL point location
  x = xmesh
  y = ymesh
  z = zmesh

  ! get elevation at point location
  x_target = real(x,kind=CUSTOM_REAL)
  y_target = real(y,kind=CUSTOM_REAL)

  elevation = 0.0_CUSTOM_REAL
  distmin = HUGEVAL

  ! get approximate topography elevation at target coordinates from free surface
  call get_topo_elevation_free_closest(x_target,y_target, &
                                       elevation,distmin, &
                                       nspec,nglob_unique,ibool,xstore_unique,ystore_unique,zstore_unique, &
                                       num_free_surface_faces,free_surface_ispec,free_surface_ijk)

  ! determines depth of GLL point
  ! note: z coordinate convention is z-axis points up
  if (distmin < HUGEVAL) then
    ! depth in m
    depth = elevation - z
  else
    ! point not found in this slice
    ! nothing to do
    return
  endif

  ! only for top 30-m layer
  ! checks GLL point depth (adds small tolerance 1.d-6
  if (depth > 30.000001_CUSTOM_REAL) then
    ! point below 30 m layer
    ! nothing to do
    return
  endif

  ! get interface value
  SUPPRESS_UTM_PROJECTION_COPY = SUPPRESS_UTM_PROJECTION

  ! project x and y in UTM back to long/lat since interface file is in long/lat
  SUPPRESS_UTM_PROJECTION = SUPPRESS_UTM_PROJECTION_VS30
  call utm_geo(long,lat,x,y,IUTM2LONGLAT)

  ! restore original flag
  SUPPRESS_UTM_PROJECTION = SUPPRESS_UTM_PROJECTION_COPY

  ! get coordinate of corner in interface model
  icornerlong = int((long - orig_x_interface_vs30) / spacing_x_interface_vs30) + 1
  icornerlat = int((lat - orig_y_interface_vs30) / spacing_y_interface_vs30) + 1

  ! avoid edge effects and extend with identical point if outside model
  if (icornerlong < 1) icornerlong = 1
  if (icornerlong > npx_interface_vs30-1) icornerlong = npx_interface_vs30-1
  if (icornerlat < 1) icornerlat = 1
  if (icornerlat > npy_interface_vs30-1) icornerlat = npy_interface_vs30-1

  ! compute coordinates of corner
  long_corner = orig_x_interface_vs30 + (icornerlong-1)*spacing_x_interface_vs30
  lat_corner = orig_y_interface_vs30 + (icornerlat-1)*spacing_y_interface_vs30

  ! compute ratio for interpolation
  ratio_xi = (long - long_corner) / spacing_x_interface_vs30
  ratio_eta = (lat - lat_corner) / spacing_y_interface_vs30

  ! avoid edge effects
  if (ratio_xi < 0.d0) ratio_xi = 0.d0
  if (ratio_xi > 1.d0) ratio_xi = 1.d0
  if (ratio_eta < 0.d0) ratio_eta = 0.d0
  if (ratio_eta > 1.d0) ratio_eta = 1.d0

  ! interpolate elevation at current point
  vs_interface = &
       interface_model_vs30(icornerlong,icornerlat)*(1.d0-ratio_xi)*(1.d0-ratio_eta) + &
       interface_model_vs30(icornerlong+1,icornerlat)*ratio_xi*(1.d0-ratio_eta) + &
       interface_model_vs30(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
       interface_model_vs30(icornerlong,icornerlat+1)*(1.d0-ratio_xi)*ratio_eta

  ! estimate vp from vs
  ! note: we're using Brocher's relationship for now in lack of better alternatives for near-surface Vp/Vs relationships.
  !       there would be for example Mayne's scaling for density, but it lacks a scaling for Vp.
  !       Also, Boore might have some more applicable scaling relations.
  !       TODO: find better scalings for near-surface Vp & density from Vs30
  vp_interface = scale_Brocher_vp_from_vs(vs_interface)
  rho_interface = scale_Brocher_rho_from_vp(vp_interface)

  ! uses vs30 value for this GLL point
  vs = vs_interface
  vp = vp_interface
  rho = rho_interface

  return

contains

    ! --------------------------------
    ! subroutine to scale missing Vp from Vs
    ! --------------------------------
    function scale_Brocher_vp_from_vs(vs_in) result(vp)

    use constants, only: CUSTOM_REAL

    implicit none

    real(kind=CUSTOM_REAL),intent(in) :: vs_in

    ! local parameters
    real(kind=CUSTOM_REAL) :: vs,vs_p2,vs_p3,vs_p4
    real(kind=CUSTOM_REAL) :: vp

    ! Brocher 2005, Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust, BSSA
    ! factors from eq. (1)
    real(kind=CUSTOM_REAL),parameter :: fac1 = 0.9409_CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: fac2 = 2.0947_CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: fac3 = -0.8206_CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: fac4 = 0.2683_CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: fac5 = -0.0251_CUSTOM_REAL

    ! scaling requires vs in km/s
    ! unit scaling factor: given in m/s -> km/s
    real(kind=CUSTOM_REAL), parameter :: unit_scale = 1.0_CUSTOM_REAL / 1000.0_CUSTOM_REAL

    ! Vs (in km/s)
    vs = vs_in * unit_scale

    vs_p2 = vs * vs
    vs_p3 = vs * vs_p2
    vs_p4 = vs * vs_p3

    ! scaling relation: eq.(1)
    vp = fac1 + fac2 * vs + fac3 * vs_p2 + fac4 * vs_p3 + fac5 * vs_p4

    ! vp in km/s, output in same unit as input: vp km/s -> m/s
    vp = vp * 1000.0_CUSTOM_REAL

    end function scale_Brocher_vp_from_vs

    ! --------------------------------
    ! subroutine to scale missing density from Vp
    ! --------------------------------
    function scale_Brocher_rho_from_vp(vp_in) result(rho)

    use constants, only: CUSTOM_REAL

    implicit none

    real(kind=CUSTOM_REAL), intent(in) :: vp_in

    ! local parameters
    real(kind=CUSTOM_REAL) :: vp,vp_p2,vp_p3,vp_p4,vp_p5
    real(kind=CUSTOM_REAL) :: rho

    ! Brocher 2005, Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust, BSSA
    ! factors from eq. (1)
    real(kind=CUSTOM_REAL),parameter :: fac1 = 1.6612_CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: fac2 = -0.4721_CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: fac3 = 0.0671_CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: fac4 = -0.0043_CUSTOM_REAL
    real(kind=CUSTOM_REAL),parameter :: fac5 = 0.000106_CUSTOM_REAL

    ! scaling requires vp in km/s
    ! unit scaling factor: given in m/s -> km/s
    real(kind=CUSTOM_REAL), parameter :: unit_scale = 1.0_CUSTOM_REAL / 1000.0_CUSTOM_REAL

    ! Vp (in km/s)
    vp = vp_in * unit_scale

    vp_p2 = vp * vp
    vp_p3 = vp * vp_p2
    vp_p4 = vp * vp_p3
    vp_p5 = vp * vp_p4

    ! scaling relation: eq.(1)
    rho = fac1 * vp + fac2 * vp_p2 + fac3 * vp_p3 + fac4 * vp_p4 + fac5 * vp_p5

    ! density rho given in g/cm^3
    ! density scaling for rho in kg/m^3: converts g/cm^3 -> kg/m^3
    ! rho [kg/m^3] = rho * 1000 [g/cm^3]
    rho = rho * 1000.0_CUSTOM_REAL

    end function scale_Brocher_rho_from_vp

! not used yet...
!
!    ! --------------------------------
!    ! subroutine to scale missing density from Vs30
!    ! --------------------------------
!    function scale_Mayne_rho_from_vs(vs_in,depth_in) result(rho)
!
!    use constants, only: CUSTOM_REAL
!
!    implicit none
!
!    real(kind=CUSTOM_REAL), intent(in) :: vs_in,depth_in
!
!    ! local parameters
!    real(kind=CUSTOM_REAL), parameter :: rho_min = 1.65 ! lower bound of density: 1.65 g/cm^3
!    real(kind=CUSTOM_REAL) :: vs,depth
!    real(kind=CUSTOM_REAL) :: rho
!
!    ! reference:
!    ! Mayne, Schneider & Martin (1999)
!    ! "Small- and large-strain soil properties from seismic flat dilatometer tests."
!    ! Pre-failure deformation characteristics of geomaterials, 1999 Balkema, Rotterdam.
!    !
!    ! see: https://www.marchetti-dmt.it/wp-content/uploads/bibliografia/mayne_1999_torino_SDMT_small_large_strain.pdf
!
!    ! vs in m/s
!    vs = vs_in
!
!    ! depth in m
!    depth = depth_in
!
!    ! avoid error of dividing by zero
!    if (depth <= 0.0_CUSTOM_REAL) depth = 0.00001
!    if (vs <= 0.0_CUSTOM_REAL) vs = 0.00001
!
!    ! note: Eq. 2 is using log(z). we assume that this is the log10 of depth.
!    !       thus, we use log10(depth) here instead of log(depth) as written in Eq. 2.
!    rho = max(rho_min, 1.0 + 1.0 / (0.614 + 58.7 * (log10(depth) + 1.095) / vs_in))
!
!    ! density rho given in g/cm^3
!    ! density scaling for rho in kg/m^3: converts g/cm^3 -> kg/m^3
!    ! rho [kg/m^3] = rho * 1000 [g/cm^3]
!    rho = rho * 1000.0_CUSTOM_REAL
!
!    end function scale_Mayne_rho_from_vs

  end subroutine model_vs30
