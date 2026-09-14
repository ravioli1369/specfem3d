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

!=====================================================================
!
!  UTM (Universal Transverse Mercator) projection from the USGS
!
!=====================================================================

! taken from the open-source package CAMx at http://www.camx.com/download
! converted by Dimitri Komatitsch to Fortran90 and slightly modified to add one parameter to the subroutine call
! and change the order of the arguments for compatibility with SPECFEM calls.
! Also converted to double precision for all calculations and all results.
! Also defined the UTM easting and northing in meters instead of kilometers for compatibility with SPECFEM calls.

! convert geodetic longitude and latitude to UTM, and back
! use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
! a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm


! original publication:
! Snyder, J.P., Map projections - a working manual, USGS professional paper 1395, 1987
! https://doi.org/10.3133/pp1395
!
! (Universal) Transverse Mercator projection - formulas for the ellipsoid (pages 61-63) are implemented
!
! note: the original version of this utm_geo(..) routine contained 2 errors in the formula
!       in factor M (rm) when converting from long/lat to UTM, and
!       in factor f3 for computing latitude (lat) when converting from UTM to long/lat.
!
!       mistakes size in original routine:
!         from UTM x/y to lon/lat and converted back again:
!           477415.5             5712313.5            ->  2.6741959317615298        51.561449486527074  (zone 31)
!           477415.50000327773   5712320.9415952200 <-  2.6741959317615298        51.561449486527074
!         has maximum UTM error ~ 7.44 m
!
!      corrected version:
!           477415.5             5712313.5            ->  2.6741959317615298        51.561449479910003  (zone 31)
!           477415.49999999994   5712313.5002520066 <-  2.6741959317615298        51.561449479910003
!         has maximum UTM error ~ 0.00025 m
!
!      error introduced in converting lon/lat to UTM:
!           2.6741959317615298        51.561449486527074 -> 477415.50000327773   5712320.9415952200  ! wrong formula
!                                                           477415.50579994149   5712313.6646418981  ! corrected


  subroutine utm_geo(rlon4,rlat4,rx4,ry4,iway)

!
!---- CAMx v6.10 2014/04/02
!
!     UTM_GEO performs UTM to geodetic (long/lat) translation, and back.
!
!     This is a Fortran version of the BASIC program "Transverse Mercator
!     Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
!     Based on algorithm taken from "Map Projections Used by the USGS"
!     by John P. Snyder, Geological Survey Bulletin 1532, USDI.
!
!     Portions Copyright 1996 - 2014
!     ENVIRON International Corporation
!
!     Modifications:
!      2012/12/02   Added logic for southern hemisphere UTM zones
!                   Use zones +1 to +60 for the Northern hemisphere, -1 to -60 for the Southern hemisphere
!                   Equator is defined as 0 km North for the Northern hemisphere, 10,000 km North for the Southern hemisphere
!
!     Input/Output arguments:
!
!        rlon4                 Longitude (degrees, negative for West)
!        rlat4                 Latitude (degrees)
!        rx4                   UTM easting (meters)
!        ry4                   UTM northing (meters)
!        UTM_PROJECTION_ZONE   UTM zone
!                              The Northern hemisphere corresponds to zones +1 to +60
!                              The Southern hemisphere corresponds to zones -1 to -60
!        iway                  Conversion type
!                              ILONGLAT2UTM = geodetic to UTM
!                              IUTM2LONGLAT = UTM to geodetic

! Some general information about UTM:
! (for more details see e.g. http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system )
!
! There are 60 longitudinal projection zones numbered 1 to 60 starting at 180 degrees W.
! Each of these zones is 6 degrees wide, apart from a few exceptions around Norway and Svalbard.
! There are 20 latitudinal zones spanning the latitudes 80 degrees S to 84 degrees N and denoted
! by the letters C to X, ommitting the letter O.
! Each of these is 8 degrees south-north, apart from zone X which is 12 degrees south-north.
!
! The UTM zone is described by the central meridian of that zone, i.e. the longitude at the midpoint of the zone,
! 3 degrees away from both zone boundary.
!
  use constants, only: PI,ILONGLAT2UTM,IUTM2LONGLAT
  use shared_input_parameters, only: UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,USE_LUNAR_PROJECTIONS

  implicit none

! input/output parameters
  double precision, intent(inout) :: rx4,ry4,rlon4,rlat4
  integer, intent(in) :: iway

! local parameters

! From http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system :
! The Universal Transverse Mercator coordinate system was developed by the United States Army Corps of Engineers in the 1940s.
! The system was based on an ellipsoidal model of Earth. For areas within the contiguous United States
! the Clarke Ellipsoid of 1866 was used. For the remaining areas of Earth, including Hawaii, the International Ellipsoid was used.
! The WGS84 ellipsoid is now generally used to model the Earth in the UTM coordinate system,
! which means that current UTM northing at a given point can be 200+ meters different from the old one.
! For different geographic regions, other datum systems (e.g.: ED50, NAD83) can be used.

! Clarke 1866
! double precision, parameter :: SEMI_MAJOR_AXIS = 6378206.4d0, SEMI_MINOR_AXIS = 6356583.8d0

! WGS84 (World Geodetic System 1984)
  double precision, parameter :: SEMI_MAJOR_AXIS = 6378137.0d0, SEMI_MINOR_AXIS = 6356752.314245d0

! Note that the UTM grids are actually Mercators which
! employ the standard UTM scale factor 0.9996 and set the Easting Origin to 500,000.
  double precision, parameter :: scfa = 0.9996d0
  double precision, parameter :: north = 0.d0, east = 500000.d0

  double precision, parameter :: DEGREES_TO_RADIANS = PI/180.d0, RADIANS_TO_DEGREES = 180.d0/PI

! local variables
  integer :: zone
  double precision :: rlon,rlat
  double precision :: e2,e4,e6,ep2,xx,yy,dlat,dlon,cm,cmr,delam
  double precision :: f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
  double precision :: rx_save,ry_save,rlon_save,rlat_save
  logical :: lsouth

  ! checks if conversion to UTM has to be done
  if (SUPPRESS_UTM_PROJECTION) then
    if (iway == ILONGLAT2UTM) then
      rx4 = rlon4
      ry4 = rlat4
    else
      rlon4 = rx4
      rlat4 = ry4
    endif
    ! all done
    return
  endif

  ! check if we need to use a lunar projection instead of UTM
  if (USE_LUNAR_PROJECTIONS) then
    ! convert using Lunar Transverse Mercator (LTM) and Lunar Polar Stereographic (LPS) projections
    call ltm_geo(rlon4,rlat4,rx4,ry4,iway)
    ! all done
    return
  endif

  ! save original parameters
  rlon_save = rlon4
  rlat_save = rlat4
  rx_save = rx4
  ry_save = ry4

  ! defines parameters of reference ellipsoid
  e2 = 1.d0 - (SEMI_MINOR_AXIS/SEMI_MAJOR_AXIS)**2
  e4 = e2*e2
  e6 = e2*e4
  ep2 = e2 / (1.d0-e2)

!
!---- Set Zone parameters
!

  lsouth = .false.
  if (UTM_PROJECTION_ZONE < 0) lsouth = .true.
  zone = abs(UTM_PROJECTION_ZONE)

  ! set central meridian for this zone
  cm = zone * 6.0d0 - 183.d0
  cmr = cm * DEGREES_TO_RADIANS

  if (iway == IUTM2LONGLAT) then
    xx = rx4
    yy = ry4
    if (lsouth) yy = yy - 1.d7
  else
    dlat = rlat4
    dlon = rlon4
  endif

!
!---- Lat/Lon to UTM conversion
!
  if (iway == ILONGLAT2UTM) then

    rlon = DEGREES_TO_RADIANS*dlon
    rlat = DEGREES_TO_RADIANS*dlat

    delam = dlon - cm
    if (delam < -180.d0) delam = delam + 360.d0
    if (delam > 180.d0) delam = delam - 360.d0

    delam = delam*DEGREES_TO_RADIANS

    ! Snyder, J.P., Map projections - a working manual, USGS professional paper 1395, 1987
    ! https://doi.org/10.3133/pp1395
    !
    ! (Universal) Transverse Mercator projection
    !
    ! page 61, eq. 3-21 for M
    f1 = (1.d0 - e2/4.d0 - 3.d0*e4/64.d0 - 5.d0*e6/256.d0)*rlat
    f2 = 3.d0*e2/8.d0 + 3.d0*e4/32.d0 + 45.d0*e6/1024.d0
    f2 = f2 * sin(2.d0*rlat)
    ! corrected: using .. + 45 e6 / 1024 instead of .. * 45 e6 / 1024
    !!wrong: f3 = 15.0 * e4 / 256.0 * 45.0 * e6 /1024.0
    f3 = 15.d0*e4/256.d0 + 45.d0*e6/1024.d0
    f3 = f3 * sin(4.d0*rlat)
    f4 = 35.d0*e6/3072.d0
    f4 = f4 * sin(6.d0*rlat)
    rm = SEMI_MAJOR_AXIS*(f1 - f2 + f3 - f4)

    if (dlat == 90.d0 .or. dlat == -90.d0) then
      xx = 0.d0
      yy = scfa*rm
    else
      ! page 61, eq. 4-20
      rn = SEMI_MAJOR_AXIS/sqrt(1.d0 - e2*sin(rlat)**2)
      ! page 61, eq. 8-13
      t = tan(rlat)**2
      ! page 61, eq. 8-14
      c = ep2*cos(rlat)**2
      a = cos(rlat)*delam

      ! page 61, eq. 8-9 for x
      f1 = (1.d0 - t + c) * a**3 / 6.d0
      f2 = 5.d0 - 18.d0*t + t**2 + 72.d0*c - 58.d0*ep2
      f2 = f2 * a**5 / 120.d0
      xx = scfa * rn * (a + f1 + f2)

      ! page 61, eq. 8-10 for y
      f1 = a**2 / 2.d0
      f2 = 5.d0 - t + 9.d0*c + 4.d0 * c**2
      f2 = f2 * a**4 / 24.d0
      f3 = 61.d0 - 58.d0*t + t**2 + 600.d0*c - 330.d0*ep2
      f3 = f3 * a**6 / 720.d0
      yy = scfa * (rm + rn * tan(rlat) * (f1 + f2 + f3))
    endif

    xx = xx + east
    yy = yy + north

!
!---- UTM to Lat/Lon conversion
!
  else

    xx = xx - east
    yy = yy - north

    ! Snyder, J.P., Map projections - a working manual, USGS professional paper 1395, 1987
    ! https://doi.org/10.3133/pp1395
    !
    ! (Universal) Transverse Mercator projection
    ! inverse formulas
    !
    ! page 63, eq. 3-24 for e_1
    e1 = sqrt(1.d0 - e2)
    e1 = (1.d0 - e1)/(1.d0 + e1)
    rm = yy/scfa
    ! page 63, eq. 7-19 for mu
    u = 1.d0 - e2/4.d0 - 3.d0*e4/64.d0 - 5.d0*e6/256.d0
    u = rm/(SEMI_MAJOR_AXIS*u)

    ! page 63, eq. 3-26 for phi_1
    f1 = 3.d0 * e1/2.d0 - 27.d0 * e1**3 / 32.d0
    f1 = f1 * sin(2.d0*u)
    f2 = 21.d0 * e1**2 / 16.d0 - 55.d0 * e1**4 / 32.d0
    f2 = f2 * sin(4.d0*u)
    f3 = 151.d0 * e1**3 / 96.d0
    f3 = f3 * sin(6.d0*u)
    rlat1 = u + f1 + f2 + f3
    dlat1 = rlat1 * RADIANS_TO_DEGREES

    if (dlat1 >= 90.d0 .or. dlat1 <= -90.d0) then
      dlat1 = dmin1(dlat1,90.d0)
      dlat1 = dmax1(dlat1,-90.d0)
      dlon = cm
    else
      c1 = ep2*cos(rlat1)**2
      t1 = tan(rlat1)**2
      f1 = 1.d0 - e2*sin(rlat1)**2
      rn1 = SEMI_MAJOR_AXIS/sqrt(f1)
      r1 = SEMI_MAJOR_AXIS*(1.d0 - e2)/sqrt(f1**3)
      d = xx/(rn1*scfa)

      ! page 63, eq. 8-17 for phi
      f1 = rn1 * tan(rlat1)/r1
      f2 = d**2 / 2.d0
      ! corrected: factor 5 + 3 T1 + .. instead of 5 * 3 T1 ..
      !!wrong: f3 = 5.d0*3.d0*t1 + 10.d0*c1 - 4.d0*c1**2 - 9.d0*ep2
      f3 = 5.d0 + 3.d0*t1 + 10.d0*c1 - 4.d0*c1**2 - 9.d0*ep2
      f3 = f3 * d**4 / 24.d0
      f4 = 61.d0 + 90.d0*t1 + 298.d0*c1 + 45.d0*t1**2 - 252.d0*ep2 - 3.d0*c1**2
      f4 = f4 * d**6 / 720.d0
      rlat = rlat1 - f1*(f2 - f3 + f4)
      dlat = rlat * RADIANS_TO_DEGREES

      ! page 63, eq. 8-18 for lambda
      f1 = 1.d0 + 2.d0*t1 + c1
      f1 = f1 * d**3 / 6.d0
      f2 = 5.d0 - 2.d0*c1 + 28.d0*t1 - 3.d0*c1**2 + 8.d0*ep2 + 24.d0*t1**2
      f2 = f2 * d**5 / 120.d0
      rlon = cmr + (d - f1 + f2)/cos(rlat1)
      dlon = rlon * RADIANS_TO_DEGREES

      if (dlon < -180.d0) dlon = dlon + 360.d0
      if (dlon > 180.d0) dlon = dlon - 360.d0
    endif
  endif

!
!----- output
!
  if (iway == IUTM2LONGLAT) then
    rlon4 = dlon
    rlat4 = dlat
    rx4 = rx_save
    ry4 = ry_save
  else
    rx4 = xx
    if (lsouth) yy = yy + 1.d7
    ry4 = yy
    rlon4 = rlon_save
    rlat4 = rlat_save
  endif

  end subroutine utm_geo


!=====================================================================
!
!  Lunar Transverse Mercator (LTM) and Lunar Polar Stereographic (LPS) projection
!
!=====================================================================

  subroutine ltm_geo(rlon4,rlat4,rx4,ry4,iway)

!     LTM_GEO performs LTM/LPS to geodetic (lon/lat) translation, and back.
!
!     This routine assumes a perfectly spherical Moon.
!     LTM implementation is a simplified Spherical Transverse Mercator, 8-degree zones.
!
!     Note: by default, the LTM implementation in the python script
!     LGRS_Coordinate_Conversion_mk7.2.py provided by USGS astrogeology site
!https://astrogeology.usgs.gov/search/map/lunar-map-projections-and-grid-reference-system-for-artemis-astronaut-surface-navigation
!     uses a spherical Moon for their LTM projection.
!
!     However, they still use the formula for a Gauss-Schreiber projection,
!     that combines a projection of an ellipsoid to a sphere followed by a
!     spherical transverse Mercator formula. Assuming the ellipsoid is perfectly
!     spherical, then the first projection is an identity transform and only the
!     second, spherical transverse Mercator formula is needed.
!     Here, we use this simplification of applying directly the spherical
!     Transverse Mercator projection, assuming a spherical Moon.
!
!     Input/Output arguments:
!
!        rlon4                 Longitude (degrees, range [-180,180])
!        rlat4                 Latitude (degrees, range [-90,90])
!        rx4                   LTM/LPS easting (meters)
!        ry4                   LTM/LPS northing (meters)
!        zone                  LTM/LPS zone
!                              LTM zones use 1 to 45, positive for Northern hemisphere,
!                              negative for Southern hemisphere
!                              LPS zones use 46 for North pole, -46 for South pole
!        iway                  Conversion type
!                              ILONGLAT2LTM = ILONGLAT2UTM = geodetic to LTM/LPS
!                              ILTM2LONGLAT = IUTM2LONGLAT = LTM/LPS to geodetic

  use constants, only: PI,ILONGLAT2UTM
  use shared_input_parameters, only: UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,USE_LUNAR_PROJECTIONS

  implicit none

  ! input/output parameters
  double precision, intent(inout) :: rx4,ry4,rlon4,rlat4
  integer, intent(in) :: iway

  ! local parameters
  double precision, parameter :: DEGREES_TO_RADIANS = PI/180.d0, RADIANS_TO_DEGREES = 180.d0/PI

  ! Lunar spherical radius
  double precision, parameter :: R = 1737400.0d0   ! mean lunar radius (meters)

  ! Lunar map projection scale factors
  double precision, parameter :: scfa = 0.999d0         ! transverse Mercator (LTM)
  double precision, parameter :: scfa_polar = 0.994d0   ! polar stereographic (LPS)

  ! False origins
  double precision, parameter :: false_east  =  250000.0d0
  double precision, parameter :: false_north = 2500000.0d0
  double precision, parameter :: false_east_polar  =  500000.0d0
  double precision, parameter :: false_north_polar =  500000.0d0

  ! local variables
  integer :: zone,z
  double precision :: rlon,rlat
  double precision :: cm,cmr,dlam
  double precision :: xx,yy
  double precision :: B,D,T,A,t_polar,rho
  logical :: use_polar

  ! checks if conversion to LTM/LPS has to be done
  if (SUPPRESS_UTM_PROJECTION) then
    if (iway == ILONGLAT2UTM) then
      rx4 = rlon4
      ry4 = rlat4
    else
      rlon4 = rx4
      rlat4 = ry4
    endif
    ! all done
    return
  endif

  ! check if lunar projection
  if (.not. USE_LUNAR_PROJECTIONS) &
    stop 'Error: USE_LUNAR_PROJECTIONS lunar projection flag is invalid in ltm_geo() routine'

  ! check zone
  if (UTM_PROJECTION_ZONE == 0 .or. abs(UTM_PROJECTION_ZONE) > 46) then
    ! user info
    print *,'Error: zone ',zone,'is invalid for Moon projections LTM/LPS; must be +/- 1-45 for LTM, +/- 46 for LPS'
    print *
    ! determine zone if it is unspecified in forward mode, compute from longitude/latitude position
    if (iway == ILONGLAT2UTM) then
      call get_ltm_zone(rlon4,rlat4,zone)
      ! user info
      print *,'       For the current lunar position lat/lon = ',sngl(rlat4),'/',sngl(rlon4)
      print *,'       the correct LTM/LPS zone would be      = ',zone
      print *
    endif
    print *,'       Please check the UTM_PROJECTION_ZONE setting in your Par_file.'
    print *
    ! returns error
    stop 'Error: zone is invalid for Moon projections LTM/LPS; must be +/- 1-45 for LTM, +/- 46 for LPS'
  endif

  ! Initialize
  use_polar = .false.

  ! zone index absolute
  zone = UTM_PROJECTION_ZONE
  if (zone == 0) zone = 1
  z = abs(zone)

  ! check if LPS or LTM
  if (z == 46) use_polar = .true.

!
!---- Lunar Polar Stereographic (LPS) projection
!
  if (use_polar) then

    if (iway == ILONGLAT2UTM) then
      ! Forward transformation: lon/lat to LPS
      ! Snyder (1987) spherical equations
      rlon = rlon4 * DEGREES_TO_RADIANS
      rlat = rlat4 * DEGREES_TO_RADIANS

      ! Calculate polar stereographic spherical scale factor
      A = 2.0d0 * R * scfa_polar

      if (zone < 0) then
        ! South pole
        t_polar = tan(PI/4.0d0 + rlat/2.0d0)
      else
        ! North pole
        t_polar = tan(PI/4.0d0 - rlat/2.0d0)
      endif

      ! Stereographic map projection, Snyder (1987)
      ! Note: here, we use +cos(lon) as done in routine spherical_stereographic_map_y()
      ! of the USGS script LGRS_Coordinate_Conversion_mk7.2.py for both, North and South poles.
      ! The standard polar stereographic projection defined by Snyder would use
      ! North Pole: x = 2 R k0 tan(pi/4 - phi/2) sin(lambda - lambda0)
      !             y = - 2 R k0 tan(pi/4 - phi/2) cos(lambda - lambda0)
      ! South Pole: x = 2 R k0 tan(pi/4 + phi/2) sin(lambda - lambda0)
      !             y = 2 R k0 tan(pi/4 + phi/2) cos(lambda - lambda0)
      ! (see Snyder 1987, page 158, chapter 21 "Stereographic Projection",
      !  section "Formulas for the Sphere", eqs. 21-5, 21-6 and 21-9,21-10.
      !  https://pubs.usgs.gov/pp/1395/report.pdf ).
      ! We take the convention from the USGS implementation to use +cos for both poles.

      ! X coordinate for a stereographic map projection
      xx = A * t_polar * sin(rlon)
      ! Y coordinate for a stereographic map projection
      yy = A * t_polar * cos(rlon)

      ! Add false Eastings and Northings
      xx = xx + false_east_polar
      yy = yy + false_north_polar

      ! Forward mode: set computed x,y values
      rx4 = xx
      ry4 = yy

    else
      ! Inverse transformation: LPS to lon/lat
      ! remove false origins (must be same values used in forward)
      xx = rx4 - false_east_polar
      yy = ry4 - false_north_polar

      ! radial distance from projection origin
      rho = sqrt(xx**2 + yy**2)

      ! at the pole: rho == 0
      if (rho == 0.0d0) then
        if (zone < 0) then
          rlat4 = -90.0d0
        else
          rlat4 = 90.0d0
        endif
        rlon4 = 0.0d0
        ! all done
        return
      endif

      A = 2.0d0 * R * scfa_polar

      ! t = tan( PI/4 ± lat/2 )  where sign depends on hemisphere in forward
      t_polar = rho / A

      if (zone < 0) then
        ! South pole
        ! forward used: tan(PI/4 + lat/2) = t  -> lat = 2*atan(t) - PI/2
        rlat = 2.0d0 * atan(t_polar) - PI/2.0d0
      else
        ! North pole
        ! forward used: tan(PI/4 - lat/2) = t  -> lat = PI/2 - 2*atan(t)
        rlat = PI/2.0d0 - 2.0d0 * atan(t_polar)
      endif

      ! longitude from sin/cos ordering used in forward: x = factor*sin(lon), y = factor*cos(lon)
      rlon = atan2(xx, yy)

      ! Inverse mode: set computed lon,lat values
      rlon4 = rlon * RADIANS_TO_DEGREES
      rlat4 = rlat * RADIANS_TO_DEGREES

      ! normalize lon to -180..180
      if (rlon4 > 180.0d0) rlon4 = rlon4 - 360.0d0
      if (rlon4 < -180.0d0) rlon4 = rlon4 + 360.0d0
    endif

  else

!
!---- Lunar Transverse Mercator (LTM) Projection
!

    ! Central meridian of 8-degree zone
    cm = -180.0d0 + (dble(z - 1) + 0.5d0) * 8.0d0
    cmr = cm * DEGREES_TO_RADIANS

    if (iway == ILONGLAT2UTM) then
      ! Forward transformation: lon/lat to LTM
      rlon = rlon4 * DEGREES_TO_RADIANS
      rlat = rlat4 * DEGREES_TO_RADIANS

      dlam = rlon - cmr

      ! Spherical Transverse Mercator (forward)
      B = cos(rlat) * sin(dlam)

      xx = 0.5d0 * R * scfa * log((1.0d0 + B) / (1.0d0 - B))
      yy = R * scfa * atan2(tan(rlat), cos(dlam))

      ! False origins
      xx = xx + false_east
      ! only for South sections
      if (zone < 0) yy = yy + false_north

      ! Forward mode: set computed x,y values
      rx4 = xx
      ry4 = yy

    else
      ! Inverse transformation: LTM to lon/lat
      ! remove false origins
      xx = rx4 - false_east
      if (zone < 0) then
        ! only for South sections
        yy = ry4 - false_north
      else
        yy = ry4
      endif

      ! Spherical TM inverse
      D = yy / (R * scfa)
      T = xx / (R * scfa)

      rlon = cmr + atan2(sinh(T), cos(D))
      rlat = asin(sin(D) / cosh(T))

      ! Inverse mode: set computed lon,lat values
      rlon4 = rlon * RADIANS_TO_DEGREES
      rlat4 = rlat * RADIANS_TO_DEGREES

      ! normalize lon to -180..180
      if (rlon4 > 180.0d0) rlon4 = rlon4 - 360.0d0
      if (rlon4 < -180.0d0) rlon4 = rlon4 + 360.0d0
    endif

  endif

  end subroutine ltm_geo

!
!------------------------------------------------------------------------------------------
!

  subroutine get_ltm_zone(lon,lat,zone)

! gets lunar LTM/LPS zone for given lon/lat (in degrees)

  implicit none

  ! input/output parameters
  double precision, intent(in) :: lon,lat
  ! local variables
  integer, intent(out) :: zone

  ! initializes
  zone = 0

  ! LTM - latitude range for LTM [-82,82], beyond is polar region
  if (lat >= -82.0d0 .and. lat <= 82.0d0) then
    ! longitudinal zone
    zone = int((lon + 180.0d0) / 8.0d0) + 1

    ! 180-degree values sometimes come up as 46. We assign it back to zone 1
    if (zone > 45) zone = zone - 45

    ! we use negative values for Southern hemisphere
    if (lat < 0.0d0) zone = -zone

  else
    ! LPS - polar region
    if (lat >= 0.0d0) then
      zone = 46    ! north pole
    else
      zone = -46   ! south pole
    endif
  endif

  end subroutine get_ltm_zone
