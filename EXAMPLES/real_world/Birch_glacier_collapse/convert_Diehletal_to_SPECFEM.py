#!/usr/bin/env python
#
# script to download and convert the Diehl et al. (2021) model to SPECFEM3D tomography file format
#
from __future__ import print_function

import os,sys
import math
import datetime
import subprocess

## numpy
import numpy as np
print("numpy version: ",np.__version__)

# UTM projection
import pyproj
from pyproj import Transformer

## VTK visualization
try:
    import vtk
except:
    print("Failed importing vtk. Please make sure to have python with vtk working properly.")
    sys.exit(1)

###########################################################################################
## PARAMETERS

# DiehlEtAl2021 repository
repository_name = 'Repository_DiehlEtAl2021_CentralAlpsTomoReloc_v07'

# model files
VP_MODEL_FILE = 'CentralAlp_VelocityModel_3D_Pg_v07.VelocityOnly.txt'
VS_MODEL_FILE = 'CentralAlp_VelocityModel_3D_Sg_v07.VelocityOnly.txt'

# VTK file of model for visualization
create_vtk_file_output = True

###########################################################################################

def download_model_data(data_dir):
    """
    downloads Diehletal repository with model data
    """
    global VP_MODEL_FILE,VS_MODEL_FILE

    # checks if we have the model files already
    if os.path.isfile(data_dir + '/' + VP_MODEL_FILE) and os.path.isfile(data_dir + '/' + VS_MODEL_FILE):
        print(f"model files found in folder: {data_dir}")
        print("")
        return

    # check data dir
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    print("downloading repository files...")

    # current working dir
    current_dir = os.getcwd()

    # change into data dir
    os.chdir(data_dir)

    # Diehl et al. 2021 repository
    repository_tarfile = repository_name + '.tar'

    # download tar file
    url = 'https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/453236/' + repository_tarfile
    cmd = 'wget ' + url
    print(f"  getting tar file {repository_tarfile} ...")
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    if status != 0:
        print("error: status returned ",status)
        sys.exit(status)
    print("")

    # check
    if not os.path.isfile(repository_tarfile):
        print(f"Error: failed to download repository file from:\n  {url}\nPlease download and extract files manually...")
        sys.exit(1)

    # untar
    if not os.path.exists(repository_name):
        cmd = 'tar -xvf ' + repository_tarfile
        print(f"  extracting tar file ...")
        print("  > ",cmd)
        status = subprocess.call(cmd, shell=True)
        if status != 0:
            print("error: status returned ",status)
            sys.exit(status)
        print("")

    # gunzip model files
    # Vp model file
    filename = VP_MODEL_FILE
    if not os.path.isfile(filename):
        gz_file = repository_name + '/' + filename + '.gz'
        # gunzip
        cmd = 'gunzip -c ' + gz_file + ' > ' + filename
        print(f"  extracting file {gz_file}...")
        print("  > ",cmd)
        status = subprocess.call(cmd, shell=True)
        if status != 0:
            print("error: status returned ",status)
            sys.exit(status)
        print("")
    # check
    if not os.path.isfile(filename):
        print(f"Error: failed to extract file {filename}\nPlease download and extract files manually...")
        sys.exit(1)

    # Vs model file
    filename = VS_MODEL_FILE
    if not os.path.isfile(filename):
        gz_file = repository_name + '/' + filename + '.gz'
        # gunzip
        cmd = 'gunzip -c ' + gz_file + ' > ' + filename
        print(f"  extracting file {gz_file}...")
        print("  > ",cmd)
        status = subprocess.call(cmd, shell=True)
        if status != 0:
            print("error: status returned ",status)
            sys.exit(status)
        print("")
    # check
    if not os.path.isfile(filename):
        print(f"Error: failed to extract file {filename}\nPlease download and extract files manually...")
        sys.exit(1)

    print("  model files ready")
    print("")

    # change back to working directory
    os.chdir(current_dir)

#----------------------------------------------------------------------------------------

def read_velocity_model_data(filename):
    """
    reads in data
    """
    print("reading file: ",filename)
    # check file exists
    if not os.path.isfile(filename):
        print(f"file {filename} doesn't exist, exiting...")
        sys.exit(1)

    with open(filename,'r') as f:
        content = f.readlines()

    nlines = len(content)
    print("  number of file lines: ",nlines)

    # check header line
    # format: 1.0 43 43 13  1       computed  20200411 134443.962
    line = content[0]
    items = line.split()
    nx = int(items[1])
    ny = int(items[2])
    nz = int(items[3])
    print(f"  number of increments: nx/ny/nz = {nx}/{ny}/{nz}")
    print("")

    # double check
    assert nx == 43, "number of X increments nx should be 43"
    assert ny == 43, "number of Y increments ny should be 43"
    assert nz == 13, "number of Z increments nz should be 13"

    # read in data
    # Skip the 9 header lines
    data_lines = content[9:]

    # Join all data lines into a single string, then split into a list of numbers.
    # This handles any number of entries per line.
    all_numbers_str = " ".join(data_lines).split()

    # Convert the list of number strings into a 1D NumPy array of floats
    data = np.array(all_numbers_str, dtype=float)

    # Check if the total number of data points is divisible by the layer size
    total_points = nx * ny * nz

    if data.size != total_points:
        print(f"Error: total data points read ({data.size}) is not matching the total {total_points}.")
        sys.exit(1)

    # return data as 3D array
    data_vel = data.reshape((nx, ny, nz),order='F')

    #debug
    #print(data_vel[:,0,0])

    return data_vel

#----------------------------------------------------------------------------------------

def save_as_GeoCSV():
    ## GeoCSV file
    header = f"""
# dataset: GeoCSV2.0
# created: {str(datetime.datetime.now())}
# delimiter: |
# global_title: CentralAlp_VelocityModel_3D
# global_model: CentralAlp_VelocityModel_3D_Pg_Sg_v07
# global_id: CentralAlp_VelocityModel_3D_Pg_Sg_v07
# global_data_revision: v07
# global_Conventions: CF-1.0
# global_Metadata_Conventions: Unidata Dataset Discovery v1.0
# global_summary: see https://www.research-collection.ethz.ch/handle/20.500.11850/453236
# global_keywords: local earthquake tomography model
# global_reference: Diehl et al. 2021
# global_author_name: Diehl, Tobias
# global_acknowledgment: model provided by SPECFEM :)
# global_repository_name: EMC
# global_repository_institution: IRIS EMC
# global_reference_ellipsoid: WGS 84
# global_geodetic_datum: UTM 32N / EPSG:32632
"""

    # needs regular lat/lon?
    # then we would need to approximate the grid, as projecting regular UTM X/Y leads to slightly different lon/lat increments
    print("GeoCSV not implemented yet...")
    sys.exit(1)

#----------------------------------------------------------------------------------------

def geo2utm(lon,lat,zone,iway=2):
    """
    from utm_geo.f90

    convert geodetic longitude and latitude to UTM, and back
    use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
    a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm
    Use zones +1 to +60 for the Northern hemisphere, -1 to -60 for the Southern hemisphere

    implicit none

    include "constants.h"

    -----CAMx v2.03

    UTM_GEO performs UTM to geodetic (long/lat) translation, and back.

    This is a Fortran version of the BASIC program "Transverse Mercator
    Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
    Based on algorithm taken from "Map Projections Used by the USGS"
    by John P. Snyder, Geological Survey Bulletin 1532, USDI.

    Input/Output arguments:

      rlon                  Longitude (deg, negative for West)
      rlat                  Latitude (deg)
      rx                    UTM easting (m)
      ry                    UTM northing (m)
      UTM_PROJECTION_ZONE  UTM zone
      iway                  Conversion type
                              ILONGLAT2UTM = geodetic to UTM
                              IUTM2LONGLAT = UTM to geodetic
    """
    PI      = math.pi
    degrad  = PI/180.0
    raddeg  = 180.0/PI

    # Clarke 1866
    #semimaj = 6378206.4
    #semimin = 6356583.8
    # WGS84 (World Geodetic System 1984) - default
    semimaj = 6378137.0
    semimin = 6356752.314245

    scfa    = 0.9996

    # some extracts about UTM:
    #
    # There are 60 longitudinal projection zones numbered 1 to 60 starting at 180 W.
    # Each of these zones is 6 degrees wide, apart from a few exceptions around Norway and Svalbard.
    # There are 20 latitudinal zones spanning the latitudes 80 S to 84 N and denoted
    # by the letters C to X, ommitting the letter O.
    # Each of these is 8 degrees south-north, apart from zone X which is 12 degrees south-north.
    #
    # To change the UTM zone and the hemisphere in which the
    # calculations are carried out, need to change the fortran code and recompile. The UTM zone is described
    # actually by the central meridian of that zone, i.e. the longitude at the midpoint of the zone, 3 degrees
    # from either zone boundary.
    # To change hemisphere need to change the "north" variable:
    #  - north=0 for northern hemisphere and
    #  - north=10000000 (10000km) for southern hemisphere. values must be in metres i.e. north=10000000.
    #
    # Note that the UTM grids are actually Mercators which
    # employ the standard UTM scale factor 0.9996 and set the
    # Easting Origin to 500,000;
    # the Northing origin in the southern
    # hemisphere is kept at 0 rather than set to 10,000,000
    # and this gives a uniform scale across the equator if the
    # normal convention of selecting the Base Latitude (origin)
    # at the equator (0 deg.) is followed.  Northings are
    # positive in the northern hemisphere and negative in the
    # southern hemisphere.
    north = 0.0
    east  = 500000.0

    IUTM2LONGLAT = 1
    ILONGLAT2UTM = 2

    #---------------------------------------------------------------
    # default uses conversion to UTM
    #iway = ILONGLAT2UTM

    # zone
    UTM_PROJECTION_ZONE = zone

    # lon/lat
    if iway == ILONGLAT2UTM:
        rlon = lon
        rlat = lat
        rx = 0.0
        ry = 0.0
    else:
        rx = lon
        ry = lat
        rlon = 0.0
        rlat = 0.0

    #---------------------------------------------------------------

    # save original parameters
    rlon_save = rlon
    rlat_save = rlat
    rx_save = rx
    ry_save = ry

    # define parameters of reference ellipsoid
    e2 = 1.0 - (semimin/semimaj)**2
    e4 = e2 * e2
    e6 = e2 * e4
    ep2 = e2/(1.0 - e2)

    #----- Set Zone parameters
    # zone
    lsouth = False
    if UTM_PROJECTION_ZONE < 0: lsouth = True
    zone = abs(UTM_PROJECTION_ZONE)

    cm = zone * 6.0 - 183.0       # set central meridian for this zone
    cmr = cm*degrad

    if iway == IUTM2LONGLAT:
        xx = rx
        yy = ry
        if lsouth: yy = yy - 1.e7
    else:
        dlon = rlon
        dlat = rlat

    #---- Lat/Lon to UTM conversion
    if iway == ILONGLAT2UTM:

        rlon = degrad*dlon
        rlat = degrad*dlat

        delam = dlon - cm
        if delam < -180.0: delam = delam + 360.0
        if delam > 180.0: delam = delam - 360.0

        delam = delam*degrad

        f1 = (1. - e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0)*rlat
        f2 = 3.0*e2/8.0 + 3.0*e4/32.0 + 45.0*e6/1024.0
        f2 = f2 * math.sin(2.0*rlat)
        f3 = 15.0*e4/256.0 + 45.0*e6/1024.0
        f3 = f3 * math.sin(4.0*rlat)
        f4 = 35.0*e6/3072.0
        f4 = f4 * math.sin(6.0*rlat)
        rm = semimaj*(f1 - f2 + f3 - f4)

        if dlat == 90.0 or dlat == -90.0:
            xx = 0.0
            yy = scfa*rm
        else:
            rn = semimaj/math.sqrt(1.0 - e2*math.sin(rlat)**2)
            t = math.tan(rlat)**2
            c = ep2 * math.cos(rlat)**2
            a = math.cos(rlat) * delam

            f1 = (1.0 - t + c) * a**3/6.0
            f2 = 5.0 - 18.0*t + t**2 + 72.0*c - 58.0*ep2
            f2 = f2 * a**5/120.0
            xx = scfa*rn*(a + f1 + f2)
            f1 = a**2/2.0
            f2 = 5.0 - t + 9.0*c + 4.0*c**2
            f2 = f2*a**4/24.0
            f3 = 61.0 - 58.0*t + t**2 + 600.0*c - 330.0*ep2
            f3 = f3 * a**6/720.0
            yy = scfa*(rm + rn*math.tan(rlat)*(f1 + f2 + f3))

        xx = xx + east
        yy = yy + north

    else:
        #---- UTM to Lat/Lon conversion
        xx = xx - east
        yy = yy - north
        e1 = math.sqrt(1.0 - e2)
        e1 = (1.0 - e1)/(1.0 + e1)
        rm = yy/scfa
        u = 1.0 - e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0
        u = rm/(semimaj*u)

        f1 = 3.0*e1/2.0 - 27.0*e1**3.0/32.0
        f1 = f1*math.sin(2.0*u)
        f2 = 21.0*e1**2/16.0 - 55.0*e1**4/32.0
        f2 = f2*math.sin(4.0*u)
        f3 = 151.0*e1**3/96.0
        f3 = f3*math.sin(6.0*u)
        rlat1 = u + f1 + f2 + f3
        dlat1 = rlat1*raddeg

        if dlat1 >= 90.0 or dlat1 <= -90.0:
            dlat1 = min(dlat1,90.0)
            dlat1 = max(dlat1,-90.0)
            dlon = cm
        else:
            c1 = ep2*math.cos(rlat1)**2
            t1 = math.tan(rlat1)**2
            f1 = 1.0 - e2*math.sin(rlat1)**2
            rn1 = semimaj/math.sqrt(f1)
            r1 = semimaj*(1.0 - e2)/math.sqrt(f1**3)
            d = xx/(rn1*scfa)

            f1 = rn1*math.tan(rlat1)/r1
            f2 = d**2/2.0
            f3 = 5.0 + 3.0*t1 + 10.0*c1 - 4.0*c1**2 - 9.0*ep2
            f3 = f3*d**2*d**2/24.0
            f4 = 61.0 + 90.0*t1 + 298.0*c1 + 45.0*t1**2 - 252.0*ep2 - 3.0*c1**2
            f4 = f4*(d**2)**3/720.0
            rlat = rlat1 - f1*(f2 - f3 + f4)
            dlat = rlat*raddeg

            f1 = 1.0 + 2.0*t1 + c1
            f1 = f1*d**2*d/6.0
            f2 = 5.0 - 2.0*c1 + 28.0*t1 - 3.0*c1**2 + 8.0*ep2 + 24.0*t1**2
            f2 = f2*(d**2)**2*d/120.0
            rlon = cmr + (d - f1 + f2)/math.cos(rlat1)
            dlon = rlon*raddeg
            if dlon < -180.0: dlon = dlon + 360.0
            if dlon > 180.0: dlon = dlon - 360.0

    if iway == IUTM2LONGLAT:
        rlon = dlon
        rlat = dlat
        rx = rx_save
        ry = ry_save
        return rlon,rlat

    else:
        rx = xx
        if lsouth: yy = yy + 1.e7
        ry = yy
        rlon = rlon_save
        rlat = rlat_save
        return rx,ry

#-----------------------------------------------------------------------------

def scale_rho_from_vp_Brocher(vp_in):
    """
    returns density scaled from vp according to Brocher(2005)'s relationship
    """
    # Brocher 2005, Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust, BSSA
    # factors from eq. (1)
    fac1 = 1.6612
    fac2 = -0.4721
    fac3 = 0.0671
    fac4 = -0.0043
    fac5 = 0.000106

    # scaling requires vp in km/s
    # (assumes input vp is given in m/s)

    # Vp (in km/s)
    vp = vp_in * 1.0/1000.0

    vp_p2 = vp * vp
    vp_p3 = vp * vp_p2
    vp_p4 = vp * vp_p3
    vp_p5 = vp * vp_p4

    # scaling relation: eq.(1)
    # (rho given in g/cm^3)
    rho = fac1 * vp + fac2 * vp_p2 + fac3 * vp_p3 + fac4 * vp_p4 + fac5 * vp_p5

    # file output needs rho in kg/m^3
    # density scaling for rho in kg/m^3: converts g/cm^3 -> kg/m^3
    # rho [kg/m^3] = rho * 1000 [g/cm^3]
    rho *= 1000.0

    return rho

#-----------------------------------------------------------------------------

def convert_2_SPECFEM3D_tomography_file(data_vp,data_vs,data_dir):
    """
    converts to SPECFEM3D format and creates tomography file
    """
    global create_vtk_file_output

    print("converting to SPECFEM3D tomography format...")

    ## velocity model
    # from: CentralAlp_VelocityModel_3D_Pg_v07.FullOutput.txt
    #origin :  latitude   longitude   rotation
    #         47  0.00    -8 30.00      0.00
    # format given in deg, min:
    #  lat 47  0.00 -> 47 deg  0.00 min -> dezimal 47.0
    #  lon -8 30.00 -> -8 deg 30.00 min -> negative? to indicate towards east/west? dezimal 8.5
    coords_origin = [ 8.5, 47.0 ]  # grid center location in CH lon/lat

    print("  velocity model:")
    print(f"  origin           : lon/lat = {coords_origin[0]} / {coords_origin[1]}")

    # converts origin lon/lat to UTM X/Y for UTM zone 32
    # (uses same geo2utm() routine as fortran version to be consistent with SPECFEM conversions)
    x,y = geo2utm(coords_origin[0],coords_origin[1],32)
    coords_origin_UTM = [ x,y ]

    print(f"  origin           : UTM X/Y = {coords_origin_UTM[0]} / {coords_origin_UTM[1]}")
    print("")

    # double-check UTM conversion with pyproj
    if 1 == 0:
        # gets UTM zone codes (uses first, single point in points[] list)
        utm_crs_list = pyproj.database.query_utm_crs_info(datum_name="WGS 84",
                                                          area_of_interest=pyproj.aoi.AreaOfInterest(west_lon_degree=coords_origin[0],
                                                                                                     south_lat_degree=coords_origin[1],
                                                                                                     east_lon_degree=coords_origin[0],
                                                                                                     north_lat_degree=coords_origin[1]))
        utm_code = utm_crs_list[0].code
        utm_epsg = "EPSG:{}".format(utm_code)
        print(f"  pyproj: detected UTM code: {utm_code}  epsg: {utm_epsg}")

        # transformer
        # transformation from WGS84 to UTM zone
        # WGS84: Transformer.from_crs("EPSG:4326", utm_epsg)
        transformer_to_utm = Transformer.from_crs("EPSG:4326", utm_epsg, always_xy=True) # input: lon/lat -> utm_x/utm_y

        utm = transformer_to_utm.transform(coords_origin[0], coords_origin[1])
        print(f"  pyproj: UTM X/Y = {utm[0]} / {utm[1]}")
        print("")

    # original Diehl gridding
    # X coordinate locations in km
    coords_X = [-400.0,-200.0,-190.0,-180.0,-170.0,-160.0,-150.0,-140.0,-130.0,-120.0,-110.0,-100.0, -90.0, -80.0, -70.0, -60.0, -50.0, -40.0, -30.0, -20.0, \
                 -10.0,   0.0,  10.0,  20.0,  30.0,  40.0,  50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, \
                190.0, 200.0, 400.0 ]

    # Y coordinate locations (same as X)
    coords_Y = coords_X

    # Z coordinate locations in km
    # convention: Z-axis pointing down (below sea-level -> positive Z, above sea-level -> negative Z)
    coords_Z = [-7.0,  -5.0,   0.0,   4.0,   9.0,  14.0,  20.0,  27.0,  35.0,  45.0,  60.0,  80.0, 300.0]

    # check gridding regular dlon/dlat?
    #coords_X_geo = []
    #for x in coords_X:
    #    utm_x = coords_origin_UTM[0] + x
    #    utm_y = coords_origin_UTM[1]
    #    lon,lat = geo2utm(utm_x,utm_y,32,iway=1)  # iway==1: utm to geo conversion
    #    coords_X_geo.append([lon,lat])
    #    #or: coords_X_geo.append(transformer_to_utm.transform(utm_x,utm_y,direction='INVERSE'))
    #for i in range(len(coords_X_geo)):
    #    print(f"{i}: {coords_X_geo[i]}")
    #    if i > 0:
    #        dlon = coords_X_geo[i][0] - coords_X_geo[i-1][0]
    #        dlat = coords_X_geo[i][1] - coords_X_geo[i-1][1]
    #        print(f"  dlon: {dlon} dlat: {dlat}")

    # note: SPECFEM tomo file needs regular grid dx/dy/dz increments
    #       thus, we remove the horizontal entries for +/- 400.0 km as they are used as additional outer points for ray tracing in SIMULPS

    # Diehletal original grid parameters
    nx = 43
    ny = 43
    nz = 13

    print("  trimming original grid model...")
    print(f"    original data shape: {data_vp.shape}")

    # Define the slicing indices for each dimension
    # Remove the first and last in nx (axis 0)
    slice_nx = slice(1, -1)
    # Remove the first and last in ny (axis 1)
    slice_ny = slice(1, -1)
    # Remove only the last in nz (axis 2)
    slice_nz = slice(None, -1)

    # Apply the slicing to the array
    model_vp = data_vp[slice_nx, slice_ny, slice_nz]
    model_vs = data_vs[slice_nx, slice_ny, slice_nz]

    print(f"    trimmed  data shape:", model_vp.shape)
    print("")

    # manual trimming
    #model_vp = np.zeros((nx-2,ny-2,nz-1),dtype='float')
    #model_vs = np.zeros((nx-2,ny-2,nz-1),dtype='float')
    #for ix in range(0,nx-2):
    #    for iy in range(0,ny-2):
    #        for iz in range(nz-1):
    #            model_vp[ix,iy,iz] = data_vp[ix+1,iy+1,iz]
    #            model_vs[ix,iy,iz] = data_vs[ix+1,iy+1,iz]

    # scale rho from Vp
    print("  scaling velocities from km/s to m/s...")
    model_vp *= 1000.0
    model_vs *= 1000.0

    # scale rho from Vp
    print("  scaling rho from Vp...")
    # scale density from vp using Brocher relationship
    model_rho = scale_rho_from_vp_Brocher(model_vp)

    # model stats min/max
    VP_MIN = model_vp.min()
    VP_MAX = model_vp.max()
    VS_MIN = model_vs.min()
    VS_MAX = model_vs.max()
    RHO_MIN = model_rho.min()
    RHO_MAX = model_rho.max()

    print("")
    print(f"  trimmed model: Vp  min/max = {VP_MIN:.2f} / {VP_MAX:.2f}")
    print(f"                 Vs  min/max = {VS_MIN:.2f} / {VS_MAX:.2f}")
    print(f"                 rho min/max = {RHO_MIN:.2f} / {RHO_MAX:.2f}")
    print("")

    # tomo grid specification
    # regular gridding for tomo file
    x_regular = np.linspace(-200000.0,200000.0,41)
    y_regular = np.linspace(-200000.0,200000.0,41)
    z_regular = np.linspace(-7000.0,80000.0,88)

    # irregular gridding of original mesh
    z_irreg = np.array(coords_Z[0:12])   # from -7.0 to 80.0
    z_irreg *= 1000.0                    # convert to m

    #print("x gridding:",x_regular)
    #print("z gridding:",z_regular)
    #print("z irreg:",z_irreg)

    # vtk file output
    if create_vtk_file_output:
        print("creating vtk file for visualization...")
        print("")

        # Create points
        points = vtk.vtkPoints()

        # Create model data
        vtk_vp = vtk.vtkFloatArray()
        vtk_vp.SetNumberOfComponents(1)  # 1 components (vp)
        vtk_vp.SetName("Vp")

        # Create model data
        vtk_vs = vtk.vtkFloatArray()
        vtk_vs.SetNumberOfComponents(1)  # 1 components (vs)
        vtk_vs.SetName("Vs")

        # Create model data
        vtk_rho = vtk.vtkFloatArray()
        vtk_rho.SetNumberOfComponents(1)  # 1 components (rho)
        vtk_rho.SetName("Density")

    # initializes header variables
    header_origin_x = sys.float_info.max
    header_origin_y = sys.float_info.max
    header_origin_z = sys.float_info.max

    header_end_x = -sys.float_info.max
    header_end_y = -sys.float_info.max
    header_end_z = -sys.float_info.max

    header_vp_min = sys.float_info.max
    header_vp_max = -sys.float_info.max

    header_vs_min = sys.float_info.max
    header_vs_max = -sys.float_info.max

    header_rho_min = sys.float_info.max
    header_rho_max = -sys.float_info.max

    print("creating tomography model...")
    print("")

    # counters for actual dimension of target model
    ndim_z = 0
    ndim_x = 0
    ndim_y = 0

    output_data = list()

    # loops over depth
    for k, z in enumerate(z_regular):
        # counter
        ndim_z += 1

        # user output
        # Show the progress.
        if k == 0:
            zero_depth = z
        elif k == 1:
            print(f'[PROCESSING] Depth range: {z_regular[0]/1000.0} to {z_regular[-1]/1000.0} (km)', flush=True)
            print(f"{zero_depth/1000.0}, {z/1000.0}", end=' ', flush=True)
        else:
            print(f", {z/1000.0}", end=' ', flush=True)

        # Go through each Y
        for j, y in enumerate(y_regular):
            # counter
            if ndim_z == 1: ndim_y += 1

            # Go through each X
            for i, x in enumerate(x_regular):
                # counter
                if ndim_z == 1 and ndim_y == 1: ndim_x += 1

                # grid point position
                # convert x/y to UTM
                x_val = coords_origin_UTM[0] - x  # original X grid coordinate seems to switch East/West direction by a negative sign
                                                  # (similar to origin lon)
                y_val = coords_origin_UTM[1] + y
                z_val = - z                       # SPECFEM convention Z-direction pointing upwards

                # find index in irregular depth array
                # Find the matching index of the depth value in the irregular z_irreg array
                k_irreg = np.max(np.where(z_irreg <= z))

                # interpolates depth
                k0 = k_irreg
                k1 = k_irreg + 1
                if k1 >= len(z_irreg): k1 = len(z_irreg) - 1

                # original model depths
                z0 = z_irreg[k0]
                z1 = z_irreg[k1]

                # depth interpolation factor
                dz_irreg = z1 - z0
                if abs(dz_irreg) > 0.0:
                    gamma = (z1 - z) / dz_irreg
                else:
                    gamma = 1.0

                #debug
                #if i == 1 and j == 1: print("interpolation: dz irreg = ",z1-z0,"z",z,"gamma",gamma)

                # bounds
                if gamma < 0.0: gamma = 0.0
                if gamma > 1.0: gamma = 1.0

                # model parameters
                vp0 = model_vp[i,j,k0]
                vp1 = model_vp[i,j,k1]

                vs0 = model_vs[i,j,k0]
                vs1 = model_vs[i,j,k1]

                rho0 = model_rho[i,j,k0]
                rho1 = model_rho[i,j,k1]

                # interpolated model parameters
                vp_val = gamma * vp0 + (1.0 - gamma) * vp1
                vs_val = gamma * vs0 + (1.0 - gamma) * vs1
                rho_val = gamma * rho0 + (1.0 - gamma) * rho1

                # data line format: #x #y #z #vp #vs #density
                output_data.append("{} {} {} {} {} {}\n".format(x_val,y_val,z_val,vp_val,vs_val,rho_val))

                # header stats
                header_origin_x = min(header_origin_x,x_val)
                header_origin_y = min(header_origin_y,y_val)
                header_origin_z = min(header_origin_z,z_val)

                header_end_x = max(header_end_x,x_val)
                header_end_y = max(header_end_y,y_val)
                header_end_z = max(header_end_z,z_val)

                header_vp_min = min(header_vp_min,vp_val)
                header_vp_max = max(header_vp_max,vp_val)

                header_vs_min = min(header_vs_min,vs_val)
                header_vs_max = max(header_vs_max,vs_val)

                header_rho_min = min(header_rho_min,rho_val)
                header_rho_max = max(header_rho_max,rho_val)

                # vtk file output
                if create_vtk_file_output:
                    # adds grid point
                    points.InsertNextPoint(x_val, y_val, z_val)
                    # adds model values
                    vtk_vp.InsertNextValue(vp_val)
                    vtk_vs.InsertNextValue(vs_val)
                    vtk_rho.InsertNextValue(rho_val)
    print("\n")

    # gridding for tomo file
    dx = -10000.0   # 10 km in m (negative sign due to sign convention of X)
    dy = 10000.0
    dz = -1000.0   # 1km in m  (we start with top layer, going downwards -> negative increments)

    # header infos
    header_dx = dx        # increment x-direction for longitudes
    header_dy = dy        # increment y-direction for latitudes
    header_dz = dz        # increment z-direction for depths (positive down)

    # convention for SPECFEM, flip origin_x to right side, end_x to left (as grid starts at top right, going left due to negative sign convention of X)
    tmp_x = header_origin_x
    header_origin_x = header_end_x
    header_end_x = tmp_x

    # convention for SPECFEM, flip origin_z to top, end_z to bottom (as grid starts at top layer, going down)
    tmp_z = header_origin_z
    header_origin_z = header_end_z
    header_end_z = tmp_z

    header_nx = ndim_x   # x-direction
    header_ny = ndim_y   # y-direction
    header_nz = ndim_z   # z-direction

    ## SPECFEM tomographic model format
    ## file header
    data_header = list()
    data_header.append("# tomography model - converted using script convert_Diehletal.py\n")
    data_header.append("#\n")
    # tomography model header infos
    #origin_x #origin_y #origin_z #end_x #end_y #end_z          - start-end dimensions
    data_header.append("#origin_x #origin_y #origin_z #end_x #end_y #end_z\n")
    data_header.append("{} {} {} {} {} {}\n".format(header_origin_x,header_origin_y,header_origin_z,header_end_x,header_end_y,header_end_z))
    #dx #dy #dz                                                 - increments
    data_header.append("#dx #dy #dz\n")
    data_header.append("{} {} {}\n".format(header_dx,header_dy,header_dz))
    #nx #ny #nz                                                 - number of models entries
    data_header.append("#nx #ny #nz\n")
    data_header.append("{} {} {}\n".format(header_nx,header_ny,header_nz))
    #vp_min #vp_max #vs_min #vs_max #density_min #density_max   - min/max stats
    data_header.append("#vp_min #vp_max #vs_min #vs_max #density_min #density_max\n")
    data_header.append("{} {} {} {} {} {}\n".format(header_vp_min,header_vp_max,header_vs_min,header_vs_max,header_rho_min,header_rho_max))
    data_header.append("#\n")

    # data record info
    data_header.append("# data records - format:\n")
    # full format: x #y #z #vp #vs #density
    data_header.append("#UTM_X #UTM_Y #Z #vp #vs #density\n")

    ## writes output file
    filename = data_dir + "/tomography_model.xyz"
    with open(filename, "w") as fp:
        fp.write(''.join(data_header))
        fp.write(''.join(output_data))
        fp.close()

    print("  tomography model written to: ",filename)
    # vtk file output
    if create_vtk_file_output:
        # finish vtk file
        # structured grid file
        # Define grid dimensions
        dimensions = (ndim_x, ndim_y, ndim_z)

        # Create a structured grid
        grid = vtk.vtkStructuredGrid()
        grid.SetDimensions(dimensions)

        grid.SetPoints(points)

        # Add model data to the grid
        grid.GetPointData().AddArray(vtk_vp)
        grid.GetPointData().AddArray(vtk_vs)
        grid.GetPointData().AddArray(vtk_rho)

        # Write grid to a VTK file
        filename = data_dir + "/tomography_model.vts"
        writer = vtk.vtkXMLStructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(grid)
        writer.Write()
        print("          vtk file written to: ",filename)

    # stats
    print("")
    print("tomographic model statistics:")
    print("  vp  min/max : {:.2f} / {:.2f} (m/s)".format(header_vp_min,header_vp_max))
    print("  vs  min/max : {:.2f} / {:.2f} (m/s)".format(header_vs_min,header_vs_max))
    print("  rho min/max : {:.2f} / {:.2f} (kg/m3)".format(header_rho_min,header_rho_max))
    print("")

    return

#----------------------------------------------------------------------------------------

def convert_Diehl_data(data_dir):
    """
    convert model files from Diehl et al. 2021 to SPECFEM3D tomography file format
    """

    print("convert Diehl et al. (2021) model:")
    print("  data dir: ",data_dir)
    print("")

    # downloads model files only if necessary
    download_model_data(data_dir)

    # reads in VP-model
    filename = data_dir + '/' + VP_MODEL_FILE
    data_vp = read_velocity_model_data(filename)

    # reads in VS-model
    filename = data_dir + '/' + VS_MODEL_FILE
    data_vs = read_velocity_model_data(filename)

    # stats
    print(f"read data: Vp min/max = {data_vp.min()} / {data_vp.max()}")
    print(f"           Vs min/max = {data_vs.min()} / {data_vs.max()}")
    print("")

    convert_2_SPECFEM3D_tomography_file(data_vp,data_vs,data_dir)

    print("done")
    print("")
    return

#----------------------------------------------------------------------------------------

def usage():
    print("./convert_Diehletal_to_SPECFEM.py data-directory")
    print("")
    sys.exit(1)

if __name__ == '__main__':
    # arguments
    if len(sys.argv) < 2: usage()

    # gets input type
    data_dir = sys.argv[1]

    convert_Diehl_data(data_dir)

