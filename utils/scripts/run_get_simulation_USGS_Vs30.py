#!/usr/bin/env python
#
# script to setup Vs30 data for a given region
#
# based on USGS global VS30:
#   see: https://earthquake.usgs.gov/data/vs30/
#
#   Metadata for GMT and geotiff downloads:
#   * Grid resolution: 30 arc-seconds (0.0083333333 degrees) ~ 920 m spacing
#   * Latitude span: -56.0 to 84.0 degrees
#   * Longitude span: -180 to 180 degrees
#   * Vs30 units: meters per second
#   * Vs30 range: 98 to 2197 m/s
#   * Vs30 in water-covered areas: 600 m/s
#
from __future__ import print_function

import os
import sys
import subprocess
import math
import datetime
import numpy

## GMT
## http://gmt.soest.hawaii.edu
## install by: > sudo apt install gmt
## setup gmt functions by:
## > source /usr/share/gmt/tools/gmt_functions.sh
try:
    cmd = 'gmt --version'
    print("> ",cmd)
    version = subprocess.check_output(cmd, shell=True)
except:
    print("Error using `gmt`")
    print("install by: > sudo apt install gmt")
    sys.exit(1)
# avoid bytes string issues with strings like b'Hello', converts to text string
if isinstance(version, (bytes, bytearray)): version = version.decode("utf-8")
version = version.strip()
print("GMT version: %s" % (version))
print("")
# get version numbers for later (grdconvert command format changes between version 5.3 and 5.4)
#elem = version.split(".")
#gmt_major = int(elem[0])
#gmt_minor = int(elem[1])

# GMT python interface
# todo: not used yet, calling gmt commands directly as shell commands...
#
#try:
#    import pygmt
#except:
#    print("Error importing module `pygmt`")
#    print("install by: > pip install -U pygmt")
#    sys.exit(1)


#########################################################################
## USER PARAMETERS

## local data directory
datadir = 'USGS_VS30'

# USGS Vs30 dataset filename
vs30_filename = 'global_vs30.grd'

#########################################################################


def check_status(status):
    if status != 0:
        print("error: status returned ",status)
        sys.exit(status)
    return


def download_Vs30():
    """
    downloads USGS Vs30 data file into current directory
    """
    global vs30_filename

    print("getting data:")

    if not os.path.isfile(vs30_filename):
        # download command
        cmd = "wget https://apps.usgs.gov/shakemap_geodata/vs30/" + vs30_filename

        print("downloading global Vs30 data:")
        print("  > ",cmd)
        status = subprocess.call(cmd, shell=True)
        check_status(status)
        print("")

    print("  ok got data")
    print("")


def extract_region(lon_min,lat_min,lon_max,lat_max):
    """
    extracts Vs30 values for given region
    """
    global vs30_filename

    # region format: #lon_min #lat_min #lon_max #lat_max (left bottom right top) in degrees
    # for example: region = (12.35, 41.8, 12.65, 42.0)
    region = (lon_min, lat_min, lon_max, lat_max)

    print("*******************************")
    print("extracting region: ",region)
    print("*******************************")
    print("")


    ## gmt
    # region format: e.g. -R123.0/132.0/31.0/40.0
    gmt_region = '-R' + str(lon_min) + '/' + str(lon_max) + '/' + str(lat_min) + '/' + str(lat_max)

    # USGS Vs30 dataset has a grid resolution of 30 arc-seconds (0.0083333333 degrees)
    # interpolation?
    #incr_dx = 0.009 # fine (~1km)
    #gmt_interval = '-I' + str(incr_dx) + '/' + str(incr_dx)

    print("  GMT:")
    print("  region  : ",gmt_region)
    #print("  interval: ",gmt_interval)
    print("")

    # cut to region
    gridfile = 'region_vs30.grd'

    cmd = 'gmt grdcut ' + vs30_filename + ' ' + gmt_region + ' -G' + gridfile

    cmd = 'gmt grdcut ' + vs30_filename + ' ' + gmt_region + ' -G' + gridfile

    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # info
    cmd = 'gmt grdinfo ' + gridfile
    print("  region file info:")
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # converts to xyz format
    xyz_file = 'region_vs30.xyz'

    cmd = 'gmt grd2xyz ' + gridfile + ' > ' + xyz_file
    print("  converting to xyz")
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("")

    # map
    plot_map(gmt_region,gridfile=gridfile)


def plot_map(gmt_region,gridfile="region_vs30.grd"):
    global datadir

    print("*******************************")
    print("plotting map ...")
    print("*******************************")

    # current directory
    dir = os.getcwd()
    print("  current directory:",dir)
    print("")

    #cmd = 'cd ' + datadir + '/' + ';'
    #status = subprocess.call(cmd, shell=True)
    #check_status(status)

    ps_file = "map.ps"
    pdf_file = "map.pdf"

    # when plotting, add a topography overlay for orientation
    topofile = 'ptopo.grd'
    # srtm 30s (30-arc seconds) ~ 1km resolution
    cmd = 'gmt grdcut @earth_relief_30s ' + gmt_region + ' -G' + topofile + ';'
    # gradient
    cmd += 'gmt grdgradient ' + topofile + ' -Nt1 -A45 -Gptopogradient.grd -V' + ';'
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)

    # gmt plotting
    cmd = 'gmt pscoast ' + gmt_region + ' -JM6i -Dh -G220 -P -K > ' + ps_file + ';'

    # topography shading
    #makecpt -Cgray -T0/1/0.01 > topo.cpt
    #cmd += 'makecpt -Cglobe -T-2500/2500/100 > topo.cpt' + ';'
    #cmd += 'makecpt -Cterra -T-2500/2500/100 > ptopo.cpt' + ';'
    #cmd += 'makecpt -Ctopo -T-2500/2500/100 > ptopo.cpt' + ';'
    #cmd += 'gmt makecpt -Crelief -T-2500/2500/100 > ptopo.cpt' + ';'
    cmd += 'gmt makecpt -Cgray -T-2500/2500/100 > ptopo.cpt' + ';'

    # seismic values
    cmd += 'gmt makecpt -Cseis -T0/2200/100 > pseis.cpt' + ';'

    # w/ topogradient
    cmd += 'gmt grdimage ' + gridfile + ' -Iptopogradient.grd -J -R -Cpseis.cpt -V -O -K >> ' + ps_file + ';'

    # only vs image
    #cmd += 'gmt grdimage ' + gridfile + ' -J -R -Cpseis.cpt -V -O -K >> ' + ps_file + ';'

    # the topography grid
    #cmd += 'gmt grdimage ' + topofile + ' -Iptopogradient.grd -J -R -Cptopo.cpt -V -O -K >> ' + ps_file + ';'

    # coast lines
    cmd += 'gmt pscoast -R -J -Di -N1/1.5p,gray40 -A1000 -W1 -O -K >> ' + ps_file + ';'

    # add a colorbar with legend
    cmd += 'gmt psscale -Dx0.5i/-0.5i+w5i/0.2i+h+e -Cpseis.cpt -Bx500+l"Vs30" -By+l"m/s" -O -K >> ' + ps_file + ';'

    cmd += 'gmt psbasemap -O -R -J -Ba1g1:"Map": -P -V  >> ' + ps_file + ';'

    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("  map plotted in file: ",ps_file)

    # imagemagick converts ps to pdf
    cmd = 'convert ' + ps_file + ' ' + pdf_file
    print("  > ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("  map plotted in file: ",pdf_file)
    print("")
    return


def create_interface_file():
    """
    reads regional .xyz file and creates interface file and information (for reading in SPECFEM3D)
    """
    global datadir

    # region file
    xyz_file = 'region_vs30.xyz'

    # reads in lon/lat/elevation
    print("*******************************")
    print("creating interface:")
    print("*******************************")
    print("  reading file " + xyz_file + " ...")
    print("")

    if not os.path.isfile(xyz_file):
        print("File not found: ",xyz_file)
        print("Please check if extracting a regional XYZ file from the global dataset has succeeded...")
        sys.exit(1)

    # data format: #lon #lat #vs
    data = numpy.loadtxt(xyz_file)

    # x/y start coordinates
    x0 = data[0,0]     # longitude
    y0 = data[0,1]     # latitude

    # get nx/ny counts
    i0 = 1
    nx = 0; ny = 0
    xmin = x0; xmax = x0
    ymin = y0; ymax = y0

    for i in range(1,len(data)):
      x = data[i,0]
      y = data[i,1]
      dx = x - x0
      x0 = x
      # stats
      if x < xmin: xmin = x
      if x > xmax: xmax = x
      if y < ymin: ymin = y
      if y > ymax: ymax = y
      # check regular intervals
      if dx < 0.0:
          ii = i + 1
          if nx > 0 and ii - i0 != nx:
              print("  non-regular nx: ",nx,ii-i0,"on line ",i+1)
          nx = ii - i0
          ny += 1
          deltay = y - y0
          y0 = y
          i0 = ii
      else:
          deltax = dx

    ii = len(data) + 1
    if nx > 0 and ii - i0 != nx:
        print("  non-regular nx: ",nx,ii-i0,"on line ",ii)
    nx = ii - i0
    ny += 1

    print("  number of points along x (NXI) and y (NETA):")
    print("  --------------------------------------------")
    print("  NXI  = ",nx)
    print("  NETA = ",ny)
    print("  xmin/xmax = ",xmin,xmax)
    print("  ymin/ymax = ",ymin,ymax)
    print("  deltax = ",deltax,"average = ",(xmax-xmin)/(nx-1))
    print("  deltay = ",deltay,"average = ",(ymax-ymin)/(ny-1))
    print("  --------------------------------------------")
    print("")

    # corresponding interface definition
    dx = deltax
    dy = deltay

    nxi = nx
    neta = ny

    dxi = dx
    deta = dy

    if dx < 0.0:
        # negative increment, starts from maximum location
        lon = xmax
    else:
        # positive increment, starts from minimum location
        lon = xmin

    if dy < 0.0:
        # negative increment, starts from maximum location
        lat = ymax
    else:
        # positive increment, starts from minimum location
        lat = ymin

    # interface definition
    # same as for DATA/meshfem3D_files/interfaces.dat
    # format:
    #   SUPPRESS_UTM_PROJECTION  NXI  NETA  LONG_MIN  LAT_MIN  SPACING_XI  SPACING_ETA
    line = '.false. ' + str(nxi) + ' ' + str(neta) + ' ' + str(lon) + ' ' + str(lat) + ' ' + str(dxi) + ' ' + str(deta)
    print ("  interface definition:")
    print("  ",line)
    print("")

    # extract only vs data
    vs = data[:,2]
    print("  Vs30 data: min/max = {} / {} (m/s)".format(vs.min(),vs.max()))
    print("")

    # writes out Vs30 interface file
    filename = 'interface_vs30.dat'

    with open(filename,'w') as f:
        # header
        f.write("# Vs30 model interface - converted by script run_get_simulation_USGS_Vs30.py\n")
        f.write("#\n")
        f.write("# providence\n")
        f.write("# created             : {}\n".format(str(datetime.datetime.now())))
        f.write("# command             : {}\n".format(" ".join(sys.argv)))
        f.write("#\n")
        f.write("# interface definition\n")
        f.write("# SUPPRESS_UTM_PROJECTION  NXI  NETA  LONG_MIN  LAT_MIN  SPACING_XI  SPACING_ETA\n")
        f.write("{} {} {} {} {} {} {}\n".format(".false.",nxi,neta,lon,lat,dxi,deta))
        f.write("# data records - Vs30 (m/s)\n")

        # vs data
        for i in range(0,len(vs)):
            f.write("%f\n" % (vs[i]) )

    print("  written: ", "./" + datadir + "/" + filename)
    print("")


def setup_Vs30(lon_min,lat_min,lon_max,lat_max):
    """
    downloads, extracts and creates a Vs30 dataset for the specified region
    """
    global datadir

    print("")
    print("*******************************")
    print("setup USGS Vs30")
    print("*******************************")
    print("")

    # current directory
    dir = os.getcwd()
    print("  current directory: ",dir)

    # creates data directory
    if not os.path.isdir(datadir):
        cmd = 'mkdir -p ' + datadir
        print("  > ",cmd)
        status = subprocess.call(cmd, shell=True)
        check_status(status)
        print("")

    # change working directory to ./USGS_VS30/
    path = dir + '/' + datadir
    os.chdir(path)
    print("  working directory: ",os.getcwd())
    print("")

    # check min/max is given in correct order
    if lat_min > lat_max:
        tmp = lat_min
        lat_min = lat_max
        lat_max = tmp

    if lon_min > lon_max:
        tmp = lon_min
        lon_min = lon_max
        lon_max = tmp

    # use longitude range in [-180,180]
    # (there could be issues with determining the correct UTM zone in routine utm.latlon_to_zone_number() below
    #  if the longitudes are > 180)
    if lon_min < -180.0 and lon_max < -180.0:
        lon_min += 360.0
        lon_max += 360.0
    if lon_min > 180.0 and lon_max > 180:
        lon_min -= 360.0
        lon_max -= 360.0

    # download global data set
    download_Vs30()

    # extract region
    extract_region(lon_min,lat_min,lon_max,lat_max)

    # print interface infos
    create_interface_file()

    print("")
    print("all done")
    print("")



def usage():
    print("usage: ./run_get_simulation_USGS_Vs30.py lon_min lat_min lon_max lat_max")
    print("   where")
    print("       lon_min lat_min lon_max lat_max - region given by points: left bottom right top")
    print("                                         for example: 12.35 42.0 12.65 41.8 (Rome)")


if __name__ == '__main__':
    # checks arguments
    if len(sys.argv) < 5:
        usage()
        sys.exit(1)

    # get arguments
    lon_min = float(sys.argv[1])
    lat_min = float(sys.argv[2])
    lon_max = float(sys.argv[3])
    lat_max = float(sys.argv[4])

    # logging
    cmd = " ".join(sys.argv)
    filename = './run_get_simulation_USGS_Vs30.log'
    with open(filename, 'a') as f:
      print("command call --- " + str(datetime.datetime.now()),file=f)
      print(cmd,file=f)
      print("command logged to file: " + filename)

    # main routine
    setup_Vs30(lon_min,lat_min,lon_max,lat_max)
