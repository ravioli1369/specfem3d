#!/usr/bin/env blender --background --python-use-system-env --python plot_with_blender.py --
#
# or for example on Mac:
#  /usr/bin/env /Applications/Blender.app/Contents/MacOS/Blender --background --python-use-system-env --python plot_blender_sphere.py --
#
#
# run with: > ./plot_with_blender.py --help
#
###############################################################################################

import sys
import os
import time
import datetime
import math

# blender
try:
    import bpy
except:
    print("Failed importing bpy. Please make sure to have Blender python working properly.")
    sys.exit(1)

print("")
print("Blender: version ",bpy.app.version_string)
print("")

# blender meshing
import bmesh

# adds path
#sys.path.append('/Applications/Blender.app/Contents/Resources/3.6/python/lib/python3.10/site-packages')

try:
    import vtk
except:
    print("Failed importing vtk. Please make sure to have Blender python with vtk working properly.")
    print("system path:")
    for path in sys.path:
        print("  ",path)
    sys.exit(1)

print("VTK: version ",vtk.vtkVersion().GetVTKVersion())
print("")

try:
    import numpy as np
except:
    print("Failed importing numpy. Please make sure to have Blender python with numpy working properly.")
    sys.exit(1)

from mathutils import Vector,Matrix
from math import radians,atan2,sin,cos

###############################################################################################
## USER parameters

## renderer
# full image size
blender_img_resolution_X = 2400
blender_img_resolution_Y = 1600

###############################################################################################

# Constants
PI = 3.141592653589793
DEGREE_TO_RAD = PI / 180.0

# Global variables
mesh_scale_factor = 1.0
mesh_origin = [0.0,0.0,0.0]

# vertical exaggeration factor
vertical_exaggeration = None

## sea-level plane
use_sea_level_plane = True

## transparent sea-level plane
use_transparent_sea_level_plane = False

# shallow coast lines can lead to z-buffering issues with the sea-level plane.
# as a work-around, one can try to add a boolean modifier to cut the plane with the main mesh object.
# this might however lead to other issues... left here to give a try.
add_sea_level_plane_modifier = False

# shifts point above/below sea-level slightly up/down to get a sharper and steeper coastline for sea-level plane
sea_level_separation = None

# baseline shift for buildings (to move above SPECFEM mesh)
shift_building_baseline = 0.0  # 0.0==no-shift default, or e.g. --shift-building-baseline=0.0005

# location labels (white or black)
location_labels_color = (0.05,0.05,0.05,1)  # white (0.9,0.9,0.9,1), black (0.05,0.05,0.05,1)

# rendering output
suppress_renderer_output = False

# cycles rendering (for better glass effect of buildings)
use_cycles_renderer = False

# animation
use_animation_dive_in = True
use_animation_rotation = True

# camera positioning
camera_elevation_offset = 0.15
camera_y_offset = -0.7
camera_angle_depth_degree = 75.0

# camera perspective
close_up_view = False
centered_view = False

# center point to focus camera view on
x_center = 0.0
y_center = 0.0
z_center = 0.0

# elevation at center point
z_elevation = 0.0

# background color
world_background_color = (1,1,1,1)  # white

# utm projections for locations
transformer_to_utm = None
utm_zone = None
use_moon_ltm = False

# timing info
time_string = None

# material appearance
use_matte_material = False

# class to avoid long stdout output by renderer
# see: https://stackoverflow.com/questions/24277488/in-python-how-to-capture-the-stdout-from-a-c-shared-library-to-a-variable/29834357
class SuppressStream(object):
    def __init__(self, stream=sys.stderr,suppress=False):
        # turns on/off suppressing stdout of renderer process
        self.SUPPRESS_STDOUT = suppress

        if self.SUPPRESS_STDOUT:
            self.orig_stream_fileno = stream.fileno()

    def __enter__(self):
        if self.SUPPRESS_STDOUT:
            self.orig_stream_dup = os.dup(self.orig_stream_fileno)
            self.devnull = open(os.devnull, 'w')
            os.dup2(self.devnull.fileno(), self.orig_stream_fileno)

    def __exit__(self, type, value, traceback):
        if self.SUPPRESS_STDOUT:
            os.close(self.orig_stream_fileno)
            os.dup2(self.orig_stream_dup, self.orig_stream_fileno)
            os.close(self.orig_stream_dup)
            self.devnull.close()


def convert_latlon_to_UTM(lat, lon):
    global transformer_to_utm
    global use_moon_ltm
    global utm_zone

    if use_moon_ltm:
        print("  converting coordinates to Moon LTM/LPS...")
        # determine zone if not specified
        if utm_zone == None:
            utm_zone = geo2ltm(lon, lat, zone=None, iway=2)
            hemisphere_auto = 'N' if utm_zone > 0 else 'S'
            print("  detected LTM/LPS zone: {}   - {} hemisphere".format(utm_zone,hemisphere_auto))
            print("")

        # convert to LTM/LPS
        ltm_x, ltm_y = geo2ltm(lon, lat, utm_zone, iway=2)
        return ltm_x, ltm_y

    # transform lat/lon to UTM
    if transformer_to_utm == None:
        # Coordinate projections
        try:
            import pyproj
        except:
            print("Failed importing pyproj.")
            sys.exit(1)
        from pyproj import Transformer

        print("  converting coordinates to UTM...")
        print("")
        # pyproj coordinate system info:
        #   WGS84                                          ==       EPSG:4326
        #   spherical mercator, google maps, openstreetmap ==       EPSG:3857
        #
        # we first need to determine the correct UTM zone to get the EPSG code.
        # for this, we take the lat/lon position and query the corresponding UTM zone for this position.
        ref_epsg = "EPSG:4326"

        # user specified UTM zone
        if not utm_zone == None:
            # Determine the hemisphere based on the zone number based on zone number 1 - 120
            #hemisphere = 'N' if utm_zone <= 60 else 'S'
            #utm_zone = utm_zone % 60  # Normalize zone number to 1-60

            # here, we will use 1-60 as input for UTM zones, positive for Northern and negative numbers for Southern hemisphere
            hemisphere = 'N' if utm_zone > 0 else 'S'

            # Construct the EPSG code
            utm_code = 32600 + abs(utm_zone) if hemisphere == 'N' else 32700 + abs(utm_zone)
            utm_epsg = "EPSG:{}".format(utm_code)

            print("  user specified UTM zone: ",utm_zone)
            print("                 UTM code: ",utm_code," epsg: ", utm_epsg)
            print("")
        else:
            # gets list of UTM codes
            utm_crs_list = pyproj.database.query_utm_crs_info(
                              datum_name="WGS 84",
                              area_of_interest=pyproj.aoi.AreaOfInterest(west_lon_degree=lon,
                                                                         south_lat_degree=lat,
                                                                         east_lon_degree=lon,
                                                                         north_lat_degree=lat))
            utm_code = utm_crs_list[0].code
            utm_epsg = "EPSG:{}".format(utm_code)

            # convert code to integer number and determine UTM zone number for info
            utm_code = int(utm_code)
            hemisphere_auto = 'N' if utm_code < 32700 else 'S'
            utm_zone_auto = utm_code - 32600 if hemisphere_auto == 'N' else utm_code - 32700

            print("  detected UTM code: ",utm_code," epsg: ", utm_epsg)
            print("           UTM zone: {}{}".format(utm_zone_auto,hemisphere_auto))
            print("")

        # transformer
        # transformation from WGS84 to UTM zone
        # WGS84: Transformer.from_crs("EPSG:4326", utm_epsg)
        #transformer_to_utm = Transformer.from_crs(ref_epsg, utm_epsg)                 # input: lat/lon -> utm_x,utm_y
        transformer_to_utm = Transformer.from_crs(ref_epsg, utm_epsg, always_xy=True) # input: lon/lat -> utm_x/utm_y

        #debug
        #print(transformer_to_utm)
        #print("debug: lon/lat ",transformer_to_utm.transform(orig_lon,orig_lat))
        #print("debug: lat/lon ",transformer_to_utm.transform(orig_lat,orig_lon))

        # user info
        utm_x,utm_y = transformer_to_utm.transform(lon,lat)
        print("       -> UTM x/y  = ",utm_x,utm_y)
        print("          backward check: orig x/y = ",transformer_to_utm.transform(utm_x,utm_y,direction='INVERSE'))
        print("")

    # converts point coordinates
    x = lon
    y = lat

    # converts to UTM location (x and y)
    x_utm,y_utm = transformer_to_utm.transform(x,y)

    return x_utm,y_utm


def geo2ltm(lon, lat, zone=None, iway=2):
    """
    Lunar Transverse Mercator (LTM) and Lunar Polar Stereographic (LPS) projection for polar regions

    this routine assumes a perfectly spherical Moon.
    LTM implementation is a simplified Spherical Transverse Mercator, 8-degree zones.

    note: by default, the LTM implementation in the python script LGRS_Coordinate_Conversion_mk7.2.py provided by USGS astrogeology site
          https://astrogeology.usgs.gov/search/map/lunar-map-projections-and-grid-reference-system-for-artemis-astronaut-surface-navigation
          uses a spherical Moon for their LTM projection.

          However, they still use the formula for a Gauss-Schreiber projection,
          that combines a projection of an ellipsoid to a sphere followed by a spherical transverse Mercator formula.
          Assuming the ellipsoid is perfectly spherical, then the first projection is an identity transform and only the second,
          spherical transverse Mercator formula is needed.
          Here, we use this simplification of applying directly the spherical transverse Mercator projection, assuming a spherical Moon.

    lon - longitude in degree range [-180,180]
    lat - latitude in degree range [-90,90]

    zone - LTM zones use 1 to 45, positive for Northern hemisphere, negative for Southern hemisphere
           LPS zones use 46 for North pole, -46 for South pole

    iway = 1  from LTM     to lon/lat
    iway = 2  from lon/lat to LTM
    """

    PI = math.pi
    degrad = PI / 180.0
    raddeg = 180.0 / PI

    # ---- Lunar spherical radius ----
    R = 1737400.0   # mean lunar radius (meters)

    # ---- Scale factor ----
    # Lunar map projection scale factors
    scfa = 0.999         # transverse Mercator (LTM)
    scfa_polar = 0.994   # polar stereographic (LPS)

    # ---- False origins ----
    false_east  =  250000.0
    false_north = 2500000.0

    false_east_polar  =  500000.0
    false_north_polar =  500000.0

    ILTM2LONGLAT = 1
    ILONGLAT2LTM = 2

    #----- Set Zone parameters
    # determine zone
    # for convenience, if zone is unspecified in forward mode, this computes it for the given longitude/latitude position
    # and returns only the zone number
    if iway == ILONGLAT2LTM and (zone is None or zone == 0):
        # checks latitute range for LTM [-82,82], beyond is polar region
        if lat >= -82.0 and lat <= 82.0:
            # LTM
            # longitudinal zone
            zone = int((lon + 180.0) // 8) + 1   # note: uses +1 because USGS's LGRS_Coordinate_Conversion_mk7.2.py has zone index starting at 1
            # 180-degree values sometimes come up as 46. We assign it back to zone 1
            if zone > 45:
                zone -= 45
            # we use negative values for Southern hemisphere
            if lat < 0.0:
                zone = -zone
        else:
            # LPS
            # polar region
            if lat >= 0.0:
                zone = 46    # north pole
            else:
                zone = -46   # south pole
        # just return zone
        return zone

    # zone is given as input, check if valid
    if zone is None or zone == 0 or int(abs(zone)) > 46:
        print(f"error: geo2ltm routine has as input zone {zone}, which is invalid. zone must be +/- 1-45 for LTM and +/- 46 for LPS")
        sys.exit(1)

    # zone index absolute
    z = int(abs(zone))

    # polar region
    use_polar = False
    # check if LPS or LTM
    if z == 46: use_polar = True

    # Lunar Polar Stereographic (LPS) projection
    if use_polar:
        if iway == ILONGLAT2LTM:
            # Forward transformation: lon/lat to LPS
            # Snyder (1987) spherical equations
            rlon = lon * degrad
            rlat = lat * degrad

            # Calculate polar stereographic spherical scale error
            A = 2.0 * R * scfa_polar

            if zone < 0:
                # South pole
                t = math.tan(PI/4.0 + rlat/2.0)
            else:
                # North pole
                t = math.tan(PI/4.0 - rlat/2.0)

            # stereographic map projection, Snyder (1987)
            # note: here, we use +cos(lon) as done in routine spherical_stereographic_map_y() of the USGS script LGRS_Coordinate_Conversion_mk7.2.py
            #       for both, North and South poles. the standard polar stereographic projection defined by Snyder would use
            #       North Pole: x = 2 R k0 tan(pi/4 - phi/2) sin(lambda - lambda0)
            #                   y = - 2 R k0 tan(pi/4 - phi/2) cos(lambda - lambda0)
            #       South Pole: x = 2 R k0 tan(pi/4 + phi/2) sin(lambda - lambda0)
            #                   y = 2 R k0 tan(pi/4 + phi/2) cos(lambda - lambda0)
            #       (see Snyder 1987, page 158, chapter 21 "Stereographic Projection", section "Formulas for the Sphere",
            #        eqs. 21-5, 21-6 and 21-9,21-10.
            #        https://pubs.usgs.gov/pp/1395/report.pdf ).
            # we'll take the convention from the USGS implementation to use +cos for both poles.
            #
            # X coordinate for a s
            x = A * t * math.sin(rlon)
            # Y coordinate for a stereographic map projection
            y = A * t * math.cos(rlon)

            # Add false Eastings and Northings
            x += false_east_polar
            y += false_north_polar

            # all done
            return x, y

        else:
            # Inverse transformation: LPS t0 lon/lat
            # remove false origins (must be same values used in forward)
            x = lon - false_east_polar
            y = lat - false_north_polar

            # radial distance from projection origin
            rho = math.hypot(x, y)

            # at the pole: rho == 0
            if rho == 0.0:
                deglat = -90.0 if zone < 0 else 90.0
                deglon = 0.0
                return deglon, deglat

            A = 2.0 * R * scfa_polar

            # t = tan( PI/4 ± lat/2 )  where sign depends on hemisphere in forward
            t = rho / A

            if zone < 0:
                # South pole
                # forward used: tan(PI/4 + lat/2) = t  -> lat = 2*atan(t) - PI/2
                latr = 2.0 * math.atan(t) - PI/2.0
            else:
                # North pole
                # forward used: tan(PI/4 - lat/2) = t  -> lat = PI/2 - 2*atan(t)
                latr = PI/2.0 - 2.0 * math.atan(t)

            # longitude from sin/cos ordering used in forward: x = factor*sin(lon), y = factor*cos(lon)
            lonr = math.atan2(x, y)

            deglon = lonr * raddeg
            deglat = latr * raddeg

            # normalize lon to -180..180
            if deglon > 180.0:
                deglon -= 360.0
            if deglon < -180.0:
                deglon += 360.0

            return deglon, deglat


    # Lunar Transverse Mercator (LTM) Projection
    # ---- Central meridian of 8-degree zone ----
    # central meridian
    cm = -180.0 + ((z - 1) + 0.5) * 8.0      # longitude of central meridian (note: z-1 because zones start at 1)
    cmr = cm * degrad

    # ----------------------------------------------------------------------
    # Forward transformation: lon/lat to LTM
    # ----------------------------------------------------------------------
    if iway == ILONGLAT2LTM:
        rlon = lon * degrad
        rlat = lat * degrad

        dlam = rlon - cmr

        # ---- Spherical Transverse Mercator (forward) ----
        B = math.cos(rlat) * math.sin(dlam)

        x = 0.5 * R * scfa * math.log((1.0 + B) / (1.0 - B))
        y = R * scfa * math.atan2(math.tan(rlat), math.cos(dlam))

        # False origins
        x += false_east
        # only for South sections
        if zone < 0:
            y += false_north

        return x, y

    # ----------------------------------------------------------------------
    # Inverse transformation: LTM to lon/lat
    # ----------------------------------------------------------------------
    else:
        # remove false origins
        x = lon - false_east
        if zone < 0:
            # only for South sections
            y = lat - false_north
        else:
            y = lat

        # ---- Spherical TM inverse ----
        D = y / (R * scfa)
        T = x / (R * scfa)

        lonr = cmr + math.atan2(math.sinh(T), math.cos(D))
        latr = math.asin(math.sin(D) / math.cosh(T))

        deglon = lonr * raddeg
        deglat = latr * raddeg

        # normalize lon to -180..180
        if deglon > 180.0:
            deglon -= 360.0
        if deglon < -180.0:
            deglon += 360.0

        return deglon, deglat


def setup_color_table(colormap, data_min, data_max):
    # creates a color table
    lut = vtk.vtkLookupTable()

    lut.SetTableRange(data_min, data_max)
    lut.SetNumberOfTableValues(256) # Set the number of table values

    #lut.SetRampToLinear()
    #lut.SetRampToSQRT()

    # VTK by default maps the value range to colors from red to blue
    # determine custom type
    if colormap == 0:
        print("  color map: default VTK")
        # nothing to special to add, let's just vtk internally do it
        # colormap is going from red to white to blue
    elif colormap == 1:
        # topo
        print("  color map: topo")
        colors_rgb = [
            [0.3,  0.3,   0.3], # gray
            [0.1,  0.1,   0.4], # blue
            [0.2,  0.5,   0.2],
            [0.25, 0.625, 0.5],
            [0.0,  0.5,   0.25],
            [0.5,  0.365, 0.0],
            [0.75, 0.625, 0.25],
            [1.0,  0.75,  0.625],
            [1.0,  0.75,  0.5],
            [1,    1,     1],   # white
        ]

    elif colormap == 2:
        # Scientific Colour Map Categorical Palette
        # https://www.fabiocrameri.ch/colourmaps/
        # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
        print("  color map: lisbon")
        # lisbon 10 Swatches
        colors_rgb255 = [
            [230, 229, 255],  #  lisbon-1 #E6E5FF
            # or start with a less white color
            #[200, 208, 237], #  lisbon-12 #C8D0ED
            [155, 175, 211],  #  lisbon-29 #9BAFD3
            [ 81, 119, 164],  #  lisbon-58 #5177A4
            [ 30,  67, 104],  #  lisbon-86 #1E4368
            [ 17,  30,  44],  #  lisbon-114 #111E2C
            [ 39,  37,  26],  #  lisbon-143 #27251A
            [ 87,  81,  52],  #  lisbon-171 #575134
            [141, 133,  86],  #  lisbon-199 #8D8556
            [201, 195, 144],  #  lisbon-228 #C9C390
            [255, 255, 217],  #  lisbon-256 #FFFFD9
        ]
        # converts the colors from 0-255 range to 0-1 range
        colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

    elif colormap == 3:
        # Scientific Colour Map Categorical Palette
        # https://www.fabiocrameri.ch/colourmaps/
        # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
        print("  color map: lajolla")
        # lajolla 10 Swatches
        colors_rgb255 = [
            [ 25,  25,   0], #  lajolla-1 #191900
            [ 51,  34,  15], # lajolla-29 #33220F
            [ 91,  48,  35], #  lajolla-58 #5B3023
            [143,  64,  61], #  lajolla-86 #8F403D
            [199,  80,  75], #  lajolla-114 #C7504B
            [224, 114,  79], #  lajolla-143 #E0724F
            [231, 148,  82], #  lajolla-171 #E79452
            [238, 181,  85], #  lajolla-199 #EEB555
            [248, 223, 124], #  lajolla-228 #F8DF7C
            [255, 254, 203], #  lajolla-256 #FFFECB
        ]
        # converts the colors from 0-255 range to 0-1 range
        colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

    elif colormap == 4:
        # Scientific Colour Map Categorical Palette
        # https://www.fabiocrameri.ch/colourmaps/
        # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
        print("  color map: lipari")
        # lipari 10 Swatches
        colors_rgb255 = [
            [  3,  19,  38], #  lipari-1 #031326
            [ 19,  56,  90], #  lipari-29 #13385A
            [ 71,  88, 122], #  lipari-58 #47587A
            [107,  95, 118], #  lipari-86 #6B5F76
            [142,  97, 108], #  lipari-114 #8E616C
            [188, 100,  97], #  lipari-143 #BC6461
            [229, 123,  98], #  lipari-171 #E57B62
            [231, 162, 121], #  lipari-199 #E7A279
            [233, 201, 159], #  lipari-228 #E9C99F
            [253, 245, 218], #  lipari-256 #FDF5DA
        ]
        # converts the colors from 0-255 range to 0-1 range
        colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

    elif colormap == 5:
        # Scientific Colour Map Categorical Palette
        # https://www.fabiocrameri.ch/colourmaps/
        # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
        print("  color map: davos")
        # davos 10 Swatches
        colors_rgb255 = [
            [  0,   5,  74], #  davos-1 #00054A
            [ 17,  44, 113], #  davos-29 #112C71
            [ 41,  82, 145], #  davos-58 #295291
            [ 67, 112, 157], #  davos-86 #43709D
            [ 94, 133, 152], #  davos-114 #5E8598
            [121, 150, 141], #  davos-143 #79968D
            [153, 173, 136], #  davos-171 #99AD88
            [201, 210, 158], #  davos-199 #C9D29E
            [243, 243, 210], #  davos-228 #F3F3D2
            [254, 254, 254], #  davos-256 #FEFEFE
        ]
        # converts the colors from 0-255 range to 0-1 range
        colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

    elif colormap == 6:
        # Scientific Colour Map Categorical Palette
        # https://www.fabiocrameri.ch/colourmaps/
        # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
        print("  color map: turku")
        # turku 10 Swatches
        colors_rgb255 = [
            [  0,   0,   0], #  turku-1 #000000
            [ 36,  36,  32], #  turku-29 #242420
            [ 66,  66,  53], #  turku-58 #424235
            [ 95,  95,  68], #  turku-86 #5F5F44
            [126, 124,  82], #  turku-114 #7E7C52
            [169, 153, 101], #  turku-143 #A99965
            [207, 166, 124], #  turku-171 #CFA67C
            [234, 173, 152], #  turku-199 #EAAD98
            [252, 199, 195], #  turku-228 #FCC7C3
            [255, 230, 230], #  turku-256 #FFE6E6
        ]
        # converts the colors from 0-255 range to 0-1 range
        colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

    elif colormap == 7:
        # Scientific Colour Map Categorical Palette
        # https://www.fabiocrameri.ch/colourmaps/
        # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
        print("  color map: berlin")
        # berlin 10 Swatches
        colors_rgb255 = [
            [158, 176, 255], #  berlin-1 #9EB0FF
            [ 91, 164, 219], #  berlin-29 #5BA4DB
            [ 45, 117, 151], #  berlin-58 #2D7597
            [ 26,  66,  86], #  berlin-86 #1A4256
            [ 17,  25,  30], #  berlin-114 #11191E
            [ 40,  13,   1], #  berlin-143 #280D01
            [ 80,  24,   3], #  berlin-171 #501803
            [138,  63,  42], #  berlin-199 #8A3F2A
            [196, 117, 106], #  berlin-228 #C4756A
            [255, 173, 173], #  berlin-256 #FFADAD
        ]
        # converts the colors from 0-255 range to 0-1 range
        colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

    elif colormap == 8:
        # Scientific Colour Map Categorical Palette
        # https://www.fabiocrameri.ch/colourmaps/
        # Crameri, F. (2018). Scientific colour maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
        print("  color map: grayC")
        colors_rgb255 = [
            [0,   0,   0],  #  grayC-1 #000000
            [35,  35,  35], #  grayC-29 #232323
            [61,  61,  61], #  grayC-58 #3D3D3D
            [86,  86,  86], #  grayC-86 #565656
            [108, 108, 108],#  grayC-114 #6C6C6C
            [130, 130, 130],#  grayC-143 #828282
            [154, 154, 154],#  grayC-171 #9A9A9A
            [182, 182, 182],#  grayC-199 #B6B6B6
            [216, 216, 216],#  grayC-228 #D8D8D8
            [255, 255, 255],#  grayC-256 #FFFFFF
        ]
        # converts the colors from 0-255 range to 0-1 range
        colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

    elif colormap == 9:
        # custom snow
        print("  color map: snow")
        colors_rgb255 = [
            [204, 204, 204], # gray
            [153, 178, 204],
            [ 71,  88, 122], #  lipari-58 #47587A
            [107,  95, 118], #  lipari-86 #6B5F76
            [142,  97, 108], #  lipari-114 #8E616C
            [188, 100,  97], #  lipari-143 #BC6461
            [229, 123,  98], #  lipari-171 #E57B62
            [231, 162, 121], #  lipari-199 #E7A279
            [255, 229, 204],
            [255, 255, 255], # white
        ]
        # converts the colors from 0-255 range to 0-1 range
        colors_rgb = [[comp / 255.0 for comp in color] for color in colors_rgb255]

    elif colormap == 10:
        # custom shakeGreen
        print("  color map: shakeGreen")
        colors_rgb = [
            [0.8,  0.8,   0.8], # gray
            [0.5,  0.5,   0.4],
            [0.5,  0.4,   0.2],
            [0.6,  0.6,   0.0], # green
            [0.72, 0.25,  0.0 ], # orange
            [0.81, 0.5,   0.0 ],
            [0.9,  0.74,  0.0 ], # yellow
            [1.0,  0.99,  0.0 ],
            [1.0,  0.99,  0.25],
            [1.0,  1.0,   1.0 ],  # white
        ]

    elif colormap == 11:
        # custom shakeRed
        print("  color map: shakeRed")
        colors_rgb = [
            [0.85, 0.85,  0.85], # gray
            [0.7,  0.7,   0.7 ],
            [0.5,  0.5,   0.5 ],
            [0.63, 0.0,   0.0 ], # red
            [0.72, 0.25,  0.0 ], # orange
            [0.81, 0.5,   0.0 ],
            [0.9,  0.74,  0.0 ], # yellow
            [1.0,  0.99,  0.0 ],
            [1.0,  0.99,  0.25],
            [1.0,  1.0,   1.0 ],  # white
        ]

    elif colormap == 12:
        # custom shakeUSGS
        # taken from a shakemap plot of the USGS
        # https://earthquake.usgs.gov/earthquakes/eventpage/us6000lqf9/shakemap/intensity
        print("  color map: shakeUSGS")
        colors_rgb = [
            [1.0,  1.0,  1.0 ], # I: white
            [0.8,  0.8,  0.8 ],
            [0.77, 0.81, 1.0 ], # II-III : light purple
            [0.5,  1.0,  0.98], # IV: turquoise
            [0.5,  1.0,  0.54], # V : green
            [1.0,  0.98, 0.0 ], # VI : yellow
            [1.0,  0.77, 0.0 ], # VII : orange
            [0.99, 0.52, 0.0 ], # VIII: darkorange
            [0.98, 0.0,  0.0 ], # IX : red
            [0.78, 0.0,  0.0 ], # X+ : darkred
        ]

    elif colormap == 13:
        # custom shakeUSGSgray
        # starts with darker gray than the default USGS
        print("  color map: shakeUSGSgray")
        colors_rgb = [
            [0.8,  0.8,  0.8 ], # gray
            [0.5,  0.5,  0.5 ],
            [0.77, 0.81, 1.0 ], # II-III : light purple
            [0.5,  1.0,  0.98], # IV: turquoise
            [0.5,  1.0,  0.54], # V : green
            [1.0,  0.98, 0.0 ], # VI : yellow
            [1.0,  0.77, 0.0 ], # VII : orange
            [0.99, 0.52, 0.0 ], # VIII: darkorange
            [0.98, 0.0,  0.0 ],  # IX : red
            [0.5, 0.0,  0.0 ],  # X+ : darkred
        ]

    elif colormap == 14:
        # custom shakeUSGSblack
        # starts with darker gray than the default USGS
        print("  color map: shakeUSGSblack")
        colors_rgb = [
            [0.15, 0.15, 0.15], # black
            [0.3,  0.3,  0.35],
            [0.77, 0.81, 1.0 ], # II-III : light purple
            [0.5,  1.0,  0.98], # IV: turquoise
            [0.5,  1.0,  0.54], # V : green
            [1.0,  0.98, 0.0 ], # VI : yellow
            [1.0,  0.77, 0.0 ], # VII : orange
            [0.99, 0.52, 0.0 ], # VIII: darkorange
            [0.98, 0.0,  0.0 ],  # IX : red
            [0.5, 0.0,  0.0 ],  # X+ : darkred
        ]

    elif colormap == 15:
        # custom shake
        # starts with dark colors than more light
        print("  color map: shakeDark")
        colors_rgb = [
            [0.15, 0.15, 0.15], # black
            [0.3,  0.3,  0.35], #
            [0.5,  0.5,  0.7],  # II-III
            [0.75, 0.75, 0.75], # IV
            [1.0,  1.0,  1.0],  # V       white
            [1.0,  0.98, 0.0 ], # VI      yellow
            [1.0,  0.77, 0.0 ], # VII     orange
            [0.99, 0.52, 0.0 ], # VIII    darkorange
            [0.98, 0.0,  0.0 ], # IX      red
            [0.7, 0.0,  0.0 ],  # X+      darkred
        ]

    elif colormap == 16:
        # custom shake
        # starts with dark colors than more light - similar to gist_earth, perceptually uniform
        print("  color map: gist_earth")
        colors_rgb = [
            [0.0,  0.06,  0.08],  # dark
            #[ 0.00, 0.13, 0.45 ],
            [ 0.00, 0.18, 0.40 ],
            [ 0.00, 0.25, 0.32 ],
            [ 0.00, 0.33, 0.25 ],
            [ 0.00, 0.41, 0.18 ],
            [ 0.02, 0.49, 0.15 ],
            [ 0.40, 0.56, 0.18 ],
            [ 0.63, 0.62, 0.28 ],
            [ 0.81, 0.68, 0.42 ],
            [ 0.96, 0.75, 0.59 ],
            [ 1.00, 0.84, 0.77 ],
            #[ 1.00, 0.95, 0.95 ],
            [ 1.00, 1.0, 1.0 ], # white
        ]

    else:
        print("Warning: colormap with type {} is not supported, exiting...".format(colormap))
        sys.exit(1)

    # sets lookup table entries
    if colormap != 0:
        # Create a vtkColorTransferFunction
        color_transfer_func = vtk.vtkColorTransferFunction()
        # add specific scalar values in the color transfer function
        for i, color in enumerate(colors_rgb):
            val = i / (len(colors_rgb) - 1.0)
            color_transfer_func.AddRGBPoint(val, color[0], color[1], color[2])
        # Calculate the color values for the lookup table by interpolating from the color transfer function
        for i in range(256):
            scalar = i / 255.0  # Normalized scalar value from 0 to 1
            color = color_transfer_func.GetColor(scalar)
            lut.SetTableValue(i, color[0], color[1], color[2], 1.0)
        print("")

    # build lookup table
    lut.Build()

    return lut


def convert_vtk_to_obj(vtk_file: str="", colormap: int=0, color_max=None) -> str:
    global mesh_scale_factor,mesh_origin
    global vertical_exaggeration,sea_level_separation

    # Path to your .vtu file
    print("converting vtk file: ",vtk_file)

    # check file
    if len(vtk_file) == 0:
        print("Error: no vtk file specified...")
        usage()
        sys.exit(1)

    if not os.path.exists(vtk_file):
        print("Error: vtk file specified not found...")
        sys.exit(1)

    # gets file extension
    extension = os.path.splitext(vtk_file)[1]
    # reads the vtk file
    if extension == '.vtk':
        # .vtk file
        reader = vtk.vtkDataSetReader()
        reader.SetFileName(vtk_file)
        reader.Update()
    elif extension == '.vtu':
        # .vtu
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtk_file)
        reader.Update()
    elif extension == '.inp':
        # AVS .inp
        #reader = vtk.vtkSimplePointsReader()
        reader = vtk.vtkAVSucdReader()
        reader.SetFileName(vtk_file)
        reader.Update()
    else:
        print("unknown vtk file extension ",extension," - exiting...")
        sys.exit(1)

    #debug
    #print(reader)

    ## scale coordinates to be in range [-1,1] for visualization
    # gets the points (vertices) from the dataset
    points = reader.GetOutput().GetPoints()

    # number of points
    num_points = points.GetNumberOfPoints()

    #print(points)
    #print(points.GetBounds())
    xmin,xmax,ymin,ymax,zmin,zmax = points.GetBounds()

    # defines the scaling factors to fit within +/- 1
    #min_coords = np.array(points.GetPoint(0))
    #max_coords = np.array(points.GetPoint(0))
    #for i in range(1, num_points):
    #    point = np.array(points.GetPoint(i))
    #    min_coords = np.minimum(min_coords, point)
    #    max_coords = np.maximum(max_coords, point)

    min_coords = np.array([xmin,ymin,zmin])
    max_coords = np.array([xmax,ymax,zmax])
    dimensions = max_coords - min_coords

    # info
    print("  minimum coordinates:", min_coords)
    print("  maximum coordinates:", max_coords)
    print("  dimensions         :", dimensions)
    print("")

    # gets data values on nodes
    data_array = None
    if reader.GetOutput().GetPointData().GetNumberOfArrays() > 0:
        data_array = reader.GetOutput().GetPointData().GetArray(0)  # Example: Accessing the first data array
        #debug
        #print(data_array)
        #info
        print("  data: ")
        print("  range = ",data_array.GetRange())
        print("")


    # determines origin of mesh to move it back to (0,0,0) and scale it between [-1,1] to better locating it in blender
    mesh_origin = min_coords + 0.5 * (max_coords - min_coords)

    # z-coordinate: leaves sea level at 0 for mesh origin
    mesh_origin[2] = 0.0

    # takes maximum size in x/y direction
    dim_max = np.maximum(dimensions[0],dimensions[1])

    if np.abs(dim_max) > 0.0:
        mesh_scale_factor = 2.0 / dim_max

    print("  mesh scaling:")
    print("  origin       :",mesh_origin)
    print("  scale factor :",mesh_scale_factor)
    print("")

    if vertical_exaggeration != None:
        print("  using vertical exaggeration factor: ",vertical_exaggeration)
        print("")

    if sea_level_separation != None:
        print("  using sea-level separation shift: ",sea_level_separation)
        print("")

    # creates an array to store scaled points
    scaled_points = vtk.vtkPoints()

    # Loop through each point, scale its coordinates, and add it to the new points array
    for i in range(num_points):
        # translation by origin
        point = np.array(points.GetPoint(i)) - mesh_origin

        # uniform scaling
        scaled_point = point * mesh_scale_factor

        # vertical exaggeration
        if vertical_exaggeration != None:
            scaled_point[2] = scaled_point[2] * vertical_exaggeration

        # shift points above/below sea-level
        if sea_level_separation != None:
            # shift points above
            if scaled_point[2] > 0.0: scaled_point[2] += sea_level_separation
            if scaled_point[2] < 0.0: scaled_point[2] -= sea_level_separation

        # stores updated points
        scaled_points.InsertNextPoint(scaled_point)

    # creates a new polydata with the scaled points
    scaled_polydata = vtk.vtkPolyData()
    scaled_polydata.SetPoints(scaled_points)

    # vertex data
    if data_array:
        num_points = data_array.GetNumberOfTuples()
        num_components = data_array.GetNumberOfComponents()
        data_min,data_max = data_array.GetRange()

        print("  data array: ")
        print("    number of points     = ",num_points)
        print("    number of components = ",num_components)
        print("    range min/max        = ",data_min,"/",data_max)
        print("")

        # checks if fixing maximum value
        if color_max:
            # Convert VTK data array to NumPy array
            array = np.zeros((num_points, num_components))
            for i in range(num_points):
                for j in range(num_components):
                    value = data_array.GetComponent(i, j)
                    array[i,j] = value

            # limit size
            if 1 == 0:
                ## determines maximum value as a multiple of 10
                # reduce first by 10%
                total_max = abs(array).max()
                total_max = 0.9 * total_max  # sets limit at 90% of the maximum

                # get maximum value in power of 10
                if total_max != 0.0:
                    total_max = 1.0 * 10**(int(np.log10(total_max)))  # example: 2.73e-11 limits to 1.e-11
                    #total_max = 1.0 * 10**(int(np.log10(total_max))-1)  # example: 2.73e-11 limits to 1.e-11
                    #total_max = 1.0 * 10**(int(np.log10(total_max))-2)  # example: 2.73e-11 limits to 1.e-12
                else:
                    total_max = 0.0
                    print("  data: color data min/max   = ",array.min(),array.max())
                    print("  data: zero color data - nothing to show")
                    # nothing left to do
                    #sys.exit(1)

            # checks if fixing maximum value
            if color_max:
                total_max = color_max

            print("  limiting color data range:")
            print("    data min/max   = ",array.min(),array.max())
            if color_max:
                print("    data total max = ",total_max," (fixed)")
            else:
                print("    data total max = ",total_max)
            print("")

            # limits range [-total_max,total_max]
            array = np.where(array < -total_max, -total_max, array)
            array = np.where(array > total_max, total_max, array)

            # in case color-max is larger than actual range, this sets an arbitrary point to the maximum value
            # to get the correct range value when plotting the data
            if abs(array).max() < total_max:
                array[0,0] = total_max

            # for shakemaps, start range with a minimum value of 0
            if 'shaking' in vtk_file or 'shakemap' in vtk_file:
                array[1,0] = 0.0

            # sets updated range back to vtk array
            print("    new data: color data min/max   = ",array.min(),array.max())
            for i in range(num_points):
                for j in range(num_components):
                    value = array[i,j]
                    data_array.SetComponent(i, j, value)

            # Inform VTK that the array has been modified
            data_array.Modified()
            data_min,data_max = data_array.GetRange()
            print("    new data: range min/max = ",data_min,"/",data_max)
            print("")

        # sets vertex data
        scaled_polydata.GetPointData().SetScalars(data_array)

    # updates the points in the original dataset with the scaled points
    reader.GetOutput().SetPoints(scaled_points)
    reader.Update()

    # output scaled bounds
    points = reader.GetOutput().GetPoints()
    xmin,xmax,ymin,ymax,zmin,zmax = points.GetBounds()
    print("  mesh dimensions after scaling: x min/max = ",xmin,xmax)
    print("                                 y min/max = ",ymin,ymax)
    print("                                 z min/max = ",zmin,zmax)
    print("")

    # Get the unstructured grid data
    unstructured_grid = reader.GetOutput()

    # Extract colors if available
    colors_array = None
    if unstructured_grid.GetPointData().GetNumberOfArrays() > 0:
        print("  grid: point data arrays = ",unstructured_grid.GetPointData().GetNumberOfArrays())
        print("")
        colors_array = unstructured_grid.GetPointData().GetArray(0)  # Assuming colors are in the first array
        print("  colors: name = ",colors_array.GetName())
        print("          range = ",colors_array.GetRange())
        print("")

    #debug
    #print("colors_array: ",colors_array)

    # convert to .ply data file
    # Path to your generated .ply file
    obj_file = 'output.ply'

    # convert the data to polydata
    geometry_filter = vtk.vtkGeometryFilter()
    geometry_filter.SetInputData(unstructured_grid)
    geometry_filter.Update()

    # creates a new polydata object
    polydata = geometry_filter.GetOutput()

    print("  polydata: initial number of points",polydata.GetNumberOfPoints())
    print("  polydata: initial number of verts",polydata.GetNumberOfVerts())
    print("")

    # checks if we have points
    # for proc****_free_surface.vtk files, only points are stored in the .vtk file
    # and the geometry_filter won't fill the polydata object.
    # here we check that we have points & verts filled, otherwise we assume to have points only in the .vtk file
    # and we will try to get a connectivity by Delauney triangulation
    if polydata.GetNumberOfPoints() == 0:
        print("  unstructured grid: number of points",unstructured_grid.GetNumberOfPoints())
        if unstructured_grid.GetNumberOfPoints() > 0:
            points = unstructured_grid.GetPoints()
            polydata.SetPoints(points)

    if polydata.GetNumberOfVerts() == 0:
        print("  unstructured grid: number of cells",unstructured_grid.GetNumberOfCells())
        if unstructured_grid.GetNumberOfCells() > 0:
            verts = unstructured_grid.GetCells()
            polydata.SetVerts(verts)
        else:
            # Perform Delaunay triangulation to generate connectivity
            print("  getting Delaunay 2D connectivity...")
            delaunay = vtk.vtkDelaunay2D()
            delaunay.SetInputData(polydata)
            delaunay.Update()

            # Get the output triangles
            triangles = delaunay.GetOutput()
            #print(triangles)
            print("  triangles: number of verts",triangles.GetNumberOfVerts())
            print("  triangles: number of cells",triangles.GetNumberOfCells())
            polydata = triangles

        print("  polydata: number of points",polydata.GetNumberOfPoints())
        print("  polydata: number of verts",polydata.GetNumberOfVerts())
        print("  polydata: number of cells",polydata.GetNumberOfCells())
        print("  polydata: number of strips",polydata.GetNumberOfStrips())
        print("  polydata: number of data arrays",polydata.GetPointData().GetNumberOfArrays())
        print("")

    # we need to set a default lookup table for the polydata set,
    # otherwise the color float values on the points won't get stored
    #
    # set lookup table for depth values to colors
    if colors_array:
        # assign color table
        lut = setup_color_table(colormap, data_min, data_max)
        colors_array.SetLookupTable(lut)

    # Write the data to PLY format
    writer = vtk.vtkPLYWriter()
    writer.SetInputData(polydata)

    # Include vertex colors if available
    if colors_array:
        writer.SetArrayName(colors_array.GetName())
        writer.SetLookupTable(lut)
        # info
        print("  writer: color mode = ",writer.GetColorMode())
        print("  writer: color array name = ",writer.GetArrayName())
        print("  writer: color component  = ",writer.GetComponent())
        print("")

    os.system('rm -f output.ply')

    writer.SetFileName(obj_file)
    writer.Write()

    if not os.path.exists(obj_file):
        print("Error writing file ",obj_file)
        sys.exit(1)

    print("")
    print("  converted to: ",obj_file)
    print("")

    # work-around for .obj files
    # however, .obj file can by default only store the mesh, not the color data on the vertices
    # we thus prefer to work with the .ply file format above.
    #
    ## save mesh as .obj file
    ## Path to your generated .obj file
    #obj_file = 'output.obj'
    #
    ## Convert the data to polydata
    #geometry_filter = vtk.vtkGeometryFilter()
    #geometry_filter.SetInputConnection(reader.GetOutputPort())
    #geometry_filter.Update()
    #
    ## Write the data to .obj format
    #writer = vtk.vtkOBJWriter()
    #writer.SetInputConnection(geometry_filter.GetOutputPort())
    #
    #os.system('rm -f output.obj')
    #
    #writer.SetFileName(obj_file)  # Output .obj file path
    #writer.Write()
    #
    #if not os.path.exists(obj_file):
    #    print("Error writing file ",obj_file)
    #    sys.exit(1)
    #
    ## appends data lines
    ##..
    ##d val
    #if data_array:
    #    min_val = data_array.GetValue(0)
    #    max_val = data_array.GetValue(0)
    #
    #    with open(obj_file, 'a') as f:
    #        f.write("# Associated data values:\n")
    #        for i in range(data_array.GetNumberOfTuples()):
    #            data_value = data_array.GetValue(i)
    #            min_val = np.minimum(min_val,data_value)
    #            max_val = np.maximum(max_val,data_value)
    #            f.write(f"# vertex {i}: {data_value}\n")
    #            #f.write(f"d {data_value}\n")
    #    print("  appended data: min/max = ",min_val,max_val)
    #    data_min = min_val
    #    data_max = max_val
    #
    #print("")
    #print("  converted to: ",obj_file)

    return obj_file

def create_blender_setup(obj_file: str="") -> None:
    global use_matte_material

    ## Blender setup
    print("blender setup:")
    print("")

    # clears all existing objects in the scene
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()

    # clears only default Cube object and any previous output & output_dfx object,
    # leaves Camera and Light object
    #print("  objects: ",bpy.data.objects.keys())
    #objs = bpy.data.objects
    #for obj in objs:
    #    print("  objects: ",obj.name)
    #    if obj.name == "Cube": objs.remove(obj, do_unlink=True)
    #    if obj.name.contains("output"): objs.remove(obj, do_unlink=True)

    print("  objects: ",bpy.data.objects.keys())
    objs = bpy.data.objects
    for obj in objs:
        print("  objects: ",obj.name)
    print("")

    # import mesh object into blender
    if len(obj_file) > 0:
        print("  importing mesh file: ",obj_file)
        # gets file extension
        extension = os.path.splitext(obj_file)[1]
        # reads the mesh object file
        if extension == '.ply':
            # Import .ply file
            # New experimental, but much faster
            bpy.ops.wm.ply_import(filepath=obj_file)
            # or standard/legacy ply importer
            #bpy.ops.import_mesh.ply(filepath=obj_file)
        elif extension == '.obj':
            # Import .obj file into Blender
            bpy.ops.import_scene.obj(filepath=obj_file)
        else:
            print("unknown mesh object file extension ",extension," - exiting...")
            sys.exit(1)

        print("  imported in blender: ",obj_file)
        print("")

    # blender info
    print("  scenes : ",bpy.data.scenes.keys())
    print("  objects: ",bpy.data.objects.keys())
    for obj in bpy.data.objects:
        print("    object: ",obj.name)
    print("")

    ## mesh object
    #obj = bpy.context.object
    #print(obj.name, ":", obj)
    #objs = bpy.context.selected_objects
    #print(", ".join(o.name for o in objs))
    # Select the imported object
    obj = bpy.data.objects['output']
    if obj == None:
        print("Error: no mesh object in blender available, exiting...")
        sys.exit(1)

    # object is a mesh
    mesh = obj.data

    #debug
    #print("  obj: ",obj)
    print("  obj type: ",obj.type)
    print("  obj data: ",mesh)
    print("  obj polygons: ",mesh.polygons)
    #print("  obj polygon 0: ",mesh.polygons[0])
    #print("  obj polygon 0 loop: ",mesh.polygons[0].loop_indices)
    #print("  obj data loops: ",mesh.loops)
    #print("  obj data loops: ",mesh.loops[0])
    print("  obj data vertex_colors: ",mesh.vertex_colors)
    #print("  obj data vertex_layers_float: ",mesh.vertex_layers_float)
    #print("  obj data vertex_layers_int: ",mesh.vertex_layers_int)
    #print("  obj data vertex_layers_string: ",mesh.vertex_layers_string)
    #print("  obj data vertex_colors: ",mesh.vertex_colors.get("a"))
    #for loop_index in mesh.polygons[0].loop_indices:
    #      index = mesh.loops[loop_index].vertex_index
    #      print("  loop: index",loop_index,obj.data.loops[loop_index].vertex_index)
    #      color = vertex_colors[loop_index].color
    #print("  obj data polygon_layers_float: ",mesh.polygon_layers_float)
    #print("  obj data polygon_layers_int: ",mesh.polygon_layers_int)
    #print("  obj data uv_layers: ",mesh.uv_layers)
    #print("  obj data vertices",mesh.vertices)
    #print("  obj data vertex 0",mesh.vertices[0])
    #print("  obj data vertex 0 coord",mesh.vertices[0].co)
    #print("  obj data vertex 0 keys",mesh.vertices.keys())
    #print("  obj data vertex 0 get",mesh.vertices.items())
    print("")

    #print("obj data vertex_colors 0: ",obj.data.vertex_colors[0])
    #print("obj data vertex_colors active data: ",obj.data.vertex_colors.active.data)

    if obj is not None:
        # Ensure the object has a mesh and vertex colors
        if obj.type == 'MESH' and not obj.data.vertex_colors is None:
            print("  Object 'output' has vertex colors and is a mesh.")
            # Access vertex color data
            #vertex_colors = obj.data.vertex_colors.active.data
            # Iterate through vertex color data
            #for poly in obj.data.polygons:
            #    for loop_index in poly.loop_indices:
            #        vertex_index = obj.data.loops[loop_index].vertex_index
            #        color = vertex_colors[loop_index].color
            #        # Print information about vertex colors
            #        #print(f"Vertex {vertex_index}: Color {color}")
        else:
            print("  Object 'output' does not have vertex colors or is not a mesh.")
    else:
        print("  Object 'output' not found.")
        sys.exit(1)

    print("")
    print("  mesh: setting up shader nodes...")
    print("")

    # assigns new material
    mat = bpy.data.materials.new(name="VertexColorMaterial")
    mat.use_nodes = True
    # node-graph
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links

    # Clear default nodes
    #for node in nodes:
    #    nodes.remove(node)

    # Create Attribute node to fetch vertex color
    color_attribute = nodes.new(type='ShaderNodeAttribute')
    color_attribute.attribute_name = "Col"  # Use "Col" as it's the default name for vertex color
    # Set the Color Attribute node to use vertex colors
    color_attribute.attribute_type = 'GEOMETRY' # 'COLOR'

    # Create Principled BSDF shader node
    # checks default node
    bsdf = nodes["Principled BSDF"]
    if bsdf == None:
        bsdf = nodes.new(type='ShaderNodeBsdfPrincipled')

    if use_matte_material:
        # more matte material appearance
        print("  using matte material appearance")
        bsdf.inputs['Metallic'].default_value = 0.1
        bsdf.inputs['Roughness'].default_value = 0.8
        bsdf.inputs['Specular'].default_value = 0.3
    else:
        # default material w/ more glossy appearance
        bsdf.inputs['Metallic'].default_value = 0.4
        bsdf.inputs['Roughness'].default_value = 0.5
        bsdf.inputs['Specular'].default_value = 0.2

    # (default) color map from vtk file
    links.new(bsdf.inputs['Base Color'],color_attribute.outputs['Color'])

    # no emission
    bsdf.inputs['Emission Strength'].default_value = 0.0
    # w/ emission (default) for brighter colors
    #links.new(bsdf.inputs['Emission'],color_attribute.outputs['Color'])
    #links.new(bsdf.inputs['Emission Strength'],color_attribute.outputs['Fac'])

    # custom mesh coloring w/ color ramp node
    if 1 == 0:
        # creates a custom color ramp node to highlight shaking regions
        ramp = nodes.new(type='ShaderNodeValToRGB')
        ramp.color_ramp.interpolation = 'LINEAR'
        #ramp.color_ramp.interpolation = 'B_SPLINE'
        # default 2 slots
        ramp.color_ramp.elements[0].position = 0.0
        ramp.color_ramp.elements[0].color = [1,1,1,1]
        ramp.color_ramp.elements[1].position = 0.4
        ramp.color_ramp.elements[1].color = [0.5,0.5,0.5,1]  # gray
        # add color slot
        ramp.color_ramp.elements.new(0.5)
        ramp.color_ramp.elements[2].color = [0.2,0.2,0.3,1]  # dark blue
        # add color slot
        ramp.color_ramp.elements.new(0.6)
        ramp.color_ramp.elements[3].color = [0.8,0.0,0.0,1]  # red
        # add color slot
        ramp.color_ramp.elements.new(0.65)
        ramp.color_ramp.elements[4].color = [1.0,0.6,0.0,1]  # yellow
        # add color slot
        ramp.color_ramp.elements.new(0.7)
        ramp.color_ramp.elements[5].color = [1,1,1,1]       # white

        # Link Math node output to Color Ramp factor input
        links.new(ramp.inputs["Fac"],color_attribute.outputs["Fac"])
        # custom Color Ramp output to Principled BSDF node
        links.new(bsdf.inputs['Base Color'],ramp.outputs['Color'])

    # adds additional emission
    if 1 == 0:
        # Create Math node to manipulate grayscale value
        math_node = nodes.new(type='ShaderNodeMath')
        math_node.operation = 'GREATER_THAN'
        math_node.inputs[1].default_value = 0.7  # Set threshold value
        links.new(math_node.inputs['Value'],color_attribute.outputs['Fac'])

        # Create RGB to BW node to convert color to black/white float value
        #rgb_to_bw = nodes.new(type='ShaderNodeRGBToBW')
        #links.new(rgb_to_bw.inputs['Color'],color_attribute.outputs['Color'])
        #links.new(math_node.inputs['Value'],rgb_to_bw.outputs['Val'])

        # emission node
        emission = nodes.new('ShaderNodeEmission')
        emission.name = "Emission"
        links.new(emission.inputs['Color'],color_attribute.outputs['Color'])
        links.new(emission.inputs['Strength'],math_node.outputs['Value'])

        # mix light emission and main image
        mix = nodes.new('ShaderNodeMixShader')
        mix.name = "Mix Shader"
        #links.new(mix.inputs['Fac'], ramp.outputs['Alpha'])
        # takes output from main BSDF node
        links.new(mix.inputs[1], bsdf.outputs[0])
        links.new(mix.inputs[2], emission.outputs[0])

        # link mixer to final material output
        material_output = nodes["Material Output"]
        links.new(material_output.inputs["Surface"], mix.outputs["Shader"])

    # Assign the material to the object
    if obj.data.materials:
        obj.data.materials[0] = mat
    else:
        obj.data.materials.append(mat)

    print("  blender mesh done")
    print("")


def add_blender_buildings(buildings_file: str="") -> None:
    """
    adds buildings given by input .ply file
    """
    global mesh_origin,mesh_scale_factor
    global use_cycles_renderer

    # checks if anything to do
    if len(buildings_file) == 0: return

    print("buildings file: ",buildings_file)
    print("")

    # check file
    if not os.path.exists(buildings_file):
        print("Error: buildings file specified not found...")
        sys.exit(1)

    # buildings need to be given as .ply (Stanford) format file
    print("  reading in .ply mesh...")

    # Enable the experimental .ply importer
    # New experimental, but much faster
    bpy.ops.wm.ply_import(filepath=buildings_file)
    # or standard/legacy ply importer
    #bpy.ops.import_mesh.ply(filepath=buildings_file)

    # Print the names of imported objects
    print("")
    for obj in bpy.context.selected_objects:
        print("  imported object:", obj.name," - type: ",obj.type)
        # select buildings mesh (name must contain dxf)
        #if obj.type == 'MESH' and 'dxf' in obj.name:
        #    bpy.context.view_layer.objects.active = obj
        #    obj_buildings = obj
        #    break
    print("")

    # gets active object
    obj = bpy.context.view_layer.objects.active

    if obj is None:
        print("  Object for buildings not found.")
        sys.exit(1)

    # note: due to different resolution of the meshes, the mesh elevation of SPECFEM mesh (e.g., AVS_shaking_map.inp)
    #       and the buildings (e.g., from SwissTopo) might be slightly off by a few meters.

    # moves and scales UTM mesh
    if 'utm' in buildings_file:
        # need to translate and scale the UTM mesh to place it within the vtk mesh
        print("  UTM mesh: moving & scaling mesh...")
        print("            mesh origin       = ",mesh_origin)
        print("            mesh scale factor = ",mesh_scale_factor)
        print("")

        # Get the mesh data
        mesh = obj.data

        # Create a BMesh object
        bm = bmesh.new()
        bm.from_mesh(mesh)

        # Translate the vertices (example: translating by (1, 0, 0))
        translation_vector = Vector(mesh_origin)
        for vert in bm.verts:
            vert.co -= translation_vector

        # Scale the vertices (example: scaling by 1.5 in all axes)
        scale_factor = mesh_scale_factor
        for vert in bm.verts:
            vert.co *= scale_factor

        # Update the mesh with the modified vertices
        bm.to_mesh(mesh)
        bm.free()

    # baseline shift
    if shift_building_baseline != 0.0:
        print("  shifting buildings up: factor = ",shift_building_baseline)
        print("                          in m  = ",shift_building_baseline / mesh_scale_factor)
        print("")
        obj.location = (0, 0, shift_building_baseline)

    # Create a new material
    mat = bpy.data.materials.new('buildingsMaterial')

    # enable node-graph edition mode
    mat.use_nodes = True

    nodes = mat.node_tree.nodes
    links = mat.node_tree.links

    # Clear default nodes
    for node in nodes:
        nodes.remove(node)

    # gets scene
    scene = bpy.context.scene

    if use_cycles_renderer:
        # Set render engine to Cycles
        scene.render.engine = 'CYCLES'

        # Create a Glass BSDF node
        mat_node = nodes.new(type='ShaderNodeBsdfGlass')
        # Set IOR (Index of Refraction) and Roughness values
        mat_node.inputs['IOR'].default_value = 1.5       # (e.g., 1.5 for typical glass)
        mat_node.inputs['Roughness'].default_value = 0.6  # Roughness value (0.0 for perfectly smooth)

    else:
        # BLENDER_EEVEE renderer
        scene.render.engine = 'BLENDER_EEVEE'

        # create a glossy node
        mat_node = nodes.new(type='ShaderNodeBsdfGlossy')
        mat_node.inputs['Roughness'].default_value = 0.5  # Roughness value (0.0 for perfectly smooth)
        mat_node.inputs['Color'].default_value = (0.8, 0.74, 0.56, 1)

        # create a glass node
        # note: transparency effect w/ eevee for buildings looks not good enough yet...
        #       glass buildings need cycles renderer.
        #       thus, i'm using glossy buildings to get some environment colors onto buildings.
        #mat_node = nodes.new(type='ShaderNodeBsdfGlass')
        # node color
        #mat_node.inputs['Color'].default_value = (0.6, 0.6, 0.6, 1)
        # enable transparency for eevee
        # material transparency
        #mat.blend_method = 'BLEND'   # 'OPAQUE', 'BLEND', ..
        #mat.shadow_method = 'HASHED'
        #mat.use_backface_culling = False
        #mat.show_transparent_back = True
        #mat.use_screen_refraction = True
        #mat.refraction_depth = 0.01

    # Create Output node
    output_node = nodes.new(type='ShaderNodeOutputMaterial')
    output_node.location = 400, 0

    # Link nodes
    links.new(mat_node.outputs['BSDF'], output_node.inputs['Surface'])

    # Assign the material to the object
    if obj.data.materials:
        obj.data.materials[0] = mat
    else:
        obj.data.materials.append(mat)

    # blender info
    print("  scenes : ",bpy.data.scenes.keys())
    print("  objects: ",bpy.data.objects.keys())
    for obj in bpy.data.objects:
        print("    object: ",obj.name)
    print("")

    print("  blender buildings mesh done")
    print("")

def add_borders(borders_file: str="") -> None:
    """
    adds AVS boundary borders given by input .inp file
    """
    global mesh_origin,mesh_scale_factor

    # checks if anything to do
    if len(borders_file) == 0: return

    print("borders file: ",borders_file)
    print("")

    # check file
    if not os.path.exists(borders_file):
        print("Error: borders file specified not found...")
        sys.exit(1)

    # borders need to be given as .inp (AVS UCD) format file
    print("  reading in .inp mesh...")

    # gets file extension
    extension = os.path.splitext(vtk_file)[1]
    # reads the vtk file
    if extension == '.inp':
        # AVS .inp
        #reader = vtk.vtkSimplePointsReader()
        reader = vtk.vtkAVSucdReader()
        reader.SetFileName(borders_file)
        reader.Update()
    else:
        print("unknown borders file extension ",extension," - exiting...")
        sys.exit(1)

    #debug
    #print(reader)

    # convert to .ply data file
    # Get the unstructured grid data
    unstructured_grid = reader.GetOutput()

    # convert the data to polydata
    geometry_filter = vtk.vtkGeometryFilter()
    geometry_filter.SetInputData(unstructured_grid)
    geometry_filter.Update()
    # creates a new polydata object
    polydata = geometry_filter.GetOutput()

    # checks if we have points
    if polydata.GetNumberOfPoints() == 0:
        print("  unstructured grid: number of points",unstructured_grid.GetNumberOfPoints())
        if unstructured_grid.GetNumberOfPoints() > 0:
            points = unstructured_grid.GetPoints()
            polydata.SetPoints(points)
        else:
            print("no points found.")
            return

    print("  polydata: number of points",polydata.GetNumberOfPoints())
    print("  polydata: number of lines",polydata.GetNumberOfLines())
    # not interested in cells, point data, etc.
    #print("  polydata: number of verts",polydata.GetNumberOfVerts())
    #print("  polydata: number of cells",polydata.GetNumberOfCells())
    #print("  polydata: number of strips",polydata.GetNumberOfStrips())
    #print("  polydata: number of data arrays",polydata.GetPointData().GetNumberOfArrays())
    print("")

    if polydata.GetNumberOfPoints() == 0:
        print("no points found in the AVS UCD file after extracting edges.")
        return

    if polydata.GetNumberOfLines() == 0:
        print("no lines found in the AVS UCD file after extracting edges.")
        return

    # moves and scales UTM point locations
    # need to translate and scale the UTM point location to place it within the vtk mesh
    print("  UTM points: moving & scaling mesh...")
    print("              mesh origin       = ",mesh_origin)
    print("              mesh scale factor = ",mesh_scale_factor)
    print("")

    # Loop through each point, scale its coordinates, and add it to the new points array
    points_vtk = polydata.GetPoints()
    num_points = polydata.GetNumberOfPoints()
    for i in range(num_points):
        # Get the current coordinates of the point
        # translation by origin
        point = np.array(points_vtk.GetPoint(i)) - mesh_origin
        # uniform scaling
        point *= mesh_scale_factor
        # Define the point of interest (X, Y, Z coordinates)
        vpoint = Vector((point[0], point[1], point[2]))
        # get elevation
        elevation = get_mesh_elevation(vpoint)
        # checks if valid
        if elevation is None: elevation = 0.1 # sets a default height
        # debug
        #print(f"    mesh elevation at point {point} = {elevation}")
        # set to vertical elevation
        point[2] = elevation
        # set modified coordinates back
        points_vtk.SetPoint(i,point)

    # Create a new curve datablock
    print("  creating Blender curve objects...")
    curve_data = bpy.data.curves.new(name='Borders', type='CURVE')
    curve_data.dimensions = '3D'
    curve_data.resolution_u = 2  # Resolution of the curve in viewport/render

    # Create a new object with the curve datablock
    curve_obj = bpy.data.objects.new('Borders', curve_data)

    # Link the object to the scene collection
    bpy.context.collection.objects.link(curve_obj)

    # Organize into a specific collection
    target_collection = bpy.data.collections.new('AVS_Borders_Lines')
    bpy.context.scene.collection.children.link(target_collection)
    print("    created new collection: 'AVS_Borders_Lines'")

    # Link the object to the target collection
    target_collection.objects.link(curve_obj)

    # Unlink from the default scene collection if it's there
    # This prevents the object from appearing in multiple collections simultaneously if it was already linked
    if curve_obj.name in bpy.context.collection.objects and bpy.context.collection != target_collection:
        bpy.context.collection.objects.unlink(curve_obj)
        #print(f"  unlinked '{curve_obj.name}' from default scene collection.")

    # Iterate through each polyline/line in the VTK data
    # VTK lines are stored as a connectivity list.
    # Each entry starts with the number of points in the polyline,
    # followed by the point indices.
    lines_vtk = polydata.GetLines()

    # Reset cursor for lines iteration
    lines_vtk.InitTraversal()

    id_list = vtk.vtkIdList()
    num_splines_created = 0

    while lines_vtk.GetNextCell(id_list):
        num_points_in_line = id_list.GetNumberOfIds()
        # A line needs at least 2 points
        if num_points_in_line < 2:
            continue

        # Create a new spline for each line/polyline
        spline = curve_data.splines.new('POLY') # Use 'POLY' for straight line segments
        spline.points.add(num_points_in_line - 1) # Add points (one is already there)

        for i in range(num_points_in_line):
            point_index = id_list.GetId(i)
            # VTK points are float[3]
            x, y, z = points_vtk.GetPoint(point_index)
            # Set the coordinates for the spline point
            # Blender spline points are (x, y, z, w) where w is weight for NURBS, not needed for POLY
            spline.points[i].co = (x, y, z, 1.0) # Set weight to 1.0 for POLY

        num_splines_created += 1

    # Make lines renderable (e.g., as tubes)
    borders_line_thickness = 0.002

    curve_data.bevel_depth = borders_line_thickness  # Thickness of the tube
    curve_data.bevel_resolution = 2 # Smoothness of the tube
    curve_data.fill_mode = 'FULL' # Make it a solid tube

    print(f"    created Blender curve 'Borders' with {num_splines_created} splines.")
    print("")
    return


def get_mesh_elevation(point_of_interest: Vector) -> float:
    """
    determines elevation of object (obj) at a given point by ray intersection
    """
    # gets vtk mesh object
    obj = bpy.data.objects['output']
    if obj == None:
        print("Info: no mesh object in blender available to determine elevation...")
        return None

    # Define the origin of the ray (above the mesh)
    ray_origin = Vector((point_of_interest.x, point_of_interest.y, 10.0))  # Adjust the Z coordinate as needed

    # Get the world matrix of the object
    matrix_world = obj.matrix_world

    # Calculate the direction of the ray (pointing downwards)
    ray_direction = Vector((0, 0, -1))

    # Transform the ray direction to world space
    ray_direction.rotate(matrix_world.to_quaternion())

    # Perform the ray casting
    success, location, _, _ = obj.ray_cast(ray_origin, ray_direction)

    if success:
        # Get the elevation (Z coordinate) of the mesh at the point of interest
        return location.z
    else:
        return None


def get_mesh_elevation_at_origin() -> float:
    """
    determines elevation of mesh at origin/center point
    """
    global x_center,y_center,z_center

    print("    origin center            = {} / {} / {}".format(x_center,y_center,z_center))

    # Define the point of interest (X, Y, Z coordinates)
    point = Vector((x_center, y_center, z_center))  # origin

    # get elevation
    elevation = get_mesh_elevation(point)

    # checks if valid
    if not elevation is None:
        # got an elevation
        print("    mesh elevation at origin = ", elevation)
    else:
        # sets a default height
        elevation = 0.1
        print("    mesh elevation at origin: no intersection found")
        print("                              setting default = ",elevation)

    return elevation

def add_camera(title: str="") -> None:
    """
    adds camera object
    """
    global x_center,y_center,z_center
    global camera_y_offset,camera_elevation_offset,camera_angle_depth_degree
    global z_elevation
    global close_up_view,centered_view

    print("  adding camera")

    # gets scene
    scene = bpy.context.scene

    # adds camera position
    bpy.ops.object.camera_add(enter_editmode=False, align='VIEW')
    # current object
    cam = bpy.context.object
    cam.name = "Camera"   #cam = bpy.data.objects["Camera"]
    scene.camera = cam

    # clip range
    # note: our objects are scaled between [-1,1].
    #       a clip end plane around 10 should be sufficient to cover the whole area.
    #       for the clip start, it is more challenging as choosing a value below 0.1 will lead to depth Z-buffer artifacts.
    #       to avoid this, for an overview scene, the default value of 0.1 is good, when zooming in into close-up views,
    #       a clip start around 0.01 is usually better, otherwise the front of the scene will be cut-off.
    scene.camera.data.clip_start = 0.1    # blender default: 0.1
    scene.camera.data.clip_end = 10.0     # blender default: 1000
    #scene.camera.data.lens = 70.0

    # determine elevation from mesh at origin point to position the camera
    z_elevation = get_mesh_elevation_at_origin()

    # elevation offset for camera positioning
    camera_elevation_offset = 0.15
    camera_y_offset = -0.7
    camera_angle_depth_degree = 75.0

    # specifics
    if title == 'Lauterbrunnen' or title == 'Zermatt':
        print("    camera setup for Lauterbrunnen/Zermatt")
        camera_elevation_offset = 0.25
        camera_angle_depth_degree = 70.0

    # center point to focus camera view on
    x_center = 0.0
    y_center = 0.0
    z_center = 0.0

    ## Rome colosseum
    if title == "Rome":
        print("    camera setup for Rome")
        camera_elevation_offset = 0.0018
        camera_y_offset = -0.013
        camera_angle_depth_degree = 82.0
        # use colosseum position as center of views
        if 1 == 1:
            # colosseum: N 41.890258 E 12.492335
            #            33 T -> utm 291957.383 4640632.710
            x_center_utm = 291957.383; y_center_utm = 4640632.710
            # need to translate and scale the UTM position to place it within the scene
            point = Vector((x_center_utm, y_center_utm, 0.0))
            translation_vector = Vector(mesh_origin)
            point -= translation_vector
            point *= mesh_scale_factor
            # set center position
            x_center = point.x
            y_center = point.y
        else:
            # directly from using cursor position in blender
            x_center = -0.01251
            y_center = -0.2263
        z_center = 0.0018

        # get elevation at x/y position
        point = Vector((x_center, y_center, 0.0))  # origin
        center_elevation = get_mesh_elevation(point)
        if not center_elevation is None:
            z_elevation = center_elevation
        else:
            z_elevation = 0.001

        #cam.location = (-0.02251, -0.24362, 0.00632)
        #cam.rotation_euler = (83.401 * DEGREE_TO_RAD, 0.000133 * DEGREE_TO_RAD, -28.0 * DEGREE_TO_RAD)

        # metallic material for buildings
        if not use_cycles_renderer:
            if 'buildingsMaterial' in bpy.data.materials.keys():
                mat = bpy.data.materials['buildingsMaterial']
                nodes = mat.node_tree.nodes
                mat_node = nodes['Glossy BSDF']
                mat_node.inputs['Roughness'].default_value = 0.25  # shinier, more metallic
                print("      using metallic buildings")

    if centered_view:
        ## centered view
        # camera location
        camera_elevation_offset = 4.0
        x_cam = x_center
        y_cam = y_center
        z_cam = z_center + camera_elevation_offset
        print("    centered scene")
        print("    center point   : x/y/z = {:.6f} / {:.6f} / {:.6f}".format(x_center,y_center,z_center))
        print("    camera location: x/y/z = {:.6f} / {:.6f} / {:.6f}".format(x_cam,y_cam,z_cam))
        print("    camera location relative to mesh elevation: ",z_elevation)
        # Set camera translation
        scene.camera.location = (x_cam,y_cam,z_cam)
        # Set camera rotation in euler angles
        scene.camera.rotation_mode = 'XYZ'
        scene.camera.rotation_euler = (0, 0, 0)  # top-down
        # Set camera fov in degrees
        scene.camera.data.angle = float(42.0 * DEGREE_TO_RAD)
    elif close_up_view:
        ## close-up view
        # adjust clip range
        scene.camera.data.clip_start = 0.01
        #scene.camera.data.clip_end = 10.0
        #scene.camera.data.lens = 70.0
        # camera location
        x_cam = x_center
        y_cam = y_center + camera_y_offset
        z_cam = z_center + z_elevation + camera_elevation_offset
        print("    close-up scene")
        print("    center point   : x/y/z = {:.6f} / {:.6f} / {:.6f}".format(x_center,y_center,z_center))
        print("    camera location: x/y/z = {:.6f} / {:.6f} / {:.6f}".format(x_cam,y_cam,z_cam))
        print("    camera location relative to mesh elevation: ",z_elevation)
        # Set camera translation
        scene.camera.location = (x_cam,y_cam,z_cam)
        # Set camera rotation in euler angles
        scene.camera.rotation_mode = 'XYZ'
        scene.camera.rotation_euler = (camera_angle_depth_degree * DEGREE_TO_RAD, 0, 0)
        # Set camera fov in degrees
        scene.camera.data.angle = float(30.0 * DEGREE_TO_RAD)
    else:
        ## overview
        # camera location
        x_cam = x_center
        y_cam = y_center - 4.0
        z_cam = z_center + 4.0
        print("    overview scene")
        print("    center point   : x/y/z = {:.6f} / {:.6f} / {:.6f}".format(x_center,y_center,z_center))
        print("    camera location: x/y/z = {:.6f} / {:.6f} / {:.6f}".format(x_cam,y_cam,z_cam))
        # Set camera translation
        scene.camera.location = (x_cam,y_cam,z_cam)
        # Set camera rotation in euler angles
        scene.camera.rotation_mode = 'XYZ'
        scene.camera.rotation_euler = (44.0 * DEGREE_TO_RAD, 0, 0)
        # Set camera fov in degrees
        scene.camera.data.angle = float(30.0 * DEGREE_TO_RAD)


def add_light() -> None:
    """
    adds a light to the scene
    """
    print("  adding light")
    # adds sun position
    bpy.ops.object.light_add(type='AREA')
    # current object
    light = bpy.context.object
    light.name = "Light"  #light = bpy.data.objects["Light"]

    light.location = (1.5, -0.5, 1.3)
    light.rotation_mode = 'XYZ'
    light.rotation_euler = (0, 40.0 * DEGREE_TO_RAD, 0)
    # sets light to rectangular (plane)
    light.data.type = 'AREA'
    light.data.shape = 'SQUARE'
    # Change the light's color
    light.data.color = (1, 1, 1)  # Set light color to white
    # intensity
    light.data.energy = 80  # W
    # Set the light's size
    light.data.size = 0.1
    light.data.use_contact_shadow = True  # Enable contact shadows


def add_plane() -> None:
    """
    adds a plane at sea-level
    """
    global use_sea_level_plane,use_transparent_sea_level_plane,add_sea_level_plane_modifier
    global close_up_view

    # checks if anything to do
    if not use_sea_level_plane: return

    print("  adding sea-level plane")
    # Create a mesh plane (to capture shadows and indicate sea level)
    bpy.ops.mesh.primitive_plane_add(size=10, enter_editmode=False, location=(0, 0, 0))

    # Get the created plane object
    plane_object = bpy.context.object

    # Set the object's material to white
    mat = bpy.data.materials.new(name="White")
    # adds transparency to the plane
    if use_transparent_sea_level_plane:
        print("    using transparent sea-level")
        # Set the material to use a Principled BSDF shader
        mat.use_nodes = True
        principled_bsdf = mat.node_tree.nodes.get('Principled BSDF')
        # color
        if close_up_view:
            principled_bsdf.inputs["Base Color"].default_value = (0.8, 0.8, 0.8, 1)
        else:
            principled_bsdf.inputs["Base Color"].default_value = (0.135, 0.135, 0.135, 1) # gray for better contrast
        # Set the shader to be transparent
        principled_bsdf.inputs['Alpha'].default_value = 0.6  # Set the alpha to control transparency
        # Set the blend mode to 'Alpha Blend'
        mat.blend_method = 'BLEND'
        # Set the shadow mode to 'Alpha Hashed'
        mat.shadow_method = 'HASHED'
    else:
        if close_up_view:
            mat.diffuse_color = (0.8, 0.8, 0.8, 1)  # similar as background color to have an infinite background
        else:
            mat.diffuse_color = (0.135, 0.135, 0.135, 1)  # gray for better contrast

    # adds the Boolean modifier
    if add_sea_level_plane_modifier:
        print("    adding sea-level plane - modifier using difference")
        boolean_modifier = plane_object.modifiers.new(name="Boolean", type='BOOLEAN')

        # Set the operation type
        boolean_modifier.operation = 'DIFFERENCE'

        # Set the operand object if provided
        mesh_object = bpy.data.objects['output']
        if mesh_object:
            boolean_modifier.object = mesh_object
        else:
            print("    Warning: mesh object 'output' not found.")

        # Set the solver type
        boolean_modifier.solver = 'FAST'   # 'FAST' or 'EXACT'

        # some further optimizations for the exact solver
        if boolean_modifier.solver == 'EXACT':
            #boolean_modifier.use_self = True           # self-intersection
            #boolean_modifier.use_hole_tolerant = True  # better results when more tolerant, but slower
            pass
        print("    modifier added")

    plane_object.data.materials.append(mat)


def add_title(title: str="") -> None:
    """
    adds a title text
    """
    global centered_view

    if len(title) > 0:
        print("  adding title text object")
        print("    title = ",title)

        # Create a new text object
        bpy.ops.object.text_add()
        text_object = bpy.context.object
        text_object.data.body = title  # Set the text content

        # Set text properties (font, size, etc.)
        text_object.data.size = 0.2  # Adjust the font size
        text_object.name = "Title_Text"
        #debug
        #print("  blender fonts available: ",bpy.data.fonts.keys())

        if 'Bfont' in bpy.data.fonts:
            text_object.data.font = bpy.data.fonts['Bfont']  # Use a specific default font
        elif 'Bfont Regular' in bpy.data.fonts:
            text_object.data.font = bpy.data.fonts['Bfont Regular']  # Use a specific default font
        elif 'Arial Regular' in bpy.data.fonts:
            text_object.data.font = bpy.data.fonts['Arial Regular']

        #text_object.data.font = bpy.data.fonts.load("/path/to/your/font.ttf")  # Replace with your font path

        # Adjust the location as needed
        x_title = 0.0
        y_title = -1.2
        z_title = 0.01

        if centered_view:
            print("    centered text")
            # get dimensions of text
            #bbox_center = text_obj.bound_box_center
            #bbox_size = text_object.dimensions
            #print("    text dimensions = ",bbox_size[0],"/",bbox_size[1],"/",bbox_size[2])
            x_title = 0.0
            y_title = -0.95

        print("    location = ",x_title,"/",y_title,"/",z_title)

        # align to center
        text_object.data.align_x = 'CENTER'
        text_object.data.align_y = 'BOTTOM_BASELINE'

        text_object.location = (x_title, y_title, z_title)

        # Set text material
        text_material = bpy.data.materials.new(name="TextMaterial")
        text_material.use_nodes = True
        text_material.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.5, 0.5, 0.5, 1)

        text_object.data.materials.append(text_material)


def add_location_labels(locations_file: str="") -> None:
    """
    adds locations defined in file
    """
    global centered_view,use_white_location_labels
    global location_labels_color

    # checks if anything to do
    if len(locations_file) == 0: return

    print("adding location labels")
    print("  locations file: ",locations_file)
    print("")

    # check file
    if not os.path.exists(locations_file):
        print("Error: locations file specified not found...")
        sys.exit(1)

    lines = []
    with open(locations_file, 'r') as file:
        lines = file.readlines()

    if len(lines) == 0:
        print("  no lines")
        return

    # moves and scales UTM positions to normalized range
    # need to translate and scale the UTM mesh position to place it within the vtk mesh
    print("  UTM mesh: moving & scaling mesh...")
    print("            mesh origin       = ",mesh_origin)
    print("            mesh scale factor = ",mesh_scale_factor)

    # Process each line in the file
    for line in lines:
        # Skip comments and empty lines
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        # Parse city name, latitude, and longitude
        # format: #city name #latitude #longitude
        parts = line.split()
        if len(parts) == 3:
            # example: Tingri 28.5762  86.6197
            loc_name = parts[0]
            lat = float(parts[1])
            lon = float(parts[2])
        elif len(parts) == 4:
            # example: Mount Everest  27.9881  86.9250
            loc_name = parts[0] + " " + parts[1]
            lat = float(parts[2])
            lon = float(parts[3])
        else:
            print(f"  Skipping malformed line: {line}")
            continue

        # Convert latitude and longitude to UTM / 3D coordinates
        x_utm, y_utm = convert_latlon_to_UTM(lat, lon)

        # sets marker type (or label type)
        if loc_name == "-":
            # marker only
            is_marker_only = True
        else:
            # label type
            is_marker_only = False

        # info
        if is_marker_only:
            # marker only
            print(f"  location label: marker at lat/lon {lat}/{lon}")
        else:
            # label text
            print(f"  location label: {loc_name} - lat/lon {lat}/{lon}")

        # Translate the vertices (example: translating by (1, 0, 0))
        x = x_utm - mesh_origin[0]
        y = y_utm - mesh_origin[1]

        # Scale the vertices (example: scaling by 1.5 in all axes)
        scale_factor = mesh_scale_factor
        x *= scale_factor
        y *= scale_factor

        # get elevation at x/y position
        point = Vector((x, y, 0.0))
        elevation = get_mesh_elevation(point)
        if not elevation is None:
            z = elevation
        else:
            z = 0.001

        # use text "-" for markers only (no text)
        if is_marker_only:
            # marker, smaller radius
            r = 0.005
        else:
            # default sphere lable radius
            r = 0.02

        # Create a small sphere (circle) to represent the location
        bpy.ops.mesh.primitive_uv_sphere_add(radius=r, location=(x, y, z))
        sphere = bpy.context.object
        sphere.name = f"{loc_name}_Marker"

        # Set marker material
        sphere_material = bpy.data.materials.new(name="TextMaterial")
        sphere_material.use_nodes = True
        if is_marker_only:
            color = (1.0, 1.0, 1.0, 1)
        else:
            color = (0.1, 0.05, 0.05, 1)
        sphere_material.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = color
        sphere_material.node_tree.nodes["Principled BSDF"].inputs["Alpha"].default_value = 0.5   # transparency
        sphere_material.node_tree.nodes["Principled BSDF"].inputs["Specular"].default_value = 0  # Reduce shine
        sphere_material.shadow_method = 'NONE'  # Blender 3.x
        sphere_material.blend_method = 'BLEND'  # 'Alpha Blend'
        sphere.data.materials.append(sphere_material)

        # Create a new text object
        if not is_marker_only:
            # location text
            bpy.ops.object.text_add()
            text_object = bpy.context.object
            text_object.data.body = loc_name  # Set the text content
            text_object.name = f"{loc_name}_Label"
            #text_object.data.use_shadow = False  # Disable text's own shadow

            # Set text properties (font, size, etc.)
            text_object.data.size = 0.05  # Adjust the font size
            #text_object.rotation_euler = (0, 0, np.radians(90))

            if 'Bfont' in bpy.data.fonts:
                text_object.data.font = bpy.data.fonts['Bfont']  # Use a specific default font
            elif 'Bfont Regular' in bpy.data.fonts:
                text_object.data.font = bpy.data.fonts['Bfont Regular']  # Use a specific default font
            elif 'Arial Regular' in bpy.data.fonts:
                text_object.data.font = bpy.data.fonts['Arial Regular']

            # Adjust the location as needed
            # align to center
            text_object.data.align_x = 'CENTER'
            text_object.data.align_y = 'BOTTOM_BASELINE'

            # get dimensions of text
            #bbox_size = text_object.dimensions
            #print("    text dimensions = ",bbox_size[0],"/",bbox_size[1],"/",bbox_size[2])

            x_text = x
            y_text = y + 0.03
            z_text = z + 0.01   # shift above mesh by tiny bit

            #debug
            #print("    mesh location = ",x_text,"/",y_text,"/",z_text)

            text_object.location = (x_text, y_text, z_text)

            # Set text material
            text_material = bpy.data.materials.new(name="TextMaterial")
            text_material.use_nodes = True
            text_material.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = location_labels_color
            text_material.node_tree.nodes["Principled BSDF"].inputs["Alpha"].default_value = 1   # transparency
            text_material.node_tree.nodes["Principled BSDF"].inputs["Specular"].default_value = 0  # Reduce shine
            text_material.shadow_method = 'NONE'  # Blender 3.x
            text_material.blend_method = 'BLEND'  # 'Alpha Blend'
            text_object.data.materials.append(text_material)

    print("")

def add_time_text():
    """
    adds a time string
    """
    global time_string

    # check if anything to do
    if time_string is None: return

    print("  adding time info")
    print("    time string: ",time_string)

    if len(time_string) == 0: return

    # time text
    bpy.ops.object.text_add()
    text_object = bpy.context.object

    text_object.data.body = time_string  # text
    #text_object.data.extrude = 0.01         # Small extrusion for visibility
    text_object.name = "Time_Label"

    # Set text properties (font, size, etc.)
    text_object.data.size = 0.05  # Adjust the font size
    #text_object.rotation_euler = (0, 0, np.radians(90))

    # load fonts
    font_dir = bpy.context.preferences.filepaths.font_directory
    for file in os.listdir(font_dir):
        if file.endswith(".ttf"):
            bpy.data.fonts.load(font_dir + file)
    # show fonts info
    #for font in bpy.data.fonts:
    #    print("    available font: ",font)
    #print("")

    # Use a specific default font
    if 'DIN-Regular Regular' in bpy.data.fonts:
        print("    using font: ",'DIN-Regular Regular')
        text_object.data.font = bpy.data.fonts['DIN-Regular Regular']
    elif 'Orbitron Regular' in bpy.data.fonts:
        print("    using font: ",'Orbitron Regular')
        text_object.data.font = bpy.data.fonts['Orbitron Regular']
    elif 'Arial Regular' in bpy.data.fonts:
        print("    using font: ",'Arial Regular')
        text_object.data.font = bpy.data.fonts['Arial Regular']
    elif 'Bfont' in bpy.data.fonts:
        print("    using font: ",'Bfont')
        text_object.data.font = bpy.data.fonts['Bfont']
    elif 'Bfont Regular' in bpy.data.fonts:
        print("    using font: ",'Bfont Regular')
        text_object.data.font = bpy.data.fonts['Bfont Regular']

    # Adjust the location as needed
    text_object.data.align_x = 'LEFT'
    text_object.data.align_y = 'BOTTOM_BASELINE'

    # get dimensions of text
    #bbox_size = text_object.dimensions
    #print("    text dimensions = ",bbox_size[0],"/",bbox_size[1],"/",bbox_size[2])

    location = (1.02, 0.0, 0.01)  # shift above sea-level by tiny bit
    print("    time text location: ",location)
    print("")
    text_object.location = location

    # Set text material
    text_material = bpy.data.materials.new(name="TimeTextMaterial")
    text_material.use_nodes = True
    text_material.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.8, 0.8, 0.8, 1)
    text_material.node_tree.nodes["Principled BSDF"].inputs["Alpha"].default_value = 1   # transparency
    text_material.node_tree.nodes["Principled BSDF"].inputs["Specular"].default_value = 0  # Reduce shine
    text_material.shadow_method = 'NONE'  # Blender 3.x
    text_material.blend_method = 'BLEND'  # 'Alpha Blend'
    text_object.data.materials.append(text_material)
    #text_object.hide_render = False # Ensure it's not hidden from render


def set_scene() -> None:
    """
    sets scene options
    """
    global blender_img_resolution_X,blender_img_resolution_Y
    global use_cycles_renderer

    print("")
    print("  setting scene options")

    # gets scene
    scene = bpy.context.scene

    # renderer options
    if use_cycles_renderer:
        # Set render engine to Cycles
        scene.render.engine = 'CYCLES'

        ## adjust settings for faster rendering w/ cycles
        #scene.cycles.device = 'GPU'          # GPU rendering
        # Set tile size
        scene.cycles.tile_size = 512 # 256
        # Lower the number of samples
        scene.cycles.samples = 128 # Adjust the number of samples as needed
        # Enable denoising
        scene.cycles.use_denoising = True

        # Adjust light path bounces
        scene.cycles.max_bounces = 4  # Adjust bounce values for diffuse, glossy, transmission, etc.

        ## additional parameters to check for better performance:
        # Use Branched Path Tracing integrator
        #scene.cycles.progressive = 'BRANCHED_PATH'
        # Adjust Branched Path Tracing settings
        #scene.cycles.use_square_samples = False  # Enable square samples for the branched path tracing
        #scene.cycles.diffuse_samples = 3  # Adjust samples per type (diffuse, glossy, etc.)
        #scene.cycles.glossy_samples = 3
        #scene.cycles.transmission_samples = 3
        #scene.cycles.ao_samples = 3
        #scene.cycles.mesh_light_samples = 3
        #scene.cycles.subsurface_samples = 3

    else:
        # Set render engine to Eevee
        scene.render.engine = 'BLENDER_EEVEE'

        # ambient occlusion
        scene.eevee.use_gtao = True
        scene.eevee.gtao_distance = 0.01
        scene.eevee.gtao_factor = 1.0

        # turns on screen space reflections
        scene.eevee.use_ssr = True
        # refraction
        scene.eevee.use_ssr_refraction = True
        scene.eevee.ssr_thickness = 0.1

        # turns on bloom
        #scene.eevee.use_bloom = True

    # render resolution
    scene.render.resolution_x = blender_img_resolution_X
    scene.render.resolution_y = blender_img_resolution_Y
    # Set the render percentage (optional, for scaling down the final output)
    #scene.render.resolution_percentage = 50  # Adjust as needed

    # sets background color
    bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = world_background_color

    print("    renderer: ",scene.render.engine)
    print("")


def render_animation() -> None:
    """
    renders a movie animation
    """
    global use_animation_dive_in,use_animation_rotation
    global suppress_renderer_output
    global camera_elevation_offset,camera_y_offset,camera_angle_depth_degree
    global z_elevation

    print("animation:")

    # gets scene
    scene = bpy.context.scene

    # initializes frames
    number_of_keyframes = 0
    total_frames = 0

    start_frame = -1
    end_frame = -1

    scene.frame_start = -1
    scene.frame_end = -1

    ## dive-in keyframes
    if use_animation_dive_in:
        # setup keyframes
        number_of_keyframes += 4
        keyframe_interval = 20

        total_frames = (number_of_keyframes-1) * keyframe_interval

        print("  dive-in:")
        print("    keyframe interval      = ",keyframe_interval)
        print("    number of keyframes    = ",number_of_keyframes)
        print("")

        # define a frame timeline
        start_frame = 0
        inter_frame_1 = keyframe_interval
        inter_frame_2 = 2*keyframe_interval
        end_frame   = total_frames

        scene.frame_start = start_frame
        scene.frame_end = end_frame

        # moves camera
        cam = bpy.data.objects["Camera"]

        # 1. start frame
        cam.keyframe_insert("location", frame=start_frame)        # (0, -4, 4)
        cam.keyframe_insert("rotation_euler", frame=start_frame)  # (44.0 * DEGREE_TO_RAD, 0, 0)

        # 2. intermediate frame
        cam.location = (0, -2.5, 1)
        cam.rotation_euler = (60.0 * DEGREE_TO_RAD, 0, 0)
        cam.keyframe_insert("location", frame=inter_frame_1)
        cam.keyframe_insert("rotation_euler", frame=inter_frame_1)
        # adjust clip range
        cam.data.clip_start = 0.1
        cam.data.keyframe_insert("clip_start", frame=inter_frame_1)
        #cam.data.clip_end = 10.0
        #cam.data.keyframe_insert("clip_end", frame=inter_frame_1)

        # 3. intermediate frame
        x_cam = 0
        y_cam = -1.2
        z_cam = z_center + z_elevation + camera_elevation_offset
        cam.location = (x_cam, y_cam, z_cam)
        cam.rotation_euler = (camera_angle_depth_degree * DEGREE_TO_RAD, 0, 0)
        cam.keyframe_insert("location", frame=inter_frame_2)
        cam.keyframe_insert("rotation_euler", frame=inter_frame_2)
        # adjust clip range
        cam.data.clip_start = 0.01
        cam.data.keyframe_insert("clip_start", frame=inter_frame_2)
        #cam.data.clip_end = 10.0
        #cam.data.keyframe_insert("clip_end", frame=inter_frame_2)

        # 4. end frame
        x_cam = x_center
        y_cam = y_center + camera_y_offset
        z_cam = z_center + z_elevation + camera_elevation_offset
        cam.location = (x_cam, y_cam, z_cam) # Ending position
        cam.rotation_euler = (camera_angle_depth_degree * DEGREE_TO_RAD, 0, 0)
        cam.keyframe_insert("location", frame=end_frame)
        cam.keyframe_insert("rotation_euler", frame=end_frame)

    ## rotation
    if use_animation_rotation:
        # appends frames for rotation
        frames = 60  # Number of frames for the animation
        keyframe_interval = 5
        number_of_keyframes += frames

        total_frames += frames * keyframe_interval

        print("  rotation:")
        print("    keyframe interval      = ",keyframe_interval)
        print("    number of keyframes    = ",frames)
        print("")

        if start_frame == -1:
            start_frame = 0
        else:
            start_frame = end_frame + 1 * keyframe_interval

        if end_frame == -1:
            end_frame = total_frames
        else:
            end_frame = total_frames + 1 * keyframe_interval

        # updates end frame count
        scene.frame_end = end_frame

        # Set the rotation parameters
        rotation_degrees = 360  # Total rotation in degrees

        # moves camera around the z-axis
        cam = bpy.data.objects['Camera']

        #if scene.frame_start == -1:
        #    # add start frame
        #    scene.frame_start = start_frame
        #    # start frame
        #    cam.keyframe_insert("location", frame=start_frame)        # (0, -4, 4)
        #    cam.keyframe_insert("rotation_euler", frame=start_frame)  # (44.0 * DEGREE_TO_RAD, 0, 0)

        # camera offset from z-axis
        x_0 = cam.location[0]
        y_0 = cam.location[1]
        z_0 = cam.location[2]

        offset_distance = ((x_0-x_center)**2 + (y_0-y_center)**2)**0.5
        angle_0 = atan2(y_0-y_center,x_0-x_center)

        # adjust clip range since rotations will be done after dive-in in more close-up views
        # to check if we zoomed in, we will use the camera's z-position which for an overview is set to +4
        if z_0 < 2.0:
            cam.data.clip_start = 0.01
            #cam.data.clip_end = 10.0

        # Set keyframes for rotation
        for frame in range(0,frames+1):
            angle = angle_0 + radians(frame * (rotation_degrees / frames))

            # Calculate camera position around the Z-axis
            x = offset_distance * cos(angle)
            y = offset_distance * sin(angle)

            # rotate around center location
            x += x_center
            y += y_center

            # new camera position
            camera_location = Vector((x, y, z_0))  # Maintain Z position

            # Set camera rotation towards the origin
            look_at = Vector((x_center, y_center, z_center + z_elevation))  # Origin point on mesh

            # gets camera rotation angle
            direction = look_at - camera_location
            rot_quat = direction.to_track_quat('-Z', 'Y')  # Point -Z axis to direction
            camera_rotation = rot_quat.to_euler()

            # there seems to be a problem with frame interpolation when the angles going past 2 PI, i.e,
            # having the angle set between [-pi,pi] from the to_euler() output and interpolating a step
            # between +pi to -pi. to avoid we always increase the angle going from [0,2pi[
            #
            # in our case here, we change & increase only the third euler angle, while the first 2 stay fixed
            # when the camera rotates around the z-axis.
            if camera_rotation[2] < 0.0: camera_rotation[2] += 2.0 * PI

            # checks difference in current camera rotation and new target rotation
            # makes sure that z-angle always increases from one frame to another,
            # otherwise the linear interpolation will between frames will lead to weird results
            if cam.rotation_euler[2] > camera_rotation[2]: camera_rotation[2] += 2.0 * PI

            # info
            if frame == 0:
                print("    initial frame: camera angles = {:.2f} / {:.2f} / {:.2f} ".format(np.rad2deg(camera_rotation.x),
                                np.rad2deg(camera_rotation.y),np.rad2deg(camera_rotation.z)))
                print("")

            frame_number = start_frame + frame * keyframe_interval
            scene.frame_set(frame_number)

            # Set camera location and keyframe for each frame
            cam.rotation_euler = camera_rotation
            cam.location = camera_location

            cam.keyframe_insert("location", frame=frame_number)
            cam.keyframe_insert("rotation_euler", frame=frame_number)

            #debug
            #print("frame: ",frame,frame_number,"angle",angle,cam.location,"rotation",cam.rotation_euler)

    print("  total number of keyframes = ",number_of_keyframes)
    print("  total number of frames    = ",total_frames)
    print("")

    # smoothing move
    #for fcurve in cam.animation_data.action.fcurves:
    #    for keyframe_point in fcurve.keyframe_points:
    #        keyframe_point.interpolation = 'LINEAR'

    # render settings
    # see: https://docs.blender.org/api/current/bpy.types.FFmpegSettings.html#bpy.types.FFmpegSettings
    # ffmpeg
    scene.render.image_settings.file_format = 'FFMPEG'   # 'AVI_JPEG', 'FFMPEG'
    # codec
    scene.render.ffmpeg.codec = 'H264'                  # compression: 'MPEG4', 'H264'
    #scene.render.ffmpeg.constant_rate_factor = 'HIGH'
    scene.render.ffmpeg.format = 'MPEG4'                # output file container format
    # frames per second (blender default is 24)
    scene.render.fps = 20

    # user output
    print("  movie fps                    = ",scene.render.fps)
    print("  movie format                 = ",scene.render.ffmpeg.format )
    print("")

    # output movie
    name = './out.anim'
    filename = name + '.mp4'
    scene.render.filepath = filename
    scene.render.use_file_extension = False

    # timing
    tic = time.perf_counter()

    # redirect stdout to null
    print("rendering animation: {} ...".format(filename))
    print("")
    # to avoid long stdout output by the renderer:
    #    Fra:1 Mem:189.14M (Peak 190.26M) | Time:00:00.68 | Syncing Sun
    #    Fra:1 Mem:189.14M (Peak 190.26M) | Time:00:00.68 | Syncing Camera
    #  ..
    suppress = suppress_renderer_output
    with SuppressStream(sys.stdout,suppress):
        # render animation
        bpy.ops.render.render(animation=True)
    print("  written to: ",filename)
    print("")

    # timing
    toc = time.perf_counter()
    if toc - tic < 100.0:
        print("elapsed time for animation render is {:0.4f} seconds\n".format(toc - tic))
    else:
        min = int((toc-tic)/60.0)
        sec = (toc - tic) - min * 60
        print("elapsed time for animation render is {} min {:0.4f} sec\n".format(min,sec))


def render_image() -> None:
    """
    renders a single image
    """
    global suppress_renderer_output

    print("image:")

    # gets scene
    scene = bpy.context.scene

    # get render settings
    render = scene.render

    ## GPU rendering
    # note: some of these settings either won't work or are only available for Cycle engine
    #       either way, it seems they slow down the rendering and the default installation is fine for now...
    #
    #scene.render.device = 'GPU'            # doesn't work w/ Blender v3.6
    #scene.render.engine = 'BLENDER_EEVEE'  # uses by default installation the OpenGL/GPU; no further GPU settings for Eevee
    #
    # for cycles - turning this on seems to slow down the rendering (?)
    #scene.render.engine = 'CYCLES'
    #scene.cycles.device = 'GPU'
    # configure GPU devices (optional, adjust based on your setup)
    #bpy.context.preferences.addons["cycles"].preferences.compute_device_type = 'METAL'  # 'CUDA' or 'METAL' for Apple Silicon
    #bpy.context.scene.cycles.device_type = 'GPU'
    # info - get_devices() to let Blender detects GPU device
    #bpy.context.preferences.addons["cycles"].preferences.get_devices()
    #print("Cycles compute device : ",bpy.context.preferences.addons["cycles"].preferences.compute_device_type)
    #for device in bpy.context.preferences.addons["cycles"].preferences.devices:
    #    device["use"] = 1 # Using all devices, include GPU and CPU
    #    print("device: ",device["name"]," - use: ",device["use"])

    # output image
    scene.render.image_settings.file_format = 'JPEG'   # 'PNG'
    name = './out'
    filename = name + '.jpg'
    scene.render.filepath = filename
    scene.render.use_file_extension = False

    # timing
    tic = time.perf_counter()

    # redirect stdout to null
    print("  rendering image: {} ...".format(filename))
    print("")
    # to avoid long stdout output by the renderer:
    #    Fra:1 Mem:462.85M (Peak 462.85M) | Time:00:00.04 | ..
    #  ..
    suppress = suppress_renderer_output
    with SuppressStream(sys.stdout,suppress):
        # Render Scene and store the scene
        bpy.ops.render.render(write_still=True)
    print("  written to: ",filename)
    print("")

    # timing
    toc = time.perf_counter()
    if toc - tic < 100.0:
        print("elapsed time for image render is {:0.4f} seconds\n".format(toc - tic))
    else:
        min = int((toc-tic)/60.0)
        sec = (toc - tic) - min * 60
        print("elapsed time for image render is {} min {:0.4f} sec\n".format(min,sec))


def render_blender_scene(title: str="", animation: bool=False) -> None:

    ## blender scene setup
    print("Setting up blender scene...")
    print("")

    ## camera
    add_camera(title)

    ## light
    add_light()

    ## background plane
    add_plane()

    # time
    add_time_text()

    # text object
    add_title(title)

    # save scene and render options
    set_scene()

    # rendering
    if animation:
      # movie rendering
      render_animation()
    else:
      # image rendering
      render_image()

    # save blend file
    dir = os.getcwd()
    name = 'out.blend'

    filename = dir + "/" + name
    bpy.ops.wm.save_as_mainfile(filepath=filename)

    print("")
    print("  saved blend file: ",filename)
    print("")


# main routine
def plot_with_blender(vtk_file: str="", image_title: str="", colormap: int=0, color_max: float=None,
                      buildings_file: str="", borders_file: str="", locations_file: str="", animation: bool=False) -> None:
    """
    renders image for (earth) sphere with textures
    """
    # set current directory, in case we need it to load files
    dir = os.getcwd()
    print("current directory: ",dir)
    print("")

    # converts .vtu to .obj file for blender to read in
    obj_file = convert_vtk_to_obj(vtk_file,colormap,color_max)

    # setup mesh node with shaders
    create_blender_setup(obj_file)

    # add buildings
    add_blender_buildings(buildings_file)

    # add borders
    add_borders(borders_file)

    # add locations
    add_location_labels(locations_file)

    # save blender scene
    render_blender_scene(image_title,animation)


def usage() -> None:
    print("usage: ./plot_with_blender.py [--vtk_file=file] [--title=my_mesh_name] [--colormap=val] [--color-max=val] [--vertical-exaggeration=val]")
    print("                              [--buildings=file] [--shift-buildings=val] [--borders=file]")
    print("                              [--locations=file] [--utm-zone=ZoneNumber] [--moon] [--time='val']")
    print("                              [--no-sea-level] [--transparent-sea-level] [--sea-level-separation=val]")
    print("                              [--centered] [--closeup] [--small] [--anim] [--matte]")
    print("                              [--with-cycles/--no-cycles] [--suppress]")
    print("                              [--help]")
    print("  with")
    print("     --vtk_file                - input mesh file (.vtk, .vtu, .inp)")
    print("     --title                   - title text (added to image rendering)")
    print("     --colormap                - color map type: 0==VTK        / 1==topo      / 2==lisbon    / 3==lajolla       / 4==lipari")
    print("                                                 5==davos      / 6==turku     / 7==berlin    / 8==grayC         / 9==snow")
    print("                                                10==shakeGreen /11==shakeRed  /12==shakeUSGS /13==shakeUSGSgray /14==shakeUSGSblack")
    print("                                                15==shakeBlack /16==gist_earth")
    print("                                                 (default is shakeUSGSgray for shakemaps)")
    print("     --color-max               - fixes maximum value of colormap for moviedata to val, e.g., 1.e-7)")
    print("     --vertical-exaggeration   - factor to scales vertical dimension")
    print("")
    print("     --buildings               - mesh file (.ply) with buildings to visualize for the area")
    print("     --shift-buildings         - moves buildings up, e.g., a factor 0.0005 (default is no shift==0.0)")
    print("     --borders                 - AVS borders file (.inp) with UTM border lines to visualize on mesh")
    print("     --locations               - file with location labels (using a format: #name #lat #lon)")
    print("     --time                    - a time string to add to the display (e.g. '13:24:00')")
    print("     --white-labels            - use white text for location labels")
    print("     --utm_zone                - use specified UTM zone number (1-60) with (+) for Northern (-) for Southern hemisphere (e.g., -58)")
    print("     --moon                    - use Moon projections (LTM/LPS) instead of UTM")
    print("     --no-sea-level            - turns off sea-level plane")
    print("     --transparent-sea-level   - turns on transparency for sea-level plane")
    print("     --sea-level-separation    - shift factor to move up/down points at sea-level for better sea-level separation")
    print("")
    print("     --centered                - centered camera view (looking from top straight down)")
    print("     --closeup                 - sets camera view closer to center of model")
    print("     --small                   - turns on small images size (400x600px) for preview")
    print("     --anim                    - turns on movie animation (dive-in and rotation by default, use --no-rotation or --no-dive-in to turn off)")
    print("     --matte                   - uses matte material appearance for mesh")
    print("")
    print("     --with-cycles/--no-cycles - turns on/off CYCLES renderer (default is off, using BLENDER_EEVEE)")
    print("     --suppress                - suppress renderer output (default is off)")
    print("     --help                    - this help for usage...")
    sys.exit(1)


if __name__ == '__main__':
    # init
    vtk_file = ""
    image_title = ""
    color_max = None
    colormap = -1
    vertical_exaggeration = None
    sea_level_separation_shift = None
    buildings_file = ""
    locations_file = ""
    borders_file = ""
    animation = False

    # reads arguments
    #print("\nnumber of arguments: " + str(len(sys.argv)))
    i = 0
    for arg in sys.argv:
        i += 1
        #print("argument "+str(i)+": " + arg)
        # get arguments
        if "--help" in arg:
            usage()
        elif "--anim" in arg:
            animation = True
        elif "--background-dark" in arg:
            world_background_color = (0.05,0.05,0.05,1)  # dark world background
        elif "--background-blue" in arg:
            world_background_color = (0.05,0.05,0.1,1)   # dark-blue world background
        elif "--background-black" in arg:
            world_background_color = (0,0,0,1)           # black world background
        elif "--borders=" in arg:
            borders_file = arg.split('=')[1]
        elif "--buildings=" in arg:
            buildings_file = arg.split('=')[1]
        elif "--centered" in arg:
            centered_view = True
        elif "--closeup" in arg:
            close_up_view = True
        elif "--color-max" in arg:
            color_max = float(arg.split('=')[1])
        elif "--colormap" in arg:
            colormap = int(arg.split('=')[1])
        elif "--locations=" in arg:
            locations_file = arg.split('=')[1]
        elif "--matte" in arg:
            use_matte_material = True
        elif "--moon" in arg:
            use_moon_ltm = True
        elif "--no-cycles" in arg:
            use_cycles_renderer = False
        elif "--no-dive-in" in arg:
            use_animation_dive_in = False
        elif "--no-rotation" in arg:
            use_animation_rotation = False
        elif "--no-sea-level" in arg:
            use_sea_level_plane = False
        elif "--sea-level-separation=" in arg:
            sea_level_separation = float(arg.split('=')[1])
        elif "--shift-buildings" in arg:
            shift_building_baseline = float(arg.split('=')[1])
        elif "--small" in arg:
            blender_img_resolution_X = 600
            blender_img_resolution_Y = 400
        elif "--suppress" in arg:
            suppress_renderer_output = True
        elif "--time=" in arg:
            time_string = arg.split('=')[1]
        elif "--title=" in arg:
            image_title = arg.split('=')[1]
        elif "--transparent-sea-level" in arg:
            use_transparent_sea_level_plane = True
        elif "--utm_zone=" in arg:
            utm_zone = int(arg.split('=')[1])
            if abs(utm_zone) < 1 or utm_zone > 60:
                print(f"Invalid UTM zone entered: {utm_zone} - Please use zones from +/- [1,60]")
                sys.exit(1)
        elif "--vertical-exaggeration=" in arg:
            vertical_exaggeration = float(arg.split('=')[1])
        elif "--vtk_file=" in arg:
            vtk_file = arg.split('=')[1]
        elif "--white-labels" in arg:
            location_labels_color = (0.9,0.9,0.9,1)      # white location labels
        elif "--with-cycles" in arg:
            use_cycles_renderer = True
        elif i >= 8:
            print("argument not recognized: ",arg)

    # sets default colormap
    if colormap == -1:
        # for shakemaps
        if 'shaking' in vtk_file or 'shakemap' in vtk_file:
            colormap = 13 # shakeUSGSgray
        else:
            colormap = 0  # VTK diverging red-blue

    # logging
    cmd = " ".join(sys.argv)
    filename = './plot_with_blender.log'
    with open(filename, 'a') as f:
      print("command call --- " + str(datetime.datetime.now()),file=f)
      print(cmd,file=f)
      print("command logged to file: " + filename)

    # main routine
    plot_with_blender(vtk_file,image_title,colormap,color_max,buildings_file,borders_file,locations_file,animation)
