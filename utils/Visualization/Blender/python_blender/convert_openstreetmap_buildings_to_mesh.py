#!/usr/bin/env python
#
# downloads openstreetmap buildings and creates 3d mesh objects saved in .ply format
#
# required python modules:
# - numpy          (numpy==1.26.4)
# - scipy          (scipy==1.13.0)
# - pygmt          (pygmt==0.11.0)
# - pyproj         (pyproj==3.6.1)
# - xarray         (xarray==2024.3.0)
# - osmnx          (osmnx==1.9.2)
# - geopandas      (geopandas==0.14.3)
# - trimesh        (trimesh==4.2.4)
# - rtree          (Rtree==1.2.0)
# - mapbox_earcut  (mapbox-earcut==1.0.1)
#
#################################################################

import sys
import os
import time
import datetime
import math

try:
    import numpy as np
except:
    print("Failed importing numpy. Please make sure to have python with numpy working properly.")
    sys.exit(1)

# elevation
# GMT python interface
try:
    import pygmt
except:
    print("Error importing module `pygmt`")
    print("install by: > pip install -U pygmt")
    sys.exit(1)

# Coordinate projections
try:
    import pyproj
except:
    print("Failed importing pyproj.")
    sys.exit(1)
from pyproj import Transformer

# for topo arrays
import xarray as xr

# OSMnx
try:
    import osmnx as ox
except:
    print("Failed importing osmnx. Please install by: > pip install osmnx")
    sys.exit(1)

# GeoJSON file import
import geopandas as gpd

# meshing
try:
    import trimesh
except:
    print("Error importing module `trimesh`")
    print("install by: > pip install -U trimesh")
    sys.exit(1)

try:
    import rtree
except:
    print("Error importing module `rtree`")
    print("install by: > pip install -U rtree")
    sys.exit(1)

try:
    import mapbox_earcut
except:
    print("Error importing module `mapbox_earcut`")
    print("install by: > pip install -U mapbox_earcut")
    sys.exit(1)

# for mesh triangulation
from scipy.spatial import Delaunay

# for triangulation filtering checks
from matplotlib.patches import Polygon

#################################################################
## USER PARAMETERS

# minimum level height for buildings (when missing further informations)
building_height_min = 3.0   # m

# level height for roofs (when roof:levels is provided)
roof_level_height = 1.5   # m

#################################################################

## globals
transformer_to_utm = None
utm_zone = None
store_scene_file = False

def download_buildings(mesh_area: np.array):
    """
    downloads OpenStreetMap (OSM) building footprints and infos
    """
    # user output
    print("downloading OSM buildings...")
    print("  mesh area: ",mesh_area)
    print("")

    # mesh area min/max
    lon_min = mesh_area[0]
    lat_min = mesh_area[1]
    lon_max = mesh_area[2]
    lat_max = mesh_area[3]

    # check
    if lon_min > lon_max or lat_min > lat_max:
        print("Wrong mesh area format, please provide in a format (lon_min,lat_min,lon_max,lat_max)")
        sys.exit(1)

    print("  lat min/max: ",lat_min,lat_max)
    print("  lon min/max: ",lon_min,lon_max)
    print("")

    # define a bounding box
    #north, south, east, west = lat_max, lat_min, lon_max, lon_min
    #target_area = (north, south, east, west)
    # new version: uses bounding box (left, bottom, right, top)
    target_area = (lon_min, lat_min, lon_max, lat_max)

    # configuration
    ox.settings.use_cache = True
    ox.settings.log_console = True
    #ox.settings.max_query_area_size = 10 * 2500000000.0    # manual increase of default query size to reduce url requests

    if 1 == 0:
        # get graph
        G = ox.graph_from_bbox(bbox=target_area, network_type="drive")

        #print("debug: G ",G)
        fig, ax = ox.plot_graph(G, node_size=0)

        # graph option to get elevations
        #raster_file = "./DEM-rome-w46075_s10/w46075_s10.tif"
        #G = ox.elevation.add_node_elevations_raster(G,raster_file)

        #print("debug: G ",G)
        #fig, ax = ox.plot_graph(G, node_size=0)

        # get one color for each node, by elevation, then plot the network
        #nc = ox.plot.get_node_colors_by_attr(G, "elevation", cmap="plasma")
        #fig, ax = ox.plot_graph(G, node_color=nc, node_size=5, edge_color="#333333", bgcolor="k")

    # download only OSM buildings feature
    tags = {"building": True, "building:part": True}

    # download buildings as a GeoDataFrame
    gdf = ox.features_from_bbox(bbox=target_area, tags=tags)

    print("")
    print("  number of all buildings: ",gdf.shape[0])

    if gdf.size == 0:
        print("Error: no buildings were found, exiting...")
        sys.exit(1)

    print("  download done.")
    print("")

    return gdf

#
#----------------------------------------------------------------------------------------
#

def prepare_osm_building_data(gdf):
    """
    prepares and cleans up buildings data array
    """
    # user output
    print("preparing OSM building data...")
    print("")

    # extracts only meaningful infos for buildings

    #bld_list = ["geometry", "building:height", "building:levels"]
    #for key in gdf.keys():
    #    if key not in bld_list:
    #        print("gdf: removing key: ",key)
    #        gdf = gdf.drop(labels=key, axis=1)

    # all possible height related infos
    cols_extended = ["building:height", "building:min_level",  \
                     "building:levels", "building:levels:roof", "building:levels:underground", \
                     "building:shape", "building:part",  \
                     "roof:levels", "roof:level", "roof:type", "roof:orientation", "roof:direction", "roof:shape", \
                     "height", "min_height", "min_level", "max_level", "maxheight", \
                     "geometry", "nodes" ]

    # main infos only
    cols = [ "geometry", "name", "height", \
             "building:height", "building:levels", "building:part", \
             "roof:levels", "roof:shape", "roof:type", "roof:direction", "roof:orientation" ]

    # adds missing keys
    # for example if downloaded gdf has no "building:height" infos
    for key in cols:
        if key not in gdf.keys():
            gdf[key] = np.nan

    # very slow...
    #for key in gdf.keys():
    #    if key not in cols:
    #        #print("gdf: removing key: ",key)
    #        gdf = gdf.drop(labels=key, axis=1)
    # faster
    gdf = gdf[cols]

    print("  GeoDataFrame keys: ",gdf.keys())
    print("")

    # view just the polygon buildings
    #gdf[gdf["geometry"].type == "Polygon"].dropna(axis=1, how="any")
    # view just the roof shapes
    #gdf[gdf['roof:shape'].notnull()]

    # keep everything (Polygon,MultiPolygon,LineString,..) except single Points geometries
    #gdf = gdf[gdf.geometry.type != "Point"]
    # keep only MultiPolygon
    #gdf = gdf[gdf.geometry.type == "MultiPolygon"]
    # keep only Polygon
    #gdf = gdf[gdf.geometry.type == "Polygon"]

    # combined mask to filter only rows with Polygon or MultiPolygon geometry types
    mask = (gdf.geometry.type == "Polygon") | (gdf.geometry.geom_type == "MultiPolygon")
    gdf = gdf[mask]

    # Drop rows with NaN values
    #gdf = gdf.dropna(axis=1,how="all")

    # checks for valid geometry
    indices_to_remove = []
    for i,geom in enumerate(list(gdf.geometry)):
        # check flags - imported building's geometry should all be valid, but just in case
        if not geom.is_valid:
            print("  building {} has invalid geometry".format(i))
            indices_to_remove.append(i)

        # check centroid info
        lon = geom.centroid.x
        lat = geom.centroid.y
        if np.isnan(lon) or np.isnan(lat):
            print("  building {} has invalid centroid lon/lat {}/{} - geometry {}".format(i,lon,lat,geometry.centroid))
            indices_to_remove.append(i)

    print("  number of usable buildings : ",gdf.shape[0])
    print("  number of invalid buildings: ",len(indices_to_remove))
    print("")

    if len(indices_to_remove) > 0:
        print("")
        gdf.iloc[indices_to_remove]
        gdf = gdf.drop(indices_to_remove)
        print("  updated number of usable buildings: ",gdf.shape[0])
        print("")

    # clear out building entries without any height infos
    print("  clearing void data entries...")
    num_entries = len(gdf)
    icount = 0
    for row in gdf.itertuples():
        found = False
        if hasattr(row, 'height'): found = True
        if hasattr(row, 'building:height'): found = True
        if hasattr(row, 'building:levels'): found = True
        if hasattr(row, 'roof:height'): found = True
        if hasattr(row, 'roof:levels'): found = True

        # remove entry if no relevant infos found
        if not found:
            gdf = gdf.drop(row.Index)
            icount += 1
    print(f"  cleared {icount} entries out of {num_entries}")
    print("")

    print("  checking for valid data entries...")

    # number to output in steps of 10/100/1000 depending on how large the number of entries is
    num_entries = len(gdf)
    if num_entries > 100000:
        noutput_info = min(10000,int(10**np.floor(np.log10(num_entries))))
    else:
        noutput_info = min(1000,int(10**np.floor(np.log10(num_entries))))

    # makes sure that the building information related to numbers can be read in as a float or as NaN.
    # info tags that are strings won't need to be changed here.

    # function to clean and convert values
    def clean_value(val):
        if isinstance(val, str):
            # Shorten string to first valid substring
            val = val.split(';')[0].split("'")[0]
            # Check if the value is a valid number
            if val.replace(".", "", 1).isdigit():
                return float(val)
            else:
                return np.nan
        else:
            return val

    # Apply the cleaning function to the 'height',.. column using apply
    gdf['height'] = gdf['height'].apply(clean_value)
    gdf['building:height'] = gdf['building:height'].apply(clean_value)
    gdf['building:levels'] = gdf['building:levels'].apply(clean_value)
    gdf['roof:levels'] = gdf['roof:levels'].apply(clean_value)
    gdf['roof:direction'] = gdf['roof:direction'].apply(clean_value)

    # function to adjust direction values
    def clean_direction(d):
        if not np.isnan(d):  # Check if the value is not NaN
            if d < 0.0:
                d += 360.0
            elif d > 360.0:
                d -= 360.0
        return d

    # check roof direction is within [0,360]
    gdf['roof:direction'] = gdf['roof:direction'].apply(clean_direction)

    if 1 == 0:
        # explicit looping
        # note: using iterrows() leads to problems as for some entries, the index gets converted
        #       and calling gdf.at[index,'height'] = h fails...
        #       instead one can use itertuples() and check if the row has the attribute, but this is quite slow.

        # iterates over each row in the GeoDataFrame
        icount = 0

        #for index, obj in gdf.iterrows():
        for row in gdf.itertuples():
            # counter
            icount += 1

            # user output
            if icount % noutput_info == 0:
                print("    entry: {} out of {}".format(icount,num_entries))

            # check height
            if hasattr(row, 'height'):
                #val = obj['height']
                val = getattr(row, 'height')
                if isinstance(val,str):
                    # shorten string, for example: "1.0;2.0" -> "1.0"
                    if val.find(";") >= 0: val = val[0:val.find(";")]
                    if val.find("'") >= 0: val = val[0:val.find("'")]
                    # check if its a valid number, for example "16.2" -> "162" -> yes, "q" -> "q" -> no
                    if val.replace(".", "", 1).isdigit():
                        h = float(val)
                    else:
                        h = np.nan
                else:
                    h = val
                #gdf.at[index,'height'] = h
                gdf.loc[row.Index, 'height'] = h

            #try:
            #except:
            #    # geodataframe entry misses column height - delete this entry
            #    #gdf = gdf.drop(index)
            #    #debug
            #    print("warning: encountered wrong index: ",index," and value h: ",h)
            #    pass
            #    #continue  # skip to next object

            # check building height
            #if hasattr(row, 'building:height'):
            if hasattr(obj, 'building:height'):
                val = obj['building:height']
                #val = getattr(row, 'building:height')

                if isinstance(val,str):
                    # shorten string, for example: "1.0;2.0" -> "1.0"
                    if val.find(";") >= 0: val = val[0:val.find(";")]
                    if val.find("'") >= 0: val = val[0:val.find("'")]
                    # check if its a valid number, for example "16.2" -> "162" -> yes, "q" -> "q" -> no
                    if val.replace(".", "", 1).isdigit():
                        h = float(val)
                    else:
                        h = np.nan
                else:
                    h = val
                gdf.at[index,'building:height'] = h
                #gdf.loc[row.Index, 'building:height'] = h

            # check levels
            #if hasattr(row, 'building:levels'):
            if hasattr(obj, 'building:levels'):
                val = obj['building:levels']
                #val = getattr(row, 'building:levels')

                if isinstance(val,str):
                    # shorten string, for example: "1;2" -> "1"
                    if val.find(";") >= 0: val = val[0:val.find(";")]
                    if val.find("'") >= 0: val = val[0:val.find("'")]
                    # check if its a valid number, for example "16.2" -> "162" -> yes, "q" -> "q" -> no
                    if val.replace(".", "", 1).isdigit():
                        lev = float(val)
                    else:
                        lev = np.nan
                else:
                    lev = val
                gdf.at[index,'building:levels'] = lev
                #gdf.loc[row.Index, 'building:levels'] = lev

            # check roof levels
            #if hasattr(row, 'roof:levels'):
            if hasattr(obj, 'roof:levels'):
                val = obj['roof:levels']
                #val = getattr(row, 'roof:levels')

                if isinstance(val,str):
                    # shorten string, for example: "1;2" -> "1"
                    if val.find(";") >= 0: val = val[0:val.find(";")]
                    if val.find("'") >= 0: val = val[0:val.find("'")]
                    # check if its a valid number, for example "16.2" -> "162" -> yes, "q" -> "q" -> no
                    if val.replace(".", "", 1).isdigit():
                        lev = float(val)
                    else:
                        lev = np.nan
                else:
                    lev = val
                gdf.at[index,'roof:levels'] = lev
                #gdf.loc[row.Index, 'roof:levels'] = lev

            # check roof direction
            #if hasattr(row, 'roof:direction'):
            if hasattr(obj, 'roof:direction'):
                val = obj['roof:direction']
                #val = getattr(row, 'roof:direction')

                if isinstance(val,str):
                    # shorten string, for example: "1;2" -> "1"
                    if val.find(";") >= 0: val = val[0:val.find(";")]
                    if val.find("'") >= 0: val = val[0:val.find("'")]
                    # check if its a valid number, for example "16.2" -> "162" -> yes, "q" -> "q" -> no
                    if val.replace(".", "", 1).isdigit():
                        d = float(val)
                    else:
                        d = np.nan
                else:
                    d = val

                # directions seem to be within the range [0,360]
                if not np.isnan(d):
                    if d < 0.0: d += 360.0
                    if d > 360.0: d -= 360.0

                # re-set direction value
                gdf.at[index,'roof:direction'] = d
                #gdf.loc[row.Index, 'roof:direction'] = d

            #debug
            #name = "" if str(obj['name']) == 'nan' else obj['name']
            #print("debug: building ",icount,"H/B:H/B:L/B:P ",obj['height'],obj['building:height'],obj['building:levels'],obj['building:part'],name)
        print("")

    # plot footprints
    if 1 == 0:
        gdf_proj = ox.project_gdf(gdf)
        fig, ax = ox.plot_footprints(gdf_proj)

    # Save the GeoDataFrame to a GeoJSON file
    if 1 == 0:
        gdf = gdf.apply(lambda c: c.astype(str) if c.name != "geometry" else c, axis=0)
        filename = "./buildings.gjson"
        gdf.to_file(filename, driver="GeoJSON")
        print("  written to file: ",filename)
        # check other formats
        #import fiona
        #fiona.supported_drivers
        # Save the GeoDataFrame to a DXF file
        #gdf.to_file("./buildings.dxf", driver="DXF")
        # direct function to GeoJSON
        #gdf.to_json("./buildings.gjson", na='drop', driver="GeoJSON")

    print("  all building setup done.")
    print("")

    return gdf

#
#----------------------------------------------------------------------------------------
#

def get_topo_DEM(mesh_area: np.array):
    """
    downloads elevation data using GMT
    """
    # mesh area min/max
    lon_min = mesh_area[0]
    lat_min = mesh_area[1]
    lon_max = mesh_area[2]
    lat_max = mesh_area[3]

    # check
    if lon_min > lon_max or lat_min > lat_max:
        print("Wrong mesh area format, please provide in a format (lon_min,lat_min,lon_max,lat_max)")
        sys.exit(1)

    print("topography:")
    print("  lat min/max: ",lat_min,lat_max)
    print("  lon min/max: ",lon_min,lon_max)
    print("")

    # GMT python interface
    # https://www.pygmt.org/dev/api/generated/pygmt.datasets.load_earth_relief.html
    #
    # GMT region format: e.g. -R123.0/132.0/31.0/40.0
    gmt_region = [lon_min,lon_max,lat_min,lat_max]

    dim_max = max((lat_max-lat_min),(lon_max-lon_min))

    # determines resolution
    if dim_max > 5.0:
        # ETOPO1 (1-arc minute)
        # (topo_grid is an xarray data object)
        res = '01m'
    elif dim_max > 1.0:
        # SRTM 3S (3-arc sec)
        res = '03s'
    else:
        # SRTM 1S (1-arc sec)
        res = '01s' # '01s'

    print("  maximum dimension extent: ",dim_max,"(deg)")
    print("  using resolution        : ",res)
    print("")

    topo_grid = pygmt.datasets.load_earth_relief(resolution=res,
                                                 region=gmt_region,
                                                 registration="gridline")
    print("")

    # checks units are in m
    if not topo_grid.units == 'meters':
        print("Error: topo grid has invalid unit ",topo_grid.units)
        sys.exit(1)

    # debug
    #lat = 40.0
    #lon = 100.0
    #elevation = get_topo_elevation(lat,lon,topo_grid)
    #print("debug: lat/lon/elevation = ",lat,lon,elevation)

    # show/write figure plot
    if 1 == 1:
        fig = pygmt.Figure()
        fig.grdimage(grid=topo_grid,
                     cmap="haxby",
                     projection="M10c",
                     frame=True,
                    )
        fig.colorbar(frame=["x+lelevation", "y+lm"])
        #fig.show()

        # save figure as jpeg image
        name = "./output_topography.jpg"
        fig.savefig(name, crop=True, dpi=720)
        print("  figure plotted to: ",name)
        print("")

    print("  topo setup done.")
    print("")

    return topo_grid

#
#----------------------------------------------------------------------------------------
#

def get_topo_elevation(lat,lon,topo_grid):
    # gets elevation interpolated by nearest neighbor method
    elevation = topo_grid.interp(lon=lon,lat=lat,method="nearest")

    # extract simple float value from returned xarray object
    elevation_val = elevation.data

    return elevation_val

#
#----------------------------------------------------------------------------------------
#

def convert_coordinates_to_UTM(gdf, mesh_area: np.array, lons, lats):
    global transformer_to_utm

    # user output
    if transformer_to_utm == None:
        print("converting coordinates to UTM...")
        print("")
        # pyproj coordinate system info:
        #   WGS84                                          ==       EPSG:4326
        #   spherical mercator, google maps, openstreetmap ==       EPSG:3857
        #
        # we first need to determine the correct UTM zone to get the EPSG code.
        # for this, we take the mesh origin point and convert it to lat/lon (GPS) position.
        # we can then query the corresponding UTM zone for this position.

        # Coordinate Reference System
        ref_epsg = str(gdf.crs)

        # check reference is WGS84
        if ref_epsg != 'epsg:4326':
            print("Error: reference coordinate system not recognized: ",ref_epsg)
            sys.exit(1)

        # mesh area min/max
        lon_min = mesh_area[0]
        lat_min = mesh_area[1]
        lon_max = mesh_area[2]
        lat_max = mesh_area[3]

        # mesh origin position
        orig_lon = 0.5 * (lon_min + lon_max)
        orig_lat = 0.5 * (lat_min + lat_max)

        # user output
        print("  origin: lon/lat  = ",orig_lon,orig_lat)
        print("")

        # gets list of UTM codes
        utm_crs_list = pyproj.database.query_utm_crs_info(
                          datum_name="WGS 84",
                          area_of_interest=pyproj.aoi.AreaOfInterest(west_lon_degree=orig_lon,
                                                                     south_lat_degree=orig_lat,
                                                                     east_lon_degree=orig_lon,
                                                                     north_lat_degree=orig_lat))
        utm_code = utm_crs_list[0].code
        utm_epsg = "EPSG:{}".format(utm_code)

        # convert code to integer number and determine UTM zone number for info
        utm_code = int(utm_code)
        hemisphere_auto = 'N' if utm_code < 32700 else 'S'
        utm_zone_auto = utm_code - 32600 if hemisphere_auto == 'N' else utm_code - 32700

        print("  detected UTM code: ",utm_code," epsg: ", utm_epsg)
        print("           UTM zone: {}{}".format(utm_zone_auto,hemisphere_auto))
        print("")

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
        utm_x,utm_y = transformer_to_utm.transform(orig_lon,orig_lat)
        print("       -> UTM x/y  = ",utm_x,utm_y)
        print("          backward check: orig x/y = ",transformer_to_utm.transform(utm_x,utm_y,direction='INVERSE'))
        print("")

    # converts point coordinates
    x = lons
    y = lats

    # converts to UTM location (x and y)
    #point = np.array([x,y])
    x_utm,y_utm = transformer_to_utm.transform(x,y)

    return x_utm,y_utm

#
#----------------------------------------------------------------------------------------
#

def create_building_mesh(gdf, mesh_area: np.array, elevation, obj, meshes, icount_added: int):
    """
    creates a 3d mesh for building specified by obj
    """
    geometry = obj['geometry']

    # gets building height
    height = get_building_height(obj)

    # gets roof height and type
    roof_type, roof_height, roof_direction, roof_orientation = get_building_roof_type(obj)

    # creates extruded building from footprint polygon
    if geometry.geom_type == 'Polygon':
        # For a single Polygon
        poly = geometry

        # creates building mesh
        b_mesh = create_extruded_building(gdf,mesh_area,height,roof_type,roof_height,roof_direction,roof_orientation,elevation,poly)

        # add to overall mesh
        if not b_mesh is None:
            meshes.append(b_mesh)
            # counter
            icount_added += 1

    elif geometry.geom_type == 'MultiPolygon':
        # For MultiPolygon
        # create "building" for each polygon
        for poly in geometry.geoms:
            # creates building mesh
            b_mesh = create_extruded_building(gdf,mesh_area,height,roof_type,roof_height,roof_direction,roof_orientation,elevation,poly)

            # add to overall mesh
            if not b_mesh is None:
                meshes.append(b_mesh)
                # counter
                icount_added += 1

    else:
        print("Error: unsupported geometry type: ",geometry.geom_type)
        print("       object location centroid : ",geometry.centroid)
        raise ValueError("Unsupported geometry type")

    return meshes,icount_added

#
#----------------------------------------------------------------------------------------
#

def get_building_height(obj):
    """
    returns the height of the building (in m)
    """
    global building_height_min

    # initialize building height
    height = building_height_min

    #debug
    #has_height_info = False

    # check building height infos
    found = False
    h = obj['building:height']
    if not np.isnan(h):
          height = h
          found = True
          #debug
          #has_height_info = True
          #name = "" if str(obj['name']) == 'nan' else " - "+obj['name']
          #print("debug: obj found building:height ",obj['building:height'],name)
    else:
        # check height instead of building:height
        h = obj['height']
        if not np.isnan(h):
              height = h
              found = True
              #debug
              #has_height_info = True
              #name = "" if str(obj['name']) == 'nan' else " - "+obj['name']
              #print("debug: obj found height ",obj['height'],name)
        else:
            # check for levels
            lev = obj['building:levels']
            if not np.isnan(lev):
                height = lev * building_height_min
                found = True
                #debug
                #has_height_info = True
                #name = "" if str(obj['name']) == 'nan' else " - "+obj['name']
                #print("debug: obj found building:levels ",obj['building:levels'],name)

    # checks
    if np.isnan(height):
        height = building_height_min

    if height == 0.0:
        height = building_height_min

    # sets a minimum building height if specified height is < 1m
    if height <= 1.0: height = 1.0

    #debug
    #print("debug: building height = ",height,has_height_info)

    return height

#
#----------------------------------------------------------------------------------------
#

def get_building_roof_type(obj):
    """
    returns the type and height of the roof (in m)
    """
    global roof_level_height

    # initialize roof height
    # there will be no roof minimum height, only if there is more roof infos we will add it
    roof_type = 'flat'
    height = 0.0
    direction = np.nan
    orientation = 'along'

    # check building height infos
    lev = obj['roof:levels']
    if not np.isnan(lev):
          height = lev * roof_level_height
          #debug
          #has_height_info = True
          #name = "" if str(obj['name']) == 'nan' else " - "+obj['name']
          #print("debug: obj found roof:levels ",obj['roof:levels'],name)

    # check building direction infos
    val = obj['roof:direction']
    if not np.isnan(val):
          direction = val
          #debug
          #has_height_info = True
          #name = "" if str(obj['name']) == 'nan' else " - "+obj['name']
          #print("debug: obj found roof:direction ",obj['roof:direction'],obj['roof:shape'],obj['roof:type'],name)

    # checks shape info
    shape = obj['roof:shape']
    if isinstance(shape, str):
        # check shape
        if shape == 'flat':
            roof_type = 'flat'

        elif shape == 'gabled':
            roof_type = 'gabled'
            # adds minimum height for gabled roof
            if height == 0.0: height = roof_level_height

        elif shape == 'hipped':
            roof_type = 'hipped'
            # adds minimum height for gabled roof
            if height == 0.0: height = roof_level_height

        elif shape == 'side_hipped':
            roof_type = 'side_hipped'
            if height == 0.0: height = roof_level_height

        elif shape == 'skillion':
            roof_type = 'skillion'
            if height == 0.0: height = roof_level_height

        elif shape == 'sawtooth':
            roof_type = 'sawtooth'
            if height == 0.0: height = roof_level_height

        elif shape == 'quadruple_saltbox':
            roof_type = 'quadruple_saltbox'
            if height == 0.0: height = roof_level_height

        elif shape == 'pyramidal':
            roof_type = 'pyramidal'
            if height == 0.0: height = 2.0 * roof_level_height

        elif shape == 'dome':
            roof_type = 'dome'
            if height == 0.0: height = 2.0 * roof_level_height

        elif shape == 'cone':
            roof_type = 'cone'
            if height == 0.0: height = 2.0 * roof_level_height

        else:
            # not implemented yet - we will set it to flat
            roof_type = 'flat'
            height = 0.0
            # debug
            #name = "" if str(obj['name']) == 'nan' else " - "+obj['name']
            #print("building has unrecognized roof shape: ", shape,name)

    # checks type info
    # this seems mostly unset, i.e., nan, with a few ones which are set to 'flat' (probably by mistake?)
    type = obj['roof:type']
    if isinstance(type, str):
        #debug
        #name = "" if str(obj['name']) == 'nan' else " - "+obj['name']
        #print("debug: obj found roof type ",roof_type,height,obj['roof:type'],obj['roof:shape'],obj['roof:levels'])

        # check type
        if type == 'flat':
            roof_type = 'flat'

        elif type == 'gabled':
            roof_type = 'gabled'
            # adds minimum height for gabled roof
            if height == 0.0: height = roof_level_height

        elif type == 'hipped':
            roof_type = 'hipped'
            # adds minimum height for gabled roof
            if height == 0.0: height = roof_level_height

        elif shape == 'dome':
            roof_type = 'dome'
            if height == 0.0: height = 2.0 * roof_level_height

        elif shape == 'cone':
            roof_type = 'cone'
            if height == 0.0: height = 2.0 * roof_level_height

        else:
            # not implemented yet - we will set it to flat
            roof_type = 'flat'
            height = 0.0
            #debug
            #name = "" if str(obj['name']) == 'nan' else " - "+obj['name']
            #print("building has unrecognized roof type: ", type,name)

    # checks orientation info
    orient = obj['roof:orientation']
    if isinstance(orient, str):
        # check orientation
        if orient == 'along':
            orientation = 'along'

        elif orient == 'across':
            orientation = 'across'

        else:
            # invalid, will use 'along'
            orientation = 'along'
        #debug
        #name = "" if str(obj['name']) == 'nan' else " - "+obj['name']
        #print("building has roof orientation: ",orient,name)

    return roof_type, height, direction, orientation

#
#----------------------------------------------------------------------------------------
#

def create_extruded_building(gdf, mesh_area: np.array, height, roof_type, roof_height, roof_direction, roof_orientation, elevation, poly):
    """
    create building mesh by extrusion
    """
    # sets building height for extrusion
    # checks if roof height needs to be added
    # only for simple 'flat' roof, more complicated roofs will be added by an extra mesh later on
    if roof_type == 'flat' and roof_height > 0.0:
        # we only add roof height to building's height
        height += roof_height

    # checks if something to do
    if abs(height) < 1.e-5:
        return None

    # polygon points
    #print("debug: exteriors ",poly.exterior)
    #print("debug: interiors ",poly.interiors,len(poly.interiors))

    # exterior points
    points = list(poly.exterior.coords)

    # convert to numpy array
    #coords_outer = np.array(points)

    # check for interior points
    points_inner = []
    for interior in poly.interiors:
        for xy in interior.coords:
            points_inner.append(xy)

    # array for inner points
    coords_inner = np.array(points_inner)

    # combine all points for triangulation
    # add inner points to point list
    points.extend(points_inner)

    # overall points as numpy array
    coords = np.array(points)

    # get all lon/lat positions as separate arrays
    lons,lats = coords.T

    # convert lon/lat to UTM
    utm_x,utm_y = convert_coordinates_to_UTM(gdf,mesh_area,lons,lats)

    coords = np.column_stack((utm_x,utm_y))

    # inner polygon
    if len(coords_inner) > 0:
        # get all lon/lat positions as separate arrays
        lons,lats = coords_inner.T
        # convert lon/lat to UTM
        utm_x,utm_y = convert_coordinates_to_UTM(gdf,mesh_area,lons,lats)
        # store utm coordinates
        coords_inner = np.column_stack((utm_x,utm_y))

    # debug timing
    #if icount % noutput_info == 0:
    #    toc2 = time.perf_counter()
    #    print("            coord {}     - elapsed time {:0.4f} s".format(coords[0],(toc2-toc)))

    #debug
    #print("new coords: ",coords)

    # triangulate and extrude
    triangle_faces = get_triangulation(coords,coords_inner)

    # check
    if len(triangle_faces) == 0:
        return None

    # creates Mesh by extrusion
    b_mesh = trimesh.creation.extrude_triangulation(vertices=coords, faces=triangle_faces,
                                                    height=height)

    # checks number of vertices in mesh
    if len(b_mesh.vertices) == 0:
        return None

    #debug
    #if not np.isnan(roof_direction): print("debug: roof direction ",roof_direction,roof_type,roof_height,"location: ",points[0])

    # adds roof
    if roof_type != 'flat' and roof_height > 0.0:
        b_mesh = add_roof_mesh(b_mesh,height,roof_type,roof_height,roof_direction,roof_orientation)

    ## topography
    # Edit vertices to have correct z values
    verts = b_mesh.vertices.view(np.ndarray)
    # adds topographic elevation to z-coordinates
    verts[:,2] += elevation

    # check
    if np.isnan(verts).any():
        #print("error: mesh verts ",verts,"height",height,"elevation",elevation)
        return None

    # updates mesh vertices
    b_mesh.vertices = verts

    # updates face normals
    b_mesh.fix_normals()

    # debug
    #print("debug: mesh is empty {} watertight {} convex {} volume {} winding_consistent {}".format(
    #      b_mesh.is_empty,b_mesh.is_watertight,b_mesh.is_convex,b_mesh.is_volume,b_mesh.is_winding_consistent))

    return b_mesh

#
#----------------------------------------------------------------------------------------
#

def get_triangulation(coords, coords_inner):

    def in_polygon(point, polygon, hole):
        # Check if a point is inside the polygon (excluding holes)
        # check if inside
        path = Polygon(polygon).get_path()
        is_inside = path.contains_point(point)
        #debug
        #print("debug: point is in polygon: ",is_inside)

        # check if point is inside a hole area
        if len(hole) > 0:
            path = Polygon(hole).get_path()
            is_inside_hole = path.contains_point(point)
            #debug
            #print("debug: point is in hole   : ",is_inside_hole)
        else:
            is_inside_hole = False

        if is_inside and not is_inside_hole:
            return True
        else:
            return False

    def filter_triangles(tris, polygon, hole):
        """
        simple check if triangle center is inside the polygon with a possible hole in it
        """
        # Filter triangles to keep only those inside the polygon (excluding holes)
        filtered_triangles = []
        for face in tris.simplices:
            # gets triangle centroid point
            positions = tris.points[face]
            triangle_center = np.mean(positions, axis=0)

            #debug
            #print("triangle: ",face,positions,"center: ",triangle_center)

            # keep triangle if inside polygon
            if in_polygon(triangle_center, polygon, hole):
                filtered_triangles.append(face)

        return np.array(filtered_triangles)

    # get triangulation using all coordinate points
    tris = Delaunay(coords)

    # check if triangles are inside outer polygon path and not inside an inner "hole" path
    polygon = coords      # note: need full coords and not coord_outer here, otherwise no triangles will get selected
    hole = coords_inner

    # Filter triangles to keep only those inside the original polygon (excluding holes)
    triangle_faces = filter_triangles(tris, polygon, hole)

    return triangle_faces


#
#----------------------------------------------------------------------------------------
#

def add_roof_mesh(b_mesh, height, roof_type, roof_height, roof_direction, roof_orientation):
    """
    adds a roof mesh to the building
    """
    # we will scale and move a simple roof primitive to match the building bounding box
    # this assumes that the initial building footprint is close to a rectangle, if not the roof would not work

    def normalize_angle(angle):
        # Ensure the angle is within [0, 360)
        normalized_angle = angle % 360.0
        # If the angle is negative, wrap it around to the positive side
        if normalized_angle < 0.0: normalized_angle += 360.0
        return normalized_angle

    def angular_difference(angle1, angle2):
        # Calculate the raw angular difference
        raw_difference = angle1 - angle2
        # Calculate the normalized angular difference in the range [-180, 180]
        normalized_difference = (raw_difference + 180) % 360 - 180
        return normalized_difference

    # determine scale and orientation of base plane of mesh
    points = b_mesh.vertices

    #debug
    #print("debug: roof: ",roof_type," roof height: ", roof_height," direction: ", roof_direction," - position: ",b_mesh.centroid)
    #print("debug: mesh original points: ",points)

    # selects only base points
    points = points[points[:,2] < 0.1]

    # takes only x/y coordinates
    base_points = points[:, :2]

    #debug
    # center of base points
    #center = np.mean(base_points, axis=0)
    #print("debug: mesh base points: ",base_points," center: ",center)

    # gets oriented 2d bounding box
    # note: the transformation T will offset and rotate to align with the origin and xyz-axis
    #       in order to move the roof to the mesh location and align it, we will have to invert the transformation T.
    T, rectangle = trimesh.bounds.oriented_bounds_2D(base_points)

    #debug
    #print("debug: 2d oriented bounds: ",T,"rectangle: ",rectangle)

    # for gabled/hipped roofs
    # check if building footprint is similar to a simple rectangle
    if roof_type != 'dome' and roof_type != 'cone':
        # we'll use the area of the footprint to see if it matches the rectangle area,
        # otherwise the shape must be quite different
        area = rectangle[0] * rectangle[1]
        # area is volume / height for our extrusion mesh
        mesh_area = b_mesh.volume / height
        #debug
        #print("debug: mesh area = ",mesh_area," - rectangle area = ",area)
        # check
        if abs(mesh_area - area) > 0.1 * area:
            #debug
            #print("roof: building footprint won't match a rectangular shape by area:",mesh_area," rectangle area:",area,"roof type: ",roof_type)
            # just returns the building mesh as is
            return b_mesh

    # constructs a roof primitive
    roof_mesh, roof_height = get_roof_primitive(roof_type,roof_height,base_points,height)

    # checks roof mesh
    if roof_mesh is None:
        # no roof created, returns building mesh as is
        return b_mesh

    # changes rectangle size for spherical roofs
    if roof_type == 'dome' or roof_type == 'cone':
        # adapts rectangle size since dome/cone uses radius, not diameters
        rectangle = rectangle / 2.0

    # mask for points making up base of roof
    if roof_type == 'sawtooth':
        # only select corner points from the base (not base points in the middle of the edges)
        mask_roof_base_points = (roof_mesh.vertices[:,2] < 0.1) & (np.abs(roof_mesh.vertices[:,0]) > 0.499)
    else:
        # shapes like gabled/hipped/.. have only 4 corners at the base
        mask_roof_base_points = roof_mesh.vertices[:,2] < 0.1

    # trimesh seems to return the larger dimension first
    if rectangle[0] > rectangle[1]:
        # rotates roof by an additional 90 degree to have longer direction align with y-direction of the roof
        # looks nicer for gabled and hipped roofs
        angle = 90.0 # in degrees
        R = trimesh.transformations.rotation_matrix(np.deg2rad(angle),
                                                    direction=[0.0, 0.0, 1.0])
        roof_mesh.apply_transform(R)

    # default orientation is 'along' which aligns with the roof primitive,
    # we only need to rotation for 'across' orientations
    if roof_orientation == 'across':
        # rotates roof by another 90 degree to have shorter direction align with y-direction of the roof
        angle = 90.0 # in degrees
        R = trimesh.transformations.rotation_matrix(np.deg2rad(angle),
                                                    direction=[0.0, 0.0, 1.0])
        roof_mesh.apply_transform(R)

    # determines roof direction
    if roof_type == 'skillion' and not np.isnan(roof_direction):
        # roof direction
        # see: https://wiki.openstreetmap.org/wiki/Key:roof:direction

        # first determine building angle
        # Compute the oriented bounding box (OBB) of the mesh
        obb = b_mesh.bounding_box_oriented
        # Extract the rotation angles from the OBB transformation matrix
        angles = trimesh.transformations.euler_from_matrix(obb.transform, 'sxyz')
        # angles variable now contains the rotation angles around the X, Y, and Z axes
        # takes angle around the Z axis (vertical) corresponds to the North direction
        north_angle = angles[2]
        # Convert the angle to degrees if needed
        b_direction = np.degrees(north_angle)
        # directions are within the range [0,360]
        b_direction = normalize_angle(b_direction)

        # angle difference between building and roof direction (in range [-180,180])
        # example: building angle =   0.0  roof direction =  80.0  -> difference =  80.0
        #          building angle = 230.0  roof direction = 210.0  -> difference = -20.0
        angle_diff = angular_difference(roof_direction, b_direction)

        # direction tolerance in degrees
        TOL_DIR = 45.0

        # checks if we need to rotate roof
        # todo - this is still unclear what the correct rotation should be:
        #        the building angle seems to be dominated by the longest length of the building,
        #        and since the skillion roof primitive aligns with the longest length, the rain run-off direction which
        #        is the roof direction for skillion roofs goes perpendicular to this longest length direction.
        #
        # for now, we only rotate if there is a +90 degree offset between building and roof - meaning that the rain run-off
        # from the roof should be in opposite direction, assuming that the skillion roof primitive has a default run-off direction
        # of -90 degree with respect to the building direction.
        #
        # checks angle difference
        if angle_diff > 90.0 - TOL_DIR and angle_diff < 90.0 + TOL_DIR:
            # 90-degree difference
            #print("debug:   90-deg building angle ",b_direction,"roof direction: ",roof_direction,"difference = ",angle_diff)
            angle = 180.0 # (in degrees) rotate by 180-degree to have rain run-off angle in opposite direction

        elif angle_diff > 180.0 - TOL_DIR and angle_diff < 180.0 + TOL_DIR:
            # 180-degree difference
            #print("debug:  180-deg building angle ",b_direction,"roof direction: ",roof_direction,"difference = ",angle_diff)
            angle = 0.0 # no rotation

        elif angle_diff > -90.0 - TOL_DIR and angle_diff < -90.0 + TOL_DIR:
            # minus 90-degree difference
            #print("debug:  -90-deg building angle ",b_direction,"roof direction: ",roof_direction,"difference = ",angle_diff)
            angle = 0.0 # no rotation

        elif angle_diff > -180.0 - TOL_DIR and angle_diff < -180.0 + TOL_DIR:
            # minus 180-degree difference
            #print("debug: -180-deg building angle ",b_direction,"roof direction: ",roof_direction,"difference = ",angle_diff)
            angle = 0.0 # no rotation

        else:
            # angle difference between +/- 45 degrees
            #print("debug:   no rot building angle ",b_direction,"roof direction: ",roof_direction,"difference = ",angle_diff)
            angle = 0.0

        # apply rotation
        if abs(angle) > 1.0:
            R = trimesh.transformations.rotation_matrix(np.deg2rad(angle),
                                                        direction=[0.0, 0.0, 1.0])
            roof_mesh.apply_transform(R)

    # scale roof in x/y dimension, scales vertically later
    scale_vec = [rectangle[0], rectangle[1], 1.0]
    # scale
    roof_mesh.apply_scale(scale_vec)

    # convert to 3D transformation (around z-axis)
    T_3d = trimesh.transformations.planar_matrix_to_3D(T)
    # takes inverse to have rotation inverse than to_origin, i.e., we want to move and rotate our roof from the origin location
    # to the actual building mesh position
    roof_mesh.apply_transform(np.linalg.inv(T_3d))

    # scale roof vertically
    scale_vec = [1.0, 1.0, roof_height]
    # scale
    roof_mesh.apply_scale(scale_vec)

    # move to top of mesh (roof sits on top of building mesh)
    translate_vec = [0.0, 0.0, height]
    roof_mesh.apply_translation(translate_vec)

    # for gabled/hipped roofs
    # let's try to move the base points of the roof to match the top corners of the building
    if roof_type != 'dome' and roof_type != 'cone':
        # tries to match corners of roof with corners of building
        roof_verts = roof_mesh.vertices[mask_roof_base_points]

        # gets closest, distance and vertex id of closest point in building mesh
        closest, distance, tid = trimesh.proximity.closest_point(b_mesh,roof_verts)

        for i,point in enumerate(roof_verts):
            #debug
            #print("debug: roof base points: ",roof_verts[i],"closest: ",closest[i],distance[i],tid[i],"height: ",height)
            # sets roof point to closest point from building
            if distance[i] < height: roof_verts[i] = closest[i]
        # updates vertex position of roof base
        roof_mesh.vertices[mask_roof_base_points] = roof_verts

    # building position
    # gets maximum height of mesh
    #bbox = b_mesh.bounding_box
    #mesh_height = bbox.bounds[:,2].max()
    #debug
    #print("debug: mesh: centroid ",b_mesh.centroid,"bounds: ",b_mesh.bounds)
    #print("debug: mesh bounding box: ", bbox.bounds)
    #print("debug: mesh max height: ", mesh_height," input height: ",height)
    # safety check
    #if abs(mesh_height - height) > 1.e-5:
    #    print("Error: mesh_height and height do not match: ",mesh_height,height)
    #    sys.exit(1)

    # combine both meshes, building and roof into a single mesh
    combined_mesh = trimesh.util.concatenate([b_mesh,roof_mesh])

    # remove degenerate faces
    combined_mesh.update_faces(combined_mesh.nondegenerate_faces())

    # remove vertices not needed
    combined_mesh.remove_unreferenced_vertices()

    # removes duplicate vertices
    combined_mesh.merge_vertices()

    return combined_mesh

#
#----------------------------------------------------------------------------------------
#

def get_roof_primitive(roof_type, roof_height, base_points, height):
    """
    creates a geometric primitive of different roof types which will need to be scaled and positioned
    """
    # Openstreetmap roof shapes
    # see: https://wiki.openstreetmap.org/wiki/Key:roof:shape
    #
    if roof_type == 'gabled':
        # gabled roof
        # based on prism
        #              _____________
        #  z          /            /\
        #  |__ y     /            /  \
        # /         /            /    \
        # x        /____________/******
        #
        tris = np.array([[ 0.5, -0.5, 0.0],
                         [ 0.0, -0.5, 1.0],
                         [-0.5, -0.5, 0.0]])
        truncation_plane_origin = [0.0, 0.5, 0.0]
        truncation_plane_normal = [0.0, 1.0, 0.0]

        roof_mesh = trimesh.creation.truncated_prisms(tris,
                                                origin=truncation_plane_origin,
                                                normal=truncation_plane_normal)

    elif roof_type == 'hipped':
        # hipped roof
        # starts like a gabled roof top and then moves top points closer to each other
        # based on prism
        #              _____________
        #  z          /            /\
        #  |__ y     /            /  \
        # /         /            /    \
        # x        /____________/******
        #
        tris = np.array([[ 0.5, -0.5, 0.0],
                         [ 0.0, -0.5, 1.0],
                         [-0.5, -0.5, 0.0]])
        truncation_plane_origin = [0.0, 0.5, 0.0]
        truncation_plane_normal = [0.0, 1.0, 0.0]

        roof_mesh = trimesh.creation.truncated_prisms(tris,
                                                origin=truncation_plane_origin,
                                                normal=truncation_plane_normal)

        # moves top points closer together
        #              _________
        #  z          /         |\
        #  |__ y     /          | \
        # /         /           |  \
        # x        /____________|****
        #
        # hipped
        verts = roof_mesh.vertices
        # moves the two top nodes closer to each other
        verts[1] = [ 0.0, -0.2, 1.0]
        verts[4] = [ 0.0,  0.2, 1.0]
        # updates roof primitive
        roof_mesh.vertices = verts

    elif roof_type == 'side_hipped':
        # side-hipped roof
        # https://wiki.openstreetmap.org/wiki/Tag:roof:shape=side%20hipped?uselang=en

        # constructs a roof primitive
        # based on prism
        # hipped on one side only
        #
        #                 ._________
        #  z            .          /\
        #  |__ y      .           /  \
        # /         .            /    \
        # x        /____________/******
        #
        tris = np.array([[ 0.5, -0.5, 0.0],
                         [ 0.0, -0.5, 1.0],
                         [-0.5, -0.5, 0.0]])
        truncation_plane_origin = [0.0, 0.5, 0.0]
        truncation_plane_normal = [0.0, 1.0, 0.0]

        roof_mesh = trimesh.creation.truncated_prisms(tris,
                                                origin=truncation_plane_origin,
                                                normal=truncation_plane_normal)
        # side-hipped
        verts = roof_mesh.vertices
        # moves only one node closer to center
        verts[1] = [ 0.0, -0.2, 1.0]
        # updates roof primitive
        roof_mesh.vertices = verts

    elif roof_type == 'skillion':
        # skillion roof
        # based on prism
        #              _____________
        #  z          /            /|
        #  |__ y     /            / |
        # /         /            /  |
        # x        /____________/***|
        #
        tris = np.array([[ 0.5, -0.5, 0.0],
                         [-0.5, -0.5, 1.0],
                         [-0.5, -0.5, 0.0]])
        truncation_plane_origin = [0.0, 0.5, 0.0]
        truncation_plane_normal = [0.0, 1.0, 0.0]

        roof_mesh = trimesh.creation.truncated_prisms(tris,
                                                origin=truncation_plane_origin,
                                                normal=truncation_plane_normal)

    elif roof_type == 'sawtooth':
        # sawtooth roof
        # based on prism
        #              _____________ ____  ____
        #  z          /            /|   /|   /|
        #  |__ y     /            / |  / |  / |
        # /         /            /  | /  | /  |
        # x        /____________/***|/***|/***|
        #
        tris = np.array([[ 0.5,    -0.5, 0.0], # 1. triangle
                         [ 0.1666, -0.5, 1.0],
                         [ 0.1666, -0.5, 0.0],
                         [ 0.1666, -0.5, 0.0], # 2. triangle
                         [-0.1666, -0.5, 1.0],
                         [-0.1666, -0.5, 0.0],
                         [-0.1666, -0.5, 0.0], # 3.triangle
                         [-0.5,    -0.5, 1.0],
                         [-0.5,    -0.5, 0.0]
                         ])
        truncation_plane_origin = [0.0, 0.5, 0.0]
        truncation_plane_normal = [0.0, 1.0, 0.0]

        roof_mesh = trimesh.creation.truncated_prisms(tris,
                                                origin=truncation_plane_origin,
                                                normal=truncation_plane_normal)


    elif roof_type == 'pyramidal':
        # creates a pyramid
        # based on prism, merging top two points
        #              _____.________
        #  z          /   .  .     /\
        #  |__ y     /  .     .   /  \
        # /         / .        . /    \
        # x        /.___________./******
        #
        tris = np.array([[ 0.5, -0.5, 0.0],
                         [ 0.0, -0.5, 1.0],
                         [-0.5, -0.5, 0.0]])
        truncation_plane_origin = [0.0, 0.5, 0.0]
        truncation_plane_normal = [0.0, 1.0, 0.0]

        roof_mesh = trimesh.creation.truncated_prisms(tris,
                                                origin=truncation_plane_origin,
                                                normal=truncation_plane_normal)

        # merges the two top nodes
        verts = roof_mesh.vertices
        verts[1] = [ 0.0, 0.0, 1.0]
        verts[4] = [ 0.0, 0.0, 1.0]
        # updates roof primitive
        roof_mesh.vertices = verts

    elif roof_type == 'quadruple_saltbox':
        # quadruple-saltbox roof
        # based on prism, merging top two points, then cutting top tip off
        #
        #  z              *___*
        #  |__ y        *____*/\
        # /           .       . \
        # x          ._________./
        #
        tris = np.array([[ 0.5, -0.5, 0.0],
                         [ 0.0, -0.5, 1.0],
                         [-0.5, -0.5, 0.0]])
        truncation_plane_origin = [0.0, 0.5, 0.0]
        truncation_plane_normal = [0.0, 1.0, 0.0]

        roof_mesh = trimesh.creation.truncated_prisms(tris,
                                                origin=truncation_plane_origin,
                                                normal=truncation_plane_normal)

        # merges the two top nodes
        verts = roof_mesh.vertices
        verts[1] = [ 0.0, 0.0, 1.0]
        verts[4] = [ 0.0, 0.0, 1.0]

        # updates roof primitive
        roof_mesh.vertices = verts

        # cuts away top tip
        roof_mesh = roof_mesh.slice_plane(plane_origin=[0.0,0.0,0.5], plane_normal=[0,0,-1], cap=True)

        # remove degenerate faces to make mesh watertight
        roof_mesh.update_faces(roof_mesh.nondegenerate_faces())
        #debug
        #print("debug: quadruple_saltbox is watertight: ",roof_mesh.is_watertight)

    elif roof_type == 'dome':
        # dome roof
        # creates a sphere and cuts away lower hemisphere
        roof_mesh = trimesh.creation.uv_sphere(radius=1)

        # cuts away lower hemisphere (cap=True to have a watertight mesh)
        roof_mesh = roof_mesh.slice_plane(plane_origin=[0.0,0.0,0.0], plane_normal=[0,0,1], cap=True)

        # determines radius to inscribe a sphere for vertical scaling of the dome
        center, radius = trimesh.nsphere.minimum_nsphere(base_points)
        # limit radius (taking the full inscribed radius becomes fairly high)
        radius *= 0.6
        # limit radius to building height
        if radius > height: radius = height
        # update roof height
        roof_height = radius

    elif roof_type == 'cone':
        # cone roof
        # creates a cone
        roof_mesh = trimesh.creation.cone(radius=1,height=1)

        # determines radius to inscribe a sphere for vertical scaling of the cone
        center, radius = trimesh.nsphere.minimum_nsphere(base_points)
        # limit radius (taking the full inscribed radius becomes fairly high)
        radius *= 0.6
        # limit radius to building height
        if radius > height: radius = height
        # update roof height
        roof_height = radius

    else:
        #debug
        #print("Roof type not supported yet: ",roof_type," - will use a flat roof...")
        return None, roof_height

    return roof_mesh, roof_height

#
#----------------------------------------------------------------------------------------
#

def determine_building_elevations(gdf, topo_grid, mesh_area: np.array):
    """
    determines for each building its topographic elevation based on the building's centroid position
    """
    print("  determining building elevations")
    print("")

    # mesh area min/max
    lon_min = mesh_area[0]
    lat_min = mesh_area[1]
    lon_max = mesh_area[2]
    lat_max = mesh_area[3]

    # buildings
    num_buildings = gdf.shape[0]

    # arrays for building positions
    lons = np.zeros(num_buildings)
    lats = np.zeros(num_buildings)
    #building_elevations = np.zeros(num_buildings)

    # timing
    tic = time.perf_counter()

    # number to output in steps of 10/100/1000 depending on how large the number of polyfaces is
    if num_buildings > 100000:
        noutput_info = min(10000,int(10**np.floor(np.log10(num_buildings))))
    else:
        noutput_info = min(1000,int(10**np.floor(np.log10(num_buildings))))

    icount = 0
    for index, obj in gdf.iterrows():
        # counter
        icount += 1

        # user output
        if icount % noutput_info == 0:
            # timing
            toc = time.perf_counter()
            print("    building: {} out of {} - elapsed time {:0.2f} s".format(icount,num_buildings,(toc-tic)))

        # Access the geometry of the building
        geometry = obj['geometry']

        # check
        if not geometry.is_valid:
            print("Error: found invalid geometry ",geometry," this should have been checked before, exiting...")
            sys.exit(1)

        # centroid position lat / lon
        lon = geometry.centroid.x
        lat = geometry.centroid.y

        # check
        if np.isnan(lon) or np.isnan(lat):
            print("Error: building {} has invalid centroid lon/lat {}/{} - geometry {}, exiting...".format(index,lon,lat,geometry.centroid))
            sys.exit(1)

        # keep inside mesh area
        if lon < lon_min: lon = lon_min
        if lon > lon_max: lon = lon_max
        if lat < lat_min: lat = lat_min
        if lat > lat_max: lat = lat_max

        # stores centroid position
        lons[icount-1] = lon
        lats[icount-1] = lat

        # single value call - too slow for many buildings...
        # gets elevations (interpolated by nearest neighbor method)
        #building_elevations[icount-1] = get_topo_elevation(lat,lon,topo_grid)
    print("")

    # check
    if icount != num_buildings:
        print("Error: number of buildings {} should match all {}".format(icount,num_buildings))
        sys.exit(1)

    # gets all base elevations (interpolated by nearest neighbor method)
    # returned array will be a 2d array[nlats,nlons]
    # this becomes too big for many buildings...
    #building_elevations = get_topo_elevation(lats,lons,topo_grid)
    # extracts 1d array needed only (array[i,i])
    #building_elevations = np.diag(building_elevations)

    # using advanced indexing to have 1d array result
    # see: https://docs.xarray.dev/en/stable/user-guide/interpolation.html#advanced-interpolation
    x = xr.DataArray(lons, dims="z")
    y = xr.DataArray(lats, dims="z")

    # gets elevation
    # default linear interpolation
    elevation = topo_grid.interp(lon=x,lat=y)
    # interpolated by nearest neighbor method - faster, but less accurate
    #elevation = topo_grid.interp(lon=x,lat=y,method="nearest")

    # extract simple float value from returned xarray object
    building_elevations = elevation.data

    print("  building elevations:")
    print("    number of elevations = {}".format(len(building_elevations)))
    print("    min/max              = {} / {}".format(building_elevations.min(),building_elevations.max()))
    print("")

    # check elevations are valid
    for i in range(len(building_elevations)):
        ele = building_elevations[i]
        # checks nan values
        if np.isnan(ele):
            print("Error: building {} has invalid elevation {} at lon/lat {}/{}".format(i,building_elevations[i],lons[i],lats[i]))
            sys.exit(1)

    return building_elevations


#
#----------------------------------------------------------------------------------------
#

def create_3d_buildings(gdf, topo_grid, mesh_area: np.array) -> None:
    """
    creates for each building a 3d mesh and stores all as .ply output file
    """
    global building_height_min,store_scene_file

    # mesh area min/max
    lon_min = mesh_area[0]
    lat_min = mesh_area[1]
    lon_max = mesh_area[2]
    lat_max = mesh_area[3]

    # makes an initial conversion to UTM to get transformer
    center_lon = 0.5 * (lon_min + lon_max)
    center_lat = 0.5 * (lat_min + lat_max)
    # convert lon/lat to UTM
    center_utm_x,center_utm_y = convert_coordinates_to_UTM(gdf,mesh_area,center_lon,center_lat)

    # buildings
    num_buildings = gdf.shape[0]

    print("creating 3D buildings...")
    print("  total number of buildings: ",num_buildings)
    print("")

    # gets all elevations for the building bases
    building_elevations = determine_building_elevations(gdf,topo_grid,mesh_area)

    # 3d buildings
    print("  creating 3d buildings")
    print("")

    # meshes for .ply output
    meshes = []

    # number to output in steps of 10/100/1000 depending on how large the number of polyfaces is
    if num_buildings > 100000:
        noutput_info = min(10000,int(10**np.floor(np.log10(num_buildings))))
    else:
        noutput_info = min(1000,int(10**np.floor(np.log10(num_buildings))))

    # timing
    tic = time.perf_counter()

    # Iterate over each row in the GeoDataFrame
    icount = 0
    icount_added = 0
    for index, obj in gdf.iterrows():
        # counter
        icount += 1

        # user output
        if icount % noutput_info == 0:
            # timing
            toc = time.perf_counter()
            print("    building: {} out of {} - elapsed time {:0.2f} s".format(icount,num_buildings,(toc-tic)))

        # gets topographic base elevation of the building
        elevation = building_elevations[icount-1]

        # creates a mesh of the building
        meshes, icount_added = create_building_mesh(gdf,mesh_area,elevation,obj,meshes,icount_added)


    print("")
    print("  done. added {} building meshes".format(icount_added))
    print("")

    # create a reference sphere for checking lat/lon positioning
    if 1 == 0:
        b_mesh = trimesh.creation.uv_sphere(radius=50.0)
        # Define the desired position as a 3D vector (x, y, z)
        lon = 12.5 # station X1 or center_lon
        lat = 42.0 #               center_lat
        # convert lon/lat to UTM
        x,y = convert_coordinates_to_UTM(gdf,mesh_area,lon,lat)
        # gets base elevation (interpolated by nearest neighbor method)
        xele = topo_grid.interp(lon=lon,lat=lat,method="nearest")
        z = float(xele.data)   # converts to single float value
        if np.isnan(z): z = 0.0
        desired_position = [x,y,z]
        # Create a translation matrix based on the desired position
        translation_matrix = trimesh.transformations.translation_matrix(desired_position)
        # Apply the translation matrix to the sphere
        b_mesh.apply_transform(translation_matrix)
        # adds to all meshes
        meshes.append(b_mesh)

    ## store meshes as geometries in a scene
    if store_scene_file:
        # Create a dictionary where keys are unique names and values are the meshes.
        geometry_dict = {f'mesh_{i}': mesh for i, mesh in enumerate(meshes)}
        # Create the trimesh.Scene object using the dictionary
        scene = trimesh.Scene(geometry_dict)
        # save scene file
        filename = './output_buildings_utm.scene.obj'
        scene.export(filename, file_type='obj')
        print("scene:")
        print(f"  contains {len(scene.geometry)} geometries (meshes).")
        print(f"  written: ",filename)
        print("")

    # combine all meshes into single mesh
    mesh = trimesh.util.concatenate(meshes)

    # remove degenerate faces
    mesh.update_faces(mesh.nondegenerate_faces())

    # remove vertices not needed
    mesh.remove_unreferenced_vertices()

    # removes duplicate vertices
    mesh.merge_vertices()

    # user output
    print("mesh:")
    print("  total number of vertices = ",len(mesh.vertices))
    print("  total number of faces    = ",len(mesh.faces))
    print("")

    # check
    if len(mesh.vertices) == 0 or len(mesh.faces) == 0:
        print("")
        print("Mesh is empty, there's nothing to save. Please check your inputs...")
        print("")
        sys.exit(1)

    # user output
    print("")
    print("mesh dimensions:")
    xmin,ymin,zmin = mesh.bounds[:][0]
    xmax,ymax,zmax = mesh.bounds[:][1]
    print("  bounding box :  x min/max = {} / {}".format(xmin,xmax))
    print("                  y min/max = {} / {}".format(ymin,ymax))
    print("                  z min/max = {} / {}".format(zmin,zmax))
    print("")

    # Save the mesh as a PLY file
    filename = './output_buildings_utm.ply'
    mesh.export(file_obj=filename);

    print("written: ",filename)
    print("")

#
#----------------------------------------------------------------------------------------
#

def convert_openstreetmap_buildings(input_file: str="", mesh_area: np.array=None) -> None:

    # user output
    print("")
    print("convert OpenStreetMap buildings")

    if len(input_file) > 0:
        print("  input file  : ",input_file)
    if not mesh_area is None:
        print("  mesh_area   : ",mesh_area)

    print("")

    # current directory
    dir = os.getcwd()
    print("  current directory: ",dir)
    print("")

    # input data
    if len(input_file) > 0:
        # reads in buildings and loads GeoJSON format file into a GeoDataFrame object
        gdf = gpd.read_file(input_file)
    else:
        # download buildings
        gdf = download_buildings(mesh_area)

    # prepare gdf data array and validate entries
    gdf = prepare_osm_building_data(gdf)

    # download topography for elevations
    topo_grid = get_topo_DEM(mesh_area)

    # make 3d buildings
    create_3d_buildings(gdf,topo_grid,mesh_area)


#
#----------------------------------------------------------------------------------------
#

def usage() -> None:
    print("usage: ./convert_openstreetmap_buildings_to_dxf.py [--mesh_area=(lon_min,lat_min,lon_max,lat_max)] [--GeoJSON_file=file] [--utm_zone=ZoneNumber]")
    print("  with")
    print("     --mesh_area             - downloads/limits buildings for area specified by (lon_min,lat_min,lon_max,lat_max)")
    print("     --GeoJSON_file          - use input mesh file (.gjson) instead of download")
    print("     --utm_zone              - use specified UTM zone number (1-60) with (+) for Northern (-) for Southern hemisphere (e.g., -58)")
    sys.exit(1)


if __name__ == '__main__':
    # init
    input_file: str = ""
    mesh_area: np.array = None

    # reads arguments
    #print("\nnumber of arguments: " + str(len(sys.argv)))
    if len(sys.argv) <= 1: usage()

    i = 0
    for arg in sys.argv:
        i += 1
        #print("argument "+str(i)+": " + arg)
        # get arguments
        if "--help" in arg:
            usage()
        elif "--mesh_area=" in arg:
            str_array = arg.split('=')[1]
            mesh_area = np.array([float(val) for val in str_array.strip('()[]').split(',')])
        elif "--GeoJSON_file=" in arg:
            input_file = arg.split('=')[1]
        elif "--utm_zone=" in arg:
            str_val = arg.split('=')[1]
            utm_zone = int(str_val)
            if abs(utm_zone) < 1 or utm_zone > 60:
                print(f"Invalid UTM zone entered: {utm_zone} - Please use zones from +/- [1,60]")
                sys.exit(1)
        elif "--store-scene" in arg:
            store_scene_file = True
        elif i >= 3:
            print("argument not recognized: ",arg)

    # logging
    cmd = " ".join(sys.argv)
    filename = './convert_openstreetmap_buildings_to_mesh.log'
    with open(filename, 'a') as f:
      print("command call --- " + str(datetime.datetime.now()),file=f)
      print(cmd,file=f)
      print("command logged to file: " + filename)

    # main routine
    convert_openstreetmap_buildings(input_file,mesh_area)
