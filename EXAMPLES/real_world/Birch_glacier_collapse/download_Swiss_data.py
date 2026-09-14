#!/usr/bin/env python
#
# script to download and process seismic data from the Swiss national seismic network
#
from __future__ import print_function
import os,sys

# uses obspy
# http://docs.obspy.org
# a python framework for seismology
import obspy
print("obspy version: ",obspy.__version__)
# for geodetic distances
from obspy.geodetics import gps2dist_azimuth

## numpy
import numpy as np
print("numpy version: ",np.__version__)

# do not show scipy warnings
import warnings
warnings.filterwarnings('ignore')

# Voronoi weighting
from scipy.spatial import Voronoi, voronoi_plot_2d

# UTM projection
import pyproj
from pyproj import Transformer

# Polygons
from shapely.geometry import Polygon

# plotting
import matplotlib as mpl
print("matplotlib version: ",mpl.__version__)

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize,LinearSegmentedColormap
from matplotlib.patches import Polygon as matplotlib_Polygon

#############################################################################
## USER parameters

## event data
# Birch glacier collapse 2025-05-28 UTC 13:24:26.162933
event_time = obspy.UTCDateTime("2025-05-28T13:24:20")
# duration in sec
event_duration = 2 * 60  # downloads 2 minutes traces

## filtering
filter_range_f1 = 1.0/50.0  # low-frequency taper range (Hz)
filter_range_f2 = 1.0/40.0
filter_range_f3 = 0.5        # high-frequency taper range (Hz)
filter_range_f4 = 1.0

## interpolation to match SEM simulation
do_interpolation = True
# to match DT from SEM simulation
DT_sem = 0.006

## amplitude correction
# amplitude correction factors
do_amplitude_correction = True
# geometric
use_correction_geom = False
correction_gamma = 0.5              # geometrical spreading factor (0.5 for surface waves, 1.0 for body waves)
# attenuation
use_correction_Q = False
correction_att_Q = 100              # attenuation factor average Q value
correction_att_v = 3.5              # attenuation factor average velocity in km/s
correction_att_f = 1.0              # attenuation factor dominant frequency (Hz)
# Voronoi cell weighting
use_correction_Voronoi = False      # station weighting according to Voronoi cell areas (suggested by Larmat et al. 2008)
# influence radius weighting
use_correction_influence_radius = True
correction_factor_influence_radius = 10.0  # 10 km radius around station
# amplitude amplification
correction_factor_amplify = 1000.0

#############################################################################

# defaults for nicer plotting
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 20, 5
plt.rcParams['lines.linewidth'] = 0.5

## file output
def save_traces_as_ascii(st,ending="ascii",show_traces=False):
    # output directory
    dir = "./network_data/"
    if not os.path.exists(dir):
        os.makedirs(dir)

    # ASCII format for gnuplot
    for i, tr in enumerate(st):
        npts = tr.stats.npts
        dt = tr.stats.delta
        net = tr.stats.network
        sta = tr.stats.station
        cha = tr.stats.channel

        # save as ascii file
        filename = dir + "{}.{}.{}.{}".format(net,sta,cha,ending)
        with open(filename,'w') as f:
            f.write("# NETWORK %s\n" % (tr.stats.network))
            f.write("# STATION %s\n" % (tr.stats.station))
            f.write("# CHANNEL %s\n" % (tr.stats.channel))
            f.write("# START_TIME %s\n" % (str(tr.stats.starttime)))
            f.write("# SAMP_FREQ %f\n" % (tr.stats.sampling_rate))
            f.write("# NDAT %d\n" % (tr.stats.npts))
            f.write("# DELTA %f\n" % (tr.stats.delta))
            for i in range(npts):
                tval = i * dt
                dat = tr[i]
                f.write("{} {}\n".format(tval,dat))

        print("written to: ",filename)

        # plotting
        if show_traces:
            duration = npts * dt
            t_s = np.arange(0, duration, dt)
            plt.clf()
            plt.plot(t_s, tr.data, color='blue', linewidth=1.5, label="{}.{}.{}".format(net,sta,cha))
            plt.show()

    print("")

#----------------------------------------------------------------------------------------

def save_STATIONS_file(inventory,stations_exclude,st):
    """
    write out STATIONS file
    """
    print("creating STATIONS file in network_data/ ...")

    # output directory
    dir = "./network_data/"
    if not os.path.exists(dir):
        os.makedirs(dir)

    # stations info
    stations_info = []
    for network in inventory.networks:
        for station in network.stations:
            # check if station gets excluded
            if station.code in stations_exclude:
                print(f"   station: *** {station.code} will be excluded ***")
            else:
                lat = station.latitude
                lon = station.longitude
                depth = 0.0
                # depth of first channel (HHE)
                if len(station) > 0:
                    channel = station[0]
                    depth = channel.depth
                    # force depth to be positive (negative depth? station in air or on podium?)
                    if depth < 0.0: depth = 0.0
                print(f"   station: {station.code} - lat/lon = {lat} / {lon} - depth = {depth}")
                # add to stations info
                # format: #station_ID #network_ID #lat #lon #elevation(ignored) #burial_depth
                line = f"{station.code}  {network.code}  {lat}  {lon}  0.0  {depth}"
                stations_info.append(line)
    print("")
    print("  found number of stations in inventory: ",len(stations_info))

    # remove stations not found in stream traces
    used_stations = set()
    used_networks = set()
    for trace in st:
        net = trace.stats.network
        sta = trace.stats.station
        used_stations.add(sta)
        used_networks.add(net)
    print("  used networks: ",used_networks)
    print("  used stations: ",used_stations)
    print("")

    stations_filtered = []
    stations_not_found = []
    for line in stations_info:
        items = line.split()
        #print("  line items: ",items)
        sta = items[0]
        net = items[1]
        if net in used_networks and sta in used_stations:
            stations_filtered.append(line)
        else:
            stations_not_found.append([sta,net])
    print("  filtered number of stations in inventory: ",len(stations_filtered))
    print("  stations not found in stream: ",stations_not_found)
    print("")

    # checks if stations used in stream are in filtered list
    stations_not_found = []
    for station in used_stations:
        found = False
        for line in stations_filtered:
            if station in line:
                found = True
        if not found:
            stations_not_found.append(station)
    if len(stations_not_found) > 0:
        print("  stations used but not found in filtered stations: ",stations_not_found)
        print("")

    # checks if anything to do
    if len(stations_filtered) == 0: return

    # file output
    filename = dir + "STATIONS"
    with open(filename,'w') as f:
        f.write("# stations file\n")
        f.write("# created by script download_Swiss_data.py\n") # providence
        f.write("# format:\n")
        f.write("#station_ID #network_ID #lat #lon #elevation(ignored) #burial_depth\n")
        for line in stations_filtered:
            f.write(line + "\n")
    print("  written to: ",filename)
    print("")

#----------------------------------------------------------------------------------------

def process_data(st,catalog,inventory,save=False):
    """
    data processing and output
    """
    global filter_range_f1,filter_range_f2,filter_range_f3,filter_range_f4
    global do_interpolation,DT_sem
    global do_amplitude_correction

    print("data processing...")

    # checks if stream has response
    attach_response = False
    trace = st[0]
    if trace:
        if 'response' in trace.stats:
            print(f"  response information available in stream...")
            print(f"    for {trace.id}:")
            print(trace.stats.response)
        else:
            attach_response = True
    if attach_response:
        print("  attaching response function from inventory to stream...")
        st.attach_response(inventory)
    print("")

    # detrending & tapering
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=0.05)

    # instrument response removal
    # (unstable procedure)
    #
    # Nanometrics Trillium
    #
    # same as STS-2?
    #   STS-2 instrument has a broadband range between 120s - 50 Hz
    #   Streckeisen STS-2.5: range between 120s - 50 Hz (https://streckeisen.swiss/en/products/sts-2.5/)
    # frequency range
    #f1 = 1.0/150.0  # low-frequency taper range (Hz)
    #f2 = 1.0/120.0
    #f3 = 50.0        # high-frequency taper range (Hz)
    #f4 = 60.0

    # landslides are mostly in frequency range 1-3 Hz (? Taiwan study)
    # select main frequency range: 0.1 - 5 Hz
    f1 = filter_range_f1  # low-frequency taper range (Hz)
    f2 = filter_range_f2
    f3 = filter_range_f3  # high-frequency taper range (Hz)
    f4 = filter_range_f4
    print("  instrument response removal:")
    print(f"    frequency range: {f1} / {f2} to {f3} / {f4} (Hz)")
    print("")

    ## acceleration
    st_accel = st.copy()
    st_accel.remove_response(output="ACC",pre_filt=(f1,f2,f3,f4))

    print("data converted to acceleration:")
    print("  vertical acceleration min/max = ",st_accel[0].data.min(),"/",st_accel[0].data.max(),"(m)")
    print("")

    # plotting
    if not save:
        st_accel.plot(size=(1800, 1600))

    # show spectrogram of first trace
    if not save:
        print("spectrogram:")
        trace = st_accel[0]
        print(trace)
        trace.spectrogram(log=True, title=f"{trace.stats.network}.{trace.stats.station} - {trace.stats.starttime}")
        print("")

    # interpolation
    if do_interpolation:
        st_accel = interpolate_traces(st_accel,DT_sem)

    # amplitude correction
    if do_amplitude_correction:
        st_accel = amplitude_correction(st_accel,catalog,inventory)

    # apply taper again to be remove spurious non-zero offsets
    st_accel.taper(max_percentage=0.05)

    # save files
    if save:
        save_traces_as_ascii(st_accel,ending="sema.ascii")

    ## velocity
    st_veloc = st.copy()
    st_veloc.remove_response(output="VEL",pre_filt=(f1,f2,f3,f4))

    print("data converted to velocity:")
    print("  vertical velocity min/max = ",st_veloc[0].data.min(),"/",st_veloc[0].data.max(),"(m)")
    print("")

    # plotting
    if not save:
        st_veloc.plot(size=(1600, 1600))

    # interpolation
    if do_interpolation:
        st_veloc = interpolate_traces(st_veloc,DT_sem)

    # amplitude correction
    if do_amplitude_correction:
        st_veloc = amplitude_correction(st_veloc,catalog,inventory)

    # apply taper again to be remove spurious non-zero offsets
    st_veloc.taper(max_percentage=0.05)

    # save files
    if save:
        save_traces_as_ascii(st_veloc,ending="semv.ascii")

    ## displacement
    st_displ = st.copy()
    st_displ.remove_response(output="DISP",pre_filt=(f1,f2,f3,f4))

    print("data converted to displacement:")
    print("  vertical displacement min/max = ",st_displ[0].data.min(),"/",st_displ[0].data.max(),"(m)")
    print("")

    # plotting
    if not save:
        st_displ.plot(size=(1600, 1600))

    # interpolation
    if do_interpolation:
        st_displ = interpolate_traces(st_displ,DT_sem)

    # amplitude correction
    if do_amplitude_correction:
        st_displ = amplitude_correction(st_displ,catalog,inventory)

    # apply taper again to be remove spurious non-zero offsets
    st_displ.taper(max_percentage=0.05)

    # save files
    if save:
        save_traces_as_ascii(st_displ,ending="semd.ascii")

#----------------------------------------------------------------------------------------

def interpolate_traces(st,DT_sem):
    """
    interpolates all traces in given stream st object using new time step DT_sem
    """
    print("interpolation:")
    print("  SEM DT                 : ",DT_sem)
    print("  total number of traces : ",len(st))
    print("")

    # first trace will determine a reference new time array
    # this avoids issues when subsequent stations have different sampling rates,
    # and the resulting new traces not the same overall length
    new_time_sampling = None

    # interpolate traces
    for i,tr in enumerate(st):
        print(f"  processing trace {tr.stats.network}.{tr.stats.station} channel {tr.stats.channel} - {i+1} out of {len(st)}...")
        # original trace data
        #print(tr.stats)
        npts = tr.stats.npts      # number of samples
        dt = tr.stats.delta       # time step
        duration = npts * dt
        # reference time
        ref_time = np.arange(0.0, duration - dt/2, dt)         # - dt/2 to avoid rounding errors which lead to an additional entry
        # new time sampling
        if new_time_sampling is None:
            new_time_sampling = np.arange(0.0, duration - DT_sem/2, DT_sem)
            # avoid having some round-off problem with an additional entry
            if (len(new_time_sampling)-1) * DT_sem > (npts-1) * dt:
                new_time_sampling = new_time_sampling[0:len(new_time)-2]    # takes 1 entry off with 0:length-2
            print("")
            print(f"  original number of samples : {npts}")
            print(f"           effective duration: {(npts-1)*dt} (s)")
            print(f"  new      number of samples : {len(new_time_sampling)}")
            print(f"           effective duration: {(len(new_time_sampling)-1)*DT_sem}")
            print("")

        # interpolate trace
        trace = tr.data.copy()
        trace_interp = np.interp(new_time_sampling, ref_time, trace)

        # info
        if i == 0:
            print(f"  original time step: {dt}")
            print(f"  original number of time steps: {npts}")
            print(f"  original trace min/max : {trace.min()} / {trace.max()}")
            print("")
            print(f"  new      time step: {DT_sem}")
            print(f"  new      number of time steps: {len(new_time_sampling)}")
            print(f"  new      trace min/max : {trace_interp.min()} / {trace_interp.max()}")
            print("")

        # replace trace
        tr.stats.delta = DT_sem
        tr.stats.npts = len(trace_interp)
        tr.data = trace_interp.copy()
    print("")

    # return modified stream object
    return st

#----------------------------------------------------------------------------------------

def amplitude_correction(st,catalog,inventory):
    """
    uses amplitude scaling for all trace in given stream to correct for geometrical spreading and attenuation effects
    """
    global correction_factor_amplify
    global use_correction_geom,correction_gamma
    global use_correction_Q,correction_att_Q,correction_att_v,correction_att_f
    global use_correction_Voronoi
    global use_correction_influence_radius,correction_factor_influence_radius

    # amplitude correction
    gamma = correction_gamma   # geometrical spreading factor (0.5 for surface waves, 1.0 for body waves)
    Q = correction_att_Q       # average Q value
    v = correction_att_v       # average velocity in km/s
    f = correction_att_f       # dominant frequency (Hz)

    print("amplitude correction:")
    if use_correction_geom:
        print(f"  using geometric weighting")
        print(f"    geometrical spreading factor : {gamma}")
    if use_correction_Q:
        print(f"  using attenuation weighting")
        print(f"    average Q value              : {Q}")
        print(f"    average velocity             : {v}")
        print(f"    dominant frequency           : {f}")
    if use_correction_Voronoi:
        print(f"  using Voronoi cell weighting")
    if use_correction_influence_radius:
        print(f"  using influence radius weighting")
        print(f"    radius                       : {correction_factor_influence_radius}")
    print("")
    print(f"  amplification constant       : {correction_factor_amplify}")
    print("")

    # check only use one or the other of the station weighting
    if use_correction_Voronoi and use_correction_influence_radius:
        print("Err: cannot weight by both station Voronoi cell and influence radius, only use one...")
        sys.exit(1)

    # determine source event location
    print("  catalog: number of events = ",len(catalog))
    if len(catalog) == 0:
        print("  no event location, no amplitude correction - returning...")
        return st
    event_locations = {}
    for i,event in enumerate(catalog):
        print(f"  event {i}: {event.resource_id}")
        event = catalog[0]
        origin = event.preferred_origin()
        if origin:
            event_lat = origin.latitude
            event_lon = origin.longitude
        else:
            print("  event has no preferred origin, skipping...")
            continue
        name = f"event{i}"
        location = { name: [ event_lat, event_lon ]}
        print("  location: ",location)
        print("")
        # add to dict
        event_locations.update(location)
    print("  event locations        : ",event_locations)

    # determine station locations
    # format like { "CH.LAUCH": [ lat, lon], "CH.FIESA": [ lat, lon], ..}
    station_locations = {}
    for network in inventory.networks:
        for station in network.stations:
            station_lat = station.latitude
            station_lon = station.longitude
            name = f"{network.code}.{station.code}"
            location = { name: [station_lat, station_lon] }
            station_locations.update(location)

    print("  station locations      : ",station_locations)
    print("  total number of traces : ",len(st))
    print("")

    # station location in UTM coordinates
    points_utm = get_stations_locations_UTM(station_locations)

    # Voronoi cell weighting
    if use_correction_Voronoi:
        # suggested by Larmat et al. 2008
        stations_voronoi_weights = determine_station_voronoi_weights(points_utm)

        # add weights to station infos
        for i,entry in enumerate(station_locations):
            station_locations[entry].append(stations_voronoi_weights[i])

    # influence radius
    if use_correction_influence_radius:
        # determine number of stations within a given radius
        stations_influence_weights = determine_station_influence_weights(points_utm)

        # add weights to station infos
        for i,entry in enumerate(station_locations):
            station_locations[entry].append(stations_influence_weights[i])

    # stats
    min_glob_org = 0.0
    max_glob_org = 0.0
    min_glob_corrected = 0.0
    max_glob_corrected = 0.0

    # interpolate traces
    for i,tr in enumerate(st):
        # get station name
        name = f"{tr.stats.network}.{tr.stats.station}"
        print(f"  processing trace {name} channel {tr.stats.channel} - {i+1} out of {len(st)}...")

        # get station location
        station_coords = station_locations.get(name)

        # distance to reference event (first in dict)
        event_coords = event_locations['event0']

        # Calculate the distance using gps2dist_azimuth
        if station_coords is not None and event_coords is not None:
            # gps2dist_azimuth(srclat,srclon,stalat,stalon)
            distance_in_m, azimuth, backazimuth = gps2dist_azimuth(event_coords[0], event_coords[1], station_coords[0], station_coords[1])
            # epicentral distance in km
            r = distance_in_m / 1000.0

            ## correct amplitude
            # geometric spreading
            # (factor 2.0 to offset forward and time-reversed backward propagation amplitude decay)
            if use_correction_geom:
                correction_factor_geom = 2.0 * (r**gamma)
            else:
                correction_factor_geom = 1.0

            # attenuation
            # (factor 2.0 to offset forward and time-reversed backward propagation amplitude decay)
            if use_correction_Q:
                correction_factor_att = 2.0 * np.exp(np.pi * f * r / (Q * v))
            else:
                correction_factor_att = 1.0

            # Voronoi station area weighting
            if use_correction_Voronoi:
                correction_factor_station_weight = station_coords[2]  # station weighting added as [lat,lon,weight] to locations
            else:
                correction_factor_station_weight = 1.0

            # station influence radius weighting
            if use_correction_influence_radius:
                correction_factor_station_weight = station_coords[2]  # station weighting added as [lat,lon,weight] to locations
            else:
                correction_factor_station_weight = 1.0

            # total correction factor
            correction_factor = correction_factor_geom * correction_factor_att * correction_factor_station_weight * correction_factor_amplify

            print(f"    distance from station ({station_coords[0]:.2f}, {station_coords[1]:.2f}) to source ({event_coords[0]:.2f}, {event_coords[1]:.2f}) = {r:.2f} km")
            print(f"    correction factor: geom {correction_factor_geom:.2f} / att {correction_factor_att:.2f} / amplify {correction_factor_amplify:.1f} / station weight {correction_factor_station_weight:.2f} - total x {correction_factor:.2f}")

            # original for vis
            trace_org = tr.data.copy()

            # correct waveform
            tr.data = tr.data * correction_factor

            # stats
            min_glob_org = min(min_glob_org,trace_org.min())
            max_glob_org = max(max_glob_org,trace_org.max())
            min_glob_corrected = min(min_glob_org,tr.data.min())
            max_glob_corrected = max(max_glob_org,tr.data.max())

            # visualize
            if 1 == 0:
                # plot waveforms
                # corrected
                plt.plot(tr.times(),tr.data,label=f"{tr.stats.station} - corrected",color='orange')
                # original
                plt.plot(tr.times(),trace_org,label=f"{tr.stats.station}",color='black',linewidth=0.5)
                plt.legend()
                plt.title(f"{tr.stats.station} - corrected waveform  x {correction_factor:.2f} (geom {correction_factor_geom:.2f} / att {correction_factor_att:.2f} / amplify {correction_factor_amplify:.1f})")
                plt.xlabel("Time (s)")
                plt.ylabel("Amplitude")
                plt.show()

    # stats
    print("")
    print(f"  original data : min/max = {min_glob_org:.4e}/{max_glob_org:.4e}")
    print(f"  corrected data: min/max = {min_glob_corrected:.4e}/{max_glob_corrected:.4e}")
    print("")

    return st

#----------------------------------------------------------------------------------------

def get_stations_locations_UTM(stations_locations):
    """
    calculates UTM coordinates of station locations
    """
    # extract stations point locations
    points_geo = []
    for entry in stations_locations:
        lat = stations_locations[entry][0]
        lon = stations_locations[entry][1]
        points_geo.append([lon,lat])

    print(f"  number of station locations: {len(points_geo)}")

    # determine centroid (center) lat/lon for UTM zone detection
    lons = np.array([point[0] for point in points_geo])
    lats = np.array([point[1] for point in points_geo])
    center_lon = np.mean(lons)
    center_lat = np.mean(lats)

    # bounding box min/max
    box_lat_min = lats.min()
    box_lat_max = lats.max()
    box_lon_min = lons.min()
    box_lon_max = lons.max()

    print(f"  centroid location          : lat/lon = {center_lat} / {center_lon}")
    print(f"  bounding box               : lat min/max = {box_lat_min} / {box_lat_max}")
    print(f"                               lon min/max = {box_lon_min} / {box_lon_max}")

    # gets UTM zone codes (uses first, single point in points[] list)
    utm_crs_list = pyproj.database.query_utm_crs_info(datum_name="WGS 84",
                                                      area_of_interest=pyproj.aoi.AreaOfInterest(west_lon_degree=center_lon,
                                                                                                 south_lat_degree=center_lat,
                                                                                                 east_lon_degree=center_lon,
                                                                                                 north_lat_degree=center_lat))
    utm_code = utm_crs_list[0].code
    utm_epsg = "EPSG:{}".format(utm_code)
    print(f"  detected UTM code          : {utm_code}  epsg: {utm_epsg}")
    print("")
    # transformer
    # transformation from WGS84 to UTM zone
    # WGS84: Transformer.from_crs("EPSG:4326", utm_epsg)
    transformer_to_utm = Transformer.from_crs("EPSG:4326", utm_epsg, always_xy=True) # input: lon/lat -> utm_x/utm_y
    points_utm = [transformer_to_utm.transform(lon, lat) for lon, lat in points_geo]

    # convert to array
    points_utm = np.array(points_utm)

    return points_utm

#----------------------------------------------------------------------------------------

def determine_station_influence_weights(points_utm):
    """
    determines weighting based on number of station close-by
    """
    global correction_factor_influence_radius
    print("  determine station influence radius weighting...")
    print(f"    influence radius         : {correction_factor_influence_radius} (km)")

    # stats
    min_distances = []
    number_of_neighbors = []

    # determines weights
    stations_influence_weights = []
    for i,point_ref in enumerate(points_utm):
        # initialize
        num_neighbors = 0
        dist_min = np.inf
        # loops over all other stations
        for j,point in enumerate(points_utm):
            # skip own reference point
            if point[0] == point_ref[0] and point[1] == point_ref[1]:
                continue

            # distance between points (in m)
            dist = np.sqrt( (point_ref[0] - point[0])**2 + (point_ref[1] - point[1])**2)
            dist_in_km = dist / 1000.0

            #debug
            #if i == 0: print(j,"point ref: ",point_ref,"point:",point,"distance",dist_in_km)

            #stats
            dist_min = min(dist_min,dist_in_km)

            # add counter for close stations
            if dist_in_km < correction_factor_influence_radius:
                num_neighbors += 1

        # weighting
        if num_neighbors == 0:
            # no influence from close-by stations
            weight = 1.0
        else:
            weight = 1.0 / num_neighbors

        # store weighting
        stations_influence_weights.append(weight)

        # stats
        min_distances.append(dist_min)
        number_of_neighbors.append(num_neighbors)

    # Print result
    #print("  station weights:")
    #for i, entry in enumerate(points_utm):
    #    weight = stations_influence_weights[i]
    #    print(f"    station {entry}: weight = {weight:.4f}")
    #print("")

    # stats
    min_distances = np.array(min_distances)
    print(f"    station distances        : closest station distance min/max = {min_distances.min():.1f} / {min_distances.max():.1f} - avg {min_distances.mean():.1f}")
    number_of_neighbors = np.array(number_of_neighbors)
    print(f"    station influence        : number of close-by stations min/max = {number_of_neighbors.min()} / {number_of_neighbors.max()} - avg {number_of_neighbors.mean():.1f}")
    weight_min = np.array(stations_influence_weights).min()
    weight_max = np.array(stations_influence_weights).max()
    print(f"    station influence weights: min/max = {weight_min} / {weight_max}")
    print("")

    return stations_influence_weights

#----------------------------------------------------------------------------------------

def determine_station_voronoi_weights(points_utm,show=False):
    """
    determine Voronoi cells for a set of stations given by UTM point locations
    """
    print("  determine station Voronoi cell area weighting...")

    # Voronoi cells
    vor = Voronoi(points_utm)

    # visualize
    if show:
        fig = voronoi_plot_2d(vor)
        #plt.plot(*zip(*points_utm), 'ro')  # Plot original UTM points
        plt.title("Voronoi Diagram - UTM Coordinates")
        plt.show()

    if show:
        # plot finite cells
        fig2, ax2 = plt.subplots()
        for region_idx in vor.point_region:
            region = vor.regions[region_idx]
            if not -1 in region:  # -1 indicates an infinite region
                polygon = vor.vertices[region]
                ax2.fill(*zip(*polygon), alpha=0.4)
        ax2.plot(vor.points[:, 0], vor.points[:, 1], 'ko')  # Plot the input points
        ax2.set_xlim(vor.min_bound[0] - 0.1, vor.max_bound[0] + 0.1)
        ax2.set_ylim(vor.min_bound[1] - 0.1, vor.max_bound[1] + 0.1)
        ax2.set_title("Finite Voronoi Regions")
        plt.show()

    # set region points outside of bounding lat/lon to infinite (-1)
    # this avoids using polygon areas that are associate with station points at the rim of our region
    for region_idx in vor.point_region:
        region = vor.regions[region_idx]
        for index,vertex in enumerate(region):
            # only finite vertex points
            if vertex >= 0:
                # get lon/lat from voronoi point location
                coord = vor.vertices[vertex]
                lon,lat = transformer_to_utm.transform(coord[0], coord[1],direction='INVERSE')
                if lon < box_lon_min or lon > box_lon_max or lat < box_lat_min or lat > box_lat_max:
                    # outside our target region, set to infinite
                    val = -1
                else:
                    # inside our target
                    val = vertex
                #print(f"  check vertex: {vertex} lon/lat {lon:.2f}/{lat:.2f} - value {val}")
                # set new vertex value
                vor.regions[region_idx][index] = val

    # Collect region areas
    region_areas = []
    max_area = 0.0
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not region or -1 in region:  # Skip infinite or invalid regions
            region_areas.append(None)
            continue

        # Build polygon from vertices
        polygon = Polygon([vor.vertices[r] for r in region])

        # Compute area (units are in square meters because we're in UTM)
        area = polygon.area
        region_areas.append(area)

        # maximum area
        max_area = max(max_area,area)
    print(f"    maximum Voronoi region area = {max_area}")

    # Print result
    #print("    Voronoi region areas:")
    #for i, area in enumerate(region_areas):
    #    if area:
    #        print(f"    region for point {i}: Area = {area:.2f}")
    #    else:
    #        print(f"    region {i} is infinite")

    # determine station weighting (normalized area)
    station_weights = []
    for area in region_areas:
        if area is None: area = max_area
        # weight by normalized area
        weight = area / max_area
        # bounds
        if weight < 0.0: weight = 0.0
        if weight > 1.0: weight = 1.0
        # add as weight
        station_weights.append(weight)

    # Print result
    #print("  station weights:")
    #for i, entry in enumerate(stations_locations):
    #    weight = station_weights[i]
    #    print(f"    station {entry}: weight = {weight:.4f}")
    #print("")

    # stats
    weight_min = np.array(station_weights).min()
    weight_max = np.array(station_weights).max()
    print(f"  station weights: min/max = {weight_min} / {weight_max}")
    print("")

    if show:
        # for plotting
        # modify polygons to avoid infinite region polygons
        # see: https://github.com/sinhrks/cesiumpy/blob/4ea2aa7d120d0a07a639aa8fc9bdc7434c63030c/cesiumpy/extension/spatial.py#L27
        new_regions = vor.regions.copy()
        new_vertices = vor.vertices.tolist()

        # determine radius
        center = np.mean(vor.points,axis=0)
        span = np.ptp(vor.points,axis=0)
        radius = span.max() * 2.0

        # Construct a map containing all ridges for a given point
        all_ridges = {}
        for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
            all_ridges.setdefault(p1, []).append((p2, v1, v2))
            all_ridges.setdefault(p2, []).append((p1, v1, v2))

        # Reconstruct infinite regions
        for point_idx, region_idx in enumerate(vor.point_region):
            vertices = vor.regions[region_idx]

            if all(v >= 0 for v in vertices):
                # doesn't contain infinite regions, can use vertices as it is
                # finite region
                #new_regions.append(vertices)
                new_regions[region_idx] = vertices
                print(f"  {point_idx}: new region {region_idx}: region org {vor.regions[region_idx]} new {new_regions[region_idx]}")
                continue

            # reconstruct a non-finite region
            ridges = all_ridges[point_idx]
            new_region = [v for v in vertices if v >= 0]

            for p2, v1, v2 in ridges:
                if v2 < 0:
                    # v2 is always finite, v1 is either finite or infinite
                    v1, v2 = v2, v1
                if v1 >= 0:
                    # finite ridge: already in the region
                    continue

                # v1 is infinite, v2 is finite
                t = vor.points[p2] - vor.points[point_idx]  # tangent
                t /= np.linalg.norm(t)
                n = np.array([-t[1], t[0]])          # normal

                midpoint = vor.points[[point_idx, p2]].mean(axis=0)
                direction = np.sign(np.dot(midpoint - center, n)) * n
                far_point = vor.vertices[v2] + direction * radius

                new_vertices.append(far_point.tolist())
                new_region.append(len(new_vertices) - 1)

            # sort region counterclockwise
            verts = np.asarray([new_vertices[v] for v in new_region])
            c = verts.mean(axis=0)
            angles = np.arctan2(verts[:,1] - c[1], verts[:,0] - c[0])
            new_region = np.array(new_region)[np.argsort(angles)]

            # finish
            #new_regions.append(new_region.tolist())
            new_regions[region_idx] = new_region.tolist()
            print(f"  {point_idx}: new region {region_idx}: region org {vor.regions[region_idx]} new {new_regions[region_idx]}")


        #print(f"  new regions: {new_regions}")
        #print(f"  new vertices: {new_vertices}")

        #print("org:")
        #for i,region in enumerate(vor.regions):
        #    print(f"  region {i}: {region}")
        #print(vor.vertices)
        #for point_idx, region_idx in enumerate(vor.point_region):
        #    print(f"  point region: {region_idx} {vor.regions[region_idx]}")

        # update voronoi data
        vor.regions = new_regions
        vor.vertices = np.array(new_vertices)

        #print("new:")
        #for i,region in enumerate(vor.regions):
        #    print(f"  region {i}: {region}")
        #print(vor.vertices)
        #for point_idx, region_idx in enumerate(vor.point_region):
        #    print(f"  point region: {region_idx} {vor.regions[region_idx]}")

        #print(vor.point_region)

        #new_regions = vor.point_region
        #polygons = [new_vertices[r].tolist() for r in new_regions]
        #polygons = [self._geometry.Polygon(p) for p in polygons]
        #polygons = [p.intersection(bbox) for p in polygons]

        #if 1 == 1:
        #    fig = voronoi_plot_2d(vor)
        #    #plt.plot(*zip(*points_utm), 'ro')  # Plot original UTM points
        #    plt.title("Voronoi Diagram - new finite regions")
        #    plt.show()

        #print("    Voronoi point region:",vor.point_region)

        # visualize
        fig = voronoi_plot_2d(vor)
        ax = fig.axes[0]
        #fig = voronoi_plot_2d(vor)
        #fig, ax = plt.subplots()
        cmap = cm.viridis
        norm = Normalize(vmin=0.0,vmax=1.0)

        for i, region_idx in enumerate(vor.point_region):
            weight = station_weights[i]
            region = vor.regions[region_idx]
            print(f"  station {i}: region {region_idx} : region {region}: weight {weight}")

            # polygon from vertices
            vertices = [vor.vertices[r] for r in region]
            polygon = matplotlib_Polygon(vertices, closed=True, facecolor=cmap(norm(weight)))
            ax.add_patch(polygon)

            # Calculate the centroid of the polygon to place the weight label (optional)
            centroid_x = np.mean([v[0] for v in vertices])
            centroid_y = np.mean([v[1] for v in vertices])
            ax.text(centroid_x, centroid_y, f"{weight:.2f}", ha='center', va='center', fontsize=8)

        #plt.plot(*zip(*points_utm), 'ro')  # Plot original UTM points
        plt.title("Voronoi Diagram - weighting")

        # Set plot limits to encompass all polygons
        all_x = [coord[0] for coord in vor.vertices]
        all_y = [coord[1] for coord in vor.vertices]
        if all_x and all_y:
          ax.set_xlim(min(all_x) - 1, max(all_x) + 1)
          ax.set_ylim(min(all_y) - 1, max(all_y) + 1)
        ax.set_aspect('equal', adjustable='box')  # Ensure polygons are drawn correctly

        # Add a colorbar to show the mapping of weights to colors
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        #cbar = fig.colorbar(sm)
        cbar = fig.colorbar(sm, ax=ax, label='Normalized Weight')
        plt.show()

    return station_weights

#----------------------------------------------------------------------------------------

def add_missing_zero_traces(st):
    """
    adds missing 'HHE' and/or 'HHN' zero-value traces to stream object
    if for example only 'HHZ' is present for a given station like CH.ROTHE
    """
    from obspy import Stream, Trace
    from obspy.core.trace import Stats

    complete_stream = Stream()

    # Get unique station names from the stream
    stations = sorted(list(set(tr.stats.station for tr in st)))

    for station_name in stations:
        # Extract traces for the current station
        station_traces = st.select(station=station_name)

        # Check for existing channels
        channels = sorted(list(set(tr.stats.channel for tr in station_traces)))

        # Define the target channels
        target_channels = ['HHN', 'HHE', 'HHZ']

        # If only 'HHZ' is present, add missing 'HHN' and/or 'HHE' as zero traces
        if len(channels) < 3:
            print(f"  station {station_name}: misses channels - got {channels}")

            # Get common stats from the existing HHZ trace to create new traces
            cha = channels[0]
            stored_trace = station_traces.select(channel=cha)[0]

            # Duration and sampling rate for zero traces
            num_samples = stored_trace.stats.npts

            # Create zero traces for missing component (e.g., HHE and/or HHN)
            for channel in target_channels:
                if channel not in channels:
                    print(f"  station {station_name}: adding zero trace for {channel}")
                    zero_data = np.zeros(num_samples, dtype=stored_trace.data.dtype)
                    zero_trace = Trace(data=zero_data)

                    # Copy relevant stats and update channel
                    zero_trace.stats = stored_trace.stats.copy()
                    #zero_trace.stats.network       = stored_trace.stats.network
                    #zero_trace.stats.station       = stored_trace.stats.station
                    #zero_trace.stats.location      = stored_trace.stats.location
                    #zero_trace.stats.starttime     = stored_trace.stats.starttime
                    #zero_trace.stats.endtime     = stored_trace.stats.endtime
                    #zero_trace.stats.sampling_rate = stored_trace.stats.sampling_rate
                    #zero_trace.stats.npts          = stored_trace.stats.npts
                    #zero_trace.stats.delta         = stored_trace.stats.delta
                    zero_trace.stats.channel       = channel

                    complete_stream.append(zero_trace)
            print("")

        # Append existing traces for this station to the completed stream
        for tr in station_traces:
            complete_stream.append(tr)

    # Sort the stream to ensure consistent ordering (optional but good practice)
    complete_stream.sort(['station', 'channel'])

    return complete_stream

#----------------------------------------------------------------------------------------

# main download
def download_and_process_data(input_type=0):
    """
    downloads data from Swiss network
    """
    global event_time,event_duration

    print("")
    print("download data:")
    print("  input type: ",input_type)
    print("")
    print("  event time: ",event_time)
    print("  duration  : ",event_duration,"(s)")
    print("")

    save = (input_type == 1)
    if save:
      print("  saving data in folder network_data/")
      print("")

    # output directory
    dir = "./network_data/"
    if not os.path.exists(dir):
        os.makedirs(dir)

    print("getting event information...")
    from obspy.clients.fdsn import Client

    # Birch glacier collapse 2025-05-28 UTC 13:24:26.162933
    #event_time = obspy.UTCDateTime("2025-05-28T13:24:20")

    # data client
    # data centers: FDSN, IRIS, ORFEUS, IPGP, ETH, GFZ, RESIF, GEONET, USGS, ..
    c = Client("ETH")

    # plus/minus 30 min
    starttime = event_time - 30 * 60
    endtime = event_time + 30 * 60
    Mw_min = 0.0
    #
    catalog = c.get_events(starttime=starttime, endtime=endtime, minmagnitude=Mw_min)
    # or
    #catalog = c.get_events(eventid=3320897)
    print("Event catalog:")
    print(catalog)
    print("")

    # map
    #catalog.plot(projection="local")

    ## stations
    print("getting station information...")

    # gets station informations
    t1 = event_time
    # adds 2 minutes
    t2 = t1 + event_duration

    ## network
    stations_exclude = []

    # Swiss networks: https://networks.seismo.ethz.ch/en/networks/

    # main broadband station network:
    #net = "CH"
    # active temporary networks:
    #net = "1I,4D,5A,8D,8R,9E,9S,C4,G2,S,XJ,XY,YZ"
    # G2: exclude LAVEY
    # S: most exclude
    #
    # best broadband networks
    net = "CH,G2"

    # stations
    # only selected
    #sta = "LAUCH,LKBD2,FIESA,GRIMS"  # station names: LAUCH, LKBD2, FIESA, GRIMS
                                      # strong ground motions(?): SLENK, STSW2, SGAK, SNTZ, SCRM, SZWD2, SFRS, SGWS
    # all available
    sta = "*"
    loc = "*"

    # band codes
    # https://epic.earthscope.org/webfm_send/2134
    # HH* - High broad band: sample rate >= 80 Hz        ; corner period >= 10s
    # BH* - Broad band     : sample rate >= 10 to < 80 Hz; corner period >= 10s
    # SH* - Short period   : sample rate >= 10 to < 80 Hz; corner period < 10s
    # EH* - Extremely short: sample rate >= 80 Hz        ; corner period < 10s
    #
    # HH* traces mostly available next to BH*, but they don't show "overlaps"
    cha = "HH*"   # or vertical component only? "BHZ"

    # after visual inspection, take out these stations
    stations_exclude = [ "MESRY","EMING","WALHA","WOLEN","ROMAN","SGT00","LAVEY" ]

    # see if station has some valuable records
    # gets inventory
    filename = dir + f"inventory.{t1}.{t2}.xml"
    if os.path.isfile(filename):
        # read in existing inventory
        from obspy import read_inventory
        inventory = read_inventory(filename, format='STATIONXML')
    else:
        # download
        #inventory = c.get_stations(network=net, station=sta, starttime=t1, endtime=t2)
        inventory = c.get_stations(network=net, station=sta, location=loc, channel=cha, starttime=t1, endtime=t2, level="response")

        print("inventory:")
        print(inventory)
        print("")
        # show map
        #inventory.plot(projection="local")

        # output as station XML file
        inventory.write(filename, format='STATIONXML')
        print("  written to: ",filename)
        print("")

    # show each single station data
    if input_type == 2:
        print("showing single station waveforms...")
        for network in inventory.networks:
            print("network: ",network.code)
            print("  number of stations: ",len(network.stations))
            for i,station in enumerate(network.stations):
                print(f"  station: {station.code} - {i+1} out of {len(network.stations)}")
                # get waveforms
                try:
                    st = c.get_waveforms(network=network.code, station=station.code, location=loc, channel=cha, starttime=t1, endtime=t2)
                except:
                    print("  no waveform data")
                    print("")
                    continue
                # show
                print(st)
                print("")
                st.plot(size=(1800, 1600))
        print("all available station waveforms shown")
        print("")
        sys.exit(0)

    # gets waveforms
    filename = dir + f"waveforms.{t1}.{t2}.mseed"
    if os.path.isfile(  filename):
        # read in existing records
        from obspy import read
        st = read(filename)
    else:
        # dowloads all seismograms as a stream
        st = c.get_waveforms(network=net, station=sta, location=loc, channel=cha, starttime=t1,endtime=t2, attach_response=True)
        print(st)
        print("")

        # output as miniseed file
        st.write(filename, format='MSEED')
        print("  written to: ",filename)
        print("")

    # exclude stations
    if len(stations_exclude) > 0:
        print("excluding stations: ",stations_exclude)
        print("  original stream : number of traces = ",len(st))
        print("")
        #print(st)
        #print("")
        # new stream
        st_filtered = st.copy()
        # Iterate through the traces in the copied stream and remove those from the selected stations
        num_excluded = 0
        for trace in list(st_filtered):  # Iterate over a copy of the list to allow removal
            if trace.stats.station in stations_exclude:
                st_filtered.remove(trace)
                num_excluded += 1
        print("  number of excluded traces = ",num_excluded)
        print("")
        # replace original with filtered stream
        if num_excluded > 0:
            print("  new stream  : number of traces = ",len(st_filtered))
            st = st_filtered
            print(st_filtered)
            print("")

    # add missing traces to have a complete set of components for each station
    st = add_missing_zero_traces(st.copy())

    # create STATIONS file
    if save:
        save_STATIONS_file(inventory,stations_exclude,st)

    # rotates to N/E/Z
    do_rotate = False
    if do_rotate and cha == "BH*":
        print("rotating to ZNE ...")
        st.rotate('->ZNE',inventory=inventory)
        st[0].stats.channel = 'BHZ'
        st[1].stats.channel = 'BHN'
        st[2].stats.channel = 'BHE'
        print(st)
        print("")
        # show
        if not save:
            st.plot()

    # trace statistics
    tr = st[0]

    #print(tr.stats)
    npts = tr.stats.npts      # number of samples
    dt = tr.stats.delta       # time step
    duration = npts * dt

    freq_Ny = 1.0/(2.0 * dt)  # Nyquist frequency
    freq_min = 1./duration    # minimal possible frequency for length of trace

    print("trace:")
    print("  dt       = ",dt)
    print("  duration = ",duration,"s"," = ",duration/60.,"min"," = ",duration/60./60.,"h")
    print("  Nyquist frequency = ",freq_Ny)
    print("  minimum frequency = ",freq_min)
    print("")

    # plotting traces
    if not save:
        st.plot(size=(1800, 1600))

    ## data processing and output
    process_data(st,catalog,inventory,save)

    print("done")
    print("")

#----------------------------------------------------------------------------------------

def usage():
    print("./download_Swiss_data.py [0==show-only/1==save/2==show-single-stations]")
    print("")
    sys.exit(1)

if __name__ == '__main__':
    # arguments
    if len(sys.argv) < 2: usage()

    # gets input type
    input_type = int(sys.argv[1])

    download_and_process_data(input_type)
