#!/bin/bash
#
# script to download topography and EMC model data
#
################################################

## User parameters

# topography resolution (topo30 == 30 arc-seconds)
topo=topo30

# shifted internal interface (in m)
toposhift=30000.0
toposcale=0.0

# EMC model
model=FWEA23.r0.0-n4.nc

################################################

# working directory
currentdir=`pwd`
echo
echo "directory: $currentdir"
echo

# check needed scripts
if [ ! -e ./run_get_simulation_topography.py ]; then ln -s ../../../utils/scripts/run_get_simulation_topography.py; fi
if [ ! -e ./run_get_simulation_USGS_Vs30.py ]; then ln -s ../../../utils/scripts/run_get_simulation_USGS_Vs30.py; fi
if [ ! -e ./run_convert_IRIS_EMC_netCDF_2_tomo.py ]; then ln -s ../../../utils/scripts/run_convert_IRIS_EMC_netCDF_2_tomo.py; fi

# region setup
if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == ""] || [ "$4" == "" ]; then
  echo
  echo "using default region"
  echo "(optional) to specify your own region, use: > ./setup_model.sh lon_min lon_max lat_min lat_max"
  echo

  ## default region
  # epicenter
  lat=21.70
  lon=95.90
  # margin +/-
  margin_lat=10.0   # includes Bangkok at around lat/lon: 13.7/100.5
  margin_lon=7.0
  shift_lat=-1.0    # to shift center location for better Bangkok distance to absorbing boundaries
  shift_lon=1.0

  # input region
  lon_min=`echo "scale=2;$lon-$margin_lon+$shift_lon" | bc`
  lon_max=`echo "scale=2;$lon+$margin_lon+$shift_lon" | bc`
  lat_min=`echo "scale=2;$lat-$margin_lat+$shift_lat" | bc`
  lat_max=`echo "scale=2;$lat+$margin_lat+$shift_lat" | bc`

  echo "epicenter: lat = $lat"
  echo "           lon = $lon"
  echo
  echo "range: lat min/max = $lat_min / $lat_max"
  echo "       lon min/max = $lon_min / $lon_max"
  echo "       using margin +/- $margin (degrees)"
  echo "       using shift center by lat/lon $shift_lat / $shift_lon"
  echo

  region="${lon_min} ${lat_min} ${lon_max} ${lat_max}"
  regionShort="${lon_min},${lat_min},${lon_max},${lat_max}"

  # extended area for velocity model and Vs30
  # from: ./run_get_simulation_topography.py 89.9 10.70 103.90 30.70 ..
  # region: (89.9, 10.7, 103.9, 30.7)
  # extended region: (88.6, 10.6, 104.61, 30.9)
  regionExtended="88.6 10.6 104.61 30.9"
  regionShortExtended="88.6,10.6,104.61,30.9"

else
  # region given as arguments
  echo
  echo "using input as region"
  echo
  # input region
  lon_min=$1
  lon_max=$2
  lat_min=$3
  lat_max=$4

  echo "range: lat min/max = $lat_min / $lat_max"
  echo "       lon min/max = $lon_min / $lon_max"
  echo

  region="${lon_min} ${lat_min} ${lon_max} ${lat_max}"
  regionShort="${lon_min},${lat_min},${lon_max},${lat_max}"

  # extended area for velocity model and Vs30
  lon_min_ext=`echo "scale=2;$lon_min-0.2" | bc`
  lon_max_ext=`echo "scale=2;$lon_max+0.2" | bc`
  lat_min_ext=`echo "scale=2;$lat_min-0.2" | bc`
  lat_max_ext=`echo "scale=2;$lat_max+0.2" | bc`
  regionExtended="${lon_min_ext} ${lat_min_ext} ${lon_max_ext} ${lat_max_ext}"
  regionShortExtended="${lon_min_ext},${lat_min_ext},${lon_max_ext},${lat_max_ext}"
fi

echo
echo "#################################################"
echo
echo "model setup..."
echo
echo "region: $region"
echo "        short $regionShort"
echo
echo "#################################################"
echo

# topo
./run_get_simulation_topography.py $region --SRTM=$topo --toposhift=$toposhift --toposcale=$toposcale

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# create additional dummy interface at 100km depth for meshing
if [ -e topo_data/ptopo.xyz.2.dat ]; then
  # dummy at 100km depth
  echo "creating dummy interface at 100km depth"
  awk '{printf("%.6f\n",-100000.0);}' topo_data/ptopo.xyz.2.dat > topo_data/ptopo.xyz.2.100km.dat
  cd DATA/meshfem3D_files/
  rm -f ptopo.xyz.2.100km.dat
  ln -s ../../topo_data/ptopo.xyz.2.100km.dat
  cd ../../
fi

# create additional dummy interface at 6km depth for meshing to accommodate topo
# topography in region has: min/max = <-3635.461182/7137.665527>
# for example: by 2 element layers
#   -> at sea level    : top elements (depth) ~ 6 / 2   = 3 km
#   -> at lowest point :                        2.4 / 2 = 1.2 km
#         highest point:                        13 / 2  = 6.5 km
#
if [ -e topo_data/ptopo.xyz.1.dat ]; then
  # dummy at 6km depth
  echo "creating dummy interface at 6km depth"
  awk '{printf("%.6f\n",-6000.0);}' topo_data/ptopo.xyz.1.dat > ./topo_data/ptopo.xyz.2.6km.dat
  cd DATA/meshfem3D_files/
  rm -f ptopo.xyz.2.6km.dat
  ln -s ../../topo_data/ptopo.xyz.2.6km.dat
  cd ../../
fi

echo
echo

##
## download IRIS model
##

# EMC model
if [ ! -e IRIS_EMC/${model} ]; then
  echo
  echo "downloading IRIS model..."
  echo
  echo "model: ${model}"
  echo

  mkdir -p IRIS_EMC
  cd IRIS_EMC/

  wget https://ds.iris.edu/ds/products/script/emcscript/meta_2.py?model=${model}
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  mv -v meta_2.py\?model=${model} ${model}.metadata.txt
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  wget https://ds.iris.edu/files/products/emc/emc-files/${model}

  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  cd ../
fi

# tomography model
echo
echo "#################################################"
echo
echo "setting up tomography model..."
echo
echo "#################################################"
echo

# get the depth from Mesh_Par_file
DEP=`grep ^DEPTH_BLOCK_KM DATA/meshfem3D_files/Mesh_Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

# Replace 'd' with '.'; for example "180.d0" -> "180.0"
DEPTH=$(echo "$DEP" | sed 's/d//')

echo "  model                     : ${model}"
echo "  extended region           : $regionShortExtended"
echo "  depth (from Mesh_Par_file): $DEPTH (km)"
echo

# EMC model
./run_convert_IRIS_EMC_netCDF_2_tomo.py --EMC_file=IRIS_EMC/${model} --mesh_area=$regionShortExtended --maximum_depth="${DEPTH}"
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mkdir -p DATA/tomo_files
mv -v tomography_model.* DATA/tomo_files/

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# USGS Vs30
echo
echo "#################################################"
echo
echo "setting up USGS Vs30"
echo
echo "#################################################"
echo
echo "  extended region           : $regionExtended"
echo

./run_get_simulation_USGS_Vs30.py $regionExtended
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# add symbolic link
if [ -e USGS_VS30/interface_vs30.dat ]; then
  echo "adding Vs30 interface to DATA/"
  cd DATA/
  ln -s ../USGS_VS30/interface_vs30.dat
  cd ../
fi

echo
echo "done"
echo
