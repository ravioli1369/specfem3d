#!/bin/bash
#
#
################################################

# smaller region - higher topo resolution (3 arc-seconds)
topo=topo15    # topo3, srtm_30m  from ./run_get_simulation_topography.py 7.709098 46.393506 8.109098 46.703506 srtm_30m
# larger region (30 arc-seconds)
#topo=topo30

# shifted topo interface
toposhift=8000.0     # shift original topography value down in m
toposcale=0.0

## region - event epicenter
lat=46.404188
lon=7.836242

# margin +/-
margin_lat=1.1
margin_lon=2.3
shift_lat=0.4    # to shift center location
shift_lon=0.4

# IRIS EMC model
use_emc=0       # 0 == not using / 1 == download IRIS model and convert to SPECFEM format
emc_model="LSP-Eucrust1.0.nc"

# USGS VS30
use_vs30=1      # 0 == not using / 1 == download USGS VS30 model and convert to SPECFEM format

################################################

# working directory
currentdir=`pwd`
echo
echo "directory: $currentdir"
echo

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

region="$lon_min $lat_min $lon_max $lat_max"
regionShort="${lon_min},${lat_min},${lon_max},${lat_max}"

# extended area for velocity model and Vs30
# extended region: 6.651392, 45.581851, 9.141411, 47.626443
regionExtended=""      # "6.651392 45.581851 9.141411 47.626443"
regionShortExtended="" # "6.651392,45.581851,9.141411,47.626443"

# set default region if no extended given
if [ "$regionExtended" == "" ]; then
  regionExtended="$region"
fi
if [ "$regionShortExtended" == "" ]; then
  regionShortExtended="$regionShort"
fi

echo
echo "#################################################"
echo
echo "model setup..."
echo
echo "region: $region"
echo "        short $regionShort"
echo

# get the depth from Mesh_Par_file
DEP=`grep ^DEPTH_BLOCK_KM DATA/meshfem3D_files/Mesh_Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

# Replace 'd' with '.'; for example "180.d0" -> "180.0"
DEPTH=$(echo "$DEP" | sed 's/d//')

echo "depth (from Mesh_Par_file): $DEPTH (km)"
echo
echo "#################################################"
echo

if [ ! -e ./run_get_simulation_topography.py ]; then
  ln -s ../../../utils/scripts/run_get_simulation_topography.py
fi

# topo
./run_get_simulation_topography.py $region --SRTM=$topo --toposhift=$toposhift --toposcale=$toposcale

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# create dummy interface at 500 m depth for USGS Vs30
if [ -e topo_data/ptopo.xyz.1.dat ]; then
  # dummy at 500m depth
  echo "creating dummy interface at 500m depth"
  awk '{printf("%.6f\n",$1-500.0);}' topo_data/ptopo.xyz.1.dat > ./topo_data/ptopo.xyz.2.500m.dat
  cd DATA/meshfem3D_files/
  rm -f ptopo.xyz.2.500m.dat
  ln -s ../../topo_data/ptopo.xyz.2.500m.dat
  cd ../../
fi

echo
echo

##
## IRIS EMC model
##
if [ "${use_emc}" == "1" ]; then
  # download IRIS EMC model
  if [ ! -e IRIS_EMC/${emc_model} ]; then
    echo
    echo "######################################"
    echo "downloading IRIS EMC model..."
    echo "######################################"
    echo
    echo "model: ${emc_model}"
    echo

    mkdir -p IRIS_EMC
    cd IRIS_EMC/

    wget https://ds.iris.edu/ds/products/script/emcscript/meta_2.py?model=${emc_model}
    # checks exit code
    if [[ $? -ne 0 ]]; then exit 1; fi

    mv -v meta_2.py\?model=${emc_model} ${emc_model}.metadata.txt
    # checks exit code
    if [[ $? -ne 0 ]]; then exit 1; fi

    wget https://ds.iris.edu/files/products/emc/emc-files/${emc_model}

    # checks exit code
    if [[ $? -ne 0 ]]; then exit 1; fi

    cd ../
  fi

  # EMC model
  echo "#######################################"
  echo "setting up IRIS EMC tomography model..."
  echo "#######################################"
  echo

  if [ ! -e ./run_convert_IRIS_EMC_netCDF_2_tomo.py ]; then
    ln -s ../../../utils/scripts/run_convert_IRIS_EMC_netCDF_2_tomo.py
  fi

  ./run_convert_IRIS_EMC_netCDF_2_tomo.py --EMC_file=IRIS_EMC/${emc_model} --mesh_area=$regionShortExtended --maximum_depth="${DEPTH}"
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  mkdir -p DATA/tomo_files
  mv -v tomography_model.* DATA/tomo_files/

  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
fi

##
## USGS Vs30
##
if [ "${use_vs30}" == "1" ]; then
  echo
  echo "#######################################"
  echo "setting up USGS VS30..."
  echo "#######################################"
  echo

  if [ ! -e ./run_get_simulation_USGS_Vs30.py ]; then
    ln -s ../../../utils/scripts/run_get_simulation_USGS_Vs30.py
  fi

  ./run_get_simulation_USGS_Vs30.py $regionExtended

  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  if [ -e USGS_VS30/interface_vs30.dat ]; then
    echo "adding Vs30 interface to DATA/"
    cd DATA/
    ln -s ../USGS_VS30/interface_vs30.dat
    cd ../
  fi
fi

echo
echo "done"
echo
