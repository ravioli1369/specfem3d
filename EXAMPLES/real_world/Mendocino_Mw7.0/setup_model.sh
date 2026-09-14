#!/bin/bash
#
# script to download topography and EMC model data
#
################################################

## User parameters

# topography resolution (topo30 == 30 arc-seconds)
topo=topo30

# shifted internal interface (in m)
toposhift=14000.0

# EMC model
model=WUS256.r0.0.nc

################################################

# input region
lon_min=$1
lon_max=$2
lat_min=$3
lat_max=$4

# usage
if [ "$1" == "" ]; then echo "./setup_model.sh lon_min lon_max lat_min lat_max"; exit 1; fi


region="${lon_min} ${lat_min} ${lon_max} ${lat_max}"
regionShort="${lon_min},${lat_min},${lon_max},${lat_max}"


# check needed scripts
if [ ! -e ./run_get_simulation_topography.py ]; then ln -s ../../../utils/scripts/run_get_simulation_topography.py; fi
if [ ! -e ./run_convert_IRIS_EMC_netCDF_2_tomo.py ]; then ln -s ../../../utils/scripts/run_convert_IRIS_EMC_netCDF_2_tomo.py; fi

echo
echo "model setup..."
echo
echo "region: $region"
echo

# topo
./run_get_simulation_topography.py $region --SRTM=$topo --toposhift=$toposhift

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

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
echo "setting up tomography model..."
echo
# gets depth from setup
DEPTH=`grep ^DEPTH_BLOCK_KM DATA/meshfem3D_files/Mesh_Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | sed -e 's/^[ \t]*//;s/[ \t]*$//'`
# convert number endings like 100.d0 -> 100.0
DEPTH=`echo $DEPTH | sed -e 's/d0/0/'`
echo "region: $regionShort"
echo "depth : $DEPTH"
echo

./run_convert_IRIS_EMC_netCDF_2_tomo.py --EMC_file=IRIS_EMC/${model} --mesh_area=$regionShort --maximum_depth=$DEPTH
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mkdir -p DATA/tomo_files
mv -v tomography_model.* DATA/tomo_files/

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo
