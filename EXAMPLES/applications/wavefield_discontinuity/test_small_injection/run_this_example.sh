#!/bin/bash

## load all necessary modules here
module load intel/2020u4 intelmpi/2020u4

## set number of processors
NPROC=4
GPU_MODE=false # set to true if using GPU

currentdir=`pwd`

## set the path to SPECFEM here
specfem_dir=~/specfem3d


echo
echo "   running a small example for wavefield injection"
echo


BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p OUTPUT_FILES/
rm -rf OUTPUT_FILES/*
mkdir -p $BASEMPIDIR
rm -rf $BASEMPIDIR/*

# links executables
mkdir -p bin
cd bin/
rm -f *
ln -s ${specfem_dir}/bin/xmeshfem3D
ln -s ${specfem_dir}/bin/xdecompose_mesh
ln -s ${specfem_dir}/bin/xgenerate_databases
#ln -s ${specfem_dir}/bin/xcombine_vol_data_vtk
ln -s ${specfem_dir}/bin/xspecfem3D
cd ../

# stores setup
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/FORCESOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# prepare to run the mesher
# at this point the meshfem3D files must be well-prepared,
# the Par_file should be set properly to turn on IS_WAVEFIELD_DISCONTINUITY,
# an empty FORCESOLUTION file should be present to indicate
# there is no external source,
# and a wavefield_discontinuity_box file to designate the interface
mkdir -p MESH/

######################################################################
## running meshfem3D, here I adopt a strategy that first run meshfem3D
## in seriel, output the mesh in CUBIT format, and then decompose
## the mesh with decompose_mesh
sed -i "/^NPROC/c\NPROC                           = 1" DATA/Par_file
sed -i "/^NPROC_XI/c\NPROC_XI                        = 1" DATA/meshfem3D_files/Mesh_Par_file
sed -i "/^NPROC_ETA/c\NPROC_ETA                        = 1" DATA/meshfem3D_files/Mesh_Par_file
sed -i "/^SAVE_MESH_AS_CUBIT/c\SAVE_MESH_AS_CUBIT              = .true." DATA/meshfem3D_files/Mesh_Par_file

echo
echo "  running mesher in seriel..."
echo
./bin/xmeshfem3D

mkdir -p MESH-default
cp -f MESH/* MESH-default

sed -i "/^NPROC/c\NPROC                           = ${NPROC}" DATA/Par_file

echo
echo "  decomposing mesh..."
echo
./bin/xdecompose_mesh $NPROC ./MESH-default $BASEMPIDIR

######################################################################
## note that this can be replaced by directly running meshfem3D
## in parallel
#sed -i "/^NPROC/c\NPROC                           = ${NPROC}" DATA/Par_file
#sed -i "/^NPROC_XI/c\NPROC_XI                        = 2" DATA/meshfem3D_files/Mesh_Par_file
#sed -i "/^NPROC_ETA/c\NPROC_ETA                        = 2" DATA/meshfem3D_files/Mesh_Par_file
#sed -i "/^SAVE_MESH_AS_CUBIT/c\SAVE_MESH_AS_CUBIT              = .false." DATA/meshfem3D_files/Mesh_Par_file
#
#echo
#echo "  running mesher in parallel..."
#echo
#mpirun -np ${NPROC} ./bin/xmeshfem3D
########################################################################

echo
echo "  generate databases..."
echo
mpirun -np $NPROC ./bin/xgenerate_databases

# generate injection wavefield by producing
# DATABASES_MPI/proc*_wavefield_discontinuity.bin files
mpirun -np $NPROC ../fk_coupling/compute_fk_injection_field

if ${GPU_MODE}; then
  sed -i "/^GPU_MODE/c\GPU_MODE                           = .true." DATA/Par_file
  echo
  echo "  launch solver on GPU..."
  echo
  CUDA_VISIBLE_DEVICES=0 mpirun -np $NPROC ./bin/xspecfem3D
else
  sed -i "/^GPU_MODE/c\GPU_MODE                           = .false." DATA/Par_file
  echo
  echo "  launch solver on CPU ..."
  echo
  mpirun -np $NPROC ./bin/xspecfem3D
fi

# compute reference seismograms using FK
../fk_coupling/compute_fk_receiver

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

