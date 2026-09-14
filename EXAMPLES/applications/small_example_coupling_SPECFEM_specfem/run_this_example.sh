#!/bin/bash

echo "running example: `date`"
currentdir=`pwd`

##
## Step 1 - create coupling points from small-scale local mesh
##
echo
echo "##############################################"
echo "STEP 1 - coupling points"
echo "##############################################"
echo

# symbolic link
rm -f DATA
ln -s DATA.local/ DATA

# clean out coupling folder
mkdir -p COUPLING_FILES
rm -rf COUPLING_FILES/*

# mesher
./run_mesher.sh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup
rm -rf OUTPUT_FILES.local
mv -v OUTPUT_FILES OUTPUT_FILES.local

rm -rf DATABASES_MPI.local
mv -v DATABASES_MPI DATABASES_MPI.local

##
## Step 2 - create wavefield solution from coarse regional simulation
##
echo
echo "##############################################"
echo "STEP 2 - regional simulation"
echo "##############################################"
echo

# symbolic link
rm -f DATA
ln -s DATA.regional/ DATA

# mesher
./run_mesher.sh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# solver
./run_solver.sh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup
rm -rf OUTPUT_FILES.regional
mv -v OUTPUT_FILES OUTPUT_FILES.regional

rm -rf DATABASES_MPI.regional
mv -v DATABASES_MPI DATABASES_MPI.regional

##
## Step 3 - run coupled small-scale local simulation
##
echo
echo "##############################################"
echo "STEP 3 - small-scale local simulation"
echo "##############################################"
echo

# symbolic link
rm -f DATA
ln -s DATA.local/ DATA

# restore local simulation folders
# (or re-run the mesher)
rm -rf OUTPUT_FILES DATABASES_MPI
mv -v OUTPUT_FILES.local OUTPUT_FILES
mv -v DATABASES_MPI.local DATABASES_MPI

# solver
./run_solver.sh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "all done"
echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo `date`



