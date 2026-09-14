#!/bin/bash

if [ "$1" == "" ]; then echo "usage: ./setup_SEM_folder.sh ending[semd/semv/sema]"; exit 1; fi

# displ traces *.semd.* or accel traces *.sema.*
ending=$1

echo
echo "using ending: $ending"
echo

currentdir=`pwd`
basedir=`basename $currentdir`

echo "current dir: $currentdir"
echo "   base dir: $basedir"
echo

# check if data ready
if [ ! network_data/ ]; then echo "check if network_data/ exists and traces are available..."; exit 1; fi

# setup SEM/ folder
mkdir -p SEM/
rm -rf SEM/*

# copy data
cp -v network_data/*.${ending}.ascii ./SEM/
echo

# rename as adjoint sources
echo "renaming endings to *.adj ..."
cd SEM/
echo "renaming endings to *.adj ..."
rename "s/.${ending}.ascii/.adj/" *.${ending}.ascii
echo "renaming code HHZ to HXZ ..."
rename 's/.HHZ./.HXZ./' *.HHZ.adj
echo "renaming code HHN to HXN ..."
rename 's/.HHN./.HXN./' *.HHN.adj
echo "renaming code HHE to HXE ..."
rename 's/.HHE./.HXE./' *.HHE.adj

# cleaning out comment lines from file
echo "removing comment lines..."
for file in *.adj; do
  sed '/^#/d' $file > __tmp__
  mv __tmp__ $file
done

cd ../
echo

# STATIONS for adjoint simulation
cp -v network_data/STATIONS DATA/STATIONS_ADJOINT

echo
echo "done"
echo
