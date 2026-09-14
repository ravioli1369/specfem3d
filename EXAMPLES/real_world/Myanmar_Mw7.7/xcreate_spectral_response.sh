#!/bin/bash

dir=$1

if [ "$1" == "" ]; then echo "usage: ./xcreate_spectral_response.sh directory[holding input files DB.BANGKOK.BX*.sema; e.g. OUTPUT_FILES/]"; exit 1; fi


# check needed scripts
if [ ! -e ./run_spectral_response.py ]; then ln -s ../../../utils/scripts/run_spectral_response.py; fi

# calculates spectral responses
./run_spectral_response.py $dir/DB.BANGKOK.BXN.sema
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

./run_spectral_response.py $dir/DB.BANGKOK.BXE.sema
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

./run_spectral_response.py $dir/DB.BANGKOK.BXZ.sema
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# plot combined figure w/ all direction (N/E/Z)
if [ -e ./OUTPUT_FILES/DB.BANGKOK.BXN.sema.spectral_response.dat ] && [ -e ./OUTPUT_FILES/DB.BANGKOK.BXE.sema.spectral_response.dat ] && [ -e ./OUTPUT_FILES/DB.BANGKOK.BXZ.sema.spectral_response.dat ] ; then
  ./plot_spectral_response_BANGKOK.py
fi

echo
echo "done"
echo

