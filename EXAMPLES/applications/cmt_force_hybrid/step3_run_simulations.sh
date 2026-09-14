#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --ntasks=8
#SBATCH --partition=debug

# function to change Par_file
change_parfile() {
    local param=$1
    local value=$2
    local file="DATA/Par_file"

    local oldstr=`grep "^$param " $file`
    local newstr="$param     =       $value"

    sed  "s?$oldstr?$newstr?g" $file  > tmp
    mv tmp $file
}

set -e

# USER DEFINED VARIABLES
GPU_MODE=.TRUE.

############## STOP HERE #######################

# run txt source
change_parfile GPU_MODE $GPU_MODE
change_parfile USE_FORCE_POINT_SOURCE .true.
change_parfile USE_CMT_AND_FORCE_SOURCE .true.
change_parfile USE_BINARY_SOURCE_FILE .false.
change_parfile USE_EXTERNAL_SOURCE_FILE  .true.
bash run.sh
mv OUTPUT_FILES OUTPUT_FILES.txt

#run binary source
change_parfile USE_FORCE_POINT_SOURCE .true.
change_parfile USE_CMT_AND_FORCE_SOURCE .true.
change_parfile USE_BINARY_SOURCE_FILE .true.
change_parfile USE_EXTERNAL_SOURCE_FILE  .true.
bash run.sh
mv OUTPUT_FILES OUTPUT_FILES.bin

#seperate run
change_parfile USE_FORCE_POINT_SOURCE .true.
change_parfile USE_CMT_AND_FORCE_SOURCE .false.
change_parfile USE_BINARY_SOURCE_FILE .false.
change_parfile USE_EXTERNAL_SOURCE_FILE  .true.
bash run.sh
mv OUTPUT_FILES OUTPUT_FILES.force

change_parfile USE_FORCE_POINT_SOURCE .false.
change_parfile USE_CMT_AND_FORCE_SOURCE .false.
change_parfile USE_BINARY_SOURCE_FILE .false.
change_parfile USE_EXTERNAL_SOURCE_FILE  .true.
bash run.sh
mv OUTPUT_FILES OUTPUT_FILES.moment
