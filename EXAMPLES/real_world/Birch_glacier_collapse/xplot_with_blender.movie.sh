#!/bin/bash

# number of image steps
N=200

##################################################

# image directory
mkdir -p movie
rm -rf movie/*

if [ ! -e ./plot_with_blender.py ]; then
  ln -s ../../../utils/Visualization/Blender/python_blender/plot_with_blender.py
fi

# update locations file
if [ -e OUTPUT_FILES/STATIONS_ADJOINT_FILTERED ]; then
  echo "updating locations file from OUTPUT_FILES/STATIONS_ADJOINT_FILTERED ..."
  awk '{print "-",$3,$4}' OUTPUT_FILES/STATIONS_ADJOINT_FILTERED > location.dat
  echo
fi

# render images
for ((i=1; i<=$N; i=i+1)); do

tt="$(printf "%6.6i" $i)"
echo "******************************************"
echo "movie image: $tt out of $N"
echo "******************************************"

# time format
# --time="time 13:24:06"
# DT = 0.006
# NTSTEP_BETWEEN_FRAMES = 100
DT=0.006
NTSTEP_BETWEEN_FRAMES=100

msec_per_file_float=$(echo "scale=0; $DT * $NTSTEP_BETWEEN_FRAMES * 1000" | bc)
msec_per_file=$(printf "%.0f" "$msec_per_file_float")
echo "msec per file: $msec_per_file"

# SEM files:
# start time     = 2025-05-28T13:24:20.000000Z
# total duration = 120.0
# adjoint simulation starts at 13:24:20.0 + 120.0 == 13:26:20.0
ORIG_HOUR=13
ORIG_MIN=26
ORIG_SEC=20

# total seconds in milliseconds
# msec_per_file becomes negative for time-reversal simulation
total_msec=$(( - msec_per_file * i + (ORIG_HOUR * 3600 + ORIG_MIN * 60 + ORIG_SEC) * 1000 ))
#echo "total msec: $total_msec"

hh=$(( (total_msec/1000) / 3600 ))
mm=$(( ((total_msec/1000) % 3600) / 60 ))
sec=$(( (total_msec/1000) % 60 ))
msec=$(( total_msec - 1000 * (hh * 3600 + mm * 60 + sec) ))
# single digit 660 -> 6
msec_int=$(( $msec / 100 ))
#echo "hh: $hh mm: $mm sec: $sec msec: $msec msec_int: $msec_int"

format_time=$(printf "%02d:%02d:%02d.%d" $hh $mm $sec $msec_int)
echo "time: $format_time"
echo

# clean
rm -f out.jpg

# render
./plot_with_blender.py --vtk_file=OUTPUT_FILES/AVS_movie_$tt.inp --color-max=8.e-17 --colormap=16 \
  --title="Birch glacier collapse" \
  --transparent-sea-level --background-dark --matte --suppress \
  --locations=location.dat --borders=AVS_boundaries_utm.CH.inp \
  --time="time $format_time"
echo

if [ -e out.jpg ]; then
  mv -v out.jpg movie/frame.$tt.jpg
  echo
fi

done

# render movie
./xcreate_ffmpeg_movie.slow.sh HD




