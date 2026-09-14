#!/bin/bash

pg=$1

if [ "$pg" == "" ]; then echo "usage: ./xplot_with_blender.shakemap.sh PG-Value[1==displ/2==vel/3==acc]"; exit 1; fi

# check if data available
if [ ! -e OUTPUT_FILES/shakingdata ]; then
  echo "no shakingdata file found in OUTPUT_FILES/"
  echo "Please check if shakemap was created before running this script."
  exit 1
fi


# set string name
case "$pg" in
1) PGVAL="PGD";;
2) PGVAL="PGV";;
3) PGVAL="PGA";;
*) echo "  peak-ground is $gp : something else not recognized";exit 1;;
esac
echo "peak-ground motion is: $PGVAL"
echo

# create shakemap file
./bin/xcreate_movie_shakemap_AVS_DX_GMT << EOF
2
-1
$pg
2
EOF
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# modify variable name
if [ -e OUTPUT_FILES/AVS_shaking_map.inp ]; then
  echo "renaming shaking map variable..."
  sed -i "s:^ a, .*: $PGVAL, b:" OUTPUT_FILES/AVS_shaking_map.inp
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  echo
fi

# update locations file
if [ -e OUTPUT_FILES/STATIONS_ADJOINT_FILTERED ]; then
  echo "updating locations file from OUTPUT_FILES/STATIONS_ADJOINT_FILTERED ..."
  awk '{print "-",$3,$4}' OUTPUT_FILES/STATIONS_ADJOINT_FILTERED > location.dat
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  echo
fi


# render image
echo
echo "******************************************"
echo "shakemap image"
echo "******************************************"

info="peak-ground motion ($PGVAL)"
echo "info string: $info"
echo

# clean
rm -f out.jpg

if [ ! -e ./plot_with_blender.py ]; then
  ln -s ../../../utils/Visualization/Blender/python_blender/plot_with_blender.py
fi

# render
./plot_with_blender.py --vtk_file=OUTPUT_FILES/AVS_shaking_map.inp --color-max=3.e-16 --colormap=16 \
  --title="Birch glacier collapse" \
  --transparent-sea-level --background-dark --matte --suppress \
  --locations=location.dat --borders=AVS_boundaries_utm.CH.inp \
  --time="$info" --centered

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# move image to OUTPUT_FILES/ folder
filename="shakemap_${PGVAL}_Birch_glacier_collapse.2025-05-28.jpg"

if [ -e out.jpg ]; then
  echo
  mv -v out.jpg OUTPUT_FILES/$filename
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
  echo
  echo "plotted to: OUTPUT_FILES/$filename"
  echo
fi

echo
echo "done"
echo
