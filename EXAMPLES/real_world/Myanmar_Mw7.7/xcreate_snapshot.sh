#!/bin/bash

# single snapshot
nstep=$1

if [ "$1" == "" ]; then echo "usage: ./xcreate_snapshot.sh moviedata-timestep [e.g. 10000]"; exit 1; fi

# check needed scripts
if [ ! -e ./plot_with_blender.py ]; then ln -s ../../../utils/Visualization/Blender/python_blender/plot_with_blender.py; fi

# check output file folder
if [ ! -e OUTPUT_FILES ]; then echo "sorry, folder OUTPUT_FILES/ not found..."; exit 1; fi

# creates movie files w/ norm of velocity
./bin/xcreate_movie_shakemap_AVS_DX_GMT << EOF
2
$nstep
$nstep
1
1
EOF

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# render images
echo "******************************************"
echo "movie image: 1 out of 1"
echo "******************************************"

# larger scale
./plot_with_blender.py --vtk_file=OUTPUT_FILES/AVS_movie_000001.inp --color-max=0.1 --colormap=2 --title="Myanmar" --transparent-sea-level --background-dark --locations=interesting_locations.dat --utm_zone=47 --sea-level-separation=0.001 --vertical-exaggeration=2.0

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mv -v out.jpg OUTPUT_FILES/image.snapshot_$nstep.jpg

echo
echo "plotted snapshot to: OUTPUT_FILES/image.snapshot_$nstep.jpg"
echo

