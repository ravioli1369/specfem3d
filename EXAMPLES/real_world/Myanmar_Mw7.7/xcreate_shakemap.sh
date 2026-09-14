#!/bin/bash

pg=$1

if [ "$pg" == "" ]; then echo "usage: ./xcreate_shakemap.sh PG-Value[1==displ/2==vel/3==acc]"; exit 1; fi

./bin/xcreate_movie_shakemap_AVS_DX_GMT << EOF
2
-1
$pg
2
EOF

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# update parameter name in shakemap file
if [ -e OUTPUT_FILES/AVS_shaking_map.inp ]; then
echo "renaming shaking map variable..."

case "$pg" in
1)
  echo "peak-ground is 1: PGD"
  sed -i "s:^ a, .*: PGD, b:" OUTPUT_FILES/AVS_shaking_map.inp
  ;;
2)
  echo "peak-ground is 2: PGV"
  sed -i "s:^ a, .*: PGV, b:" OUTPUT_FILES/AVS_shaking_map.inp
  ;;
3)
  echo "peak-ground is 3: PGA"
  sed -i "s:^ a, .*: PGA, b:" OUTPUT_FILES/AVS_shaking_map.inp
  ;;
*)
  echo "peak-ground is $pg: something else"
  ;;
esac
fi

echo
echo "done"
echo

