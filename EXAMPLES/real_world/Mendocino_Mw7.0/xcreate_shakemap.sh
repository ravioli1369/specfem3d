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

echo
echo "done"
echo


