#!/bin/bash

# simulation steps
NSTEP=`grep ^NSTEP DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

echo "Par_file:"
echo "  NSTEP = $NSTEP"
echo

# creates movie files w/ norm of velocity

./bin/xcreate_movie_shakemap_AVS_DX_GMT << EOF
2
1
$NSTEP
1
1
EOF

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo

