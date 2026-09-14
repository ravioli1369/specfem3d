#!/usr/bin/env python
from __future__ import print_function

import os
import sys

SEMoutput='MESH'
CUBIToutput='MESH_GEOCUBIT'

# checks for path for modules
found_lib = False
for path in sys.path:
    if "geocubitlib" in path:
        found_lib = True
        break
if not found_lib:
    sys.path.append('../../../CUBIT_GEOCUBIT/geocubitlib')
    sys.path.append('../../../CUBIT_GEOCUBIT/')
#print("path:")
#for path in sys.path: print("  ",path)
#print("")

import cubit
try:
    cubit.init([""])
except:
    pass

cubit.cmd('reset')
cubit.cmd('brick x 67000 y 134000 z 60000')
cubit.cmd('volume 1 move x 33500 y 67000 z -30000')
cubit.cmd('brick x 67000 y 134000 z 60000')
cubit.cmd('volume 2 move x 100500 y 67000 z -30000')
cubit.cmd('merge all')

# Meshing the volumes
elementsize = 3000.0

cubit.cmd('volume 1 size '+str(elementsize))
cubit.cmd('volume 2 size '+str(elementsize))
cubit.cmd('mesh volume 1 2')


from geocubitlib import boundary_definition,exportlib

# boundary
boundary_definition.define_bc(parallel=True)

# file export
exportlib.collect(outdir=CUBIToutput)
exportlib.e2SEM(outdir=SEMoutput)

cubit.cmd('save as "meshing.cub" overwrite')


