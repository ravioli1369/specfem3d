#!/usr/bin/env python
#############################################################
#
# script uses ACIS surface formats
#
# note: this script seems to work with CUBIT version > 12.2
#          meshing takes about 15 minutes (without refinement)
#
#
#############################################################
from __future__ import print_function

import os
import sys
import os.path
import time

# checks path for GEOCUBIT modules
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

#
# CUBIT
#
try:
    import cubit
except ImportError:
    print("Error: Importing cubit as python module failed")
    print("       Please check your PYTHONPATH settings...")
    print("")
    print("current path: ")
    print(sys.path)
    print("")
    sys.exit("Import cubit failed")

#cubit.init([""])
cubit.init(["-noecho","-nojournal"])

#
# GEOCUBIT
#
try:
    from geocubitlib import boundary_definition
    from geocubitlib import cubit2specfem3d
except:
    import boundary_definition
    import cubit2specfem3d

# time stamp
print("#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))


# working directory
cwd = os.getcwd()
print("# current working directory: " + str(cwd))
if cwd[len(cwd)-14:len(cwd)] != "Mount_StHelens":
  print("")
  print("# Please run this script from example directory: SPECFEM3D/EXAMPLES/applications/Mount_StHelens/")
  print("")

print("#")
print("## cubit version:")
print("#")
cubit.cmd('version')

print("#")
print("# running meshing script...")
print("#")
print("# note: this script uses topography surface in ACIS format")
print("#       meshing will take around 15 min")
print("#")

# clean workspace
cubit.cmd('reset')

# uses developer commands
cubit.cmd('set developer commands on')
cubit.cmd('set import mesh tolerance 1')

#############################################################
#
# 0. step: loading topography surface
#
#############################################################
print("#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
print("# loading topo surface...")

# topography surface
if os.path.exists("topo.cub"):
    print("# opening existing topography surface")
    # topography surface
    # previously run, just reopen the cubit file
    cubit.cmd('open "topo.cub"')
else:
    # topo surface doesn't exist yet, this creates it:
    print("# reading in topography surface")
    # reads in topography points and creates sheet surface
    # old: execfile("./read_topo.py")
    exec(open("./read_topo.py").read())
    # clear
    cubit.cmd('reset')
    # now reopen the cubit file
    cubit.cmd('open "topo.cub"')

# healing the surface...
cubit.cmd('Auto_clean volume 1 small_surfaces small_curve_size 10')
cubit.cmd('regularize volume 1')

#############################################################
#
# 1. step: creates temporary brick volume
#
#############################################################
print("#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
print("# creating brick...")

# creates temporary brick volume
if os.path.exists("topo_1.cub"):
    print("# opening existing volume 1")
    # clears workspace
    cubit.cmd('reset')
    # previously run, just reopen the cubit file
    cubit.cmd('open "topo_1.cub"')
else:
    # new brick volume (depth will become 1/2 * 20,000 = 10,000 m)
    cubit.cmd('brick x 15000 y 22000 z 20000')

    # moves volume to UTM coordinates of topography surface
    cubit.cmd('volume 2 move x 561738. y 5116370. z 0 ')

    # temporary backup
    cubit.cmd('save as "topo_1.cub" overwrite')

#############################################################
#
# 2. step: creates volume with topography surface
#
#############################################################
print("#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
print("# creating volume with topography...")

# topography surface
if os.path.exists("topo_2.cub"):
    print("# opening existing volume 2")
    # clears workspace
    cubit.cmd('reset')
    # previously run, just reopen the cubit file
    cubit.cmd('open "topo_2.cub"')
else:
    print("#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
    print("# imprinting volume, this will take around 1 min, please be patience...")
    cubit.cmd('imprint all')
    cubit.cmd('merge all')

    # exports only surfaces which will create single volume
    cubit.cmd('export acis "topo_2.acis" surface 3 10 12 14 15 9 ascii overwrite')

    # backup
    cubit.cmd('save as "topo_2.cub" overwrite')

#############################################################
#
# 3. step: manipulate ACIS file to create a single volume
#
#############################################################
# checks if new file available
if not os.path.exists("topo_2.acis"):
  print("")
  print("# error creating new volume, please check manually...")
  print("")
  cubit.cmd('pause')
# clears workspace
cubit.cmd('reset')

# single volume
if os.path.exists("topo_3.cub"):
    print("# opening existing volume 3")
    # clears workspace
    cubit.cmd('reset')
    # previously run, just reopen the cubit file
    cubit.cmd('open "topo_3.cub"')
else:
    # imports surfaces and merges to single volume
    cubit.cmd('import acis "topo_2.acis" ascii merge_globally')

    print("#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
    print("# creating new volume, this will take another 2 min...")
    cubit.cmd('create volume surface all heal')

    # backup
    cubit.cmd('save as "topo_3.cub" overwrite')

#############################################################
#
# 4. step: create mesh
#
#############################################################
print("#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))

# topography surface
if os.path.exists("topo_4.cub"):
    print("# opening existing volume 4")
    # clears workspace
    cubit.cmd('reset')
    # previously run, just reopen the cubit file
    cubit.cmd('open "topo_4.cub"')
else:
    print("# initial meshing...")
    print("# (will take around 7 min)")

    # optional: refining mesh at surface
    #
    # note: refining using ACIS surface format takes a long time ... (up to 3 hours)
    DO_TOPO_REFINEMENT = False

    # Meshing the volumes
    if DO_TOPO_REFINEMENT == False:
        elementsize = 500.0
        cubit.cmd('volume all size '+str(elementsize))
        # note: we will mesh first the topography surface, then sweep down the mesh
        # topography surface
        #cubit.cmd('control skew surface 12')
        cubit.cmd('surface 12 submap smooth off')
        cubit.cmd('surface 12 scheme submap')
        cubit.cmd('mesh surface 12')
        # propagates mesh down for whole volume
        cubit.cmd('volume 1  redistribute nodes off')
        cubit.cmd('volume 1 scheme sweep source surface 12 target surface 7')
        cubit.cmd('mesh volume all')
        # draw/update mesh lines for visualization
        # this will draw also the tripling layer mesh lines in case
        cubit.cmd('draw volume all')
    else:
        # optional surface refinement
        # time stamp
        print("#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
        print("# refining surface mesh...")
        print("# (will take around 3 hours)")
        # starts with a crude mesh
        elementsize = 2000.0
        cubit.cmd('volume all size '+str(elementsize))
        # sets meshing type
        # explicitly sets scheme for topography surface
        cubit.cmd('surface 12 submap smooth off')
        cubit.cmd('surface 12 scheme submap')
        # uses a sweep algorithm in vertical (Z) direction
        cubit.cmd('volume all scheme sweep Vector 0 0 1')
        # initial coarse mesh
        cubit.cmd('mesh volume all')
        # optional smoothing to improve mesh quality (takes up to 5 min)
        cubit.cmd('volume all smooth scheme condition number beta 1.0 cpu 5')
        cubit.cmd('smooth volume all')
        # refines global mesh
        cubit.cmd('refine volume 1 numsplit 1')
        cubit.cmd('draw volume all')
        # optional smoothing
        cubit.cmd('smooth volume all')
        # refines elements at topography surface
        cubit.cmd('refine surface 12 numsplit 1 bias 1.0 depth 1')
        cubit.cmd('draw volume all')
        # optional smoothing to improve mesh quality (takes up to 10 min)
        cubit.cmd('volume all smooth scheme condition number beta 2.3 cpu 10')
        cubit.cmd('smooth volume all')
        # displays final mesh
        cubit.cmd('draw volume all')


    # time stamp
    print("#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
    print("# done meshing...")

    # backup
    cubit.cmd('save as "topo_4.cub" overwrite')

#### End of meshing

# re-indexes
cubit.cmd('compress all')

# avoids assigning empty blocks
cubit.cmd('set duplicate block elements on')

###### This is boundary_definition.py of GEOCUBIT
#..... which extracts the bounding faces and defines them into blocks
print("#### DEFINE BC #######################")
boundary_definition.entities=['face']
boundary_definition.define_bc(boundary_definition.entities,parallel=True)

# given the tolerance on the bounding box surfaces is too narrow, the xmin/xmax/ymin/ymax boundaries might not be found.
# we can create the missing absorbing blocks manually here:
if 1 == 0:
    print("## surfaces")
    cubit.cmd('list surface')
    print("#")
    #box = cubit.get_bounding_box('surface', 5)
    #print("# surface 5: bounding box = ",box)

    cubit.cmd('block 1003 face in surface 5')
    cubit.cmd('block 1003 name "face_abs_xmin"')

    cubit.cmd('block 1004 face in surface 4')
    cubit.cmd('block 1004 name "face_abs_ymin"')

    cubit.cmd('block 1005 face in surface 3')
    cubit.cmd('block 1005 name "face_abs_xmax"')

    cubit.cmd('block 1006 face in surface 2')
    cubit.cmd('block 1006 name "face_abs_ymax"')


#### Define material properties for the 3 volumes ################
print("#### DEFINE MATERIAL PROPERTIES #######################")
cubit.cmd('block 1 name "elastic 1"')        # elastic material region
cubit.cmd('block 1 attribute count 7')
cubit.cmd('block 1 attribute index 1 1')      # flag for material: 1 for 1. material
cubit.cmd('block 1 attribute index 2 2800')   # vp
cubit.cmd('block 1 attribute index 3 1500')   # vs
cubit.cmd('block 1 attribute index 4 2300')   # rho
cubit.cmd('block 1 attribute index 5 9999.0') # Qkappa
cubit.cmd('block 1 attribute index 6 150.0')  # Qmu
cubit.cmd('block 1 attribute index 7 0')      # anisotropy_flag

# optional saves backups
cubit.cmd('export mesh "top.e" dimension 3 overwrite')
cubit.cmd('save as "meshing.cub" overwrite')

#### Export to SPECFEM3D format using cubit2specfem3d.py of GEOCUBIT

os.system('mkdir -p MESH')
cubit2specfem3d.export2SPECFEM3D('MESH')

# all files needed by SCOTCH are now in directory MESH

# time stamp
print("# all done")
print("#" + time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))


