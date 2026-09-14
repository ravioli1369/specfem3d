# Couple with injection - example of hybrid method run Specfem/Specfem

This is a very small example to test the coupling between SPECFEM3D and itself.
In other words, it runs a coarse regional simulation to create an initial wavefield for a distant event,
and then re-injects that wavefield at the coupling boundaries of a small-scale local simulation.


## How to run

This coupling example uses the following steps:
1. create the coupling boundary points of the small-scale local mesh by running the mesher
   with the small-scale local setup in folder `DATA.local/`

2. run the coarse regional simulation using DATA.regional/ to store the wavefield at the boundary points
   into a separate folder `COUPLING_FILES/`

3. run the small-scale local simulation with injecting the wavefield at the boundary points

The default setup uses 4 MPI processes.

To launch the full coupling example with the above 3 steps, type:
```
./run_this_example.sh
```


## Notes

* The coupling example here between coarse regional and small-scale local mesh uses a flat top surface.
  This allows the coarse regional simulation to locate all coupling points at the exact locations.
  The resulting seismograms for stations inside the small-scale region should therefore become (almost) identical.

  You can compare the seismograms in the output folder `OUTPUT_FILES/` from the small-scale local simulation, with the
  seismograms in `OUTPUT_FILES.regional/` from the coarse regional simulation.

* In case coarse regional and small-scale local simulation would use realistic topographic interfaces, there might be a difference
  in topography resolution and thus coupling points from the small-scale coupling boundary could be located slightly
  outside the coarser regional mesh. Locating those coupling points would move them closer to the nearest coarse mesh elements.
  Storing the wavefield and re-injecting it would then lead to smaller artefacts.

  To work-around this, one would have to make sure that the top of the small-scale mesh surfaces match exactly
  the top of the coarse mesh surfaces at the injection boundaries.
  This would either require to tweak the topo interfaces for the in-house mesher or use an external mesher to resolve this.
  For "smoother" incoming wavefields and small discrepencies between surface points, the injected wavefield however
  might still be good enough.

* To save storage space, we use `NTSTEP_BETWEEN_OUTPUT_SAMPLE == 2` in the `Par_file` of the coarse regional simulation.
  This seems to provide enough time samplings for the injection. One would have to test this out a bit for other examples.
  In general, the storage of the wavefield components (velocity & traction) over all coupling boundary points and over all time steps
  can lead to significant file sizes.

  In any case, the stored wavefield will be (time) interpolated during the preparation stage of the small-scale local simulation,
  as the time steps `DT` of the coarse regional and small-scale local simulations are likely different. Time interpolation is done by a
  Catmull-Rom scheme that provides cubic interpolation between time samples. This scheme seems to exhibit a good balance between accuracy
  and speed for these kind of simulations. Thus, storing the wavefield at lower sampling rates can be fine for many applications.


## Recipes for more complex modeling

In the following, we give some short remarks how this example setup could be tweaked a bit to become more realistic.
All the mentioned scripts below will use python and require some additional installations on your system, in particular they
all need [GMT](http://gmt.soest.hawaii.edu) functionalities. Please check out the header descriptions in the scripts themself
to see what requirements they have.

* Topography:

  You can use the script in `utils/scripts/run_get_simulation_topography.py` to download and setup the topography surface
  (including one "downshifted" surface to optimize the mesh with the in-house mesher).

  For example, to download a 3-arc seconds topography for the specified coarse regional region, you would type:
  ```
  ./run_get_simulation_topography.py 7.2 47.5 8.6 48.5 --SRTM=topo3 --toposhift=10000.0 --toposcale=0.0
  ```
  The script will modify the `DATA/meshfem3D_files/Mesh_Par_file` and `DATA/Par_file`,
  as well as the `DATA/meshfem3D_files/interfaces.dat` files accordingly.

  A higher resolution topography can be used for the small-scale mesh:
  ```
  ./run_get_simulation_topography.py 7.76 47.95 7.90 48.05 --SRTM=high --dx=0.0003 --toposhift=1000.0 --toposcale=0.0
  ```

* Tomography models from IRIS EMC:

  To use a tomographic model for the region, you can check out available [IRIS EMC models](https://ds.iris.edu/ds/products/emc-earthmodels/).
  For this particular example, the region is covered for example by the
  [LSP_Eucrust1.0 model](https://ds.iris.edu/ds/products/emc-lsp_eucrust10/) by Lu et al. (2018).

  To convert this model to a tomography model usable in SPECFEM3D, you would download the data by:
  ```
  # setup a data folder
  mkdir -p IRIS_EMC; cd IRIS_EMC/

  # NetCDF data file
  wget https://ds.iris.edu/files/products/emc/emc-files/LSP-Eucrust1.0.nc

  # meta data description
  wget https://ds.iris.edu/ds/products/script/emcscript/meta_2.py?model=LSP-Eucrust1.0.nc
  ```
  and create the tomography file for SPECFEM3D by the `utils/scripts/run_convert_IRIS_EMC_netCDF_2_tomo.py` script:
  ```
  ./run_convert_IRIS_EMC_netCDF_2_tomo.py --EMC_file=IRIS_EMC/LSP-Eucrust1.0.nc --mesh_area=7.2,47.5,8.6,48.5 --maximum_depth=80.0
  ```
  using the area and depth of the coarse regional mesh used above. This script will create the needed `tomography_model.xyz` file, which you
  can for example copy into a `DATA/tomo_files/` folder.

  Once you have a tomography file, you would adjust the `TOMOGRAPHY_PATH` parameter in `DATA/Par_file`, as well
  as modify the material section accordingly in `DATA/meshfem3D_files/Mesh_Par_file`:
  ```
  # number of materials
  NMATERIALS                      = 1
  # define the different materials in the model as :
  # #material_id  #rho  #vp  #vs  #Qkappa #Qmu  #anisotropy_flag #domain_id
  #     Q                : quality factor
  #     anisotropy_flag  : 0=no anisotropy/ 1,2,.. check with implementation in aniso_model.f90
  #     domain_id        : 1=acoustic / 2=elastic
  # for tomography:
  #   #material_id #type-keyword #domain-name #tomo-filename #tomo_id #domain_id
  #   example:
  #   -1 tomography elastic tomography_model.xyz 0 2
  #
  -1 tomography elastic tomography_model.xyz 0 2
  ```

* USGS VS30 near-surface model:

  For small-scale local simulations, the mesh can be constructed such that the top elements have an enough fine grid sampling
  that a near-surface layer can be setup. For spectral-element meshes, you would want to have at least 2 grid points
  in vertical direction within the top 30-m depth. This avoids further smearing out the surface velocities to much bigger depths.

  To impose velocities from the [USGS VS30 model](https://earthquake.usgs.gov/data/vs30/), you can download and
  create a Vs30 interface file into a data folder `USGS_VS30/` by using the script `utils/scripts/run_get_simulation_USGS_Vs30.py`:
  ```
  # setup data folder
  mkdir -p USGS_VS30

  # download and create interface file
  ./run_get_simulation_USGS_Vs30.py 7.6 47.8 8.2 48.2
  ```
  This will download the USGS VS30 model, and create a file `interface_vs30.dat` that will be needed for SPECFEM3D.

  The interface file should be copied into the small-scale local `DATA/` folder if one wants to used it during the meshing stage.
  Whenever the `xgenerate_databases` tool finds this file, it will overimpose the Vs velocities as well as Vp and density derived by
  a crustal rock relationship from Brocher (2005) to the top GLL points within the 30-m depth layer.

The rest is up to your imagination! We can't wait to see what new and exciting research you come up with ;)
