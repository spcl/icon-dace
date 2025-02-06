# ICON Verification Pipeline

## Setting Up The Container / Docker
If you want to use the docker container (I would suggest WSE on Windows, or directly working on Ubuntu, if Mac use Docker otherwise it won't work):
```sh
# RM if you have it
sudo docker rm icon_gpu_container
# Clone this repo
git clone --recursive git@spclgitlab.ethz.ch:benweber/icon-dace-integration.git
cd icon-dace
git checkout yakup/ICON_24_10_merge_v2
cd ..
# More convenient to use in Docker container this way
git clone --recursive git@github.com:spcl/dace.git
# Checkout f2dace/dev
cd dace
git checkout f2dace/dev
cd ..
# Create the Docker container (assuming you have GPUs, otherwise without --gpu all)
# Replace path <path>/<to>/<icon-dace>
sudo docker run --name icon_gpu_container --gpus all -it -v /home/primrose/Work/IconGrounds/icon-dace:/home/icon/icon icon_nvhpc_rev1:latest
apt update -y && apt upgrade -y
# You might need to do chmod -R 777 on all elements
# Go to icon user (sudo + ICON sometimes have problems)
su icon
cd ~
# Create pip env and install dace
python3 -m venv .venv
source .venv/bin/activate
cd icon
pip install --upgrade setuptools wheel pip
pip install -e ./dace[testing]
```
To re-run again:
```sh
sudo docker start -ai icon_gpu_container
```

To attach the 2nd terminal:
```sh
sudo docker exec -it icon_gpu_container /bin/bash
```

To build ICON go to the icon-dace repo now:
If you are in Docker, do not forget:
```sh
su icon
cd /home/icon/icon
source ~/.venv/bin/activate
```

To download and load the container:
```sh
wget -O icon_demo_docker_image.tar.gz https://polybox.ethz.ch/index.php/s/uwFpDFFTBMmMZrN/download
docker load --input icon_demo_docker_image.tar.gz
```

## Adding a New SDFG to the Verification Pipeline

### What Happens?
We want to run ICON by replacing a Fortran call with a call to an SDFG.
Ben's existing scripts generate Fortran bindings, copy in and copy out functions.
There are manual steps if we want to achieve these, what does need to happen:

1. Generate SDFG (assuming we have a working SDFG)
2. Compile SDFG
3. Generate bindings, copy_in and copy_out and verify functions
4. Insert them into correct places
5. Build ICON

To get all of these steps we need to pass the following information to `dace-genfi` and `dace-substf` scripts:
1. Provide locations to where we need to cut the function calls and write in the new call to the SDFG
2. Provide where SDFG inputs are defined (the types to copy in and out)

### How to Add the SDFG?

Let's say you want to add a new SDFG, let's go through all the steps you need to do:
1. Copy the SDFG you want to validate to `icon-model/sdfgs`
    1. For the tutorial I have already copied `radiation.sdfgz`
2. You need to define the file `optimize_${SDFG_NAME}` that optimized the SDFG.
    1. If you have already optimize it oustide, then you can just load and save it.
    2. I have the file defined that does nothing.
    3. If there are members in structs that are not supposed to be copied over then you need define, one name at a time in `unused_names.txt`
    4. If there are structs that are not supposed to be copied (the struct becomes empty) then you need to write that in `meta_data.yaml`
3. You need to define the integrations:
   1. The yaml for ingerations `integrations.yaml` can look like this. It is a list of the following format.
   ```yaml
   ./src/atm_phy_nwp/mo_nwp_ecrad_interface.f90:
      radiation:
        - [432, 442]
   ```
   ```
    <file where we have the function call that we want replace or verify>
      <name of the SDFG>:
        - [<lines we want to cut start>:<lines we want to cut end>]
   ```
   If we need to replace 1 call per SDFG it works fine.
4. You need to define the `create_sdfg` method in `generate.py`
    1. Just check for yout SDFG name and copy the template from a previous function.
    ```python
       if args.sdfg_name == "radiation":
            reate_radiation_sdfg(args.output_folde required_associations)
    ```
5. You need to add the dictionary/mapping of SDFG function input arguments to Fortran objects/stuff being passed. The map `<key>:<value>` is `<C++ input argument name>:<Fortran input obj name>`. Example:
   ```python
    sub_dict3 = dict(
        aerosol="ecrad_aerosol",
        cloud="ecrad_cloud",
        config="ecrad_conf",
        flux="ecrad_flux",
        gas="ecrad_gas",
        single_level="ecrad_single_level",
        thermodynamics="ecrad_thermodynamics",
        iendcol="i_endidx_rad",
        istartcol="i_startidx_rad",
        ncol="nproma_sub",
        nlev="nlev",
        nulout="0",
        sym_iendcol="i_endidx_rad",
        sym_istartcol="i_startidx_rad",
    )
    RADIATION_ASSOCIATIONS = {
        "mo_nwp_ecrad_interface.f90": {
            432: {442: sub_dict3},
        }
    }
    ```
6. You need to define the modules needed for the used types in the function. It is map of `<key>:<value>` equals `<Fortran type name>:<Fortran module where it is defined>`. Example:
   ```python
    RADIATION_MODULE_DEFINITIONS = {
        "config_type": "radiation_config",
        "aerosol_type": "radiation_aerosol",
        "single_level_type": "radiation_single_level",
        "thermodynamics_type": "radiation_thermodynamics",
        "gas_type": "radiation_gas",
        "flux_type": "radiation_flux",
        "cloud_type": "radiation_cloud",
        "t_opt_ptrs": "mo_ecrad",
        "ecrad_hyd_list": "mo_ecrad",
        "ecrad_iqr": "mo_ecrad",
        "ecrad_iqs": "mo_ecrad",
        "ecrad_iqg": "mo_ecrad",
        "nulout": "radiation_io",
        "ncol": "mo_parallel_config",
        "nproma_sub": "mo_parallel_config",
        "c_null_ptr": "iso_c_binding",
    }
   ```

Once these steps are done you can call the preprocessing script and then the configure command to compile. But some remarks about the command:

The `SDFG_LIB_PATHS` variable needs to include all compiled library locations. The `SDFG_LIBS` variable needs to have all library link commands (it should be the `-l<sdfg name>`).
```sh
export SDFG_LIB_PATHS="-L<path/to/sdfg_lib.so> -L..."
export SDFG_LIBS="-l<sdfg_name>"
```

Depending on the mode you need to define the some compiler definitions.

1. If you want substite Fortran function with SDFG (replace):
```sh
-DDACE_SUBST_ENABLE
```
2. If you want to just call Fortran functions:
```sh
# No definitions
```
3. If you want to call both and numerically verify call:
```sh
-DDACE_SUBST_ENABLE -DDACE_SUBST_VERIFY
```

The complete command then looks like this:
```sh
./pipeline.sh
cd icon-scratchpad/icon-model
mkdir -p build/verification/
cd build/verification/
export SDFG_LIB_PATHS="-L/home/primrose/Work/IconGrounds/icon-dace/.dacecache/radiation/build"
export SDFG_LIBS="-lradiation"
../../config/generic/gcc \
    CC=mpicc \
    FC=mpif90 \
    FCFLAGS="-g -O2 -I/usr/include/ -Wall -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-missing-include-dirs -DDACE_SUBST_VERIFY -DDACE_SUBST_ENABLE" \
    LDFLAGS="-L/usr/lib/x86_64-linux-gnu/ ${SDFG_LIB_PATHS}" \
    LIBS="-leccodes -lnetcdff -lnetcdf -lopenblas ${SDFG_LIBS}" \
    --enable-acm-license \
    --disable-mixed-precision \
    --disable-edmf \
    --disable-les \
    --disable-ocean \
    --disable-jsbach \
    --disable-coupling \
    --disable-aes \
    --disable-rte-rrtmgp \
    --enable-ecrad \
    --disable-mpi \
    --disable-mpi-checks \
    --disable-openmp \
    --disable-loop-exchange \
    --enable-dace-subst=no \
    --enable-explicit-fpp
make -j16
```

After this you need to make run scripts and run by:
```sh
./make_runscripts --all
cd run
./exp.exclaim_R2B09.run
```

Some small netcdf files and grids I have already hardcoded into the repository.

Final question:

## How Does the Preprocessing Pipeline Works?

We need to replace existing source codes with new ones to make things work. For this purpose the preprocessing does the following:

1. Copy the whole `icon-model` ICON sources as symlinks to `icon-scratchpad` directory.
2. Removes the symlinks from `.run` files and copies them normally because ICON also symlinks them and then `make_runscripts` dont work for a weird reason.
3. It compiles the SDFGs defined in `integrations.yaml` and generates bindings for them. What it does:
4. It calls `generate_sdfgs.py`, that scripts processess `integrations.yaml`
5. It loads, optimizes and compiles each SDFG
6. It generates the binding files as `${SDFG_NAME}_bindings.f90` to call the extern C SDFG functions. Using `dace-genfi.py`
7. It calls `dace-substf.py` to replace the calls to the function in the file specified by the `integrations.yaml` to insert the `copy_in`, `copy_out`, `call sdfg()`, `verify` functions. The symlink of the original file is removed and the new file is pasted instead in the scratchpad folder.
8. This way the code can be compiled with the new sources.