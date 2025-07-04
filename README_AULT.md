# ICON Verification Pipeline

```bash
# To be able to compile netcdf-c
spack load nvhpc
spack load openblas%nvhpc
spack load hdf5%gcc
unset CC
unset CXX
unset F90
unset CXXFLAGS
unset LDFLAGS
unset LD_LIBRARY_PATH
unset CPATH
unset FCFLAGS
export CC=nvc
export CXX=nvc++
export F90=nvfortran
export CXXFLAGS=" -I/scratch/ybudanaz/local/include -L/scratch/ybudanaz/local/lib -I$(spack location -i libxml2)/include -I$(spack location -i libxml2)/include/libxml2 ${CXXFLAGS} -fPIC"
export LDFLAGS="-L/scratch/ybudanaz/local/lib -L$(spack location -i libxml2)/lib -L$(spack location -i nvhpc)/Linux_x86_64/25.1/cuda/lib64 ${LDFLAGS}"
export FCFLAGS="-fPIC -L/scratch/ybudanaz/local/lib -L$(spack location -i nvhpc)/Linux_x86_64/25.1/cuda/lib64 ${FCFLAGS}"
export LD_LIBRARY_PATH="$(spack location -i libxml2)/lib:$(spack location -i nvhpc)/Linux_x86_64/25.1/cuda/lib64:/scratch/ybudanaz/local/lib:${LD_LIBRARY_PATH}"
export CPATH="/scratch/ybudanaz/local/include:$(spack location -i libxml2)/include/libxml2:$(spack location -i libxml2)/include:${CPATH}"
export FCFLAGS=" -L/scratch/ybudanaz/local/lib -I/scratch/ybudanaz/local/include ${FCFLAGS} -fPIC ${FCFLAGS}"
../../config/cscs/ault_yakup_dace.gpu.a100.nvidia # Ignore rsync error I guess
make -j
./make_runscripts --all
cp -R /scratch/ybudanaz/icon-dace/icon-model/grids /scratch/ybudanaz/icon-dace/icon-model/build/verification/grids
chmod a+w /scratch/ybudanaz/icon-dace/icon-model/build/verification/grids/icon_grid_0013_R02B04_G.nc
cd run
./exp.exclaim_ape_R2B09.run
```