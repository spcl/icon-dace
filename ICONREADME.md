# Steps

Clone DaCe, checkout `f2dace/staging`
Clone icon-artifacts, checkout `new_sched_solve_nh`
Clone icon-dace, checkout `yakup/ICON_24_10_merge_v2`

Go to icon-artifacts/velocity
Download data by calling 
```bash
./download_data_nproma20480.sh
```

Build CPU library by
```sh
spack load gcc@14.2
spack load netcdf-fortran@4.6.1%gcc@14.2
spack load cuda%gcc@14.2
spack load openmpi%gcc@14.2
export _USE_CUDA_EVENTS=0 # Use CUDA events for timing if set, else Host C++ timers
export _RELEASE=1 # Rel flags
export _TILE=0 # Experimental flag, that causes numerical invalidation sometimes, runs at stage 8
export _REDUCE_BITWIDTH_TRANSFORMATION=0 # Runs on stage8, irrelevant
export _BUILD_LIB_FOR_SOLVE_NH=0 # Only when the integration to solve nh is necessary
export GENCODE_NUMBER=80 # Ampere
python -m utils.stages.compile_gpu_stage5 --no-optimize --compile
```

Build ICON on AULT, CPU/Only:
```sh
spack load gcc@14.2
spack load netcdf-c@4.9.2%gcc@14.2
spack load netcdf-fortran@4.6.1%gcc@14.2
spack load cuda
spack load openmpi%gcc@14.2
spack load libxml2%gcc@14.2
spack load openblas%gcc@14.2
export _USE_CUDA_EVENTS=0 # Use CUDA events for timing if set, else Host C++ timers
export _RELEASE=1 # Rel flags
export _TILE=0 # Experimental flag, that causes numerical invalidation sometimes, runs at stage 8
export _REDUCE_BITWIDTH_TRANSFORMATION=0 # Runs on stage8, irrelevant
export _BUILD_LIB_FOR_SOLVE_NH=0 # Only when the integration to solve nh is necessary
export GENCODE_NUMBER=80 # Ampere
export _RELEASE=TRUE
export _USE_NVHPC=0
export _CUDA_ARCH=80
export SPATH=$(pwd)
export VELOCITY_PATH=${SCRATCH}/icon-artifacts/velocity
export DACE_LDFLAGS="-L${VELOCITY_PATH} -DDACE_SUBST_VERIFY=1 -DDACE_SUBST_ENABLE=1"
export DACE_LIBS="-lvelocity_gpu_stage5_standalone_release"
export DACE_FCFLAGS="-DDACE_SUBST_VERIFY=1 -DDACE_SUBST_ENABLE=1 -L${VELOCITY_PATH} -Wl,-rpath,${VELOCITY_PATH}"
export LD_LIBRARY_PATH="${VELOCITY_PATH}:${LD_LIBRARY_PATH}"

export F90=mpifort
export FC=mpifort
export F77=mpifort
export CC=mpicxx
export CXX=mpicxx
export OMPI_CC=gcc
export OMPI_CXX=g++
export OMPI_FC=gfortran

export OPENBLAS_PREFIX=$(spack location -i openblas%gcc@14.2)
export LIBXML2_PREFIX=$(spack location -i libxml2%gcc@14.2)
# Do not use local netcdf
export NETCDF_FORTRAN_PREFIX=$(spack location -i netcdf-fortran@4.6.1%gcc@14.2)
export NETCDF_C_PREFIX=$(spack location -i netcdf-c%gcc@14.2)

export CXXFLAGS=" \
  -I${NETCDF_C_PREFIX}/include \
  -I${LIBXML2_PREFIX}/include \
  -I${LIBXML2_PREFIX}/include/libxml2 \
  -I${OPENBLAS_PREFIX}/include \
  -fPIC  -Wno-implicit-function-declaration \
  ${CXXFLAGS}  \
"

export CFLAGS=" \
  -I${NETCDF_C_PREFIX}/include \
  -I${LIBXML2_PREFIX}/include \
  -I${LIBXML2_PREFIX}/include/libxml2 \
  -I${OPENBLAS_PREFIX}/include \
  -fPIC -Wno-implicit-function-declaration \
  ${CFLAGS}  \
"

export CPATH=" \
  ${NETCDF_C_PREFIX}/include:\
  ${LIBXML2_PREFIX}/include/libxml2:\
  ${LIBXML2_PREFIX}/include:\
  ${OPENBLAS_PREFIX}/include:\
  ${CPATH} \
"

export CPLUS_PATH=" \
  ${NETCDF_C_PREFIX}/include:\
  ${LIBXML2_PREFIX}/include/libxml2:\
  ${LIBXML2_PREFIX}/include:\
  ${OPENBLAS_PREFIX}/include:\
  ${CPATH} \
"

export FCFLAGS=" \
  -I${NETCDF_FORTRAN_PREFIX}/include \
  -I${NETCDF_C_PREFIX}/include \
  -I${OPENBLAS_PREFIX}/include \
  -fPIC \
  ${FCFLAGS} \
"

export LDFLAGS=" \
  -L${NETCDF_FORTRAN_PREFIX}/lib \
  -L${NETCDF_C_PREFIX}/lib \
  -Wl,-rpath,${NETCDF_FORTRAN_PREFIX}/lib \
  -Wl,-rpath,${NETCDF_C_PREFIX}/lib \
  -L${LIBXML2_PREFIX}/lib \
  -L${OPENBLAS_PREFIX}/lib \
  -Wl,-rpath,${OPENBLAS_PREFIX}/lib \
  -Wl,-rpath,${LIBXML2_PREFIX}/lib \
  ${DACE_LDFLAGS} -L${VELOCITY_PATH} -Wl,-rpath,${VELOCITY_PATH} \
  ${LDFLAGS} \
"

export LD_LIBRARY_PATH=" \
  ${NETCDF_FORTRAN_PREFIX}/lib:\
  ${NETCDF_C_PREFIX}/lib:\
  ${LIBXML2_PREFIX}/lib:\
  ${OPENBLAS_PREFIX}/lib:\
  ${VELOCITY_PATH} \
  ${LD_LIBRARY_PATH} \
"

mkdir -p icon-model/build/verification
cd icon-model/build/verification
../../config/cscs/ault_yakup_dace_cpu
make -j
cp -R ../../grids .
./make_runscripts --all
cd run
```