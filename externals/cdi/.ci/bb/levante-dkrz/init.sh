# Initialize utility functions:
. "$(cd "$(dirname "${BASH_SOURCE}")"; pwd)/../common_utils.sh"

#
# Initializes the environment.
#
init_env ()
{
  switch_for_module ruby/3.0.2-gcc-11.2.0
}

#
# Sets variables for tests with GCC.
#
init_gcc ()
{
  init_env
  switch_for_module gcc/11.2.0-gcc-11.2.0 openmpi/4.1.2-gcc-11.2.0

  CC=gcc
  CXX=g++
  FC=gfortran
  MPICC=mpicc
  MPIFC=mpif90
  MPI_LAUNCH="$(which mpirun)"

  ECCODES_ROOT='/sw/spack-levante/eccodes-2.21.0-4ywkk4'
  NETCDF_ROOT='/sw/spack-levante/netcdf-c-4.8.1-6qheqr'
  PPM_ROOT='/work/mh0287/m300488/libcdi-ci-sw/install/ppm-1.0.8-gcc-11.2.0'
  YAXT_ROOT='/work/mh0287/m300488/libcdi-ci-sw/install/yaxt-0.9.3-gcc-11.2.0'

  # The installations of NetCDF and ecCodes libraries do not provide '*.la'
  # files, which would trigger the RPATH injection:
  export LD_LIBRARY_PATH="${ECCODES_ROOT}/lib64:${NETCDF_ROOT}/lib:${LD_LIBRARY_PATH-}"
}

#
# Sets variables for tests with Intel Classic.
#
init_intelclassic ()
{
  init_env
  # Try to make sure that the compiler works with the system gcc:
  module unload gcc
  # For whatever reason, Intel Classic 2021.5.0 is enabled with
  # intel-oneapi-compilers/2022.0.1-gcc-11.2.0:
  switch_for_module intel-oneapi-compilers/2022.0.1-gcc-11.2.0 openmpi/4.1.2-intel-2021.5.0

  AR=xiar
  CC=icc
  CXX=icpc
  FC=ifort
  MPICC=mpicc
  MPIFC=mpif90
  MPI_LAUNCH="$(which mpirun)"

  ECCODES_ROOT='/sw/spack-levante/eccodes-2.21.0-3ehkbb'
  NETCDF_ROOT='/sw/spack-levante/netcdf-c-4.8.1-2k3cmu'
  PPM_ROOT='/work/mh0287/m300488/libcdi-ci-sw/install/ppm-1.0.8-intel-classic-2021.5.0'
  YAXT_ROOT='/work/mh0287/m300488/libcdi-ci-sw/install/yaxt-0.9.3-intel-classic-2021.5.0'

  # The installations of NetCDF and ecCodes libraries do not provide '*.la'
  # files, which would trigger the RPATH injection:
  export LD_LIBRARY_PATH="${ECCODES_ROOT}/lib64:${NETCDF_ROOT}/lib:${LD_LIBRARY_PATH-}"
}
