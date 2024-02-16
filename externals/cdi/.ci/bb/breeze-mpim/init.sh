# Initialize utility functions:
. "$(cd "$(dirname "${BASH_SOURCE}")"; pwd)/../common_utils.sh"

#
# Initializes the environment.
#
init_env ()
{
  switch_for_module automake autoconf libtool swig ruby python
}

#
# Initializes the Debian libtool
init_debian_libtool ()
{
  export PATH="/usr/bin:${PATH-}"
  export ACLOCAL_PATH="/usr/share/aclocal:${ACLOCAL_PATH-}"

  # Prepend ACLOCAL_PATH with the path to Automake:
  module unload automake
  switch_for_module automake
}

#
# Sets variables for tests with GCC.
#
init_gcc ()
{
  init_env

  switch_for_module gcc/12.1.0 mpich/3.4.3-gcc-12.1.0

  CC=gcc
  CXX=g++
  FC=gfortran
  MPICC=mpicc
  MPIFC=mpif90
  MPI_LAUNCH="$(which mpirun)"

  ECCODES_ROOT='/sw/bullseye-x64/packages/gcc-12.1.0/eccodes-2.26.0'
  NETCDF_ROOT='/sw/bullseye-x64/packages/gcc-12.1.0/netcdf-c-4.9.0'
  PPM_ROOT='/data/mpi/sclab/sip/m300488/libcdi-ci-sw/ppm-1.0.8.1-gcc-12.1.0'
  YAXT_ROOT='/data/mpi/sclab/sip/m300488/libcdi-ci-sw/yaxt-0.9.3.1-gcc-12.1.0'

  # The installations of NetCDF and ecCodes libraries do not provide '*.la'
  # files, which would trigger the RPATH injection:
  export LD_LIBRARY_PATH="${ECCODES_ROOT}/lib:${NETCDF_ROOT}/lib:${LD_LIBRARY_PATH-}"
}

#
# Sets variables for tests with NVHPC.
#
init_nvhpc ()
{
  init_env

  switch_for_module nvhpc/22.3 mpich/3.4.3-nvhpc-22.3

  CC=nvc
  CXX=nvc++
  FC=nvfortran
  MPICC=mpicc
  MPIFC=mpif90
  MPI_LAUNCH="$(which mpirun)"

  ECCODES_ROOT='/sw/bullseye-x64/packages/nvhpc-22.3/eccodes-2.26.0'
  NETCDF_ROOT='/sw/bullseye-x64/packages/nvhpc-22.3/netcdf-c-4.9.0'
  PPM_ROOT='/data/mpi/sclab/sip/m300488/libcdi-ci-sw/ppm-1.0.8.1-nvhpc-22.3'
  YAXT_ROOT='/data/mpi/sclab/sip/m300488/libcdi-ci-sw/yaxt-0.9.3.1-nvhpc-22.3'

  # The installations of NetCDF and ecCodes libraries do not provide '*.la'
  # files, which would trigger the RPATH injection:
  export LD_LIBRARY_PATH="${ECCODES_ROOT}/lib:${NETCDF_ROOT}/lib:${LD_LIBRARY_PATH-}"
}

#
# Sets variables for tests with Clang.
#
init_clang ()
{
  # Make the system gfortran the first gfortran in the PATH:
  export PATH="/usr/bin:${PATH-}"

  init_env

  switch_for_module clang/14.0.6 openmpi/4.1.3-clang-14.0.6

  CC=clang
  CXX=clang++
  FC=gfortran
  MPICC=mpicc
  MPIFC=mpif90
  MPI_LAUNCH="$(which mpirun) --oversubscribe"

  ECCODES_ROOT='/sw/bullseye-x64/packages/clang-14.0.6/eccodes-2.26.0'
  NETCDF_ROOT='/sw/bullseye-x64/packages/clang-14.0.6/netcdf-c-4.9.0-openmpi-4.1.3'
  PPM_ROOT='/data/mpi/sclab/sip/m300488/libcdi-ci-sw/ppm-1.0.8.1-clang-14.0.6'
  YAXT_ROOT='/data/mpi/sclab/sip/m300488/libcdi-ci-sw/yaxt-0.9.3.1-clang-14.0.6'

  # The installations of NetCDF and ecCodes libraries do not provide '*.la'
  # files, which would trigger the RPATH injection:
  export LD_LIBRARY_PATH="${ECCODES_ROOT}/lib:${NETCDF_ROOT}/lib:${LD_LIBRARY_PATH-}"
}

#
# Sets variables for tests with NAG.
#
init_nag ()
{
  # Make the system gcc the first gcc in the PATH:
  export PATH="/usr/bin:${PATH-}"

  init_env

  switch_for_module nag/7.1 mpich/3.4.3-nag-7.1

  CC=gcc
  CXX=g++
  FC=nagfor
  MPICC=mpicc
  MPIFC=mpif90
  MPI_LAUNCH="$(which mpirun)"

  ECCODES_ROOT='/sw/bullseye-x64/packages/nag-7.1/eccodes-2.26.0'
  NETCDF_ROOT='/sw/bullseye-x64/packages/nag-7.1/netcdf-c-4.9.0-mpich-3.4.3'
  PPM_ROOT='/data/mpi/sclab/sip/m300488/libcdi-ci-sw/ppm-1.0.8.1-nag-7.1'
  YAXT_ROOT='/data/mpi/sclab/sip/m300488/libcdi-ci-sw/yaxt-0.9.3.1-nag-7.1'

  # The installations of NetCDF and ecCodes libraries do not provide '*.la'
  # files, which would trigger the RPATH injection:
  export LD_LIBRARY_PATH="${ECCODES_ROOT}/lib:${NETCDF_ROOT}/lib:${LD_LIBRARY_PATH-}"
}
