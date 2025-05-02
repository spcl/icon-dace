# Initialize utility functions:
. "$(cd "$(dirname "${BASH_SOURCE}")"; pwd)/../common_utils.sh"

#
# Initializes the environment.
#
init_env ()
{
  AUTOMAKE_ROOT='/project/d56/libcdi-ci-sw/gcc-11.2.0-haswell/automake-1.16.3-46gthis'
  export PATH="${AUTOMAKE_ROOT}/bin:${PATH-}"
  # Tell the custom installation of Automake where the libtool macros are:
  export ACLOCAL_PATH="/usr/share/aclocal:${ACLOCAL_PATH-}"
}

#
# Sets variables for tests with Cray.
#
init_cray ()
{
  init_env
  switch_for_module craype PrgEnv-cray cce/12.0.3 cray-mpich

  # Build and test against NetCDF that does not support MPI parallel invocations
  # (parallel NetCDF4 tests are known to fail in this case when run from the
  # user home directory on Daint):
  switch_for_module cray-netcdf

  # Uncomment the following line to test against MPI-capable NetCDF:
  # switch_for_module cray-netcdf-hdf5parallel cray-hdf5-parallel

  CC=cc
  CXX=CC
  FC=ftn
  MPI_LAUNCH="$(which srun) -p cscsci -C gpu -A d56 -t 05:00"

  ECCODES_ROOT='/project/d56/libcdi-ci-sw/cce-12.0.3-haswell/eccodes-2.24.2-o2a4fw3'
  PPM_ROOT='/project/d56/libcdi-ci-sw/cce-12.0.3-haswell/scales-ppm-1.0.8-44zlrlu'
  YAXT_ROOT='/project/d56/libcdi-ci-sw/cce-12.0.3-haswell/yaxt-0.9.2.1-enz3pcz'

  # The installations of ecCodes, PPM and YAXT do not provide '*.la' files,
  # which would trigger the RPATH injection:
  export LD_LIBRARY_PATH="${ECCODES_ROOT}/lib64:${PPM_ROOT}/lib:${YAXT_ROOT}/lib:${LD_LIBRARY_PATH-}"
}

#
# Sets variables for tests with PGI.
#
init_pgi ()
{
  init_env
  # We use deprecated versions (the most recent compatible with PGI though) of
  # the Cray packages and have to make sure that the default versions are
  # unloaded (otherwise, we get various warnings and errors):
  module unload cray-mpich cray-netcdf cray-netcdf-hdf5parallel cray-hdf5 cray-hdf5-parallel
  switch_for_module craype PrgEnv-pgi/6.0.8 pgi/20.1.1 cray-mpich/7.7.15

  # Build and test against NetCDF that does not support MPI parallel invocations
  # (parallel NetCDF4 tests are known to fail in this case when run from the
  # user home directory on Daint):
  switch_for_module cray-netcdf/4.7.4.0 cray-hdf5/1.12.0.0

  # Uncomment the following line to test against MPI-capable NetCDF:
  # switch_for_module cray-netcdf-hdf5parallel/4.7.4.0 cray-hdf5-parallel/1.12.0.0

  CC=cc
  CXX=CC
  FC=ftn
  MPI_LAUNCH="$(which srun) -p cscsci -C gpu -A d56 -t 05:00"

  ECCODES_ROOT='/project/d56/libcdi-ci-sw/pgi-20.1.1-haswell/eccodes-2.24.2-hwtl5nr'
  PPM_ROOT='/project/d56/libcdi-ci-sw/pgi-20.1.1-haswell/scales-ppm-1.0.8-z2nxqya'
  YAXT_ROOT='/project/d56/libcdi-ci-sw/pgi-20.1.1-haswell/yaxt-0.9.2.1-3orop7g'

  # The deprecated versions of the Cray packages are not in the default linker
  # search path:
  export LD_LIBRARY_PATH="${MPICH_DIR}/lib:${NETCDF_DIR}/lib:${HDF5_DIR}/lib:${LD_LIBRARY_PATH-}"

  # The installations of ecCodes, PPM and YAXT do not provide '*.la' files,
  # which would trigger the RPATH injection:
  export LD_LIBRARY_PATH="${ECCODES_ROOT}/lib64:${PPM_ROOT}/lib:${YAXT_ROOT}/lib:${LD_LIBRARY_PATH}"
}
