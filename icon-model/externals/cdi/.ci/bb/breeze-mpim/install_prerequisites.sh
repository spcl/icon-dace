#!/bin/bash

set -eu
unset CDPATH

module purge
module load git

install_dir='/data/mpi/sclab/sip/m300488/libcdi-ci-sw'
work_dir="$(pwd)/build"
make_cmd='make -j22'

mkdir -p "${work_dir}" && cd "${work_dir}"

# Get PPM:
wget https://swprojects.dkrz.de/redmine/attachments/download/525/ppm-1.0.8.2.tar.gz
tar xvf ppm-1.0.8.2.tar.gz
ppm_src_dir="${work_dir}/ppm-1.0.8.2"
ppm_name_tag='ppm-1.0.8.2'
ppm_config_args='--enable-MPI --disable-netcdf --disable-hdf5 --disable-parmetis --disable-metis --disable-crypto'

# Get YAXT:
git clone -b release-0.10.0 https://gitlab.dkrz.de/dkrz-sw/yaxt.git
git -C yaxt cherry-pick 602493aad8c6e817f32c9a4889fc2a271573f896
yaxt_src_dir="${work_dir}/yaxt"
yaxt_name_tag='yaxt-0.10.0'
yaxt_config_args=''

export CC='mpicc'
export FC='mpif90'

# Install for GCC 12.1.0:
module load mpich/4.1.2-gcc-12.1.0
mpi_name_tag='mpich-4.1.2'
compiler_name_tag='gcc-12.1.0'

# Install PPM:
build_dir="${ppm_name_tag}-${mpi_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${ppm_src_dir}/configure" ${ppm_config_args} --prefix="${install_dir}/${ppm_name_tag}-${mpi_name_tag}-${compiler_name_tag}" FCFLAGS='-g -O2 -fallow-argument-mismatch'
  $make_cmd
  $make_cmd check
  $make_cmd install )

# Install YAXT:
build_dir="${yaxt_name_tag}-${mpi_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${yaxt_src_dir}/configure" ${yaxt_config_args} --prefix="${install_dir}/${yaxt_name_tag}-${mpi_name_tag}-${compiler_name_tag}"
  $make_cmd
  $make_cmd check
  $make_cmd install )

module unload mpich/4.1.2-gcc-12.1.0

# Install for NVHPC 22.3:
module load mpich/4.1.2-nvhpc-23.7
mpi_name_tag='mpich-4.1.2'
compiler_name_tag='nvhpc-23.7'

# Install PPM:
build_dir="${ppm_name_tag}-${mpi_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${ppm_src_dir}/configure" ${ppm_config_args} --prefix="${install_dir}/${ppm_name_tag}-${mpi_name_tag}-${compiler_name_tag}" FCFLAGS='-g -O2 -tp sandybridge' CFLAGS='-g -O2 -tp sandybridge'
  $make_cmd
  $make_cmd check
  $make_cmd install )

# Install YAXT:
build_dir="${yaxt_name_tag}-${mpi_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${yaxt_src_dir}/configure" ${yaxt_config_args} --prefix="${install_dir}/${yaxt_name_tag}-${mpi_name_tag}-${compiler_name_tag}" FCFLAGS='-g -O2 -tp sandybridge' CFLAGS='-g -O2 -tp sandybridge'
  $make_cmd
  $make_cmd check
  $make_cmd install )

module unload mpich/4.1.2-nvhpc-23.7

# Install for Clang 14.0.6:
module load openmpi/4.1.3-clang-14.0.6
mpi_name_tag='openmpi-4.1.3'
compiler_name_tag='clang-14.0.6'

# Install PPM:
build_dir="${ppm_name_tag}-${mpi_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${ppm_src_dir}/configure" ${ppm_config_args} --prefix="${install_dir}/${ppm_name_tag}-${mpi_name_tag}-${compiler_name_tag}" FCFLAGS='-g -O2 -no-pie'
  $make_cmd
  $make_cmd check
  $make_cmd install )

# Install YAXT:
build_dir="${yaxt_name_tag}-${mpi_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${yaxt_src_dir}/configure" ${yaxt_config_args} --prefix="${install_dir}/${yaxt_name_tag}-${mpi_name_tag}-${compiler_name_tag}" FCFLAGS='-g -O2 -no-pie'
  $make_cmd
  $make_cmd check
  $make_cmd install )

module unload openmpi/4.1.3-clang-14.0.6

# Install for NAG 7.1.7125:
module load mpich/4.1.2-nag-7.1.7125
mpi_name_tag='mpich-4.1.2'
compiler_name_tag='nag-7.1.7125'

# Install PPM:
build_dir="${ppm_name_tag}-${mpi_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${ppm_src_dir}/configure" ${ppm_config_args} --prefix="${install_dir}/${ppm_name_tag}-${mpi_name_tag}-${compiler_name_tag}" FC=no
  $make_cmd
  $make_cmd check
  $make_cmd install )

# Install YAXT:
build_dir="${yaxt_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${yaxt_src_dir}/configure" ${yaxt_config_args} --prefix="${install_dir}/${yaxt_name_tag}-${mpi_name_tag}-${compiler_name_tag}"
  $make_cmd
  $make_cmd check
  $make_cmd install )

module unload mpich/4.1.2-nag-7.1.7125
