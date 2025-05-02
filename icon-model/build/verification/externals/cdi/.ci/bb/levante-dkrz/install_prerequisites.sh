#!/bin/bash

set -eu
unset CDPATH

module purge
module load git

install_dir='/work/mh0287/m300488/libcdi-ci-sw/install'
work_dir="$(pwd)/build"
make_cmd='make -j22'

mkdir -p "${work_dir}" && cd "${work_dir}"

# Get PPM 1.0.8:
wget https://swprojects.dkrz.de/redmine/attachments/download/517/ppm-1.0.8.tar.gz
tar xvf ppm-1.0.8.tar.gz
ppm_src_dir="${work_dir}/ppm-1.0.8"
ppm_name_tag='ppm-1.0.8'
ppm_config_args='--enable-MPI --disable-netcdf --disable-hdf5 --disable-parmetis --disable-metis --disable-crypto'

# Get YAXT 0.9.3:
git clone --depth=1 -b release-0.9.3 https://gitlab.dkrz.de/dkrz-sw/yaxt.git
yaxt_src_dir="${work_dir}/yaxt"
yaxt_name_tag='yaxt-0.9.3'
yaxt_config_args=''

export CC='mpicc'
export FC='mpif90'

# Install for GCC 11.2.0:
compiler_name_tag='gcc-11.2.0'
module load openmpi/4.1.2-gcc-11.2.0

# Install PPM:
build_dir="${ppm_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${ppm_src_dir}/configure" ${ppm_config_args} --prefix="${install_dir}/${ppm_name_tag}-${compiler_name_tag}"
  $make_cmd
  $make_cmd check
  $make_cmd install )

# Install YAXT:
build_dir="${yaxt_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${yaxt_src_dir}/configure" ${yaxt_config_args} --prefix="${install_dir}/${yaxt_name_tag}-${compiler_name_tag}"
  $make_cmd
  $make_cmd check
  $make_cmd install )

module unload openmpi/4.1.2-gcc-11.2.0

# Install for Intel Classic 2021.5.0:
compiler_name_tag='intel-classic-2021.5.0'
module load openmpi/4.1.2-intel-2021.5.0

# Install PPM:
build_dir="${ppm_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${ppm_src_dir}/configure" ${ppm_config_args} --prefix="${install_dir}/${ppm_name_tag}-${compiler_name_tag}"
  $make_cmd
  $make_cmd check
  $make_cmd install )

# Install YAXT:
build_dir="${yaxt_name_tag}-${compiler_name_tag}"
mkdir "${build_dir}"
( cd "${build_dir}"
  "${yaxt_src_dir}/configure" ${yaxt_config_args} --prefix="${install_dir}/${yaxt_name_tag}-${compiler_name_tag}"
  $make_cmd
  $make_cmd check
  $make_cmd install )

module unload openmpi/4.1.2-intel-2021.5.0
