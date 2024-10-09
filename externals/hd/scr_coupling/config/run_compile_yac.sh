#!/bin/sh
#
# run_compile_yac.sh - Compile YAC on DKRZ HLRE Levante
# 
# Copyright (C) 2022 DKRZ, MPI-M
# SPDX-License-Identifier: BSD-3-Clause
# See ./LICENSES/ for license information
#
# Authors: Moritz Hanke (DKRZ)
# Contact: <hanke@dkrz.de>
#_________________________________________
#
#module load python3                 # must be loaded before mpiifort
#module load openmpi/4.1.2-intel-2021.5.0               # for mpiifort
#module load netcdf-fortran/4.5.3-intel-oneapi-mpi-2021.5.0-intel-2021.5.0
#module load netcdf-c/4.8.1-intel-oneapi-mpi-2021.5.0-intel-2021.5.0
module load intel-oneapi-mkl/2022.0.1-gcc-11.2.0

### spack find -dp /k6xq5g
autoreconf -i

configure \
--with-yaxt-root=/work/gg0302/g260062/GCOAST_oas5/yaxt \
--with-netcdf-root=/sw/spack-levante/netcdf-c-4.8.1-f7hh57 \
--with-fyaml-include=/usr/include \
--with-fyaml-lib=/usr/lib \
--with-external-lapack=fortran --with-external-mtime \
--prefix=/work/gg0302/g260062/GCOAST_oas5/yac \
CC=/sw/spack-levante/openmpi-4.1.2-yfwe6t/bin/mpicc \
FC=/sw/spack-levante/openmpi-4.1.2-yfwe6t/bin/mpifort \
CFLAGS="-O2 -g"  \
FCFLAGS="-O2 -g"  \
MTIME_CFLAGS="-I/work/gg0302/g260062/GCOAST_oas5/mtime/include"
FORTRAN_LAPACK_CLIBS= \
MPI_LAUNCH="$(which mpirun)"

make install

#make check
