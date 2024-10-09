# compile_toy_models.sh - Script to compile the toy models for a coupled system with YAC
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Author: Moritz Hanke (DKRZ) and Ha Ho-Hagemann (Hereon, Germany)
# Contact: ha.hagemann@hereon.de
#_________________________________________
#
#
# First, please adapt WORK_DIR & SRC_DIR if neccessary
WORK_DIR=$(pwd)
SRC_DIR=/work/gg0302/g260062/GCOAST_oas5

YAXT_DIR=${SRC_DIR}/yaxt
YAC_DIR=${SRC_DIR}/yac

export PKG_CONFIG_PATH=${YAC_DIR}/lib/pkgconfig:$PKG_CONFIG_PATH

LPATH="/sw/spack-levante/netcdf-c-4.8.1-7dq6g2/lib"
NETCDF_LIB="-L/sw/spack-levante/netcdf-c-4.8.1-7dq6g2/lib -lnetcdf -Wl,-rpath,$LPATH"

NC_INC1="-I/sw/spack-levante/netcdf-c-4.8.1-7dq6g2/include" 
NC_INC2="-I/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/include" 
NETCDF_INC="$NC_INC2 $NC_INC1"

echo ${NETCDF_LIB}
echo ${NETCDF_INC}

#=================== CHOOSE HERE WHAT TO COMPILE !!! ==============
compile_ocean="yes"
#compile_ocean="no"

compile_land="yes"
#compile_land="no"

#================================================================
if [ $compile_ocean == "yes" ];then
 mpicc -o BUILD/toy_ocean.exe -O0 -g ../code/src/toy_ocean.c  `pkg-config --libs yac`  `pkg-config --cflags yac` `pkg-config --libs-only-L yac | sed 's/-L/-Wl,-rpath,/g'`

ls -l BUILD/toy_ocean.exe
fi

#================================================================
if [ $compile_land == "yes" ];then
 mpicc -o BUILD/toy_land.exe -O0 -g ../code/src/toy_land.c  `pkg-config --libs yac`  `pkg-config --cflags yac` `pkg-config --libs-only-L yac | sed 's/-L/-Wl,-rpath,/g'`

ls -l BUILD/toy_land.exe
fi

#================================================================
cd ${WORK_DIR}

