#!/usr/bin/bash

set -e -o pipefail

mkdir build
cd build

export MPI_ROOT="/opt/nvidia/hpc_sdk/Linux_x86_64/23.11/comm_libs/openmpi4/bin"

cmake -DCMAKE_C_COMPILER="${MPI_ROOT}/mpicc" -DCMAKE_CXX_COMPILER="${MPI_ROOT}/mpic++" -DCMAKE_Fortran_COMPILER="${MPI_ROOT}/mpif90" ..

make -j 8
