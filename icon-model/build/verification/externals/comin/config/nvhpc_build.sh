#!/usr/bin/bash

set -e -o pipefail

mkdir build
cd build

cmake \
    -DCMAKE_VERBOSE_MAKEFILE=ON \
    -DBUILD_TESTING=ON \
    -DCOMIN_ENABLE_PYTHON_ADAPTER=ON \
    -DCOMIN_ENABLE_REPLAY_TOOL=ON \
    -DCMAKE_C_COMPILER="${MPI_ROOT}/mpicc" \
    -DCMAKE_CXX_COMPILER="${MPI_ROOT}/mpic++" \
    -DCMAKE_Fortran_COMPILER="${MPI_ROOT}/mpif90" ..

make -j 8
