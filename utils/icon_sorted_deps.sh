#!/bin/sh

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

my_dir=`dirname "$0"`

for python in python3 python; do
  "${python}" "${my_dir}/mkhelper/deplist.py" --help >/dev/null 2>&1 && break
done

"${python}" "${my_dir}/mkhelper/deplist.py" --reverse -t icon -f - <<_EOF
icon: stdc++ cuda rocm mpi netcdf-fortran rte-rrtmgp ecrad rttov sct yaxt cdi serialbox2 mtime blas lapack yac tixi eccodes hdf5 zlib comin fortran-support
netcdf-fortran: netcdf
netcdf: hdf5 zlib
cdi: eccodes netcdf aec mpi yaxt ppm
yac: lapack mtime fyaml netcdf yaxt mpi
lapack: blas
sct: hdf5 mpi
hdf5: mpi aec zlib
yaxt: mpi
tixi: xml2
rttov: netcdf-fortran hdf5 lapack
ecrad: netcdf-fortran
serialbox2: netcdf stdc++
eccodes: aec
cuda: stdc++
rocm: stdc++ cuda
xml2: zlib
ppm: mpi netcdf
_EOF

