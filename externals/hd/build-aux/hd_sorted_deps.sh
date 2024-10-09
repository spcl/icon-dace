#!/bin/sh

MY_DIR=`dirname "$0"`

"${MY_DIR}/mkhelper/deplist.py" --reverse -t hd -f - <<_EOF
hd: mpi netcdf-fortran oasis yac
netcdf-fortran: netcdf
yac: netcdf mpi
oasis: netcdf mpi
_EOF

