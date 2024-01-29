#/bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# Checks that the ICON executable does not depend on the MPI library.
# Must be run in the root build directory after a successful build of a non-MPI
# configuration (i.e.
# --disable-mpi --disable-parallel-netcdf --disable-coupling --disable-yaxt).

set -eu
set -o pipefail

executable='bin/icon'

needed=$(readelf -d "${executable}" | grep '(NEEDED)')
if test 0 -ne $? || test -z "${needed}"; then
  echo "ERROR: failed to get a list of direct library dependencies of '${executable}'" >&2
  exit 1
fi

# In the following we assume that if the executable depends on a shared MPI
# library, it's soname contains the 'libmpi' substring:
if echo "${needed}" | grep 'libmpi'; then
  cat >&2 <<_EOF
ERROR: detected a direct dependency of '${executable}' on an MPI library (see above)
_EOF
  exit 1
fi
