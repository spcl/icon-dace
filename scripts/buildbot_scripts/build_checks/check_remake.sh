#!/bin/bash

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

# Checks that:
#   1) the second call to make checks the version but does not rebuild;
#   2) the output of the third call to make is identical to the output of the
#      second one.
# Must be run in the root build directory after a successful build (i.e. after
# the first call to make).

set -eu
set -o pipefail

log=$(mktemp)
trap 'rm -f -- "${log}"' EXIT

make V=0 >"${log}" 2>&1 || {
  echo "ERROR: make failed in $(pwd)" >&2
  exit 1
}

grep '^  GEN      version\.c' "${log}" >/dev/null || {
  test $? -eq 1 && {
  cat "${log}" >&2
  cat >&2 <<_EOF
---
ERROR: make did not try to regenerate 'version.c' (see the full log above)
_EOF
} || {
  echo "ERROR: failed to grep '${log}'" >&2
}
  exit 1
}

grep '^  FCLD     bin/icon' "${log}" >/dev/null && {
  cat "${log}" >&2
  cat >&2 <<_EOF
---
ERROR: make regenerated ICON executable (see the full log above)
_EOF
  exit 1
} || {
  test $? -eq 1 || {
    echo "ERROR: failed to grep '${log}'" >&2
    exit 1
  }
}

make V=0 2>&1 | diff -u "${log}" - || {
  cat >&2 <<_EOF
ERROR: the output of the third call to make is not identical to the output of
       the second call to make (see above)
_EOF
  exit 1
}
