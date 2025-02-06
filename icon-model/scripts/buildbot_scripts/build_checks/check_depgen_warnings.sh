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

# Checks whether the dependency generator produces any warnings.
# Must be run in the root build directory after a successful build.

set -eu
set -o pipefail

log=$(mktemp)
trap 'rm -f -- "${log}"' EXIT

make >"${log}" 2>&1 || {
  echo "ERROR: make failed in $(pwd)" >&2
  exit 1
}

warnings=$(grep '^deplist\.py: ' "${log}") && {
  cat "${log}" >&2
  cat >&2 <<_EOF
---
ERROR: the dependency generator reported warnings (see the full log above):
${warnings}
_EOF
  exit 1
} || {
  test $? -eq 1 || {
    echo "ERROR: failed to grep '${log}'" >&2
    exit 1
  }
}
