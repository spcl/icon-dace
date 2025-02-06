#! /bin/bash

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

#
# Check that python3 is actually loaded by module
#

set -eu
set -o pipefail

die () { echo "ERROR: $*" >&2; exit 1; }

EXPECTED_PYTHON=$(module show python3 2>&1 | 
    awk '/^prepend-path\tPATH / {print $3 "/python"}')

module load python3
REAL_PYTHON=$(type -p python)

if [ "$EXPECTED_PYTHON" != "$REAL_PYTHON" ]
then
    die "python is '$REAL_PYTHON' after 'module load python3'" \
        "(should be '$EXPECTED_PYTHON')"
fi

