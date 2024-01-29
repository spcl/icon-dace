#!/bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)

echo "Removing all spack environement in $SCRIPT_DIR"
rm  -f "$SCRIPT_DIR"/v*/*/spack.lock
rm -rf "$SCRIPT_DIR"/v*/*/.spack-env
