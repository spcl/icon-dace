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

set -x

cp iconplot.zonal.parallel.s iconplot.zonal.parallel.list \
   iconplot.zonal-vert.ncl iconplot.zonal-time.ncl \
   iconplot.maps.ncl $TMPDIR
cd $TMPDIR
pwd

/e/uhome/for1han/bin/pshell -p7 -f iconplot.zonal.parallel.list
iconplot.zonal.parallel.s
ncl iconplot.maps.ncl
