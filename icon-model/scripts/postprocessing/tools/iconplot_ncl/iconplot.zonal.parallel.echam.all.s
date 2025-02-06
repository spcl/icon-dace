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

cp iconplot.zonal.parallel.echam.s iconplot.zonal.parallel.echam.list \
   iconplot.zonal-vert.echam.ncl iconplot.zonal-time.echam.ncl \
   iconplot.maps.echam.ncl $TMPDIR
cd $TMPDIR
pwd

/e/uhome/for1han/bin/pshell -p7 -f iconplot.zonal.parallel.echam.list
iconplot.zonal.parallel.echam.s
ncl iconplot.maps.echam.ncl
