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

# ----------------------------------------------------------------------------
# prepare SCM input data for multiple points
#
# attention: lon between -180 and 180!!
#
# work flow SCM from ICON input:
#  - ICON ini:   read_icon_ana_oper_mem1_40km.s
#  - ICON run:   run_ICON_4_SCMini
#  - SCM ini:    get_SCM_data_ICON.py
#  - SCM extpar: create_SCM_extpar_ICON.py
#  - SCM run:    run_SCM_ICONini
#  - plot SCM:   plot-scm-*.py
#
# 2021-02-16 Martin Koehler
# ----------------------------------------------------------------------------

for ((ii=1; ii<=17; ii+=1)) ; do

  lat_scm='-60'
  lon_scm=$((ii*20-180))
  echo '--- processing lat/lon: ' $lat_scm $lon_scm ' ---'

  python3 get_SCM_data_ICON.py $lat_scm $lon_scm

done
