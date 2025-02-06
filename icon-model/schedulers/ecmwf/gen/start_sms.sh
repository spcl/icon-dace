#!/bin/ksh

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# Define suite base directory
basedir=/home/ms/de/${USER}/ICON_r2B06_30d

# change working directory
cd ${basedir}/sms

# play suite:
cdp << EOF
   myalias
   play icon.def
   sms
   begin icon
EOF



