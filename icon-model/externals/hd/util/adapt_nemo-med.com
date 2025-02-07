#!/bin/ksh
#
# adapt_nemo-med.com - Reads land sea mask from NEMO-input and adapts it for use by convert_inflow.com
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
#
# *** Script that reads Land sea mask from NEMO-input file and adapts it for the input
#     to convert_inflow.com
#
set -ex
#
# *********** Data for Program *************************
#dnam=$1
dnam=runoff_med7km.nc
dnout=lsm_med7km.nc
#
cdo selvar,lsm,nav_lat,nav_lon $dnam $dnout
ncrename -h -O -v nav_lat,lat $dnout
ncrename -h -O -v nav_lon,lon $dnout
ncrename -h -O -v lsm,mask $dnout
ncatted -O -h -a coordinates,mask,o,c,"lon lat" $dnout
#

