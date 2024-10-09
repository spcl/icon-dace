#!/bin/bash -ex

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

YEAR_START=2930
YEAR_END=2959
INTERVAL=30

EXPID=slo1016

export EXPID
export DATADIR=/work/mh0287/users/stephan/Icon/Git_repos/icon.oes.mergecpl1710/experiments


export CDO='cdo -s -P 4'
export CHUNK=2  #number of years in each output file  
export GRID=R2B6
export LEV=L64
#export PLOTGRID_2D='r3600x1800'
export POOL=/mnt/lustre01/work/mh0033/m211054/projects/icon/FX
OQSDIR=$(pwd)  #add here path to the OQs scripts 

# Create plots
YEAR=$YEAR_START
while [ $YEAR -le $(( YEAR_END - INTERVAL + 1)) ]
do

    Y1=$YEAR Y2=$(( YEAR + INTERVAL - 1 )) ${OQSDIR}/plot_ocean.sh
    YEAR=$(( YEAR + INTERVAL ))
    
done

exit
