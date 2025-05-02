#!/bin/bash

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2024, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

#----------------------------------------------------------------------
# JN 21.02.17, script to flip the latdata resulting from the application of the programm 
# Calculate_relative_transitions_remap_and_flipLats.f90 
# -- i.e. get input in a format required for JSBACH simulations
# (julia.nabel@mpimet.mpg.de)
# SF 24.05.22, flip dimensions lat,lon,time -> time,lat,lon
# (stefanie.falk@lmu.de)

echo "--> flip_latdata_LUH2_transitions.bash"

workingPath=${1}
inFilePrefix=${2}
outFilePrefix=${3}
startYear=${4}
endYear=${5}

year=${startYear}
while [ ${year} -le ${endYear} ]; do

    year_char=${year}
    if [ ${year} -le 999 ]   ; then   year_char=0${year} ; fi

    #ncpdq -a time,lat,lon ${workingPath}/${inFilePrefix}${year_char}.nc ${workingPath}/${outFilePrefix}${year_char}_tmp.nc
    cdo -s -invertlatdata ${workingPath}/${inFilePrefix}${year_char}.nc ${workingPath}/${outFilePrefix}${year_char}.nc
    #rm ${workingPath}/${outFilePrefix}${year_char}_tmp.nc

    (( year = year + 1 ))
done
