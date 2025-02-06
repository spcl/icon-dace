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
# JN 06.12.16, script to set negative LUH2 harvest to zero
# -- i.e. get input required for JSBACH simulations
# (julia.nabel@mpimet.mpg.de)

echo "--> prevent_negative_harvest_LUH2_and_add_attributes.bash"

workingPath=${1}
filePrefix=${2}
startYear=${3}
endYear=${4}
calledFrom=${5}
luhRelease=${6}

thisDate=`date +%Y-%m-%d`

year=${startYear}
while [ ${year} -le ${endYear} ]; do
    year_char=${year}
    if [ ${year} -le 999 ]   ; then   year_char=0${year} ; fi

    file=${workingPath}/${filePrefix}${year_char}.nc
    cdo -s -ifthenelse -gec,0 ${file} ${file} -gec,0 ${file} ${workingPath}/tmp.nc

    mv ${workingPath}/tmp.nc ${file}
#    ncatted -h -O -a '.',global,d,, ${file}
    ncatted -h -O -a history,global,o,c,"${thisDate} \n  created with: ${calledFrom} \n from: ${luhRelease}" ${file}

    (( year = year + 1 ))
done
