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
# JN 13.01.17, script to remap and scale LUH2 states to LUH1 states ("CMIP6 way")
# (logic based on rewrite_lu_states.csh by SW/TR (?))
# --> creates two files, a mapped file and then a mapped file scaled such that cover fractions sum up to 1
# (julia.nabel@mpimet.mpg.de)

echo "--> scale_and_remap_states.bash"

workPath=${1}/
inFilePrefix=${2}
prefixMapped=${3}
prefixMappedCF1=${4}

startYear=${5}
endYear=${6}

cdoGrid=${7}

year=$startYear
while [ $year -le $endYear ]; do

    year_char=${year}
    if [ ${year} -le 999 ]   ; then   year_char=0${year} ; fi

    echo ${year_char}

    # calculate weights file
    [[ ! -f ${workPath}tmp_weights.nc ]] && cdo -s gencon,${cdoGrid} ${workPath}${inFilePrefix}${year_char}.nc ${workPath}tmp_weights.nc
    
    # map
    echo "... output to ${workPath}${prefixMapped}${year_char}.nc"
    cdo -s -remap,${cdoGrid},${workPath}tmp_weights.nc ${workPath}${inFilePrefix}${year_char}.nc ${workPath}${prefixMapped}${year_char}.nc
    # Remove lat_bnd from output file
    ncks -C -O -x -v lat_bnds ${workPath}${prefixMapped}${year_char}.nc ${workPath}${prefixMapped}${year_char}.nc
    # scale: the area fractions in several 0.25deg cells in LUH do not sum up to 1, 
    #        because water is excluded (sum(area fractions) = 1 - ice - water) [beside accuracy issues this is a static mask in LUH2]
    #           -> this is particularly true for coastline gridcells!
    #        Here the vegetation fractions are scaled up proportionally such that sum(area fractions) = 1 
    #        -- a sum if 1 is required for JSBACH
    # --> according to TR (08.02.17) this scaling statistically averages out for costline gridcells
    cdo -s -add -selvar,gothr ${workPath}${prefixMapped}${year_char}.nc -add -selvar,gpast ${workPath}${prefixMapped}${year_char}.nc \
           -add -selvar,gsecd ${workPath}${prefixMapped}${year_char}.nc -selvar,gcrop ${workPath}${prefixMapped}${year_char}.nc ${workPath}tmp_added_${year_char}.nc

    cdo -s setctomiss,0. -gtc,0.01 ${workPath}tmp_added_${year_char}.nc ${workPath}tmp_miss_${year_char}.nc
    cdo -s -mul ${workPath}tmp_added_${year_char}.nc ${workPath}tmp_miss_${year_char}.nc ${workPath}tmp_added_miss_${year_char}.nc 

    echo "... output to ${workPath}${prefixMappedCF1}${year_char}.nc"
    cdo -s -div ${workPath}${prefixMapped}${year_char}.nc ${workPath}tmp_added_miss_${year_char}.nc ${workPath}${prefixMappedCF1}${year_char}.nc
   
    # Remove temporary files
    rm ${workPath}tmp_added_${year_char}.nc ${workPath}tmp_miss_${year_char}.nc ${workPath}tmp_added_miss_${year_char}.nc 

    (( year = year + 1 ))
done

rm ${workPath}tmp_*
