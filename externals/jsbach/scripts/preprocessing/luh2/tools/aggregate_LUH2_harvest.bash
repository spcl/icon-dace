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
# JN 06.12.16, script to aggregate and prepre LUH2 mass harvest variables for JSBACH
# (julia.nabel@mpimet.mpg.de)
# SF 18.05.22, add config file to load nco and cdo versions in MAIN
# (stefanie.falk@lmu.de)

echo "--> aggregate_LUH2_harvest.bash"

transitionsFileWithPath=${1}
outfolder=${2}

variablesToAggregate=( primf_bioh primn_bioh secmf_bioh secyf_bioh secnf_bioh )

if [ ! -d ${outfolder} ]; then
    mkdir ${outfolder}
fi

cd ${outfolder}

#aggregate all mass harvest variables and rename variable
cdoCommand=""
index=0
(( stopIndex = ${#variablesToAggregate[@]} - 1 ))
while [ $index -lt ${stopIndex} ]; do
    cdoCommand="${cdoCommand} -add -selvar,${variablesToAggregate[$index]} ${transitionsFileWithPath}"
   (( index = index + 1 ))
done
cdoCommand="${cdoCommand} -selvar,${variablesToAggregate[$index]} ${transitionsFileWithPath}"
cdo -s -chname,primf_bioh,harvest ${cdoCommand} tmp_aggregated_harvest_selyear.nc

# invert latdata (required for expected T63 mapping resulting from the f90 programm)
cdo -s -invertlatdata tmp_aggregated_harvest_selyear.nc tmp_aggregated_harvest_selyear_flippedLatsAndData.nc
#cdo -s -invertlat -invertlatdata tmp_aggregated_harvest_selyear.nc tmp_aggregated_harvest_selyear_flippedLatsAndData.nc

# remove wrong attributes
ncatted -a long_name,harvest,d,, -a standard_name,harvest,d,, tmp_aggregated_harvest_selyear_flippedLatsAndData.nc

# convert harvest variable to double (currently it is float but the f90 mapping programm is only implemented for double)
${ncap} -s 'harvest=double(harvest)' tmp_aggregated_harvest_selyear_flippedLatsAndData.nc tmp_aggregated_harvest_selyear_double_flippedLatsAndData.nc

# JSBACH and the JSBACH scripting expect one file per year
cdo -s -splityear tmp_aggregated_harvest_selyear_double_flippedLatsAndData.nc LUH_harvest_

rm tmp_aggregated_harvest_selyear.nc tmp_aggregated_harvest_selyear_flippedLatsAndData.nc tmp_aggregated_harvest_selyear_double_flippedLatsAndData.nc
