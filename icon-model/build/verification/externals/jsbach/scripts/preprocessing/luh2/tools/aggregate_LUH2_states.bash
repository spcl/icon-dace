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
# JN 23.07.18, script to aggregate LUH2 type states to LUH1 type states - here with rangelands treated depending
# on the fstnf value in the staticData_quarterdeg.nc (if 1, i.e. forests in LUH2, there is a conversion i.e. to pasture; if 0: no conversion i.e. natural)
# 
# called as ./aggregate_LUH2_states.bash statesFile outfolder fstnfFile
# (julia.nabel@mpimet.mpg.de)

echo "--> aggregate_LUH2_states.bash"

statesFile=${1}
outfolder=${2}
fstnfFile=${3}

variablesToAggregateFor_gcrop=( c3ann c3per c4ann c4per c3nfx )

if [ ! -d ${outfolder} ]; then
    mkdir ${outfolder}
fi

cd ${outfolder}

#-- get the pasture and natural vegetation share from rangelands
cdo -s -mul -selvar,range ${statesFile} ${fstnfFile} tmp_rangeland_to_pasture.nc
cdo -s -mul -selvar,range ${statesFile} -mulc,-1 -subc,1 ${fstnfFile} tmp_rangeland_to_nat_veg.nc

#--- aggregate gother: primf and primn
cdo -s -chname,primf,gothr -add -selvar,primf ${statesFile} -selvar,primn ${statesFile} tmp_aggregated_gothr_selyear.nc
ncatted -a long_name,gothr,d,, tmp_aggregated_gothr_selyear.nc

#--- aggregate gsecd: secdf secdn range*
cdo -s -chname,secdf,gsecd -add -selvar,secdf ${statesFile} -selvar,secdn ${statesFile} tmp2_aggregated_gsecd_selyear.nc
cdo -s -add tmp2_aggregated_gsecd_selyear.nc tmp_rangeland_to_nat_veg.nc tmp_aggregated_gsecd_selyear.nc
ncatted -a long_name,gsecd,d,, tmp_aggregated_gsecd_selyear.nc

#--- aggregate gpast: pastr urban range*
cdo -s -chname,pastr,gpast -add -selvar,pastr ${statesFile} -selvar,urban ${statesFile} tmp2_aggregated_gpast_selyear.nc
cdo -s -add tmp2_aggregated_gpast_selyear.nc tmp_rangeland_to_pasture.nc tmp_aggregated_gpast_selyear.nc
ncatted -a long_name,gpast,d,, tmp_aggregated_gpast_selyear.nc

#--- aggregate gcrop
cdoCommand=""
index=0
(( stopIndex = ${#variablesToAggregateFor_gcrop[@]} - 1 ))
while [ $index -lt ${stopIndex} ]; do
    cdoCommand="${cdoCommand} -add -selvar,${variablesToAggregateFor_gcrop[$index]} ${statesFile}"
   (( index = index + 1 ))
done
cdoCommand="${cdoCommand} -selvar,${variablesToAggregateFor_gcrop[$index]} ${statesFile}"
cdo -s -chname,${variablesToAggregateFor_gcrop[0]},gcrop ${cdoCommand} tmp_aggregated_gcrop_selyear.nc
ncatted -a long_name,gcrop,d,, tmp_aggregated_gcrop_selyear.nc

# join the 4 files
cdo -s merge tmp_aggregated_?????_selyear.nc tmp_all_states.nc
rm tmp*_aggregated_*_selyear.nc

# JSBACH and the JSBACH scripting expect one file per year
cdo -s -splityear tmp_all_states.nc LUH_states_
rm tmp_*