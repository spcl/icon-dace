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
# JN 23.07.18, script to aggregate LUH2 type transitions to LUH1 type transitions - here with rangelands treated depending
# on the fstnf value in the staticData_quarterdeg.nc (if 1, i.e. forests in LUH2, there is a conversion i.e. to pasture; if 0: no conversion, i.e. natural vegetation)
#
# called as ./aggregate_LUH2_transitions_treat_rangelands_depending_on_fstnf_value.bash transFileWithPath outfolder fstnfFile
# (julia.nabel@mpimet.mpg.de)
#
# SF 18.05.22, add config file to load nco and cdo versions in MAIN
# (stefanie.falk@lmu.de)

echo "--> aggregate_LUH2_transitions.bash"

transFileWithPath=${1}
outfolder=${2}
fstnfFile=${3}

# NOTE: range is listed twice. Rangelands are aggregated to pasture or secondary vegetation depending on the fstnf value
org_states=( primf primn urban secdf secdn pastr range range c3ann c3per c4ann c4per c3nfx )
aggregated_states=( primary primary pasture secondary secondary pasture secondary pasture crop crop crop crop crop )

nonExistingTransitions=( primf_to_secdf primn_to_secdn range_to_range )
# and there are obviously no transitions to primf or primn

if [ ! -d ${outfolder} ]; then
    mkdir ${outfolder}
fi

cd ${outfolder}

state_from_index=0
fileList=""
while [ ${state_from_index} -lt ${#org_states[@]} ]; do
  state_to_index=0
  state_from=${org_states[$state_from_index]}
  aggregated_state_from=${aggregated_states[$state_from_index]}
  while [ ${state_to_index} -lt ${#org_states[@]} ]; do
    state_to=${org_states[$state_to_index]}
    aggregated_state_to=${aggregated_states[$state_to_index]}

    thisVar="${state_from}_to_${state_to}"
    echo "current var: ${thisVar}"

    if [ ! ${aggregated_state_from} = ${aggregated_state_to} ]; then
      thisTargetVar="${aggregated_state_from}2${aggregated_state_to}"

      # there are no transitions to primf or primn
      if [ ${state_to} = primf ] || [ ${state_to} = primn ]; then
        (( state_to_index = state_to_index + 1 ))
        continue
      fi

      doesNotExist=false
      for nonTrans in ${nonExistingTransitions[@]}; do
        if [ ${nonTrans} = ${thisVar} ]; then
          doesNotExist=true
          break
        fi
      done
      if [ $doesNotExist = true ]; then
        (( state_to_index = state_to_index + 1 ))
        continue
      fi

      echo "... adding ${thisVar} to ${thisTargetVar}"

      # get the data
      cdo -s -selvar,${thisVar} ${transFileWithPath} tmp.nc

      # convert to double (currently they are float but the f90 mapping program is only implemented for double)
      ${ncap} -s "${thisVar}=double(${thisVar})" tmp.nc tmp1.nc

      # treat rangelands
      if [ ${state_from} = range ] && [ ${aggregated_state_from} = secondary ]; then
          # only take those transitions where LUH2 does not assume forest (fstnf = 0); no land use change
          cdo -s -mul tmp1.nc -mulc,-1 -subc,1 ${fstnfFile} tmp2.nc
          mv tmp2.nc tmp1.nc
      elif [ ${state_to} = range ] && [ ${aggregated_state_to} = secondary ]; then
          cdo -s -mul tmp1.nc -mulc,-1 -subc,1 ${fstnfFile} tmp2.nc
          mv tmp2.nc tmp1.nc
      elif [ ${state_to} = range ] && [ ${aggregated_state_to} = pasture ]; then
          # only take those transitions where LUH2 assumes forest (fstnf = 1); land use change
          cdo -s -mul tmp1.nc ${fstnfFile} tmp2.nc
          mv tmp2.nc tmp1.nc
      elif [ ${state_from} = range ] && [ ${aggregated_state_from} = pasture ]; then
          cdo -s -mul tmp1.nc ${fstnfFile} tmp2.nc
          mv tmp2.nc tmp1.nc
      fi

      thisAggregationFile=tmp_${thisTargetVar}.nc
      if [ -f ${thisAggregationFile} ]; then
        cdo -s -add ${thisAggregationFile} tmp1.nc tmp2.nc
        mv tmp2.nc ${thisAggregationFile}
      else
        cdo -s -chname,${thisVar},${thisTargetVar} tmp1.nc ${thisAggregationFile}
        # remove old not suitable attributes
        ncatted -a standard_name,${thisTargetVar},d,, ${thisAggregationFile}
        fileList="${fileList} ${thisAggregationFile}"
      fi

      rm -f tmp.nc tmp1.nc tmp2.nc
    fi
    (( state_to_index = state_to_index + 1 ))
  done

  (( state_from_index = state_from_index + 1 ))
done

# join the files
cdo -s -merge ${fileList} tmp_all_transitions.nc
rm ${fileList}

# JSBACH and the JSBACH scripting expect one file per year
cdo -s -splityear tmp_all_transitions.nc LUH_transitions_
rm tmp_all_transitions.nc

