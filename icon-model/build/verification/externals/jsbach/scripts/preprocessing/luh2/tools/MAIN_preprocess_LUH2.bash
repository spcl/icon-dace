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
# Script to preprocess LUH2 data for JSBACH CMIP6-like simulations
# (updated according to methodology used for recent GCB simulations)
#----------------------------------------------------------------------
# JN 2019-2022 (julia.nabel@mpimet.mpg.de)
# SF July-September 2022 (stefanie.falk@lmu.de)
#--------------------------------------------------------------------
### Batch Queuing System is SLURM
#SBATCH --output=YOUR-TAG_preprocess_LUH2_forcing.o%j
#SBATCH --error=YOUR-TAG_preprocess_LUH2_forcing.o%j
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=yourname@YOUR-INSTITUTE-DOMAIN
#SBATCH --account=YOUR-PROJECT
#SBATCH --partition=interactive
#SBATCH --mem-per-cpu=3000
#SBATCH --exclusive    #-- exclusive required when aggregate_LUH2_transitions.bash or aggregate_LUH2_harvest.bash work on a large transitions_selectedYears.nc file
#SBATCH --ntasks=1

#############################################################################
# settings in examples/user_config to be adapted to user needs

ln -sf YOUR-CONFIG-FILE user_config
source ./user_config

#############################################################################
# settings probably not to be adapted any more

source ./config

selEndYearTransitions=$(expr $selEndYearStates - 1)
cdoGrid="t${ires}grid"
LUH_statesFilePrefix='LUH_states_'
LUH_mappedStatesFilePrefix="LUH_states_T${ires}_"
LUH_mappedStatesCF1FilePrefix="LUH_states_T${ires}_CF1_"
stateTagScaled='_scaled'
LUH_scaledMappedStatesFilePrefix="LUH_states_T${ires}${stateTagScaled}_"
LUH_transitionsPrefix="LUH_transitions_"
LUH_mappedTransitionsPrefix="LUH_transitions_T${ires}_"
LUH_mappedFlippedTransitionsPrefix="LUH_transitions_T${ires}_flipped_"
LUH_harvestFilePrefix='LUH_harvest_'
LUH_mappedHarvestPrefix="LUH_harvest_T${ires}_"

selTrans=transitions_selectedYears.nc
selStates=states_selectedYears.nc

stateTagUnscaled='_CF1' # i.e. use "LUH_states_T${ires}_CF1_" files
transitionsTag='_flipped'

#############################################################################
#############################################################################
# For standard LUH2 data processing
if [ $user_defined_states_trans = false ]; then
  statesFile="states.nc"
  transitionFile="transitions.nc"
  echo "Looking for files:" $statesFile $transitionFile in ${pathToOrgData}
fi

# prepare paths
scaledAndMappedStatesPath="${mainWorkingPath}/updated_states_out"
mappedHarvestPath="${mainWorkingPath}/lu_out_harvest"
forMappingTransPath="${mainWorkingPath}/lu_out"

if [ ${orgDataHasAlreadyBeenAggregated} = true ] && [ ${resolutionHasAlreadyBeenUsedForMapping} = false ]; then
  aggregatedStatesPath="${mainWorkingPathForAggregatedData}/updated_states_out"
  if [ ! -d "${aggregatedStatesPath}" ]; then
    echo "ERROR: flag indicates completed aggregation, but path to aggregated scenario states data '${aggregatedStatesPath}' does not exist, please check!"
    exit
  fi
  aggregatedHarvestPath="${mainWorkingPathForAggregatedData}/lu_out_harvest"
  if [ ! -d "${aggregatedHarvestPath}" ]; then
    echo "ERROR: flag indicates completed aggregation, but path to aggregated scenario harvest data '${aggregatedHarvestPath}' does not exist, please check!"
    exit
  fi
  if [ ${processTransitions} = true ]; then
    aggregatedTransPath="${mainWorkingPathForAggregatedData}/lu_out"
    if [ ! -d "${aggregatedTransPath}" ]; then
      echo "ERROR: flag indicates completed aggregation, but path to aggregated scenario trans data '${aggregatedTransPath}' does not exist, please check!"
      exit
    fi
  fi
else
  # thus either aggregation is required, or resolution has already been used for mapping - in both cases the folders are/will be below the current resolution folder
  aggregatedStatesPath="${scaledAndMappedStatesPath}"
  aggregatedHarvestPath="${mappedHarvestPath}"
  aggregatedTransPath="${forMappingTransPath}"
fi

if [ ${resolutionHasAlreadyBeenUsedForMapping} = true ]; then
  if [ ! -d "${aggregatedStatesPath}" ]; then
    echo "ERROR: flag indicates completed CMIP5 processing, but path to aggregated scenario states data '${aggregatedStatesPath}' does not exist, please check!"
    exit
  fi
  if [ ! -d "${aggregatedHarvestPath}" ]; then
    echo "ERROR: flag indicates completed CMIP5 processing, but path to aggregated scenario harvest data '${aggregatedHarvestPath}' does not exist, please check!"
    exit
  fi
  if [ ${processTransitions} = true ]; then
    if [ ! -d "${aggregatedTransPath}" ]; then
      echo "ERROR: flag indicates completed CMIP5 processing, but path to aggregated scenario trans data '${aggregatedTransPath}' does not exist, please check!"
      exit
    fi
  fi
else
  [ ! -d ${mainWorkingPath} ] && mkdir ${mainWorkingPath}
fi

cd ${mainWorkingPath}

#############################################################################
# main

if [ ! ${resolutionHasAlreadyBeenUsedForMapping} = true ]; then
  #############################################################################
  #-- extract the required years from the original data
    if [ ! ${orgDataHasAlreadyBeenAggregated} = true ]; then
        cdo -s -setmissval,"${missval}" -selyear,${selStartYear}/${selEndYearStates} ${pathToOrgData}/${statesFile} ${mainWorkingPath}/${selStates}
        # the transitions file is also required if [ ! ${processTransitions} = true ] because of the data on harvest
        cdo -s -setmissval,"${missval}" -selyear,${selStartYear}/${selEndYearTransitions} ${pathToOrgData}/${transitionFile} ${mainWorkingPath}/${selTrans}
    fi

  #############################################################################
  #################################################### first step: "CMIP5 way"
  # write timeinfo.nml
cat > timeinfo.nml <<EOF
&timeinfo
  ires      =${ires}
  year_start=${selStartYear}
  year_end  =${selEndYearTransitions}
  workingPath='${mainWorkingPath}'
/
EOF

  #-- process states
  echo "########## process states"
  if [ ! ${orgDataHasAlreadyBeenAggregated} = true ]; then 
    ${scriptsFolder}/${tools_dir}/aggregate_LUH2_states.bash ${mainWorkingPath}/${selStates} "${aggregatedStatesPath}" ${fstnfFile}
  else
    [ ! -d ${scaledAndMappedStatesPath} ] && mkdir ${scaledAndMappedStatesPath}
    cp -s ${aggregatedStatesPath}/${LUH_statesFilePrefix}????.nc ${scaledAndMappedStatesPath}
  fi
  ${scriptsFolder}/${tools_dir}/scale_and_remap_states.bash "${scaledAndMappedStatesPath}" ${LUH_statesFilePrefix} \
      ${LUH_mappedStatesFilePrefix} ${LUH_mappedStatesCF1FilePrefix} ${selStartYear} ${selEndYearStates} ${cdoGrid}

  #-- process transitions
  if [ ${processTransitions} = true ]; then
    echo "########## process transitions"
    if [ ! ${orgDataHasAlreadyBeenAggregated} = true ]; then 
      ${scriptsFolder}/${tools_dir}/aggregate_LUH2_transitions.bash ${mainWorkingPath}/${selTrans} "${aggregatedTransPath}" ${fstnfFile}
    else
      [ ! -d ${forMappingTransPath} ] && mkdir ${forMappingTransPath}
      cp -s ${aggregatedTransPath}/${LUH_transitionsPrefix}????.nc ${forMappingTransPath}
    fi
    if [ ! -d "${mainWorkingPath}/lu_out_gauss" ]; then
      mkdir "${mainWorkingPath}/lu_out_gauss"
    fi
    ${scriptsFolder}/${bin_dir}/Calculate_relative_transitions_remap_and_flipLats.x
    ${scriptsFolder}/${tools_dir}/flip_latdata_LUH2_transitions.bash "${mainWorkingPath}/lu_out_gauss" ${LUH_mappedTransitionsPrefix} ${LUH_mappedFlippedTransitionsPrefix} \
        ${selStartYear} ${selEndYearTransitions} 
  fi

  #-- process harvest
  echo "########## process harvest"
  if [ ! ${orgDataHasAlreadyBeenAggregated} = true ]; then 
    ${scriptsFolder}/${tools_dir}/aggregate_LUH2_harvest.bash ${mainWorkingPath}/${selTrans} "${aggregatedHarvestPath}"
  else
    [ ! -d ${mappedHarvestPath} ] && mkdir ${mappedHarvestPath}
    cp -s ${aggregatedHarvestPath}/${LUH_harvestFilePrefix}????.nc ${mappedHarvestPath}
  fi
  ${scriptsFolder}/${bin_dir}/Remap_harvest_and_flipLats.x
  ${scriptsFolder}/${tools_dir}/prevent_negative_harvest_LUH2_and_add_attributes.bash "${mappedHarvestPath}" ${LUH_mappedHarvestPrefix} \
      ${selStartYear} ${selEndYearTransitions} "${calledFrom}" "${luhRelease}"
fi

#############################################################################
#################################################### second step: "CMIP6 way"
#-> scale states and transitions with given veg_ratio_max to better achieve LUH crop and pasture area

#-- process states
echo "########## scale states"
if [ ! -d "${mainWorkingPath}/${scaledStatesOutputFolder}" ]; then
  mkdir "${mainWorkingPath}/${scaledStatesOutputFolder}"
fi
${scriptsFolder}/${tools_dir}/scale_states_with_veg-ratio-max_v2.1-check-0.bash  "${scaledAndMappedStatesPath}" "${mainWorkingPath}/${scaledStatesOutputFolder}" ${LUH_mappedStatesCF1FilePrefix} \
    ${LUH_scaledMappedStatesFilePrefix} ${vegRatioMaxFile} ${vrmVarName} ${selStartYear} ${selEndYearStates} "${calledFrom}" "${luhRelease}"

#-- process transitions
if [ ${processTransitions} = true ]; then
  # write scaling_of_transition.nml
cat > scaling_of_transition.nml <<EOF
&scaling_of_transition
  resolution=${ires}
  year_start=${selStartYear}
  year_end  =${selEndYearTransitions}     !Note: states required for the next year also!
  workingPath='${mainWorkingPath}'
  vegRatioMaxFileNameWithPath='${vegRatioMaxFile}'
  vrmVarName='${vrmVarName}'
  nrOfVarsvegMaxFile=${nrOfVarsvegMaxFile}
  nrOfDimsVegMaxFile=${nrOfDimsVegMaxFile}
  folderStateFiles='/${scaledStatesOutputFolder}/'
  folderUnscaledStateFiles='/updated_states_out/'
  folderTransitionFiles='/lu_out_gauss/'
  outputFolder='/${scaledTransitionsOutputFolder}/'
  stateTagScaled='${stateTagScaled}'
  stateTagUnscaled='${stateTagUnscaled}'
  transitionsTag='${transitionsTag}'
  calledFrom='${calledFrom}'
  luhRelease='${luhRelease}'
/
EOF

  echo "########## scale transitions"
  if [ ! -d "${mainWorkingPath}/${scaledTransitionsOutputFolder}" ]; then
    mkdir "${mainWorkingPath}/${scaledTransitionsOutputFolder}"
  fi
  ${scriptsFolder}/${bin_dir}/Adapting_transitions_after_scaling_of_states.x
fi 

echo "***Finished***"
