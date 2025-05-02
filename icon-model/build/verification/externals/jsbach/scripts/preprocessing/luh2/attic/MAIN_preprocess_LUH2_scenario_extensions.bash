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
# Script to preprocess LUH2 extensions data for JSBACH CMIP6 simulations
#----------------------------------------------------------------------
# JN 11.03.19 (julia.nabel@mpimet.mpg.de)
# SF 15.06.22 (stefanie.falk@lmu.de)
#--------------------------------------------------------------------
### Batch Queuing System is SLURM
#SBATCH --output=T63_preprocess_LUH2_extensions_dynveg.o%j
#SBATCH --error=T63_preprocess_LUH2_extensions_dynveg.o%j
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=yourname@mpimet.mpg.de
#SBATCH --account=bm1241
#SBATCH --partition=interactive
#SBATCH --mem-per-cpu=3000
#SBATCH --exclusive    #-- exclusive required when aggregate_LUH2_transitions.bash or aggregate_LUH2_harvest.bash work on a large transitions_selectedYears.nc file
#SBATCH --ntasks=1

#############################################################################
# settings to be adapted to user needs

ln -sf examples/user_config_sfa_lamaclima_test user_config
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
# iterate over scenarios
if [ ! ${#scenPart[@]} -eq ${#scenario[@]} ]; then
  echo "ERROR: scenPart needs to have same number of elements as scenario, please check!"
  exit
fi

scenIndex=0
while [ ${scenIndex} -lt ${#scenPart[@]} ]; do
  thisScenario=${scenario[${scenIndex}]}
  echo "============================================ ${thisScenario}"
  thisScenPart=${scenPart[${scenIndex}]}

  filesPrefix="multiple-" 
  filesCommon="_extension_input4MIPs_landState_ScenarioMIP_UofMD-"
  filesSuffix="-2-1-e_gn_2015-2100.nc"

  # For standard LUH2 data processing
  if [ $user_defined_states_trans = false ]; then
    statesFile="${filesPrefix}states${filesCommon}${thisScenPart}${filesSuffix}"
    transitionFile="${filesPrefix}transitions${filesCommon}${thisScenPart}${filesSuffix}"
    echo $statesFile $transitionFile
    exit
  fi

  # prepare paths
  scenarioWorkPath=${mainWorkingPath}/${thisScenario}
  scaledAndMappedStatesPath="${scenarioWorkPath}/updated_states_out"
  mappedHarvestPath="${scenarioWorkPath}/lu_out_harvest"
  forMappingTransPath="${scenarioWorkPath}/lu_out"

  if [ ${orgDataHasAlreadyBeenAggregated} = true ] && [ ${resolutionHasAlreadyBeenUsedForMapping} = false ]; then
    aggregatedStatesPath="${mainWorkingPathForAggregatedData}/${thisScenario}/updated_states_out"
    if [ ! -d "${aggregatedStatesPath}" ]; then
      echo "ERROR: flag indicates completed aggregation, but path to aggregated scenario states data '${aggregatedStatesPath}' does not exist, please check!"
      exit
    fi
    aggregatedHarvestPath="${mainWorkingPathForAggregatedData}/${thisScenario}/lu_out_harvest"
    if [ ! -d "${aggregatedHarvestPath}" ]; then
      echo "ERROR: flag indicates completed aggregation, but path to aggregated scenario harvest data '${aggregatedHarvestPath}' does not exist, please check!"
      exit
    fi
    if [ ${processTransitions} = true ]; then
      aggregatedTransPath="${mainWorkingPathForAggregatedData}/${thisScenario}/lu_out"
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
    [ ! -d ${scenarioWorkPath} ] && mkdir ${scenarioWorkPath}
  fi

  cd ${scenarioWorkPath}

  #############################################################################
  # main

  if [ ! ${resolutionHasAlreadyBeenUsedForMapping} = true ]; then
    #############################################################################
    #-- extract the required years from the original data
      if [ ! ${orgDataHasAlreadyBeenAggregated} = true ]; then
          cdo -s -selyear,${selStartYear}/${selEndYearStates} ${pathToOrgData}/${statesFile} ${scenarioWorkPath}/${selStates}
          # the transitions file is also required if [ ! ${processTransitions} = true ] because of the data on harvest
          cdo -s -selyear,${selStartYear}/${selEndYearTransitions} ${pathToOrgData}/${transitionFile} ${scenarioWorkPath}/${selTrans}
      fi

    #############################################################################
    #################################################### first step: "CMIP5 way"
    # write timeinfo.nml
cat > timeinfo.nml <<EOF
&timeinfo
  ires      =${ires}
  year_start=${selStartYear}
  year_end  =${selEndYearTransitions}
  workingPath='${scenarioWorkPath}'
/
EOF

    #-- process states
    echo "########## process states"
    if [ ! ${orgDataHasAlreadyBeenAggregated} = true ]; then 
      ${scriptsFolder}/${tools_dir}/aggregate_LUH2_states.bash ${scenarioWorkPath}/${selStates} "${aggregatedStatesPath}"
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
        ${scriptsFolder}/${tools_dir}/aggregate_LUH2_transitions.bash ${scenarioWorkPath}/${selTrans} "${aggregatedTransPath}"
      else
        [ ! -d ${forMappingTransPath} ] && mkdir ${forMappingTransPath}
        cp -s ${aggregatedTransPath}/${LUH_transitionsPrefix}????.nc ${forMappingTransPath}
      fi
      if [ ! -d "${scenarioWorkPath}/lu_out_gauss" ]; then
        mkdir "${scenarioWorkPath}/lu_out_gauss"
      fi
      ${scriptsFolder}/${bin_dir}/Calculate_relative_transitions_remap_and_flipLats.x
      ${scriptsFolder}/${tools_dir}/flip_latdata_LUH2_transitions.bash "${scenarioWorkPath}/lu_out_gauss" ${LUH_mappedTransitionsPrefix} ${LUH_mappedFlippedTransitionsPrefix} \
          ${selStartYear} ${selEndYearTransitions} 
    fi

    #-- process harvest
    echo "########## process harvest"
    if [ ! ${orgDataHasAlreadyBeenAggregated} = true ]; then 
      ${scriptsFolder}/${tools_dir}/aggregate_LUH2_harvest.bash ${scenarioWorkPath}/${selTrans} "${aggregatedHarvestPath}"
    else
      [ ! -d ${mappedHarvestPath} ] && mkdir ${mappedHarvestPath}
      cp -s ${aggregatedHarvestPath}/${LUH_harvestFilePrefix}????.nc ${mappedHarvestPath}
    fi
    ${scriptsFolder}/${bin_dir}/Remap_harvest_and_flipLats.x
    ${scriptsFolder}/${tools_dir}/prevent_negative_harvest_LUH2_and_add_attributes.bash "${mappedHarvestPath}" ${LUH_mappedHarvestPrefix} \
        ${selStartYear} ${selEndYearTransitions} "${calledFrom}" "${luhRelease} - scenario ${thisScenario}"
  fi

  #############################################################################
  #################################################### second step: "CMIP6 way"
  #-> scale states and transitions with given veg_ratio_max to better achieve LUH crop and pasture area

  #-- process states
  echo "########## scale states"
  if [ ! -d "${scenarioWorkPath}/${scaledStatesOutputFolder}" ]; then
    mkdir "${scenarioWorkPath}/${scaledStatesOutputFolder}"
  fi
  ${scriptsFolder}/${tools_dir}/scale_states_with_veg-ratio-max_v2.1-check-0.bash  "${scaledAndMappedStatesPath}" "${scenarioWorkPath}/${scaledStatesOutputFolder}" ${LUH_mappedStatesCF1FilePrefix} \
      ${LUH_scaledMappedStatesFilePrefix} ${vegRatioMaxFile} ${vrmVarName} ${selStartYear} ${selEndYearStates} "${calledFrom}" "${luhRelease} - scenario ${thisScenario}"

  #-- process transitions
  if [ ${processTransitions} = true ]; then
    # write scaling_of_transition.nml
cat > scaling_of_transition.nml <<EOF
&scaling_of_transition
  resolution=${ires}
  year_start=${selStartYear}
  year_end  =${selEndYearTransitions}     !Note: states required for the next year also!
  workingPath='${scenarioWorkPath}'
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
  luhRelease='${luhRelease} - scenario ${thisScenario}'
/
EOF

    echo "########## scale transitions"
    if [ ! -d "${scenarioWorkPath}/${scaledTransitionsOutputFolder}" ]; then
      mkdir "${scenarioWorkPath}/${scaledTransitionsOutputFolder}"
    fi
    ${scriptsFolder}/${bin_dir}/Adapting_transitions_after_scaling_of_states.x
  fi 

  (( scenIndex = scenIndex + 1 ))

done

echo "***Finished successfully***"
