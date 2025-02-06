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

#-----------------------------------------------------------------------------
# submit as a serial job to mistral -> adapt wall_clock_limit to your needs!
# ==> NOTE: unfortunately wget seems not to be installed on the nodes - i.e. execute interactively...
# ===> hmmm, wget does work when reserving a node for interactive work ... 
#       ===> salloc --partition=compute --nodes=1 --time=03:00:00 --account bm0891
# -------> THUS DO IT ON THE LOGIN NODE :/
#----------------------------------------------------------------------------
#
###############################################################################
### Batch Queuing System is SLURM
#SBATCH --job-name=get_raw_data_LUH2_scenarios
#SBATCH --output=get_raw_data_LUH2_scenarios.o%j
#SBATCH --error=get_raw_data_LUH2_scenarios.o%j
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=5120 
#SBATCH --ntasks=1
#SBATCH --account mj0060
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --mail-user=julia.nabel@mpimet.mpg.de # Set your e-mail address

#----------------------------------------------------------------------------

#example: http://gsweb1vh2.umd.edu/LUH2/LUH2_v2f/IMAGE/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-f_gn_2015-2100.nc

urlprefix="http://gsweb1vh2.umd.edu/LUH2/LUH2_v2f/"

scenSource=( "IMAGE/" "GCAM34/" "MESSAGE/" "GCAM60/" "AIM/" "MAGPIE/" )
scenPart=( "IMAGE-ssp126" "GCAM-ssp434" "MESSAGE-ssp245" "GCAM-ssp460" "AIM-ssp370" "MAGPIE-ssp585" )

filesPrefix="multiple-" 
filekinds=( "states" "transitions" "management" )
filesCommon="_input4MIPs_landState_ScenarioMIP_UofMD-"
filesSuffix="-2-1-f_gn_2015-2100.nc"


wget --output-document LUH2_v2f_README.pdf http://gsweb1vh2.umd.edu/LUH2/LUH2_v2f_README.pdf

index=0
while [ ${index} -lt ${#scenSource[@]} ]; do 
    for kind in ${filekinds[@]}; do
        thisSource=${scenSource[${index}]}
        thisScenPart=${scenPart[${index}]}
        thisFile=${filesPrefix}${kind}${filesCommon}${thisScenPart}${filesSuffix}
        toGet=${urlprefix}${thisSource}${thisFile}
        echo ${toGet}
        wget --output-document ${thisFile} ${toGet}
    done
    (( index = index + 1 ))
done



