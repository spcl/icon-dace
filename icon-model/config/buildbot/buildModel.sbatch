#!/bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#=============================================================================
#SBATCH --job-name=build-icon
#SBATCH --nodes=1

#set -x
# abort on error
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)

if test "${SLURM_SUBMIT_DIR:-0}" != '0'; then
  ICON_DIR="${SLURM_SUBMIT_DIR}"
else
  ICON_DIR=$(cd "${SCRIPT_DIR}/../.."; pwd)
fi
echo "SCRIPT_DIR:        |${SCRIPT_DIR}|"
echo "ICON_DIR:          |${ICON_DIR}|"
echo "BB_NAME:           |${BB_NAME}|"
echo "BB_SYSTEM:         |${BB_SYSTEM}|"
echo "BB_experimentList: |${BB_EXPERIMENTLIST}|"
echo "BB_experiment:     |${BB_EXPERIMENT}|"


defaultExperimentList='merge2rc'

experimentList=${experimentList:-${BB_EXPERIMENTLIST:-${defaultExperimentList}}}
experimentLogfile="${ICON_DIR}/run/buildbotsExperiments.log"

echo "experimentList:    |${experimentList}|"
echo "experimentLogfile: |${experimentLogfile}|"


# plan: simplify the buildbot server config and take all decisions within the
# icon repo in a transparent way

# In case an experiment is given and the list is `tolerance`, let's generate new tolerance data for that experiment! {{{
if [[ -n "${BB_EXPERIMENT-}" ]] && [[ "${experimentList}" = "tolerance" ]]; then
  cd "${ICON_DIR}/scripts/buildbot_scripts";
  ./create_list_tolerance "${BB_EXPERIMENT}"
# }}}

# In case an experiment is given and the list is `select-members`, select the members for that experiment! {{{
elif [[ -n "${BB_EXPERIMENT-}" ]] && [[ "${experimentList}" = "select-members" ]]; then
  cd "${ICON_DIR}/scripts/buildbot_scripts";
  ./create_list_select-members "${BB_EXPERIMENT}"
# }}}

# In case only an experiment is given {{{
elif [ -n "${BB_EXPERIMENT}" ]; then

  # first check if the given expriment file/template exists
  fullExperimentPath="${ICON_DIR}/run/${BB_EXPERIMENT}"
  if [ ! -f "${fullExperimentPath}" ]; then
    echo "ERROR: Cannot find experiment ${fullExperimentPath}"
    exit 1;
  else
    echo "Found single tests: ${fullExperimentPath}"
  fi

  # create a temporary list with the given experiments
  experimentList="buildbot_explist.tmp"
  # creata a new list with the given experiment
  cd ${ICON_DIR}/scripts/buildbot_scripts; ./create_all_builders "${experimentList}";
  cd ${ICON_DIR} ; ./scripts/buildbot_scripts/addexp --list "${experimentList}" "${BB_EXPERIMENT}"
# }}}
else # use the given/default experiment list  {{{

  # check if the given list is part of the repository. if not, then run the corresponding script to create the list
  cd ${ICON_DIR}/scripts/buildbot_scripts;
  if test ! -f experiment_lists/${experimentList}; then ./create_list_${experimentList} ; fi
fi #}}}

# list the expriments for this builder and list| this will be used by buildbot later on
cd ${ICON_DIR}; ./scripts/buildbot_scripts/lsexperiments --list ${experimentList} ${BB_NAME} > ${experimentLogfile}

cat ${experimentLogfile}

# build the model
cd ${ICON_DIR}; ./scripts/buildbot_scripts/build --list ${experimentList} ${BB_NAME}
