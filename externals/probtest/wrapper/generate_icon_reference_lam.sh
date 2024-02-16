#!/bin/bash
# script that generates tolerances and references for the LAM experiment family
# EXP      : the name of the experiment to be run (within the LAM family)
# BUILDER  : the name of the builder used by buildbot (needed to find/store references)

# check EXP and BUILDER env variables
if [[ -z "${EXP}" ]]; then
  echo "environment variable EXP unset. must contain the name of the experiment"
  exit 1
fi
if [[ -z "${BUILDER}" ]]; then
  echo "environment variable BUILDER unset. must contain the name of the buildbot builder"
  exit 1
fi

# probtest directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/.. >/dev/null 2>&1 && pwd )"

# set common variables
REFERENCE_DIR=$(pwd)/probtest-reference/${BUILDER}
ICON_DIR=$(pwd)
export PROBTEST_CONFIG=${EXP}.json

# initialize configuration file
python ${SCRIPT_DIR}/probtest.py init --codebase-install ${ICON_DIR} --experiment-name ${EXP} --file-id NetCDF "*atm_3d*" --reference ${REFERENCE_DIR} --config ${PROBTEST_CONFIG}

python ${SCRIPT_DIR}/probtest.py run-ensemble --submit-command "sbatch --wait --account=g110" || exit 1
python ${SCRIPT_DIR}/probtest.py stats --ensemble || exit 1
python ${SCRIPT_DIR}/probtest.py tolerance || exit 1
python ${SCRIPT_DIR}/probtest.py check --input-file-ref stats_ref.csv --input-file-cur stats_1.csv
python ${SCRIPT_DIR}/probtest.py check-plot --input-file-ref stats_ref.csv --input-file-cur stats_1.csv

echo "copying reference from stats_ref.csv to ${REFERENCE_DIR}/reference/${EXP}.csv"
if [[ ! -d ${REFERENCE_DIR}/reference ]]; then
  mkdir -p ${REFERENCE_DIR}/reference
fi
cp stats_ref.csv ${REFERENCE_DIR}/reference/${EXP}.csv
