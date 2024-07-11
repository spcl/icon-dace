#! /bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

set +eu

nvsmi_logger_PID=0
function kill_nvsmi()
{
    set +x
    if (( nvsmi_logger_PID != 0 ))
    then
        kill $nvsmi_logger_PID
    fi
}
trap kill_nvsmi ERR
trap kill_nvsmi EXIT

lrank=${OMPI_COMM_WORLD_LOCAL_RANK:-$SLURM_LOCALID}

# Local Task 0 runs always the nvidia-smi logger
# To enable logging of nvidia-smi by setting the following in your run script
# (either directly or via create_target_header)
#     export ENABLE_NVIDIA_SMI_LOGGER=yes
if [[ "$lrank" == 0 ]] && [[ ${ENABLE_NVIDIA_SMI_LOGGER:-"no"} == "yes" ]]
then
    set +x
    basedir=${basedir:=.}
    # Start logger in background. It will be killed by the ERR trap or kill_nvsmi.
    loop_repetition_time=500000000 # in nano seconds
    while sleep 0.$(( ( 1999999999 - 1$(date +%N) ) % loop_repetition_time ))
    do
        LC_TIME=en_US date -Ins
        nvidia-smi --format=csv --query-gpu=index,power.draw,utilization.gpu,temperature.gpu,memory.used
    done > nvsmi.log.${lrank} &
    nvsmi_logger_PID=$!
    set -x
fi

gpus=(0 1 2 3)

export CUDA_VISIBLE_DEVICES=${gpus[$((lrank % ${#gpus[@]} ))]}

export KMP_AFFINITY=scatter

echo "lrank: $lrank (${OMPI_COMM_WORLD_LOCAL_RANK}:-${SLURM_LOCALID}), gpu: $CUDA_VISIBLE_DEVICES"

"$@"
return=$?

kill_nvsmi
exit $return

