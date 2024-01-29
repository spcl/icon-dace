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

lrank=$OMPI_COMM_WORLD_RANK
compute_tasks=$OMPI_COMM_WORLD_SIZE

echo compute_tasks: $compute_tasks, lrank: $lrank

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

# Put GPU 0 last to use it primarily for IO-Procs as GPU0 is some times used by other users.
gpus=(1 2 3 4 5 6 7 0)

# split the blocks for later use of gpnl nodes as I/O server
if (( lrank < compute_tasks ))
then
    echo Compute process $OMPI_COMM_WORLD_RANK on $(hostname)
    export CUDA_VISIBLE_DEVICES=${gpus[$((lrank % ${#gpus[@]} ))]}
else
    echo IO process $OMPI_COMM_WORLD_RANK on $(hostname)
fi

export KMP_AFFINITY=scatter

$@
return=$?

echo "nvsmi_logger_PID $nvsmi_logger_PID"
kill_nvsmi
exit $return

