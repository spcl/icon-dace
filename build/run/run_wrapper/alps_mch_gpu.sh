#!/usr/local/bin/bash -l

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

export LOCAL_RANK=$SLURM_LOCALID
export GLOBAL_RANK=$SLURM_PROCID
# The first 4 ranks get assigned to GPUs
# Subsequent ranks get no GPU
export GPUS=(3 2 1 0)
export NUMA=(0 1 2 3 0 1 2 3)
export NUMA_NODE=${NUMA[$LOCAL_RANK]}

export CUDA_VISIBLE_DEVICES=${GPUS[$SLURM_LOCALID%8]}

sanitizer_wrapper="../../run/run_wrapper/sanitizer_wrapper.sh"

if [[ $# -lt 2 ]]; then
    # only ICON binary as argument
    numactl --cpunodebind=$NUMA_NODE --membind=$NUMA_NODE bash -c "$@"
elif [[ "$1" = "--compute-sanitizer" ]]; then
    module use $USER_ENV_ROOT/modules
    module load nvhpc
    shift
    numactl --cpunodebind=$NUMA_NODE --membind=$NUMA_NODE bash -c "${sanitizer_wrapper} $@"
else
    echo "Error: Unrecognised command line option: $1"
    exit 1
fi
