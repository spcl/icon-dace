#!/bin/bash
# Submit an SC2026 ICON integration run with libvelocity.so LD_PRELOAD'd.
#
# Usage:
#   ./sbatch_sc2026.sh <ATM_TIMESTEP> [VT_PREC=fp64] [GRID=0010_R02B04]
#
# ATM_TIMESTEP : atmosphere time step in seconds (e.g. 2, 4, 8)
# VT_PREC      : VT precision: fp16 | fp32 | fp64 (default fp64)
# GRID         : <gridID>_<refinement>, e.g. 0010_R02B04 (default),
#                0008_R02B05, 0002_R02B06, 0050_R02B03
#
# VT_DIR may be overridden via the environment. It must contain the
# libvelocity_gpu_stage8_solve_nh_integration_release.${VT_PREC}.so files.
set -euo pipefail

export ATM_TIMESTEP=${1:? "Usage: $0 <ATM_TIMESTEP> [VT_PREC=fp64] [GRID=0010_R02B04]"}
export VT_PREC=${2:-fp64}
export GRID=${3:-0010_R02B04}

VT_DIR="${VT_DIR:-/capstor/scratch/cscs/pmazumde/sc2026-ad-test/icon-vt-dace}"
VT_SO="${VT_DIR}/libvelocity_gpu_stage8_solve_nh_integration_release.${VT_PREC}.so"

if [[ ! -f "$VT_SO" ]]; then
    echo "Error: $VT_SO not found" >&2
    exit 1
fi

# Per-precision symlink dir so the linker's 'libvelocity.so' resolves to the
# chosen precision's .so.
LINK_DIR="${VT_DIR}/lib_${VT_PREC}"
mkdir -p "$LINK_DIR"
ln -sf "$VT_SO" "${LINK_DIR}/libvelocity.so"

# Locate build/verification/ regardless of invocation cwd
SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)
BUILD_VERIF=${SCRIPT_DIR%/run}

export EXPNAME="sc2026_dt${ATM_TIMESTEP}_${VT_PREC}_${GRID}"

sbatch --job-name="${EXPNAME}" \
       --output="${BUILD_VERIF}/run/LOG.SAVEME-SER-${VT_PREC^^}.${EXPNAME}.%j.o" \
       --error="${BUILD_VERIF}/run/LOG.SAVEME-SER-${VT_PREC^^}.${EXPNAME}.%j.o" \
       --export=ALL,ATM_TIMESTEP=${ATM_TIMESTEP},GRID=${GRID},LD_PRELOAD=${VT_SO} \
       "${BUILD_VERIF}/run/exp.sc2026.run"
