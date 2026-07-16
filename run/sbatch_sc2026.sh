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
#
# Serialization is gated by (env, passed through to the job):
#   SERDE_GEN_START/END    : dump generations in [START, END). Default 0..51.
#   SERDE_GEN_STRIDE       : within that window, dump every S-th generation.
#                            Default 1.
#   SERDE_GEN_LIST         : explicit generations, comma/space separated, max 64.
#                            When set, replaces the window/stride gate entirely.
set -euo pipefail

export ATM_TIMESTEP=${1:? "Usage: $0 <ATM_TIMESTEP> [VT_PREC=fp64] [GRID=0010_R02B04]"}
export VT_PREC=${2:-fp64}
export GRID=${3:-0010_R02B04}

# vt_serde knob passthroughs (default values preserved if env unset)
export USE_VT_GPU=${USE_VT_GPU:-1}
export SERDE_GEN_START=${SERDE_GEN_START:-0}
export SERDE_GEN_END=${SERDE_GEN_END:-51}
export NDYN_SUBSTEPS_OVERRIDE=${NDYN_SUBSTEPS_OVERRIDE:-0}

# EXPNAME: tags the experiment dir so parallel runs don't overwrite each other
case "${USE_VT_GPU}" in
  0|F|f|false|FALSE|.false.|.FALSE.) VARIANT="vanilla" ;;
  *)                                  VARIANT="gpu${VT_PREC}" ;;
esac
# ss tag = effective dyn substeps. 0 means "no override" → namelist default (5).
SS="${NDYN_SUBSTEPS_OVERRIDE}"
[[ "$SS" == "0" ]] && SS=5
export EXPNAME="sc2026_dt${ATM_TIMESTEP}_ss${SS}_${VARIANT}_${GRID}"

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

sbatch --job-name="${EXPNAME}" \
       --output="${BUILD_VERIF}/run/LOG.SAVEME-SER-${VT_PREC^^}.${EXPNAME}.%j.o" \
       --error="${BUILD_VERIF}/run/LOG.SAVEME-SER-${VT_PREC^^}.${EXPNAME}.%j.o" \
       --export=ALL,ATM_TIMESTEP=${ATM_TIMESTEP},GRID=${GRID},LD_PRELOAD=${VT_SO} \
       "${BUILD_VERIF}/run/exp.sc2026.run"
