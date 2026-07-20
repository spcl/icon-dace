#!/bin/bash
# Submit all 6 SC2026 SNR variants (§8.1–8.6) for one grid.
#
# Usage:
#   [SERDE_GEN_END=<N>] [SERDE_GEN_STRIDE=<S>] [SERDE_GEN_LIST=<g1,g2,...>] \
#     ./sbatch_all_sc2026.sh <GRID>
#     GRID             : <gridID>_<refinement>, e.g. 0010_R02B04
#
#   The serde knobs below are env passthroughs to sbatch_sc2026.sh.
#     SERDE_GEN_END    : serialized-window size (physics_generation 0..N-1).
#                        Default 51.
#     SERDE_GEN_STRIDE : within the window, dump every S-th generation.
#                        Default 1 (every generation).
#     SERDE_GEN_LIST   : explicit generations to dump, comma/space separated,
#                        max 64 (e.g. logspace points over a long run). When
#                        set, it replaces the window/stride gate entirely.
set -euo pipefail

GRID=${1:? "Usage: $0 <GRID>  (e.g. 0010_R02B04)"}

SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)
SBATCH="$SCRIPT_DIR/sbatch_sc2026.sh"

# 8.1 vanilla FP64, ss5 (reference)
USE_VT_GPU=0 NDYN_SUBSTEPS_OVERRIDE=5  "$SBATCH" 8 fp64 "$GRID"
# 8.2 vanilla FP64, ss10 (temporally-refined reference)
USE_VT_GPU=0 NDYN_SUBSTEPS_OVERRIDE=10 "$SBATCH" 8 fp64 "$GRID"
# 8.3 VT FP64
USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  "$SBATCH" 8 fp64 "$GRID"
# 8.4 VT FP32
USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  "$SBATCH" 8 fp32 "$GRID"
# 8.5 VT FP16
USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  "$SBATCH" 8 fp16 "$GRID"
# 8.6 VT BF16
USE_VT_GPU=1 NDYN_SUBSTEPS_OVERRIDE=5  "$SBATCH" 8 bf16 "$GRID"
