#!/bin/bash
# Regenerate run/exp.sc2026.run (via make_runscripts) and apply the SBATCH
# patches that ICON's default template generator doesn't emit for our daint_gpu
# target.
#
# Usage:   ./run/make_sc2026_runscript.sh          (from build/verification/)
#   or:    bash make_sc2026_runscript.sh           (from build/verification/run/)
set -euo pipefail

# Locate build/verification/ regardless of where this was invoked from
SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)
BUILD_VERIF=${SCRIPT_DIR%/run}
cd "$BUILD_VERIF"

./make_runscripts sc2026

RUN="run/exp.sc2026.run"

# --- SBATCH header: --nodes=1, --partition=normal, --time=04:00:00 ---
sed -i \
  -e 's|^#SBATCH --nodes=2$|#SBATCH --nodes=1\n#SBATCH --partition=normal\n#SBATCH --time=04:00:00|' \
  -e 's|^: \${no_of_nodes:=2}|: ${no_of_nodes:=1}|' \
  "$RUN"

# --- uenv: the job must load it itself, since sbatch cannot be called from
#     inside a uenv session (libslurm-uenv-mount rc=-3000) ---
sed -i '/^#SBATCH --account=/ a\
#SBATCH --uenv=icon/25.2:v1@santis\
#SBATCH --view=default' "$RUN"

# --- ulimit -s unlimited right before ${START_MODEL} invocation ---
sed -i '/^\${START_MODEL}/ i\
ulimit -s unlimited' "$RUN"

# --- ATM_TIMESTEP guard, inserted before the EXPNAME line ---
sed -i '/^export EXPNAME="sc2026"$/ i\
: ${ATM_TIMESTEP:? "Error: ATM_TIMESTEP is not set. Please provide a value."}\
' "$RUN"

# --- Drop the hardcoded EXPNAME="sc2026" so the submitter-exported value survives ---
sed -i '/^export EXPNAME="sc2026"$/d' "$RUN"

echo "[patched] $RUN"
