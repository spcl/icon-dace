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

# make_runscripts reads the target from run/set-up.info. configure records
# use_target='default', which emits no SLURM batch header at all; force the
# daint_gpu target so the header (account, node/queue directives) is generated.
sed -i "s/^use_target=.*/use_target='daint_gpu'/" run/set-up.info

./make_runscripts sc2026

RUN="run/exp.sc2026.run"

# --- SBATCH header: --nodes=1, --partition=normal, --time=04:00:00 ---
sed -i \
  -e 's|^#SBATCH --nodes=2$|#SBATCH --nodes=1\n#SBATCH --partition=normal\n#SBATCH --time=04:00:00|' \
  -e 's|^: \${no_of_nodes:=2}|: ${no_of_nodes:=1}|' \
  "$RUN"

# --- uenv: the job must load it itself, since sbatch cannot be called from
#     inside a uenv session (libslurm-uenv-mount rc=-3000). The daint_gpu
#     target already emits the account line to anchor after. ---
sed -i '/^#SBATCH --account=/ a\
#SBATCH --uenv=icon/25.2:v1@santis\
#SBATCH --view=default' "$RUN"

# --- launcher: the daint_gpu target defaults to mpiexec, which is not on PATH
#     under the uenv; use SLURM's srun. ---
sed -i 's|^export START="mpiexec -n \$mpi_total_procs"$|export START="srun -n $mpi_total_procs --ntasks-per-node $mpi_procs_pernode --threads-per-core=1 --cpus-per-task $OMP_NUM_THREADS"|' "$RUN"

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
