#!/bin/bash
# ault (A100) counterpart of make_sc2026_runscript.sh.
#
# make_runscripts has no ault machine target (create_target_header ships
# daint_* and balfrin_gpu only), so the runscript comes out for the "default"
# target with no SLURM batch header at all. This injects an ault header
# directly. Unlike daint there is no uenv, but the launcher is swapped to srun
# to match the daint script.
#
# The launcher is srun, same as the daint script. The job inherits the
# submitter's environment via sbatch --export=ALL, so submit from a shell that
# has the nvhpc runtime, netcdf and the VT .so on LD_LIBRARY_PATH (see
# SC2026_HOWTO.ault.md).
#
# Usage:   ./run/make_sc2026_runscript.ault.sh     (from build/verification/)
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)
BUILD_VERIF=${SCRIPT_DIR%/run}
cd "$BUILD_VERIF"

./make_runscripts sc2026

RUN="run/exp.sc2026.run"

# --- SLURM header: no ault target emits one, so inject after the shebang.
#     Set SBATCH_ACCOUNT before running, or edit the account below. ---
sed -i "1a\\
#SBATCH --account=${SBATCH_ACCOUNT:-g34}\\
#SBATCH --partition=amda100\\
#SBATCH --nodelist=ault25\\
#SBATCH --gres=gpu:a100:1\\
#SBATCH --nodes=1\\
#SBATCH --time=04:00:00" "$RUN"

# --- single node ---
sed -i 's|^: \${no_of_nodes:=2}|: ${no_of_nodes:=1}|' "$RUN"

# --- launcher: default target records the nvhpc mpiexec in START; swap to srun
#     to match the daint script. ---
sed -i 's|^export START=.*mpiexec.*|export START="srun -n $mpi_total_procs --ntasks-per-node $mpi_procs_pernode --threads-per-core=1 --cpus-per-task $OMP_NUM_THREADS"|' "$RUN"

# --- ulimit -s unlimited right before ${START_MODEL} invocation ---
sed -i '/^\${START_MODEL}/ i\
ulimit -s unlimited' "$RUN"

# --- ATM_TIMESTEP guard, inserted before the EXPNAME line ---
sed -i '/^export EXPNAME="sc2026"$/ i\
: ${ATM_TIMESTEP:? "Error: ATM_TIMESTEP is not set. Please provide a value."}\
' "$RUN"

# --- Drop the hardcoded EXPNAME="sc2026" so the submitter-exported value survives ---
sed -i '/^export EXPNAME="sc2026"$/d' "$RUN"

echo "[patched ault] $RUN"
