#!/bin/bash
# ault (A100) counterpart of make_sc2026_runscript.sh.
#
# make_runscripts has no ault machine target (create_target_header ships
# daint_*, balfrin_gpu, lumi_*, etc. but nothing for ault), so the runscript
# comes out for the "default" target: no SLURM header, CPU cache-blocking
# (nproma=48, nblocks_c=0), and multiple MPI ranks. This rewrites those to the
# single-GPU values the daint_gpu target uses (nproma=0, nblocks_c=1, one rank)
# and injects an ault SLURM header.
#
# The grid directory is taken from $grids_folder at run time (it must hold
# icon_grid_<id>_<ref>_G.nc); the same variable also drives the grid symlink.
# The job inherits the submitter's environment via sbatch --export=ALL, so submit
# from a shell that has the nvhpc runtime and netcdf on LD_LIBRARY_PATH and
# grids_folder / VT_DIR set (see SC2026_HOWTO.ault.md).
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
#SBATCH --ntasks=1\\
#SBATCH --time=04:00:00" "$RUN"

# --- single MPI rank (default target sets 4) ---
sed -i 's|mpi_procs_pernode:=[0-9]\+|mpi_procs_pernode:=1|' "$RUN"
sed -i 's|^: \${no_of_nodes:=[0-9]\+}|: ${no_of_nodes:=1}|' "$RUN"

# --- GPU blocking: whole cell dimension in one block, as the daint_gpu target
#     does. The VT no_nproma library requires nblk_c == 1. ---
sed -i 's|^nproma=[0-9]\+|nproma=0|'       "$RUN"
sed -i 's|^nproma_sub=[0-9]\+|nproma_sub=800|' "$RUN"
sed -i 's|^nblocks_c=[0-9]\+|nblocks_c=1|' "$RUN"

# --- launcher: the binary is MPI-enabled, so it needs an mpirun launcher even
#     for one rank; a direct exec hangs in MPI_Init (no PMIx bootstrap). ---
sed -i 's|^export START=.*|export START="mpirun -n 1"|' "$RUN"

# --- helper include: the template sources add_required_* relative to ${thisdir},
#     but that helper lives in run/; point the include there. ---
sed -i 's|\${thisdir}/add_required_atmo_non-hydrostatic_files|${thisdir}/run/add_required_atmo_non-hydrostatic_files|' "$RUN"

# --- grid directory: the experiment template hardcodes an absolute path; take
#     it from $grids_folder instead (same variable used for the grid symlink). ---
sed -i 's|^atmo_grid_folder=.*|atmo_grid_folder=${grids_folder:?set grids_folder to the directory holding icon_grid_<id>_<ref>_G.nc}|' "$RUN"

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
