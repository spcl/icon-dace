#!/bin/bash

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2024, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

#-----------------------------------------------------------------------------
# Master script to generate ICON-Land initial files
#
# It makes use of the existing scripts:
#  1) jsbach4_ini_files_from_gauss.sh
#     - Remapping of Gaussian grid JSBACH3 initial files to the ICON grid
#     - Assigning the data to the different jsbach4 ic and bc files
#
#  2) extpar4jsbach_mpim_icon.sh
#     - run extpar to generate soil texture data as well as albedo, roughness
#       length, forest fraction and LAI and vegetation fraction climatologies.
#
#  3) jsbach4_ini_files_from_extpar.ksh
#     Generate ic/bc files containing the new expar data:
#     - Glacier and lake masks replaced (-> bc_land_frac)
#     - Albedo, roughness length, forest fraction, lai_clim and veg_fract
#       replaced. Note: albedo_veg_vis, -veg_nir, -soil_vis and -soil_nir are
#       still interpolated from Gaussian grid
#     - Rooting depth (and maxmoist) replaced, additional variables:
#       FR_SAND, FR_SILT, FR_CLAY, FR_OC, SUB_FR_SAND, SUB_FR_SILT,
#       SUB_FR_CLAY and SUB_FR_OC
#
#  4) adapt_extpar_file.bash
#     The extpar data file read by the NWP atmosphere contains several
#     variables, that are also included in the bc_land files. For consistency,
#     these variables are replaced by the respective bc_land file variables.
#
# This approach is meant to be preliminary. It documents the current process
# of initial data generation. The aim is however, to generate all initial data
# from extpar in the not so far future.
#-----------------------------------------------------------------------------
set -e

scripts_dir=$(dirname $0)
cd $scripts_dir
scripts_dir=$(pwd)

# Run all sections of this script?
run_jsbach4_ini_files_from_gauss=true
run_extpar4jsbach_mpim_icon=true
run_jsbach4_ini_files_from_extpar=true
rm_bc_files_from_gauss=true
run_adapt_extpar_file=true

#
# Settings - also needed in sub-scripts
# -------------------------------------
# ICON grid used
export icon_grid_rootdir=/pool/data/ICON/grids/public/mpim
export atmGridID=0012
export refinement=R02B04
export icon_grid_file=icon_grid_${atmGridID}_${refinement}_G.nc
export icon_grid=${icon_grid_rootdir}/${atmGridID}/${icon_grid_file}

export coupled=false       # land sea mask for coupled experiment?
if [[ ${coupled} == true ]]; then
  export oceGridID=0035      # Required for coupled configurations
  export grid_label=$atmGridID-$oceGridID
else
  export grid_label=$atmGridID
fi
# Initial Extpar file (needed by jsbach4_ini_files_from_gauss.sh)
export initial_extpar_file=${icon_grid_rootdir}/${grid_label}/icon_extpar_0012_R02B04_G_20161124_tiles.nc
# Extpar file read by the NWP atmosphere at run time, only used in adapt_extpar_file below
export coupled_extpar_file=/pool/data/ICON/grids/public/mpim/0012-0035/icon_extpar_oceLSM_a0012_R02B04_o0035_R02B06_20161124_tiles.nc

if [[ ${coupled} == true ]]; then
  export fractional_mask=${icon_grid_rootdir}/${grid_label}/fractional_mask/fractional_lsm_${atmGridID}_${oceGridID}.nc
else
  export fractional_mask=${initial_extpar_file}
  # export fractional_mask=${icon_grid}          # Leads to integer LSM!
fi

# Definitions for the new extpar, ic and bc files
export revision=r00xx      # ic/bc revision directory that will be generated

export start_year=1850
export end_year=1850

# Directory for extpar source code and the extpar file (to be generated)
export extpar_dir=/work/mj0060/m220053/extpar/extpar4jsbach
today=$(date  +%Y%m%d)
export extpar_file=icon_extpar4jsbach_${grid_label}_${today}.nc

export work_dir=/scratch/m/$USER/ini_files/work/${grid_label}/${revision}  # temporary working directory
export output_root_dir=/work/mj0060/m220053/ini_files     # root directory for the new ic/bc files
export bc_file_dir=${output_root_dir}/${grid_label}/land/${revision}

# Shared parameters
export min_fract=0.000001 # minimum land grid cell fraction if not complete ocean
export max_fract=0.999999 # maximum land grid cell fraction if not complete land

# 1. Run jsbach4_ini_files_from_gauss
if [[ ${run_jsbach4_ini_files_from_gauss} == true ]]; then
  ./jsbach4_ini_files_from_gauss.sh
fi

# 2. Run extpar4jsbach
if [[ ${run_extpar4jsbach_mpim_icon} == true ]]; then

  . $MODULESHOME/init/bash
  # Clone and make extpar4jsbach
  [[ ! -d ${extpar_dir%/*} ]] && mkdir -p ${extpar_dir%/*}
  if [[ ! -d ${extpar_dir}/bin ]]; then
    cd ${extpar_dir%/*}
    # setup for levante (DKRZ)
    module load git/2.31.1-gcc-11.2.0
    module load gcc/11.2.0-gcc-11.2.0
    git clone --recursive git@gitlab.dkrz.de:m212005/extpar4jsbach.git ${extpar_dir}
    cd ${extpar_dir}
    git submodule update
    ./configure.levante.gcc
    source modules.env
    make -j 4
  else
    cd ${extpar_dir}
    source modules.env
  fi
  module load nco

  # Generate extpar data for Jsbach
  cd ${scripts_dir}
  ./extpar4jsbach_mpim_icon.sh
fi

# 3. Run jsbach4_ini_files_from_extpar
if [[ ${run_jsbach4_ini_files_from_extpar} == true ]]; then
  . $MODULESHOME/init/bash
  module load nco

  ./jsbach4_ini_files_from_extpar.ksh
fi

# 4. Move new bc/ic files to bc_file_dir and remove preliminary directories
if [[ ${rm_bc_files_from_gauss} == true ]]; then
  [[ -d ${bc_file_dir} ]] || mkdir ${bc_file_dir}
  mv ${bc_file_dir}_with_extpar/*.nc                ${bc_file_dir}
  rmdir ${bc_file_dir}_with_extpar
  mv ${bc_file_dir}_from_gauss/bc_land_sso_????.nc  ${bc_file_dir}
  mv ${bc_file_dir}_from_gauss/ic_land_soil_????.nc ${bc_file_dir}
  rm -fr ${bc_file_dir}_from_gauss
  echo ""
  echo "$(basename $0): Output moved to"
  echo "     ${bc_file_dir}"
fi

# 5. Run jsbach4_ini_files_from_extpar
if [[ ${run_adapt_extpar_file} == true ]]; then
  ./adapt_extpar_file.bash
fi

echo "====  Initial and boundary file generation completed. ===="
echo ""

exit 0
