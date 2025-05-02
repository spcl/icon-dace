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
#  0) generate_fractional_mask.sh
#     - Generation of a fractional land sea mask for coupled setups:
#       Remapping of the ocean mask to the atmosphere grid
#
#  1) extpar4jsbach_mpim_icon.sh
#     - run extpar to generate soil texture data as well as albedo, roughness
#       length, forest fraction and LAI and vegetation fraction climatologies.
#
#  2) jsbach4_ini_files_from_gauss.sh
#     - Remapping of Gaussian grid JSBACH3 initial files to the ICON grid
#     - Read additional soil parameters from 0.5 deg grid and do the remapping
#       to the ICON grid.
#     - Get glacier, land sea mask and orographic data based on extpar data
#     - Assigning the data to the different jsbach4 ic and bc files
#
#  3) jsbach4_ini_files_from_extpar.sh
#     Generate ic/bc files containing the new expar data:
#     - Lake mask replaced (-> bc_land_frac)
#     - Albedo, roughness length, forest fraction, lai_clim and veg_fract
#       replaced. Note: albedo_veg_vis, -veg_nir, -soil_vis and -soil_nir are
#       still interpolated from Gaussian grid
#     - Rooting depth (and maxmoist) replaced, additional variables:
#       FR_SAND, FR_SILT, FR_CLAY, SUB_FR_SAND, SUB_FR_SILT and
#       SUB_FR_CLAY
#
#  4) adapt_nwp_extpar_file.sh
#     The extpar data file read by the NWP atmosphere contains several
#     variables, that are also included in the bc_land files. For consistency,
#     these variables are replaced by the respective bc_land file variables.
#
# This approach is meant to be preliminary. It documents the current process
# of initial data generation. The aim is however, to generate all initial data
# from extpar in the not so far future.
#-----------------------------------------------------------------------------
#
# Notes on the extpar files
#
# In this script we use three different variables defining extpar files:
# - initial_extpar_file: An extpar file generally available from the pool,
#      that contains e.g. orographic parameters
# - nwp_extpar_file: The extpar file read by the NWP-atmosphere at runtime.
#      It will be adapted to the newly generated ic/bc files in
#      adapt_nwp_extpar_file.sh.
# - extpar_file: The extpar file we generate here using extpar4jsbach_mpim_icon.sh
#      (unless it is already available from /pool/data/JSBACH/icon/extpar4jsbach/).
#      This extpar file contains e.g. soil textures and other soil parameters.
#
# Generally, you should only need one or two different extpar files:
#   In setups with AES atmosphere
#      you should use the "extpar_file" also as "initial_extpar_file".
#      (An nwp_extpar_file is not used.)
#   In setups with NWP atmosphere
#      you should use the "nwp_extpar_file" as "initial_extpar_file", if
#      provided from DWD. If you do not have this file available, it is
#      possible to use the "extpar_file" instead - and to also use this file
#      as nwp_extpar_file.
#
#-----------------------------------------------------------------------------
set -e

# Which sections of this script should be run?
#  It is possible to run one task after the other, which might be useful in case of problems

run_extpar4jsbach_mpim_icon=false  # Generally available in /pool/data/JSBACH/icon/extpar4jsbach/
run_generate_fractional_mask=false # Only needed for coupled setups with new grid combination
run_jsbach4_ini_files_from_gauss=true
run_jsbach4_ini_files_from_extpar=true
rm_bc_files_from_gauss=true        # true: move ic/bc files into output directory
run_adapt_nwp_extpar_file=false    # Only needed for configurations with NWP atmosphere,
                                   # and only possible after "rm_bc_files_from_gauss"

#-----------------------------------------------------------------------------
#
# Settings - also needed in sub-scripts
# -------------------------------------
# Variables which typically need to be adapted to the target grid(s) and directory
#                   and file paths (from here up to "### End of settings ###" line:
# - work_dir_base and scratch_dir_base
# - output_root_dir
# - atmGridID and refinement
# - coupled
# - oceGridID (if coupled=true)
# - list of years
# - extpar_dir (needs to be writable by user if generating new extpar file)
# - extpar_source_dir (parent dir needs to be writeable by user if extpar source is to be compiled)
# - initial_extpar_file
# - nwp_extpar_file (for simulations with NWP atmosphere)
# - minimum/maximum fractions for grid cells that are not completely ocean/land (if coupled=false)
#
# --------------------------------------
# Base directories for constructing work and output directories (see below)
work_dir_base=/work/mh0287/${USER}    # levante
scratch_dir_base=/scratch/m/$USER     # levante
#scratch_dir_base=/work/mh0287/${USER} # breeze4

# --------------------------------------
# ICON grids used
icon_grid_rootdir=/pool/data/ICON/grids/public   # for seamless also check /pool/data/edzw-shadow/

export atmGridID=0049
#export atmGridID=0012 # seamless
export atmRes=R02B04
export icon_grid=${icon_grid_rootdir}/mpim/${atmGridID}/icon_grid_${atmGridID}_${atmRes}_G.nc

export coupled=false       # land sea mask for coupled experiment?
if [[ ${coupled} == true ]]; then
  export oceGridID=0035      # Required for coupled configurations
  export oceRes=R02B06
                             # We need the global ocean grid/mask file: *_G.nc
  export icon_grid_oce=${icon_grid_rootdir}/mpim/${oceGridID}/icon_mask_${oceGridID}_${oceRes}_G.nc
  export grid_label=$atmGridID-$oceGridID
else
  export grid_label=$atmGridID
fi

# --------------------------------------
# Output and working directories
export revision=r00xx      # ic/bc revision directory that will be generated
export work_dir=${scratch_dir_base}/ini_files/work/${grid_label}/${revision}  # temporary working directory
export output_root_dir=${work_dir_base}/ini_files     # root directory for the new ic/bc files
export bc_file_dir=${output_root_dir}/${grid_label}/land/${revision}

# --------------------------------------
# List of years:
# a) First and last year for a series of bc_land_frac files, e.g. for historical simulations
start_year=1850
end_year=1850
year_list=""
for ((yr=${start_year}; yr<=${end_year}; yr++)); do
  export year_list="${year_list} $yr"
done
#---
# b) Define list of years
#export year_list="1850 1979 1992 2005"

# --------------------------------------
# extpar4jsbach file ("extpar_file"; compare above "Notes on the extpar files")
#
today=$(date +%Y%m%d) # time stamp for generated extpar file; set to a fixed YYYYMMDD to re-use file from different date
if [[ ${run_extpar4jsbach_mpim_icon} == false ]]; then
  # Use an already existing file and don't run extpar
  # -----
  export extpar_dir=/pool/data/JSBACH/icon/extpar4jsbach/mpim           # Check in 'mpim' and 'dwd' directories
  export extpar_file=icon_extpar4jsbach_${atmGridID}_20240827_tiles.nc  #   for the grid ID and corresponding date tag
  # -----
  #export extpar_dir=/pool/data/JSBACH/icon/extpar4jsbach/dwd            # seamless
  #export extpar_file=icon_extpar4jsbach_${atmGridID}_20240801_tiles.nc  # seamless
  # -----
  #export extpar_dir=${output_root_dir}                                  # Use extpar file from previous run
  #export extpar_file=icon_extpar4jsbach_${atmGridID}_${today}_tiles.nc  #   with run_extpar4jsbach_mpim_icon=true
else
  # Run extpar and generate a new extpar4jsbach file
  export extpar_version=v5.14 # extpar version if compiling extpar from scratch
  # Directory for extpar source code
  export extpar_source_dir=/work/mh0287/m212005/extpar/extpar4jsbach # on levante or breeze4
  # Directory and file name for generated extpar file
  export extpar_dir=${output_root_dir}
  export extpar_file=icon_extpar4jsbach_${atmGridID}_${today}_tiles.nc

  # Directory for extpar input data
  # If not running on levante or breeze4, add a case for your machine here
  case $(hostname) in
    levante*|*.lvt.dkrz.de)
      export extpar_data_dir=/work/pd1167/extpar-input-data
      ;;
    breeze4)
      export extpar_data_dir=/work/mh0287/icon-preprocessing/extpar-input-data
      ;;
    *)
      echo ERRO: Unknown host: $(hostname)
      exit 1
      ;;
  esac
fi

# --------------------------------------
# NWP Extpar file:  Only used in setups with NWP atmosphere  - compare above "Notes on the extpar files"
#export nwp_extpar_file=${icon_grid_rootdir}/${grid_label}/icon_extpar_oceLSM_a0012_R02B04_o0035_R02B06_20161124_tiles.nc
export nwp_extpar_file=${extpar_dir}/${extpar_file}       # Use extpar file defind above (to be generated or existing)

# --------------------------------------
# Initial Extpar file:  Used in jsbach4_ini_files_from_gauss.sh  - compare above "Notes on the extpar files"
# ----- Setups with AES atmosphere
export initial_extpar_file=${extpar_dir}/${extpar_file}   # Use extpar file defind above (to be generated or existing)
# ----- Setups with NWP atmosphere
#export initial_extpar_file=${nwp_extpar_file}            # Use NWP extpar file defined above    # seamless

# --------------------------------------
# Fractional mask file
if [[ ${coupled} == true ]]; then
  if [[ ${run_generate_fractional_mask} == true ]]; then
    export fractional_mask=${output_root_dir}/${grid_label}/fractional_mask/fractional_lsm_${atmGridID}_${oceGridID}.nc
  else
    # Needs to exist
    export fractional_mask=${icon_grid_rootdir}/mpim/${grid_label}/fractional_mask/fractional_lsm_${atmGridID}_${oceGridID}.nc
  fi
else
  export fractional_mask=${initial_extpar_file}
  # export fractional_mask=${icon_grid}          # Leads to integer LSM!
fi

# --------------------------------------
# Shared parameters
#
# minimum/maximum fractions for grid cells that are not completely ocean/land
#   0.000001 / 0.999999: Currently used; seams rather suitable for coupled setups, where water conservation is an issue.
#   0.001 / 0.999: Used in coupled and uncoupled setups until August 2023
#   0.25 / 0.75: Corresponding to coupled setup with matching grids, e.g. R2B5/R2B6
#   0.5 /0.5: No fractional grid cells
if [[ ${coupled} == true ]]; then
  export min_fract=0.000001 # minimum land grid cell fraction if not complete ocean
  export max_fract=0.999999 # maximum land grid cell fraction if not complete land
else
  #--- AMIP setup corresponding to coupled setup with matching grids, e.g. R2B5/R2B6
  export min_fract=0.25     # minimum land grid cell fraction if not complete ocean
  export max_fract=0.75     # maximum land grid cell fraction if not complete land
  #---  AMIP setup with integer land sea mask
  #export min_fract=0.5      # minimum land grid cell fraction if not complete ocean
  #export max_fract=0.5      # maximum land grid cell fraction if not complete land
fi

### End of settings ###

#-----------------------------------------------------------------------------
# Some preparations
scripts_dir=$(dirname $0)
cd $scripts_dir
scripts_dir=$(pwd)

function finish {
  if [[ -d ${work_dir} ]]; then
    echo "Cleaning up work dir ${work_dir}"
    rm -rf ${work_dir} >& /dev/null
  fi
}
trap finish EXIT

if [[ ! -d ${output_root_dir} ]]; then
  echo "Directory 'output_root_dir' for output needs to exist (now: ${output_root_dir}). "
  exit 1
fi
if [[ ! -w ${output_root_dir} ]]; then
  echo "No write permission for output directory ${output_root_dir}."
  exit 1
fi

# Create directory to save all scripts used for this initial file revision
[[ -d ${bc_file_dir}/scripts ]] || mkdir -p ${bc_file_dir}/scripts


#-----------------------------------------------------------------------------
# 0. Generate fractional land sea mask
if [[ ${run_generate_fractional_mask} == true ]]; then
  if [[ ${coupled} != true ]]; then
    echo " Skipping the generation of fractional masks as this is only needed in coupled setups."
  else
    ./generate_fractional_mask.sh
    cp ./generate_fractional_mask.sh ${bc_file_dir}/scripts
  fi
fi

#-----------------------------------------------------------------------------
# 1. Run extpar4jsbach
if [[ ${run_extpar4jsbach_mpim_icon} == true ]]; then

  [[ ! -d ${extpar_dir} ]] && [[ ${extpar_dir} != ${output_root_dir} ]] && mkdir -p ${extpar_dir}

  if [[ ! -d ${extpar_data_dir} ]]; then
    echo "Directory with extpar input data doesn't exist:"
    echo "    ${extpar_data_dir}"
    echo "Preferably, since this is a very large data volume, this should be a clone of"
    echo "    https://gitlab.dkrz.de/extpar-data/extpar-input-data"
    echo "already existing in your file system. If you need to download this yourself"
    echo "contact jonas.jucker@c2sm.ethz.ch for access to the git lsf repository."
    exit 1
  fi

  . $MODULESHOME/init/bash

  # Clone extpar source from repository if necessary
  if [[ ! -d ${extpar_source_dir} ]]; then
    module purge
    module load git
    [[ ! -d ${extpar_source_dir%/*} ]] && mkdir -p ${extpar_source_dir%/*}
    cd ${extpar_source_dir%/*}
    if ! git clone --recursive git@github.com:C2SM-RCM/extpar.git ${extpar_source_dir##*/} >&/dev/null ;then
      if ! git clone --recursive git@gitlab.dkrz.de:m212005/extpar4jsbach.git ${extpar_source_dir##*/} >&/dev/null ;then
        echo "No valid git repository to clone extpar"
        echo "Make sure that you have password-less SSH read access to the extpar repository at"
        echo "   https://github.com/C2SM-RCM/extpar"
        echo "To get access, contact jonas.jucker@c2sm.ethz.ch"
        echo "Or check if there is already a compiled version of extpar on your computer by"
        echo "someone else and set variable 'extpar_source_dir' to point to it."
        exit 1
      else
        echo "Cloned extpar from git@gitlab.dkrz.de:m212005/extpar4jsbach.git"
      fi
    else
      echo "Cloned extpar from git@github.com:C2SM-RCM/extpar.git"
    fi
    cd ${extpar_source_dir}
    git checkout ${extpar_version}
    git submodule update
  fi

  # Configure and compile/make extpar if necessary assuming that the existence of the
  # modules.env file indicates a complete installation of extpar with the binaries
  # and python scripts already in the bin directory
  if [[ ! -f ${extpar_source_dir}/modules.env ]]; then
    if [[ ! -w ${extpar_source_dir} ]]; then
      echo "No write permission for extpar source directory."
      exit 1
    fi
    cd ${extpar_source_dir}
    # If not running on levante or breeze4, add a case for your machine here
    case $(hostname) in
      levante*|*.lvt.dkrz.de)
        ./configure.levante.gcc
        ;;
      breeze4)
        ./configure.breeze4.gcc
        ;;
      *)
        echo ERRO: Unknown host: $(hostname)
        exit 1
        ;;
    esac

    source modules.env
    make -j 4
  else
    source ${extpar_source_dir}/modules.env
  fi

  # Generate extpar data for Jsbach
  cd ${scripts_dir}
  module load nco
  ./extpar4jsbach_mpim_icon.sh
  cp ./extpar4jsbach_mpim_icon.sh ${bc_file_dir}/scripts
fi

#-----------------------------------------------------------------------------
# 2. Run jsbach4_ini_files_from_gauss
if [[ ${run_jsbach4_ini_files_from_gauss} == true ]]; then
  ./jsbach4_ini_files_from_gauss.sh
  cp ./jsbach4_ini_files_from_gauss.sh ${bc_file_dir}/scripts
fi

#-----------------------------------------------------------------------------
# 3. Run jsbach4_ini_files_from_extpar
if [[ ${run_jsbach4_ini_files_from_extpar} == true ]]; then
  . $MODULESHOME/init/bash
  module load nco

  if [[ ! -f ${extpar_dir}/${extpar_file} ]]; then
    echo "${extpar_dir}/${extpar_file} does not exist."
    exit 1
  fi

  ./jsbach4_ini_files_from_extpar.sh
  cp ./jsbach4_ini_files_from_extpar.sh ${bc_file_dir}/scripts
fi

#-----------------------------------------------------------------------------
# 4. Move new bc/ic files to bc_file_dir and remove preliminary directories
if [[ ${rm_bc_files_from_gauss} == true ]]; then
  [[ -d ${bc_file_dir} ]] || mkdir ${bc_file_dir}
  mv ${bc_file_dir}_with_extpar/*.nc                ${bc_file_dir}
  rmdir ${bc_file_dir}_with_extpar
  mv ${bc_file_dir}_from_gauss/bc_land_sso.nc  ${bc_file_dir}
  mv ${bc_file_dir}_from_gauss/ic_land_soil.nc ${bc_file_dir}
  rm -fr ${bc_file_dir}_from_gauss
  echo ""
  echo "$(basename $0): Output moved to"
  echo "     ${bc_file_dir}"
fi

#-----------------------------------------------------------------------------
# 5. Run adapt_nwp_extpar_file
if [[ ${run_adapt_nwp_extpar_file} == true ]]; then
  ./adapt_nwp_extpar_file.sh
  cp ./adapt_nwp_extpar_file.sh ${bc_file_dir}/scripts
fi

cp $0 ${bc_file_dir}/scripts

echo "====  Initial and boundary file generation completed. ===="
echo ""

#-----------------------------------------------------------------------------
exit 0
