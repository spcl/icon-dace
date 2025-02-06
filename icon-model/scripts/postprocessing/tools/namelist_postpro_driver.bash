#!/bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#
#==========================================================================
#      Driver script for postprocessing and visualization of
#          ICON test case variables and namelist parameters
#==========================================================================
#
# History: 
# Original Version: Kristina Froehlich (DWD)
# Modified (Namelist parameters...) by Constantin Junk (MPI-M) (2010-01-21)
#
#==============================================================================
#
#
#==============================================================================

#==========================================================================
#                          USER'S SPECIFICATIONS 
#--------------------------------------------------------------------------
# 1. About the model output
#--------------------------------------------------------------------------
# 1.1 The directory in which the model output can be found. 
# Don't forget the trailing "/".

# for automatic testing
dir=`pwd -P`
icon_path=${dir%%scripts/postprocessing/tools}

# 1.2 The data file name is constructed in the same way as in the model.
# The name is composed of experiment name, hor. and vert. resolution,
# and a file index. The experiment name is the first part of the 
# file name of your model output. In this case: test_SBR=solid body rotation 
# test case

export EXP="hat_jww_echam_cld-conv-rad"

# 1.3 Spatial resolution
# (This script assumes that a string indicating the spatial resolution,
#  e.g., "_R2B04L31", has been added after ${EXP} as part of the name
#  of the model output.)

horizontal_resolution="R2B04"
vertical_resolution="L31"

resolution="${horizontal_resolution}${vertical_resolution}"

# 1.4 define grid optimization and configuration string which appears on the plot
gridopt="spr0.90"
export CONFIG_STRING="${EXP} ${resolution} ${gridopt}"


# 1.5 Define location of input grid-file
export GRIDFILE="${icon_path}grids/icon${horizontal_resolution}-grid_${gridopt}.nc"

#--------------------------------------------------------------------------
# 2. Decide what to do
#--------------------------------------------------------------------------

#decide which LEVEL and TIMESTEP should be plotted
export LEVEL=10
export TIMESTEP=20 

#Plot file format, pdf, ps...
export PFMT="pdf"

#Output filename
export PNAME="${EXP}_${resolution}_level-${LEVEL}_step-${TIMESTEP}_varplot"

#Directory of Namelist and experiment data
export DIRI="${icon_path}experiments/$EXP/"

# Where should the plot files be located? Don't forget the trailing "/".

#Plot file path
export DIRO="${DIRI}plots/"

#Input filename
export FNAM="${DIRI}${EXP}_${horizontal_resolution}${vertical_resolution}_0001.nc"


#--------------------------------------------------------------------------
# 3. specifiy namelist parameters
#--------------------------------------------------------------------------

# the user can specify a string which contains the namelist parameters 
# to be printed into the title section of the plot. The namelist parameters are stored 
# in the header file ${FNAM}. Use `ncdump -c ${FNAM}` (section "// global attributes") 
# to see the namelist parameters.

# Two examples how to specify a namelist string:

# namstr=("nroot" "n_dom" "run_ctl:i_cell_type" "dynamics_ctl:itime_scheme")

# namstr=("nroot" "n_dom" "run_ctl:i_cell_type" "dynamics_ctl:itime_scheme"\
#        "run_ctl:iequations" "io_ctl:lwrite_z3" "diffusion_ctl:hdiff_multfac"\
#        "io_ctl:lwrite_omega" "io_ctl:lwrite_precip" "io_ctl\:nsteps_diag"\
#        "run_ctl:dtime")

# if the namelist string is empty or missing, a default namelist parameters will be printed
# which contains dtime, iequation, ldynamics, iforcing 

#namstr=("run_ctl:dtime" "")

#--------------------------------------------------------------------------
#                    END OF USER'S SPECIFICATIONS 
#==========================================================================

#export DAY=`expr ${TIMESTEP} \* ${DTIME} / 86400`

echo
echo "**********************************************************"
echo "******************  ICON Postprocessing ******************"
echo "**********************************************************"
echo 
echo "=== Postprocessing started..."


#===================== Namelist ===============================

#check if namelist string was specified by user
if [ ! $namstr ]; then
   echo "==="
   echo "=== Warning: no namelist string specified. Use default namelist string."
   echo "==="
   namstr=("run_ctl:dtime" "run_ctl:iequations" "run_ctl:ldynamics" "run_ctl:iforcing")
fi

#get number namestring parameters specified by user
n_namstr=${#namstr[@]}
i=0

#write namelist string in ascii file ${DIRO}list.txt
while [ ${i} -lt ${n_namstr} ]; do
  echo ${namstr[${i}]} >> ${DIRO}list.txt
  let  i=${i}+1
done


#===================== Info ==================================
# Temporary variables

#export NCL_SCRIPT_DIR=`pwd`"/namelist_postpro_scripts/"

# The directories for intermediate data and plots will be created, if 
# not already there

if [ ! -d ${DIRO} ]; then
   mkdir -p ${DIRO} 
fi

# We have to call NCL from the directory where our own colormaps and resource 
# files are located, so that they can be loaded by the NCL scripts correctly. 

#cd ${NCL_SCRIPT_DIR}

# Make sure that NCL finds the color maps defined in this directory.

export NCARG_COLORMAP_PATH=$NCARG_ROOT/lib/ncarg/colormaps:./



#------------------------------------------------------------------------
# plot (see the ncl script for details)
#------------------------------------------------------------------------     

echo "=== start plotting."

ncl namelist_postpro.ncl

echo 
echo "=== Postprocessing finished."
echo "=== The plots can be found in "${DIRO}

rm ${DIRO}list.txt
unset namstr

exit
#