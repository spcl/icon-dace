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
#      Driver script for postprocessing and visulization of the 
#           solid body rotation test case 
#==========================================================================
#
# History: 
# Initial version by Daniel Rainert
# Modified version by Constantin Junk (MPI-M) (2010-10-18)
#
#==============================================================================
#
#
#==============================================================================

#==========================================================================
#                          USER'S SPECIFICATIONS 
#--------------------------------------------------------------------------
set -x
if [ "x$1" != "x" ]
then
  set_env=$1 
else
  set_env=/null
fi

# 1. About the model output
#--------------------------------------------------------------------------
# 1.1 The directory in which the model output can be found. 
# Don't forget the trailing "/".

# for automatic testing
dir=`pwd -P`
icon_path=${dir%%scripts/postprocessing/testcases}
icon_path=${basedir}/

# 1.2 The data file name is constructed in the same way as in the model.
# The name is composed of experiment name, hor. and vert. resolution,
# and a file index. The experiment name is the first part of the 
# file name of your model output. In this case: SBR=solid body rotation

EXP="test_SBR"


# 1.3 Shape of control volume (3 = triangle, 6 = hexagon/pentagon) 

export cell_type=3

# 1.4 Spatial resolution
# (This script assumes that a string indicating the spatial resolution,
#  e.g., "_R2B04L31", has been added after ${EXP} as part of the name
#  of the model output.)

horizontal_resolution="R2B05"

VERTICAL_RES=60  #needed for NCL script
vertical_resolution="L${VERTICAL_RES}"

# 1.5 The experiment identifier that will appear in the plots, e.g.,
# "SBR R2B05L31 spr0.90" and model name

CONFIG_STRING="SBR tri ${horizontal_resolution}${vertical_resolution} spr0.90"
MODEL="ICOHAM"

# 1.6 Define location of input grid-file
export GRIDFILE="${icon_path}grids/icon${horizontal_resolution}-grid_spr0.90.nc"

#--------------------------------------------------------------------------
# 2. Decide what to do
#--------------------------------------------------------------------------

#Plot file format
export PFMT="pdf"
#Output filename
export PNAME="SBR_error"

#The ncl scripts can read/plot the velocity field, and the tracers Q4 
#Set if we want to plot the velocity (0=not, 1=yes)
ivel=1

#--------------------------------------------------------------------------
#                    END OF USER'S SPECIFICATIONS 
#==========================================================================

echo
echo "**********************************************************"
echo "*  ICON Postprocessing for Solid Body Rotation test case *"
echo "**********************************************************"
echo 
echo "=== Postprocessing started..."

if [ -f ${set_env} ] 
then 
  echo " "
  echo " !!!!! Use setting from ./${set_env}"
  echo " "
  source ./${set_env}
fi


#===================== Info ==================================
# Temporary variables

export NCL_SCRIPT_DIR=`pwd`"/SBR_postpro_scripts/"


# We have to call NCL from the directory where our own colormaps and resource 
# files are located, so that they can be loaded by the NCL scripts correctly. 

cd ${NCL_SCRIPT_DIR} 

# Make sure that NCL finds the color maps defined in this directory.

export NCARG_COLORMAP_PATH=$NCARG_ROOT/lib/ncarg/colormaps:./

# For the NCL plotting scripts

export MODEL

# Where should the plot files be located? Don't forget the trailing "/".
export DIRI="${icon_path}experiments/$EXP/"
#Plot file path
export DIRO="${DIRI}plots/"
# The directories for intermediate data and plots will be created, if 
# not already there
if [ ! -d ${DIRO} ]; then
   mkdir -p ${DIRO} 
fi
export VERTICAL_RES
vertical_resolution="L${VERTICAL_RES}"
#Input filename
export CONFIG_STRING="SBR tri ${horizontal_resolution}${vertical_resolution} spr0.90"
export RESOLUTION="${horizontal_resolution}${vertical_resolution}"
export ExpName=${ExpName}
export FNAM="${ExpName}_DOM01_ML_0001.nc"
# 1.6 Define location of input grid-file
export GRIDFILE="${DIRI}/${grid_name}.nc"
export PFMT
export PNAME=${EXP}_${PNAME}
export ivel

#------------------------------------------------------------------------
# plot (see the ncl script for details)
#------------------------------------------------------------------------     

echo "=== start plotting."

ncl solid_body_init.ncl

#------------------------------------------------------------------------
# Finish...
#------------------------------------------------------------------------

echo 
echo "=== Postprocessing finished for the SBR test case."
echo "=== The plots can be found in "${DIRO}

unset DIRI
unset DIRO
unset FNAM
unset GRIDFILE
unset PFMT
unset PNAME
unset CONFIG_STRING
unset RESOLUTION
unset MODEL
unset VERTICAL_RES
