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
#==========================================================================================
#      Driver script which calls ncl-scripts that calculate the PS l1-,l2- and l_inf- error,
#                     (norm_err_std.ncl) each day and plot the time evolution 
#                      of the PS standardized error norms (plot_err_evol.ncl).
#                      
#==========================================================================================
#
# History: 
# Pilar Ripodas, DWD
# Modified by Constantin Junk, MPI-M, 2010-12-07
#
#==============================================================================
#
# packages needed:
#    - cdo
#    - ncl
#
#==============================================================================
#==============================================================================
#
#
#==============================================================================

check_error()
{

# Check if the first parameter (return status) is not OK

  if [ $1 -ne 0 ]
  then

# Stop running this script and return the error status
    echo "ERROR: $2"
    exit $1
  fi
}

#==========================================================================
#                          USER'S SPECIFICATIONS 
#--------------------------------------------------------------------------
#
# Check if a parameter is given

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
set -ex

# for automatic testing
dir=`pwd -P`
icon_path=${dir%%scripts/postprocessing/testcases}
icon_path=${basedir}/

# 1.2 The data file name is constructed in the same way as in the model.
# The name is composed of experiment name, hor. and vert. resolution,
# and a file index. The experiment name is the first part of the 
# file name of your model output. In this case: test_SBR=solid body rotation 
# test case

EXPN="test_hat_jww"

# 1.3 Spatial resolution
# Specifiy the root division 'Root', the number of bisections 'Bis' and the 
# vertical resolution 'ver_res'

Root=2
Bis=4
ver_res="L31"


# 1.4 define grid optimization and configuration string which appears on the plot
gridopt="spr0.90_M4"

if [ -f ${set_env} ] 
then 
  echo " "
  echo " !!!!! Use setting from ./${set_env}"
  echo " "
  source ./${set_env}
fi


export EXPN
export Root
export Bis
export ver_res

export hor_res="R${Root}B0${Bis}"
export CONFIG_STRING="${EXPN} ${hor_res}${ver_res} ${gridopt}"

# 1.5 Define location of input grid-file

export DIRI="${icon_path}experiments/$EXPN/"
export GridFileN="${DIRI}/${grid_name}.nc"

#--------------------------------------------------------------------------
# 2. Decide what to do
#--------------------------------------------------------------------------

#plot variable 'VarN' with dimension 'VarD' (momentarily only VarN="pres_sfc" and VarD=2 possible)
VarN=${VarN:-pres_sfc}
export VarN
export VarD=2

#Plot file format, pdf, ps...
export PFMT="eps"

#Output filename
export PNAME="${EXPN}_${hor_res}${var_res}_PS_error_norms"

#Directory of Namelist and experiment data, assumed that input file looks like
# ${DIRI}${EXPN}_${hor_res}${ver_res}_000*.nc

# Where should the plot files be located? Don't forget the trailing "/".
#Plot file path
export DIRO="${DIRI}plots/"

#ncl script dir
export NCL_SCRIPT_DIR=`pwd`"/error_postpro_scripts/"

#--------------------------------------------------------------------------
#                    END OF USER'S SPECIFICATIONS 
#==========================================================================

echo
echo "**********************************************************"
echo "******************  ICON Postprocessing ******************"
echo "**********************************************************"
echo 
echo "=== Postprocessing started..."

export ErrorFile="${DIRI}PS_error_${EXPN}_${hor_res}.txt"
export IconFileN="${DIRI}${VarN}.nc"

#get the data

cdo -s selname,${VarN} ${DIRI}${ExpName}_0001.nc ${IconFileN}
check_error $? "Calling of cdo" 

#   for file in `ls ${DIRI}${EXPN}_${hor_res}${ver_res}_000[2-9].nc`
#   do
#    cdo selname,${VarN} $file varout.nc
#    cdo cat ${IconFileN} ${DIRI}varout.nc ${DIRI}${VarN}_i.nc
#    rm ${IconFileN}
#    mv ${DIRI}${VarN}_i.nc ${IconFileN}
#    rm ${DIRI}varout.nc
#   done

# The directories for intermediate data and plots will be created, if 
# not already there

if [ ! -d ${DIRO} ]; then
   mkdir -p ${DIRO} 
fi

# We have to call NCL from the directory where our own colormaps and resource 
# files are located, so that they can be loaded by the NCL scripts correctly. 

cd ${NCL_SCRIPT_DIR}

# Make sure that NCL finds the color maps defined in this directory.

export NCARG_COLORMAP_PATH=$NCARG_ROOT/lib/ncarg/colormaps:./

#------------------------------------------------------------------------
# plot (see the ncl script for details)
#------------------------------------------------------------------------     

echo " "
echo "=== CALCULATE error norms of surface pressure."
echo " "
ncl norm_err_std.ncl
check_error $? "calling of ncl norm_err_std.ncl"

echo " "
echo "=== PLOT error norms of surface pressure."
echo " "
ncl plot_err_evol.ncl
check_error $? "calling of ncl ncl plot_err_evol.ncl"

echo "  "
echo "=== Postprocessing finished."
echo "=== The plots can be found in "${DIRO}

#rm ${DIRI}${VarN}.nc
#rm ${ErrorFile}

exit 0
#
