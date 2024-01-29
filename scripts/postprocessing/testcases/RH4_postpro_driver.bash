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
#      Driver script for postprocessing and visulization for the
#                Rossby-Haurwitz (wavenumber 4) test
#==========================================================================
#
# History:
# Initial version by Hui Wan (MPI-M, 2009-05)
#
# Short description:
# This script interpolates the orginal ICOHDC output to selected 
# pressure level, then makes contour plots. 
#
# Software needed:
# - CDO (Climate Data Operators, www.mpimet.mpg.de/cdo)
#   for interplation and for the spectral transform;
# - NCL (NCAR Command Language, www.ncl.ucar.edu)
#   for visualization.
#
#==============================================================================
# Functions

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
if [ "x$1" != "x" ]
then
  set_env=$1 
else
  set_env=/null
fi

#--------------------------------------------------------------------------
# 1. About the model output
#--------------------------------------------------------------------------
# 1.1 The directory in which the model output can be found.
# Don't forget the trailing "/".

# for automatic testing
dir=`pwd -P`
icon_path=${dir%%scripts/postprocessing/testcases}
icon_path=${basedir}/

# Output frequency expressed as number of time slices per day.
# (The value you give here must be ASCII representations of integer values,
# such as "1" for daily output, and "4" for 6-hourly data.)

export output_frequency=1

# 1.2 The data file name is constructed in the same way as in the model.
# The name is composed of experiment name, hor. and vert. resolution (cf. 1.5),
# and a file index (cf. 1.6), if the output has been split on several files.
# The experiment name is the first part of the file name of your model output.
# You may have used something like "A0001" just for simplicity.

EXP="Rossby_Haurwitz"
#DOMAIN=""

# 1.3 The experiment identifier that will appear in the plots, e.g.,
# "brk R2B5L31"

config_string="RH4 tri R2B04L31 spr0.90"

# 1.4 Shape of control volume (3 = triangle, 6 = hexagon/pentagon)

export cell_type=3

# 1.5 Spatial resolution
# (This script assumes that a string indicating the spatial resolution,
#  e.g., "_R2B04L31", has been added after ${EXP} as part of the name
#  of the model output.)

horizontal_resolution="R2B04"
vertical_resolution="L31"

# 1.6 Whether the model output is in a single NetCDF file or
# has been split into separate files. (Here we assume that 
# in case of a single file, the string "_0001" has been added to 
# the name of the model output after ${EXP}"_RxBxxLxx".
# In case of output splitting, the first file is labeled "_0001",
# the second one "_0002", and so on.)

data_file_split=1     # 0 = single file; 1 = multiple files

# If the output is split, how often is a new file started? 

dt_trigger_file_hour=720

#--------------------------------------------------------------------------
# 2. Decide what to do
#--------------------------------------------------------------------------
# Which model day should be processed?

export diag_day=15

# This postprocessing tools draws the following panels on the same page:
# - surface pressure
# - 500 hPa geopotential height
# - 850 hPa zonal wind, meridional wind, relative vorticity, divergence,
#   temperature and vertical velocity (omega).
#
# Before plotting, vertical interpolation is needed. Set the following 
# variables to let the script do interplation and/or visualization. 

cn_plot_option=1   # (1=ON,0=OFF)
interp_option=1    # (1=ON,0=OFF)


# In which format (pdf, eps or ps)?

export plot_file_format="pdf"


# Remove these files after finishing visualization?
# (1=REMOVE,0=SAVE FOR LATER USE)

rm_tmp_files=1

# Check if there is a local file 'set_env' whith includes a special setting of the values

if [ -f ${set_env} ] 
then 
  echo " "
  echo " !!!!! Use setting from ./${set_env}"
  echo " "
  source ./${set_env}
fi

model_data_path="${icon_path}experiments/$EXP/"

# Where should the plot files be located? Don't forget the trailing "/".

plot_file_path="${model_data_path}plots/"

# Now specify the directory in which the pressure level data and other
# intermediate files should be placed. Don't forget the trailing "/".

tmp_data_path="${model_data_path}tmp/"

#--------------------------------------------------------------------------
#                    END OF USER'S SPECIFICATIONS
#==========================================================================

echo
echo "**********************************************************"
echo "***       ICON Tool Kit for Idealized Test Cases       ***"
echo "**********************************************************"
echo
echo "=== Postprocessing started for the Rossby-Haurwitz wave test."

# Temporary variables

export script_path=`pwd`
resolution=${horizontal_resolution}${vertical_resolution}
fori=${model_data_path}${EXP}_${resolution}
ftmp=${tmp_data_path}${EXP}_${resolution}

# The directories for intermediate data and plots will be created, if
# not already there

if [ ! -d ${plot_file_path} ]; then
   mkdir -p ${plot_file_path}
fi
if [ ! -d ${tmp_data_path} ]; then
   mkdir -p ${tmp_data_path}
fi

#---------------------------------------------------------------
# Select the single time step of data that will be processed
#---------------------------------------------------------------
# It is easy if model output is in a single file ...
#
if [ ${data_file_split} -eq 0 ]; then 

   fid=1
   tid=`expr 1 + ${diag_day} \* ${output_frequency}`

# Otherweise we need to compute ...
else
   skip=`expr ${diag_day} \* 24 / ${dt_trigger_file_hour}`
   rmdr=`expr ${diag_day} \* 24 % ${dt_trigger_file_hour}`

   # If the time step we want is the last step of a certain file ...
   if [ $rmdr -eq 0 ]; then
      fid=$skip
      tid=`expr ${dt_trigger_file_hour} \* ${output_frequency} / 24`

   # if it is before the last
   else
      fid=`expr $skip + 1`
      tid=`expr $rmdr \* ${output_frequency} / 24`
   fi

   # if in the first file, take into account the initial time step
   if [ $fid -eq 1 ]; then
      tid=`expr $tid + 1`
   fi
fi
  
   fid_str=$(printf "%04d" $fid) 

#==========================================================================
# Perform vertical interpolation and make contour plots
#==========================================================================

if [ $interp_option -eq 1 ];then
   echo
   echo "=== Performing vertical interpolation ..."

   if [ ! -f ${fori}_${fid_str}.nc ]; then
      echo "!!! Need model output ${fori}_${fid_str}.nc but cannot find it!"
      exit 1
   fi

   cdo selname,PHIS,PS,ZF3,U,V,VOR,DIV,T,OMEGA -seltimestep,$tid  \
       ${fori}_${fid_str}.nc ${ftmp}_day${diag_day}.nc

   check_error $? "In scripte RH4_postpro_driver.bash: cdo call 1"

   cdo ml2pl,85000,50000 ${ftmp}_day${diag_day}.nc     \
                         ${ftmp}_day${diag_day}_pres.nc

   check_error $? "In scripte RH4_postpro_driver.bash: cdo call 2"

   echo "=== Done."
fi

if [ $cn_plot_option -eq 1 ];then

   export Model="ICOHDC"
   export Data=${ftmp}_day${diag_day}_pres.nc
   export PlotPath=${plot_file_path}
   export Resolution=${resolution}
   export ConfigStr=${config_string}
   export ExpName=${EXP}_${resolution}
   echo
   echo "=== Making countour plots..."
   echo

   if [ ! -f ${Data} ]; then
      echo "!!! Need data file ${Data} for plotting but cannot find it!"
      exit 1
   fi

   ncl ${script_path}/RH4_postpro.ncl

   check_error $? "In scripte RH4_postpro_driver.bash: calling of RH4_postpro.ncl"

   echo "=== Done."
fi

#------------------------------------------------------------------------
# Clean up
#------------------------------------------------------------------------ 
echo
echo "=== Postprocessing finished for the RH wave test."

if [ $cn_plot_option -eq 1 ];then
   echo "=== The plots can be found in "${plot_file_path}
fi

if [ $rm_tmp_files -eq 1 ]; then

   rm ${tmp_data_path}/${EXP}_${resolution}*.nc

   if [ `ls ${tmp_data_path} |wc -l` -eq 0 ]; then
      rm -rf ${tmp_data_path}
   fi
   echo "=== The temporary data have been removed." 
else
   echo "=== The temporary data can be found in "${tmp_data_path}
fi

exit 0
