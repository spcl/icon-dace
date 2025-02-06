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
#           Jablonowski-Williamson baroclinic wave test 
#==========================================================================
#
# History: 
# Initial version by Hui Wan (MPI-M, 2009-04)
#
# Short description:
# This script interpolates the orginal ICOHDC output on hybrid 
# vertical levels to the 850 hPa pressure level, then makes
# contour plots. If the user wishes, the kinetic energy spectra
# at various pressure levels can also be diagnosed and plotted.
#
# Software needed:
# - ${cdo} (Climate Data Operators, www.mpimet.mpg.de/${cdo}) 
#   for interplation and for the spectral transform;
# - NCL (NCAR Command Language, www.ncl.ucar.edu)
#   for visualization.
#
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

set -ex

#==========================================================================
#                          USER'S SPECIFICATIONS 
# Check if a parameter is given

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


echo "JWw_postpro_driver.bash \$\$ $$"

# for automatic testing
dir=`pwd -P`
icon_path=${dir%%scripts/postprocessing/testcases}
icon_path=${basedir}/

# Output frequency expressed as number of time slices per day. 
# (The value you give here must be ASCII representations of integer values, 
# such as "1" for daily output, and "4" for 6-hourly data.)

export output_frequency=4

# 1.2 The data file name is constructed in the same way as in the model.
# The name is composed of experiment name, hor. and vert. resolution (cf. 1.5),
# and a file index (cf. 1.6), if the output has been split on several files.
# The experiment name is the first part of the file name of your model output. 
# You may have used something like "A0001" just for simplicity.

EXP="test_hat_jww"
DOMAIN=""

# 1.3 The experiment identifier that will appear in the plots, e.g.,
# "brk R2B5L31" 

config_string="JWw tri R2B04L31 spr0.90"

# 1.4 Shape of control volume (3 = triangle, 6 = hexagon/pentagon) 

export cell_type=3

# At which location (cell centers or corners) is the vorticity field provided

export vorticity_at_corners=1  # 0 = center, 1 = corner

# 1.5 Spatial resolution
# (This script assumes that a string indicating the spatial resolution,
#  e.g., "_R2B04L31", has been added after ${EXP} as part of the name
#  of the model output.)

horizontal_resolution="R2B04"
vertical_resolution="L47"

# 1.6 Whether the model output is in a single NetCDF file or 
# has been split into separate files. (Here we assume that
# in case of a single file, the string "_0001" has been added to
# the name of the model output after ${EXP}"_RxBxxLxx".
# In case of output splitting, the first file is labeled "_0001",
# the second one "_0002", and so on.)

data_file_split=0     # 0 = single file; 1 = multiple files

# Note that for the 850 hPa contour plots, only the first 10 days 
# of model output are needed. For the kinetic energy spectra
# and for calculation of the difference norms, we need all the 30 days.

#--------------------------------------------------------------------------
# 2. Decide what to do
#--------------------------------------------------------------------------
# A complete set of the contour plots includes
#
# - evolution of surface pressure and 850 hPa temperature (days 4-10);
# - 850 hPa vorticity, divergence and omega at days 7 and 9, focusing 
#   on the main signals;
# - 850 hPa vorticity, divergence and omega at day 9, showing the 
#   region 25N-75N, 0-360 degrees longitude. 
#
# You can let the script produce the whole set by setting the following 
# variable to 1, or switch off the contour plotting completely (e.g., when
# you only want to diagnose the kinetic energy spectrum) by setting the 
# variable to zero,

cn_plot_option=2

# or set it to any other integer, and switch on or off each individual 
# variable separately (1=ON,0=OFF):

plot_ps=1
plot_temp=1
plot_vor=1
plot_div=1
plot_omega=1

# in which format (pdf, eps or ps)?

export plot_file_format="eps"

#--------------------------------------------------------------------------
# Before plotting, vertical interpolation is needed. If you have already
# done this and stored the data, turn this step off.
# Like for plotting, you can choose to interpolate all the variables
# listed above (except for surface pressure) by setting the following 
# variable to 1, to turn off interpolation completely by setting it to 0, 

interp_option=1    # (1=ON,0=OFF)

# or switch on or off individual variables (1=ON,0=OFF):

interp_temp=1
interp_vor=1
interp_div=1
interp_omega=1

#--------------------------------------------------------------------------
# 850 hPa kinetic energy spectra at days 15 and 30 will be 
# calculated and plotted if the next variable is set to 1.

ke_spectrum_diag=0

# For this diagnosis the ICON model output (divergence and vorticity) 
# on the geodesic grid needs to be first interpolated to a Gaussian grid, 
# and then transformed into spectral coefficients. 

# Choose the Gaussian grid by specifying the triangular truncation:

trunc=85

# For the horizontal interpolation from the ICON grid to Gaussian grid,
# a large part of the time will be spent on calculating the remapping 
# weights. This may take very long when the resolution is high.
# Therefore we suggest calculating the weights only once and store them 
# for later use. 

compute_remap_weights=0   # (1=ON,0=OFF)

# For the remapping, we also need to know which optimization was used 
# to generate the ICON grid. For example,
#   "ori"     : original icosahedral grid without optimization
#   "hro"     : Heikes-Randall optimization
#   "spr0.90" : spring dynamics, with spring coefficient 0.90

grid_optimization="spr0.90"

# Remove these files after finishing the diagnoses? 
# (1=REMOVE,0=SAVE FOR LATER USE)

rm_tmp_files=1

cdo="cdo"
if [ -f ${set_env} ] 
then 

  echo " "
  echo " !!!!! Use setting from ./${set_env}"
  echo " "
  source ./${set_env}
fi


# Where should the plot files be located? Don't forget the trailing "/".

model_data_path="${icon_path}experiments/$EXP/"

plot_file_path="${model_data_path}plots/"

#--------------------------------------------------------------------------
# Now specify the directory in which the pressure level data and other
# intermediate files (excluding the remapping weights) should be placed.
# Don't forget the trailing "/".

tmp_data_path="${model_data_path}tmp/"

resolution=${horizontal_resolution}${vertical_resolution}
# fori=${model_data_path}${EXP}${DOMAIN}_${resolution}
# ftmp=${tmp_data_path}${EXP}_${resolution}

# For the NCL plotting scripts
export script_path=`pwd`"/JWw_postpro_scripts/"
export Model="ICOHDC"
export DataPath=${tmp_data_path}
export PlotPath=${plot_file_path}
export Resolution=${resolution}
export ConfigStr=${config_string}
if [ x${ExpName} == "x" ]; then
  echo "Error in ExpName:", $ExpName
  exit
fi
export ExpName=${ExpName}
export plotBaseName=${plotBaseName}
fori=${model_data_path}${ExpName}
ftmp=${tmp_data_path}${ExpName}



# Location of the file that contains /will contain the weights. 
# Don't forget the trailing "/".

remap_weights_path="${model_data_path}remap_weights/"

# The remapping weights file generated by this script will be named, 
# e.g., icon_R2B04_spr0.90_cell_to_T159.nc 

#--------------------------------------------------------------------------
#                    END OF USER'S SPECIFICATIONS 
#==========================================================================

echo
echo "**********************************************************"
echo "***       ICON Tool Kit for Idealized Test Cases       ***"
echo "**********************************************************"
echo 
echo "=== Postprocessing started for the JW wave test."
#===================== Info ==================================
# echo "===================== Info =================================="
#  echo "EXP=$EXP"
#  echo "output_frequency=${output_frequency}"
#  echo "config_string=${config_string}"
#  echo "cell_type=${cell_type}"
#  echo "horizontal_resolution=${horizontal_resolution}"
#  echo "vertical_resolution=${vertical_resolution}"
#  echo "data_file_split=${data_file_split}"
#  echo "interp_option=${interp_option}"
#  echo "export plot_file_format=${plot_file_format}"
#  echo "rm_tmp_files=${rm_tmp_files}"
# 
# echo "cn_plot_option=${cn_plot_option}"
# echo "plot_ps=${plot_ps}"
# echo "plot_temp=${plot_temp}"
# echo "info_plot_vor=${plot_vor}"
# echo "plot_div=${plot_div}"
# echo "plot_omega=${plot_omega}"
# echo "interp_temp=${interp_temp}"
# echo "interp_vor=${interp_vor}"
# echo "interp_div=${interp_div}"
# echo "interp_omega=${interp_omega}"
# echo "ke_spectrum_diag=${ke_spectrum_diag}"
# echo "trunc=${trunc}"
# echo "compute_remap_weights=${compute_remap_weights}"
# echo "grid_optimization=${grid_optimization}"
# 
# echo "Path:"
# echo "model_data_path=${model_data_path}"
# echo "plot_file_path=${plot_file_path}"
# echo "tmp_data_path=${tmp_data_path}"
# echo "remap_weights_path=${remap_weights_path}"
# echo "===================== Info =================================="

#===================== Info ==================================
# Temporary variables


# The directories for intermediate data and plots will be created, if 
# not already there

if [ ! -d ${plot_file_path} ]; then
   mkdir -p ${plot_file_path} 
fi
if [ ! -d ${tmp_data_path} ]; then
   mkdir -p ${tmp_data_path} 
fi

# If necessary, combine the daily output into a single file

if [ ${data_file_split} -eq 1 ]; then
   echo
   echo "=== Merging model output ..."
   ${cdo} copy ${fori}*.nc ${ftmp}.nc

   check_error $? "In script JWw_postpro_driver.bash: part 'Merging model output'"

   echo "=== Done."
else
   echo
   echo "=== Copying model output ..."
   cp ${fori}*.nc ${ftmp}.nc

   check_error $? "In script JWw_postpro_driver.bash: part 'Copying model output'"

   echo "=== Done."
fi

# We have to call NCL from the directory where our own colormaps and resource 
# files are located, so that they can be loaded by the NCL scripts correctly. 

cd ${script_path} 

# Make sure that NCL finds the color maps defined in this directory.

export NCARG_COLORMAP_PATH=$NCARG_ROOT/lib/ncarg/colormaps:./


#==========================================================================
# Perform vertical interpolation and make contour plots
#==========================================================================
if [ ${interp_option} -eq 1 ]; then
   interp_temp=1
   interp_vor=1
   interp_div=1
   interp_omega=1
elif [ ${interp_option} -eq 0 ]; then
   interp_temp=0
   interp_vor=0
   interp_div=0
   interp_omega=0
fi

if [ ${cn_plot_option} -eq 1 ]; then
   plot_ps=1
   plot_temp=1
   plot_vor=1
   plot_div=1
   plot_omega=1
elif [ ${cn_plot_option} -eq 0 ]; then
   plot_ps=0
   plot_temp=0
   plot_vor=0
   plot_div=0
   plot_omega=0
fi


#------------------------------------------------------------------------
# plot ps at selected time steps (see the ncl script)
#------------------------------------------------------------------------     
if [ $plot_ps -eq 1 ]; then
   export JWw_VarName="PS"
   echo "=================================="
   echo "=== $(date)"
   echo "=== Plotting surface pressure ..."
   ncl JWw_plot_PS-T_evol.ncl 

   check_error $? "In script JWw_postpro_driver.bash: part 'Plotting surface pressure'"

   echo "=== $(date)"
   echo "=== Done."
fi

T_var=""

#------------------------------------------------------------------------
# interpolate and plot temperature
#------------------------------------------------------------------------
if [ $interp_temp -eq 1 ];then
   echo "=================================="
   echo "=== $(date)"
   echo "=== Interpolating temperature ..."

   ${cdo} selname,ta    ${ftmp}".nc" ${ftmp}"_T_hyb.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating temperature' ${cdo} call 1"

   ${cdo} selname,ps   ${ftmp}".nc" ${ftmp}"_PS.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating temperature' ${cdo} call 2"

   ${cdo} -chname,orog,geosp -mulc,9.81 -selname,orog ${ftmp}".nc" ${ftmp}"_PHIS.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating temperature' ${cdo} call 3"

   ${cdo} merge ${ftmp}"_T_hyb.nc" ${ftmp}"_PS.nc" ${ftmp}"_PHIS.nc" ${ftmp}"_T-PS-PHIS.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating temperature' ${cdo} call 4"

   ${cdo} -v ml2pl,85000 -chname,ps,aps ${ftmp}"_T-PS-PHIS.nc" ${ftmp}"_T850.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating temperature' ${cdo} call 5"

   rm  ${ftmp}"_T_hyb.nc" ${ftmp}"_PS.nc" ${ftmp}"_PHIS.nc" ${ftmp}"_T-PS-PHIS.nc"

   echo "=== $(date)"
   echo "=== Done."
fi

if [ $plot_temp -eq 1 ]; then
   export JWw_VarName="T"
   echo "=================================="
   echo "=== $(date)"
   echo "=== Plotting 850 hPa temperature ..."
   ncl JWw_plot_PS-T_evol.ncl 

   check_error $? "In script JWw_postpro_driver.bash: part 'Plotting 850 hPa temperature'"

   echo "=== $(date)"
   echo "=== Done."
fi

#------------------------------------------------------------------------
# interpolate and plot divergence 
#------------------------------------------------------------------------
if [ $interp_div -eq 1 ]; then
   echo "=================================="
   echo "=== $(date)"
   echo "=== Interpolating divergence ..."
   ${cdo} selname,sd  ${ftmp}".nc" ${ftmp}"_DIV_hyb.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating divergence' ${cdo} call 1"

   ${cdo} selname,ps   ${ftmp}".nc" ${ftmp}"_PS.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating divergence' ${cdo} call 2"

   ${cdo} merge ${ftmp}"_DIV_hyb.nc" ${ftmp}"_PS.nc" ${ftmp}"_DIV-PS.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating divergence' ${cdo} call 3"

   ${cdo} -v ml2pl,85000 -chname,ps,aps ${ftmp}"_DIV-PS.nc" ${ftmp}"_DIV850.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating divergence' ${cdo} call 4"

   echo "=== $(date)"
   rm  ${ftmp}"_DIV_hyb.nc" ${ftmp}"_PS.nc" ${ftmp}"_DIV-PS.nc"
   echo "=== Done."
fi

if [ $plot_div -eq 1 ]; then
   export JWw_VarName="DIV"
   echo "=================================="
   echo "=== $(date)"
   echo "=== Plotting 850 hPa divergence ..."
   ncl JWw_plot_DIV-VOR_day0709.ncl 

   check_error $? "In script JWw_postpro_driver.bash: part 'Plotting 850 hPa divergence' (JWw_plot_DIV-VOR_day0709.ncl)"

   ncl JWw_plot_DIV-VOR_day09.ncl

   check_error $? "In script JWw_postpro_driver.bash: part 'Plotting 850 hPa divergence' (JWw_plot_DIV-VOR_day09.ncl)"

   echo "=== $(date)"
   echo "=== Done."
fi

#------------------------------------------------------------------------
# interpolate and plot omega 
#------------------------------------------------------------------------
if [ $interp_omega -eq 1 ]; then
   echo "=================================="
   echo "=== $(date)"
   echo "=== Interpolating omega ..."
   ${cdo} selname,wap  ${ftmp}".nc" ${ftmp}"_OMEGA_hyb.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating omega' ${cdo} call 1"

   ${cdo} selname,ps     ${ftmp}".nc" ${ftmp}"_PS.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating omega' ${cdo} call 2"

   ${cdo} merge ${ftmp}"_OMEGA_hyb.nc" ${ftmp}"_PS.nc" ${ftmp}"_OMEGA-PS.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating omega' ${cdo} call 3"

   ${cdo} -v ml2pl,85000 -chname,ps,aps ${ftmp}"_OMEGA-PS.nc" ${ftmp}"_OMEGA850.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating omega' ${cdo} call 4"

   echo "=== $(date)"
   rm  ${ftmp}"_OMEGA_hyb.nc" ${ftmp}"_PS.nc" ${ftmp}"_OMEGA-PS.nc"
   echo "=== Done."
fi

if [ $plot_omega -eq 1 ]; then
   export JWw_VarName="OMEGA"
   echo "=================================="
   echo "=== $(date)"
   echo "=== Plotting 850 hPa omega ..."
   ncl JWw_plot_DIV-VOR_day0709.ncl 

   check_error $? "In script JWw_postpro_driver.bash: part 'Plotting 850 hPa omega' (JWw_plot_DIV-VOR_day0709.ncl)"

   ncl JWw_plot_DIV-VOR_day09.ncl

   check_error $? "In script JWw_postpro_driver.bash: part 'Plotting 850 hPa omega' (JWw_plot_DIV-VOR_day09.ncl)"

   echo "=== $(date)"
   echo "=== Done."
fi

#------------------------------------------------------------------------
# interpolate and plot vorticity 
#------------------------------------------------------------------------
if [ $interp_vor -eq 1 ]; then
   echo "=================================="
   echo "=== $(date)"
   echo "=== Interpolating vorticity ..."
   ${cdo} selname,svo  ${ftmp}".nc" ${ftmp}"_VOR_hyb.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating vorticity' ${cdo} call 1"

   if [ $vorticity_at_corners -eq 1 ]; then
     ${cdo} remapdis,${ftmp}"_VOR_hyb.nc" -selname,ps ${ftmp}".nc" ${ftmp}"_PS.nc"
   else
     ${cdo} selname,ps  ${ftmp}".nc" ${ftmp}"_PS.nc"
   fi

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating vorticity' ${cdo} call 2"

   ${cdo} merge ${ftmp}"_VOR_hyb.nc" ${ftmp}"_PS.nc" ${ftmp}"_VOR-PS.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating vorticity' ${cdo} call 3"

   ${cdo} -v ml2pl,85000 -chname,ps,aps ${ftmp}"_VOR-PS.nc" ${ftmp}"_VOR850.nc"

   check_error $? "In script JWw_postpro_driver.bash: part 'Interpolating vorticity' ${cdo} call 4"

   rm  ${ftmp}"_VOR_hyb.nc" ${ftmp}"_PS.nc" ${ftmp}"_VOR-PS.nc"
   echo "=== $(date)"
   echo "=== Done."
fi

if [ $plot_vor -eq 1 ]; then
   export JWw_VarName="VOR"
   echo "=================================="
   echo "=== $(date)"
   echo "=== Plotting 850 hPa vorticity ..."
   ncl JWw_plot_DIV-VOR_day0709.ncl 

   check_error $? "In script JWw_postpro_driver.bash: part 'Plotting 850 hPa vorticity' (JWw_plot_DIV-VOR_day0709.ncl)"

   ncl JWw_plot_DIV-VOR_day09.ncl

   check_error $? "In script JWw_postpro_driver.bash: part 'Plotting 850 hPa vorticity' (JWw_plot_DIV-VOR_day09.ncl)"

   echo "=== $(date)"
   echo "=== Done."
fi

#------------------------------------------------------------------------
# Diagnose kinetic energy spectra
#------------------------------------------------------------------------

weights=${remap_weights_path}"icon_"${horizontal_resolution}

if [ ${cell_type} -eq 3 ]; then
   weights=${weights}"_"${grid_optimization}"_cell_to_T"${trunc}".nc"
elif [ ${cell_type} -eq 6 ]; then
   weights=${weights}"_"${grid_optimization}"_vert_to_T"${trunc}".nc"
else
   echo "Wrong choice of cell_type. Should be 3 or 6 !"
   exit 1
fi

#-----------------------------
# compute remapping weights

if [ ${compute_remap_weights} -eq 1 ]; then

 if [ ! -d ${remap_weights_path} ]; then
    mkdir -p ${remap_weights_path} 
 fi

 if [ -e ${ftmp}"_DIV850.nc" ]; then
    echo
    echo "=== Computing remapping weights (ICON to Gaussian) ..."
    ${cdo} gendis,t${trunc}grid -selname,sd ${ftmp}"_DIV850.nc" ${weights}

    check_error $? "In script JWw_postpro_driver.bash: part 'Computing remapping weights (ICON to Gaussian)'"

    echo "=== Done."
 else
    echo "    Error: File not found ("${ftmp}"_DIV850.nc)" 
    exit 1
 fi

fi

#-------------------------------------------------------
# convert the grid point data to spectral coefficients 

if [ ${ke_spectrum_diag} -eq 1 ]; then

# Check if the divergence and vorticity data are available

  if [ ! -e ${ftmp}"_DIV850.nc" ]; then
     echo "    Error: File not found ("${ftmp}"_DIV850.nc)" ; exit 1
  fi
  if [ ! -e ${ftmp}"_VOR850.nc" ]; then
     echo "    Error: File not found ("${ftmp}"_VOR850.nc)" ; exit 1
  fi
  if [ ! -e ${weights} ]; then
     echo "    Error: File not found ("${weights}")" ; exit 1
  fi

  echo
  echo "=== Computing kinetic energy spectrum ..."

# Interpolate to Gaussian grid

  ${cdo} -r remap,t${trunc}grid,${weights} -selname,sd  \
         ${ftmp}"_DIV850.nc" ${ftmp}"_DIV850_T"${trunc}"_gp.nc"

  check_error $? "In script JWw_postpro_driver.bash: part 'Interpolate to Gaussian grid' ${cdo} call 1"

  ${cdo} -r remap,t${trunc}grid,${weights} -selname,svo  \
         ${ftmp}"_VOR850.nc" ${ftmp}"_VOR850_T"${trunc}"_gp.nc"

  check_error $? "In script JWw_postpro_driver.bash: part 'Interpolate to Gaussian grid' ${cdo} call 2"

  ${cdo} merge ${ftmp}"_DIV850_T"${trunc}"_gp.nc"  \
            ${ftmp}"_VOR850_T"${trunc}"_gp.nc"  \
            ${ftmp}"_DIV-VOR850_T"${trunc}"_gp.nc"

  check_error $? "In script JWw_postpro_driver.bash: part 'Interpolate to Gaussian grid' ${cdo} call 3"

           
# Transform to spectral space 

  ${cdo} gp2sp ${ftmp}"_DIV-VOR850_T"${trunc}"_gp.nc"   \
            ${ftmp}"_DIV-VOR850_T"${trunc}"_spec.nc"

  check_error $? "In script JWw_postpro_driver.bash: part 'Transform to spectral space'"

  rm ${ftmp}"_DIV850_T"${trunc}"_gp.nc" ${ftmp}"_VOR850_T"${trunc}"_gp.nc"

# Compute and plot the kinetic energy spectrum

  export SpecData=${ftmp}"_DIV-VOR850_T"${trunc}
  export vor_coeff_name="VOR"
  export div_coeff_name="DIV"

  echo "=================================="
  echo "=== $(date)"
  ncl ke_spectrum_cal.ncl

   check_error $? "In script JWw_postpro_driver.bash: part 'Compute and plot the kinetic energy spectrum' (ke_spectrum_cal.ncl)"

  ncl ke_spectrum_plot.ncl

   check_error $? "In script JWw_postpro_driver.bash: part 'Compute and plot the kinetic energy spectrum' (ke_spectrum_plot.ncl)"

  echo "=== $(date)"
  echo "=== Done."
fi

#------------------------------------------------------------------------
# Clean up
#------------------------------------------------------------------------

echo 
echo "=== Postprocessing finished for the JW wave test."
echo "=== The plots can be found in "${plot_file_path}

if [ $rm_tmp_files -eq 1 ]; then

   rm ${tmp_data_path}/${ExpName}*.nc

   if [ `ls ${tmp_data_path} |wc -l` -eq 0 ]; then
      rm -rf ${tmp_data_path}
   fi 
   echo "=== The temporary data have been removed."

else
   echo "=== The pressure level data and spectral coefficients can be found in "${tmp_data_path}
fi
