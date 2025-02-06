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
#      Driver script for postprocessing and visulization of a  
#                   vertical cross sesction for ICONAM/ICOHDC   
#
#==========================================================================
#
# History: 
# Initial version by Pilar Ripodas (DWD, 2010-11) (driver based on a driver 
# from Hui, ncl script based on a script from Guenther)
# Software needed:
# - CDO (Climate Data Operators, www.mpimet.mpg.de/cdo) 
#   for interplation and for the spectral transform;
# - NCL (NCAR Command Language, www.ncl.ucar.edu)
#   for visualization.
# CDO is used to interpolate the icon file to a gaussian grid, NCL is used 
# to generate and plot the vertical cross section
#
#==============================================================================
#
#==========================================================================
#                          USER'S SPECIFICATIONS 
#--------------------------------------------------------------------------
# 0. About the dynamical core ( Hydrostatic or Non-Hydrostatic )
#--------------------------------------------------------------------------

set -x

export Model="ICOHDC"
#export Model="ICONAM"
#EDIR="exp39"
#EDIR="exp_echam"
EDIR="hat_icoham_ape"
#--------------------------------------------------------------------------
# 1. About the model output
#--------------------------------------------------------------------------


# 1.1 In which directory is the data file located?
# (Don't forget the trailing "/")

#export model_data_path="/e/uscratch/mripodas/icon-dev/experiments/hat_heldsuarez_2tl/"
#export model_data_path="/e/uscratch/mripodas/icon-dev/experiments/nh35_heldsuarezB4/"
#export model_data_path="/e/uscratch/mripodas/icon-dev/experiments/Guehnter/"

export model_data_path="/scratch/mpi/CC/mh0287/users/kristina/icon-dev/experiments/${EDIR}/"
#export model_data_path="/e/uhome/gzaengl/icon-dev/experiments/${EDIR}/"

# temporary dir if files are from 'permission denied' dir's
export model_data_path2="/scratch/mpi/CC/mh0287/users/kristina/icon-dev/experiments/${EDIR}/"

#1.2 Name of file containing the data
export file_name="hat_icoham_ape_iconR2B04-grid_0001.nc"
#export file_name="JWmoistturbnh_DOM01_R2B05L40_0010.nc"
#export file_name="JW_R2B5_nh_DOM01_R2B05L60_0008.nc"
#export file_name="JW_echam_DOM01_R2B05L31_0001.nc"

# 1.3 Shape of control volume (3 = triangle, 6 = hexagon/pentagon) 
export cell_type=3

# 1.4 Spatial resolution
horizontal_resolution="R2B04"

# 1.5 Specify if the file contains the global model domain or a refined area. 
# In case it is a refined area, then give a name that you asociate to it 
# (it will be appended to the weights file for the interpolation to a gaussian grid)

i_ref=0        # (0 for the global domain, any other value for a refined area domain)

#ref_name="gzaengl"    #
ref_name=""


#--------------------------------------------------------------------------
# 2. Set the desired variable and time step for the vcs
#    Set also the (lat,lon) values for the transect
#    If the range of values is already know, the desired contour levels can be 
#    set manually, if not, they are set by NCL
#--------------------------------------------------------------------------

export var="ozone"

timestep="1"

export i_manual_levels=F # (T, the contour levels are set manually 
                         #  any other value are automatically set by NCL,
                         
if [ $i_manual_levels == "T" ]; then
 export minlev=170
 export maxlev=330
 export interval=2
fi

#
# if the transect has no constant longitude or not constant latitud, then 
# the size of the transect can not be larger than 180. degress (the 
# NCL function gc_latlon would consider the shortest arc between the 2 points)
# if leftlat=rightlat  or leftlon=rightlon, there is no restriction
#

lela="-90"
rila="90"
lelo="0"
rilo="0"

export leftlat="${lela}.0"
export rightlat="${rila}.0"
export leftlon="${lelo}.0"
export rightlon="${rilo}.0"

#leftlon=`echo $leftlon | sed 's/\./ /' | awk '{ print $1 }'`

#--------------------------------------------------------------------------


# To do the vcs, fisrt an interpolation to a gaussian grid is done, 
# Choose the Gaussian grid by specifying the triangular truncation:

export trunc=85      #106

# For the horizontal interpolation from the ICON grid to Gaussian grid,
# a large part of the time will be spent on calculating the remapping 
# weights. This may take very long when the resolution is high.
# Therefore we suggest calculating the weights only once and store them 
# for later use. 

compute_remap_weights=1   # (1=ON,0=OFF)

# The remapping weights file generated by this script will be named, 
# e.g., icon_R2B04_spr0.90_tri_cell_to_T159.nc 
# If the weights are already available, specify the location:
# (Don't forget the trailing "/")

remap_weights_path='/scratch/mpi/CC/mh0287/users/kristina/icon-dev/experiments/weights/'

# For the remapping, we also need to know which optimization was used 
# to generate the ICON grid. For example,
#   "ori"     : original icosahedral grid without optimization
#   "hro"     : Heikes-Randall optimization
#   "spr0.90" : spring dynamics, with spring coefficient 0.90

grid_optimization="spr0.90"


# Where should the plot files be located? Don't forget the trailing "/".

#export plot_file_path="${model_data_path}/"

export plot_file_path="/scratch/mpi/CC/mh0287/users/kristina/icon-dev/experiments/plots/"


# Name of the output plot file

#export plot_file="${var}_step${timestep}_${file_name%.nc}_transect"
export plot_file="${var}_step${timestep}_${horizontal_resolution}_transect${lela}-${lelo}_${rila}-${rilo}"

# What format do you want, pdf, eps or ps?

export plot_file_format="ps"

#--------------------------------------------------------------------------
# Now specify the directory in which the intermediate files 
# (excluding the remapping weights) should be placed. 
# Don't forget the trailing "/".

# model_data path differs if dir is not own
#tmp_data_path="${model_data_path}"
tmp_data_path="${model_data_path2}"

# Remove these files after finishing the diagnoses? 
# (1=REMOVE,0=SAVE FOR LATER USE)

rm_tmp_files=0


#--------------------------------------------------------------------------
# Do you want CDO to run in silence mode, or to report everything
# it is doing?

cdo_silence=0   #( 1 = silence mode; 0 = detailed report )

#--------------------------------------------------------------------------
#                    END OF USER'S SPECIFICATIONS 
#==========================================================================

# Temporary variables
export script_path='/scratch/mpi/CC/mh0287/users/kristina/icon-dev/scripts/postprocessing/tools/'

#
#directory with the colormap files
export NCARG_COLORMAPS=${script_path}color_map:${NCARG_COLORMAPS}
#

fori=${model_data_path}${file_name}
ftmp=${tmp_data_path}${file_name%.nc}

if [ $cdo_silence -eq 1 ]; then
   silence='-s'
else
   silence=''
fi

# The directories for intermediate data and plots will be created, if 
# not already there

if [ ! -d ${plot_file_path} ]; then
   mkdir -p ${plot_file_path} 
fi
if [ ! -d ${tmp_data_path} ]; then
   mkdir -p ${tmp_data_path} 
fi

#==========================================================================
# Prepare remappping weights
#==========================================================================

if [ ${i_ref} -eq 0 ]; then
  weights=${remap_weights_path}"icon_"${horizontal_resolution}
else 
  weights=${remap_weights_path}"icon_"${horizontal_resolution}"_"${ref_name}
fi

if [ ${cell_type} -eq 3 ]; then
   weights=${weights}"_"${grid_optimization}"_cell_to_T"${trunc}".nc"
elif [ ${cell_type} -eq 6 ]; then
   weights=${weights}"_"${grid_optimization}"_vert_to_T"${trunc}".nc"
else
   echo "Wrong choice of cell_type. Should be 3 or 6 !"
   exit 1
fi

# compute weights

if [ ${compute_remap_weights} -eq 1 ]; then

 if [ ! -d ${remap_weights_path} ]; then
    mkdir -p ${remap_weights_path} 
 fi

 echo
 echo "=== Computing remapping weights (ICON to Gaussian) ..."

 label=$(printf "%04d" ${evol_istart})
# cdo $silence gendis,t${trunc}grid -seltimestep,1 -selname,PS ${fori} ${weights}
  cdo  $silence gendis,t${trunc}grid -seltimestep,1 -selname,PS ${fori} ${weights}

 check_error $? "In script plot_vcs_driver.bash: part 'Computing remapping weights (ICON to Gaussian)'"
 echo "=== Done."
fi

#======================================================================================
# Select the time step and the variable to plot, interpolate to the gaussian grid
#======================================================================================


 if [ $Model == "ICONAM" ]; then
   cdo $silence  selname,${var},ZH3,ZF3  -seltimestep,${timestep} ${fori} ${ftmp}_${var}_height.nc

   cdo $silence remap,t${trunc}grid,${weights} \
                       ${ftmp}_${var}_height.nc  ${ftmp}_${var}_height_T${trunc}.nc


 fi
 if [ $Model == "ICOHDC" ]; then
   cdo $silence  selname,${var},ZF3,PHIS  -seltimestep,${timestep} ${fori} ${ftmp}_${var}_height.nc

   cdo $silence remap,t${trunc}grid,${weights} \
                       ${ftmp}_${var}_height.nc  ${ftmp}_${var}_height_T${trunc}.nc

 fi

 export file_int=${ftmp}_${var}_height_T${trunc}.nc


#======================================================================================
# run the ncl script to do the plot
#======================================================================================

 ncl  ${script_path}/plot_vcs.ncl

#------------------------------------------------------------------------
# Clean up
#------------------------------------------------------------------------

if [ $rm_tmp_files -eq 1 ]; then

  rm ${ftmp}_${var}_height.nc 
  rm ${ftmp}_${var}_height_T${trunc}.nc

fi


exit
