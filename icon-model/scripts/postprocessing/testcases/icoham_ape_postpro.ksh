#! /bin/ksh

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#--- General specifications ---

Model=ICOHAM
EXP=hat_ape_qobs_moist_RH00-clip  # String will show up in the plot
exp_config="$Model APE w/o O3"

cell_type=3
horizontal_resolution=R2B04
vertical_resolution=L31

model_data_path=/scratch/work/mh0287/hwan/exps/icon-dev/output_201007_port_phy/$EXP/
remap_weights_path=$model_data_path/remap_weights/
plot_file_path=$model_data_path/plots/
tmp_data_path=$model_data_path/tmp/

compute_remap_weights=1
grid_optimization="spr0.90"
trunc=63

do_computation=1
make_plot=1
plot_file_format=pdf

#--- Vertical cross section of zonal mean climate ---

diag_climate=1
clim_istart=171
clim_iend=172
TopHeight=35
plev_clim=""

ref_config="same run, days 151-220"  # String will show up in the difference plot
ref_exp_name=$EXP
ref_resolution=${horizontal_resolution}${vertical_resolution}
ref_timerange="151-220"
ref_data_path=$tmp_data_path

#--- Evolution of zonal mean/variance on individual vertical levels ---

evol_istart=1
evol_iend=220
plev_evol="85000"
hlev_evol=""
mlev_evol=""

#==========
# set env

set_env="Post.${EXP}.set_env."`date +%F-%H%M%S-%N`

cd ./APE_postpro

# Remove old ${set_env} file
 
if [ -f ${set_env} ]; then
  rm -f ${set_env}
fi

# Write new ${set_env} file
 
echo Model=$Model                                 >  ${set_env}
echo EXP=$EXP                                     >> ${set_env}
echo export exp_config=\"$exp_config\"            >> ${set_env}

echo cell_type=$cell_type                         >> ${set_env}
echo horizontal_resolution=$horizontal_resolution >> ${set_env}
echo vertical_resolution=$vertical_resolution     >> ${set_env}

echo model_data_path=$model_data_path             >> ${set_env}
echo remap_weights_path=$remap_weights_path       >> ${set_env}
echo plot_file_path=$plot_file_path               >> ${set_env}
echo tmp_data_path=$tmp_data_path                 >> ${set_env}

echo compute_remap_weights=$compute_remap_weights >> ${set_env}
echo grid_optimization=$grid_optimization         >> ${set_env}
echo trunc=$trunc                                 >> ${set_env}

echo do_computation=$do_computation               >> ${set_env}
echo make_plot=$make_plot                         >> ${set_env}
echo export plot_file_format=$plot_file_format    >> ${set_env}

echo diag_climate=$diag_climate                   >> ${set_env}
echo clim_istart=$clim_istart                     >> ${set_env}
echo clim_iend=$clim_iend                         >> ${set_env}
echo export TopHeight=$TopHeight                  >> ${set_env}
echo plev_clim=\"${plev_clim}\"                   >> ${set_env}
echo export ref_config=\"$ref_config\"            >> ${set_env}
echo ref_exp_name=$ref_exp_name                   >> ${set_env}
echo ref_resolution=$ref_resolution               >> ${set_env}
echo ref_data_path=$ref_data_path                 >> ${set_env}
echo ref_timerange=$ref_timerange                 >> ${set_env}

echo evol_istart=$evol_istart                     >> ${set_env}
echo evol_iend=$evol_iend                         >> ${set_env}
echo plev_evol=\"${plev_evol}\"                   >> ${set_env}
echo hlev_evol=\"${hlev_evol}\"                   >> ${set_env}
echo mlev_evol=\"${mlev_evol}\"                   >> ${set_env}

#=============
# Run scripts
#
# Note: Make sure that
#   1. all variable listed below exist in model output and 
#      have been registered in "lookup_variable.ksh";
#   2. PS and PHIS are processed before any other variable
#      (because they are needed for the vertical 
#      interpolation from model levels to pressure levels).
#
# --- Zonal mean climate (vertical cross section)---
 
for var in PS PHIS T U V OMEGA Qv Qw Qi ACLC ; do
  ./zonal_clim.ksh -v ${var} -e ${set_env}
done

# --- Evolution plots ---

for var in PS PHIS T ; do
  ./lat-time.ksh -v ${var} -e ${set_env}
done

#==========
# Clean up

if [ -f ${set_env} ]; then
  rm -f ${set_env}
fi
exit
