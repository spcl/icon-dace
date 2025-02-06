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

expnum="exp03"
expnum2="exp03"   #for comparison in plots, data must pre-exist

Model=ICONAM
EXP=APE
exp_config="$Model APE exp $expnum"   # string will show up in the plot

cell_type=3
horizontal_resolution=R2B04
vertical_resolution=L90

model_data_path="/e/uwork/mkoehler/icon/experiments/"${expnum}"/"
remap_weights_path="/e/uwork/mkoehler/icon/experiments/remap_weights/"
plot_file_path=$model_data_path/plots/
tmp_data_path=$model_data_path/tmp/

compute_remap_weights=0
grid_optimization="spr0.90"
trunc=85

do_computation=1
make_plot=1
plot_file_format=ps
wkOrientation=landscape

#--- Vertical cross section of zonal mean climate ---

diag_climate=1
clim_istart=1
clim_iend=2
TopHeight=50

ref_config="$expnum2"   # string will show up in the difference plot
ref_exp_name=$EXP
ref_resolution=${horizontal_resolution}${vertical_resolution}
ref_timerange="1-1"
ref_data_path="/e/uwork/mkoehler/icon/experiments/"${expnum2}"/tmp/"

#--- Evolution of zonal mean/variance on individual vertical levels ---

evol_istart=1
evol_iend=2
plev_evol=""
hlev_evol=""
mlev_evol="80,69,52,38,11"

#==========
# set env

set_env="Post.${EXP}.set_env."`date +%F-%H%M%S-%N`

cd ./APE_postpro

# Remove old ${set_env} file
#
if [ -f ${set_env} ]; then
  rm -f ${set_env}
fi

# Write new ${set_env} file
#
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
echo export wkOrientation=$wkOrientation          >> ${set_env}

echo diag_climate=$diag_climate                   >> ${set_env}
echo clim_istart=$clim_istart                     >> ${set_env}
echo clim_iend=$clim_iend                         >> ${set_env}
echo export TopHeight=$TopHeight                  >> ${set_env}
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

# NOTE: Any variable listed below must exist in model output and 
# have been registered in "lookup_variable.ksh".

# --- Zonal mean climate (vertical cross section) ---

#for var in T U V W QV Q1 QC Q2 QI Q3 Q4 Q5 CC ; do
#  ./zonal_clim.ksh -v $var -e ${set_env}
#done

# DWD parallel execution (alternative to loop)
cat > zonal.list << EOF_LIST
  zonal_clim.ksh -v T  -e ${set_env}
  zonal_clim.ksh -v U  -e ${set_env}
  zonal_clim.ksh -v V  -e ${set_env}
  zonal_clim.ksh -v W  -e ${set_env}
  zonal_clim.ksh -v CC -e ${set_env}
  zonal_clim.ksh -v QV -e ${set_env}
  zonal_clim.ksh -v QC -e ${set_env}
  zonal_clim.ksh -v QI -e ${set_env}
  zonal_clim.ksh -v Q1 -e ${set_env}
  zonal_clim.ksh -v Q2 -e ${set_env}
  zonal_clim.ksh -v Q3 -e ${set_env}
  zonal_clim.ksh -v Q4 -e ${set_env}
  zonal_clim.ksh -v Q5 -e ${set_env}
EOF_LIST
/e/uhome/for1han/bin/pshell -p7 -f zonal.list

# --- Evolution plots ---                                                 
                                                                          
for var in T ; do                                                 
  ./lat-time.ksh -v ${var} -e ${set_env}                                  
done                                                                      
                                                                          
#==========                                                               
# Clean up                                                                
                                                                          
if [ -f ${set_env} ]; then                                                
  rm -f ${set_env}                                                        
fi                               
exit
