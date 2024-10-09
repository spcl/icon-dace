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
# Script to generate initial soil data from model output.
# The intention is to reduce equilibration time after model initialization.
#                                             Veronika Gayler, September 2023
#-----------------------------------------------------------------------------
set -e

output=/work/mj0060/m220053/seamless/icon-nwp_jsbach-main-vga/build_intel/experiments/vga0437/vga0437_lnd_2d_10900101T000000Z.nc
outdate=1100-01-01
outtime=00:00:00
ic_soil_from_gauss=/work/mj0060/m220053/ini_files/0012-0035/land/r0002/ic_land_soil_1850.nc
ic_soil_new=ic_land_soil_1850_vga0437.nc

CDO="cdo -s -f nc4 -b F64"  # tested with cdo-2.0.6

for var in hydro_wtr_soil_sl_box hydro_ice_soil_sl_box hydro_weq_snow_box; do
  $CDO selvar,${var} -seltime,${outtime} -seldate,${outdate} ${output} ${var}.nc
done

# Initial snow cover
$CDO setname,"snow" hydro_weq_snow_box.nc snow.nc

# Initial soil moisture - for each soil levels
$CDO setname,"layer_moist" -add hydro_wtr_soil_sl_box.nc hydro_ice_soil_sl_box.nc layer_moist.nc
lev_rst=($($CDO showlevel layer_moist.nc))
lev_ic=($($CDO showlevel -selvar,layer_moist ${ic_soil_from_gauss}))
$CDO splitlevel layer_moist.nc layer_moist_
for i in 0 1 2 3 4; do
  while [[ $(echo ${lev_rst[$i]} | wc -c ) < 7 ]]; do  lev_rst[$i]=0${lev_rst[$i]}; done
  $CDO setlevel,${lev_ic[$i]} layer_moist_${lev_rst[$i]}.nc layer_moist_$i.nc
done
$CDO -O merge layer_moist_[0-4].nc layer_moist.level.nc

# Initial snow moisture - vertically integrated, not used in ICON-Land
$CDO setname,"init_moist" -vertsum layer_moist.nc init_moist.nc
$CDO -O merge snow.nc layer_moist.level.nc init_moist.nc new_ic_soil_vars.nc

# Relace variables extracted from restart file
$CDO replace ${ic_soil_from_gauss} new_ic_soil_vars.nc ${ic_soil_new}

echo "New initial conditions file: ${ic_soil_new}"

# Clean up
rm hydro_weq_snow_box.nc hydro_ice_soil_sl_box.nc hydro_wtr_soil_sl_box.nc
rm layer_moist_*.nc
rm snow.nc layer_moist.nc layer_moist.level.nc init_moist.nc
rm new_ic_soil_vars.nc
