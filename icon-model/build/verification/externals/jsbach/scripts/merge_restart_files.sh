#!/bin/sh

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

#-------------------------------------------------------------------------------
# Script to manipulate an ICON restart file, to import the land carbon and
# nlcc states of another experiment.
# All variables of the CARBON, DISTURB, ASSIMI and PPLCC processes, as well as
# all fraction variables will be replaced.
#
#  (Preliminary method to replace read_cpools and read_fpc of jsbach3)
#
# Based on proceeding of Susanne Fuchs, MPI-M    August 2022
#-------------------------------------------------------------------------------
set -e

# Restart file with equilibrated nlcc and carbon pools
src_restart=/work/mh0287/m300879/restart_files_for_exp_with_les/tra0308_restart_atm_26500101_renamed_dynveg.nc
# Restart file you want to modify
dst_restart=/work/mj0060/m220053/icon-esm/esm-merge-icon-2.4.6.1-dev/experiments/vga0403/restart/vga0403_restart_atm_18710101.nc
# Name of the restart file that will be generated
new_restart=vga0403-tra0308-26500101_restart_atm_18710101.nc

#-------------------------------------------------------------------------------
src_restart_name=${src_restart##*/} ; src_restart_name=${src_restart_name%%.nc}
dst_restart_name=${dst_restart##*/} ; dst_restart_name=${dst_restart_name%%.nc}

TRAVARS=$(cdo -showname ${src_restart} | tr -s [:blank:] '\n'  | grep '^carbon_\|^disturb_\|^assimi_\|^pplcc_\|^nlcc_\|^fract_' | sort -n | tr '\n' ',' | sed 's/,$//')

cdo selname,$TRAVARS ${src_restart} ${src_restart_name}_selected_vars.nc

# Adapt the source restart file variables to the destination grid sea land and
# glacier mask. This is needed, if src and dst experiments did not use identical
# ic and bc file revisions.
cdo -gtc,0.001 -sub -selname,fract_box  ${dst_restart} \
                    -selname,fract_glac ${dst_restart}   non_glacier_land.nc
cdo selname,$TRAVARS ${dst_restart} ${dst_restart_name}_selected_vars.nc
cdo ifthenelse non_glacier_land.nc \
    ${src_restart_name}_selected_vars.nc \
    ${dst_restart_name}_selected_vars.nc \
    ${src_restart_name}_selected_vars_on_dst_land.nc

cdo delete,name=$TRAVARS ${dst_restart} ${dst_restart_name}_without_vars_to_replace.nc

cdo merge ${dst_restart_name}_without_vars_to_replace.nc \
    ${src_restart_name}_selected_vars_on_dst_land.nc ${new_restart}

echo ""
echo " Manipulated restart file: ${new_restart}"
echo ""

# Clean up
rm -f non_glacier_land.nc
rm -f ${src_restart_name}_selected_vars.nc
rm -f ${dst_restart_name}_selected_vars.nc
rm -f ${src_restart_name}_selected_vars_on_dst_land.nc
rm -f ${dst_restart_name}_without_vars_to_replace.nc

# Test printings
echo ""
echo "-------- Test ----------"
testvar=carbon_c_sum_natural_ta_box
echo ${src_restart}:
cdo infov -selvar,${testvar} ${src_restart}
echo ${dst_restart}:
cdo infov -selvar,${testvar} ${dst_restart}
echo ${new_restart}:
cdo infov -selvar,${testvar} ${new_restart}
