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

#_____________________________________________________________________________
# Script to construct an initial condition file for the HD model from a
# restart file of an existing experiment.
#
# The restart can be either a single NetCDF file or a multi-file restart directory.
# No remapping between different grids is performed.
#
# Requirements: cdo and nco
#_____________________________________________________________________________

set -eu

src_restart=/work/bm1235/k203123/nextgems_prefinal/experiments/ngc4008/work/ngc4008_restart_atm_20300101T000000Z.mfr
src_grid_label=r2b8_0033_0016
src_exp_id=ngc4008
src_date=20300101
output_hdstart_file=hdstart_${src_grid_label}_from_restart_${src_exp_id}_${src_date}.nc

CDO="cdo -s"

function finish {
  rm temp?_$$.nc *ncap* >& /dev/null
}
trap finish EXIT

if [[ -d $src_restart ]]; then
  cdo collgrid,hd_overlflow_res_box,hd_riverflow_res_box,hd_baseflow_res_box $src_restart/patch1_*.nc $(basename $src_restart .mfr).nc
  src_restart=$(basename $src_restart .mfr).nc
fi
ncrename -d cells,cell -d layers_1,bresnum -d layers_5,rresnum ${src_restart} temp1_$$.nc
ncwa -O -a y,time temp1_$$.nc temp2_$$.nc
ncatted -O -a cell_methods,,d,, temp2_$$.nc
${CDO} chname,hd_riverflow_res_box,FRFMEM,hd_baseflow_res_box,FGMEM,hd_overlflow_res_box,FLFMEM temp2_$$.nc temp3_$$.nc
ncap2 -O -s 'defdim("oresnum",$bresnum.size); FLFMEM_new[$oresnum,$cell]=0.0;FLFMEM_new(:,:)=FLFMEM(:,:);FLFMEM_new@long_name="content of the overflow reservoir"' \
            temp3_$$.nc temp4_$$.nc
ncks -x -v FLFMEM temp4_$$.nc temp5_$$.nc
ncrename -v FLFMEM_new,FLFMEM temp5_$$.nc ${output_hdstart_file}
