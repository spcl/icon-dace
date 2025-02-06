#!/bin/bash

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

#########################################################
# Part 1 in MJO-Analysis Suite
#   * Uses CDO to compute seasonal mean for 
#     total precipitation and
#     low-level zonal wind (U at 850hPa)
#   * Calls NCL-Script 
#     "mjo_analysis_clivar_mean_state.ncl"
#     with required input information 
#     dimensions etc.
#   * NCL-Script plots seasonal means
#
#-------------------------------------------------------
# DWD, FE 13, Julia Keller, 02/2016
########################################################


for season in summer winter
do
  if [ ${season} == "summer" ]
  then
    cdo timmean -selmon,5,6,7,8,9,10 ${filepath}${dataset}'_TOT_PREC_daily.grb2' ${filepath}${dataset}'_TOT_PREC_nh_'${season}'.grb2'
    cdo timmean -selmon,5,6,7,8,9,10 ${filepath}${dataset}'_U850_'${dataext}'mean.grb2' ${filepath}${dataset}'_U850_nh_'${season}'.grb2'
  fi
  if [ ${season} == "winter" ]
  then
    cdo timmean -selmon,1,2,3,4,11,12 ${filepath}${dataset}'_TOT_PREC_daily.grb2' ${filepath}${dataset}'_TOT_PREC_nh_'${season}'.grb2'
    cdo timmean -selmon,1,2,3,4,11,12 ${filepath}${dataset}'_U850_'${dataext}'mean.grb2' ${filepath}${dataset}'_U850_nh_'${season}'.grb2'
  fi


  echo "Process "${season} 
  ncl 'infile1="'${filepath}${dataset}'_TOT_PREC_nh_'${season}'.grb2"' \
      'infile2="'${filepath}${dataset}'_U850_nh_'${season}'.grb2"' \
      'plotdir="'${plotpath}'"' \
      'datainfo="'${dataset}'"'  \
      'season="'${season}'"'    \
      lats=${nooflats}          \
      lons=${nooflons}          \
      mjo_analysis_clivar_mean_state.ncl

done



