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

#########################################################
# Part 5 in MJO-Analysis Suite
#   * Loops through all variables
#   * Calls NCL-Skript 
#     "mjo_analysis_clivar_correlations.ncl"
#     and provides required information on
#     dimensions etc.
#   * NCL Script
#     - Opens FULL data (not anomalies) 
#  
#-------------------------------------------------------
# DWD, FE 13, Julia Keller, 02/2016
########################################################



for var in  TOT_PREC OLR U200 U850 V850 
do
  if [ "${var}" == "U200" ]
  then
    variable=${var}
    varname="UWND"
    fileinpart=${var}"_daymean"
  elif [ "${var}" == "U850" ]
  then
    variable=${var}
    varname="UWND"
    fileinpart=${var}"_daymean"
  elif [ "${var}" == "V850" ]
  then
    variable=${var}
    varname="VWND"
    fileinpart=${var}"_daymean"
  elif [ "${var}" == "TOT_PREC" ]
  then
    variable=${var}
    varname="TOT_PREC"
    fileinpart=${var}"_daily"
  elif [ "${var}" == "RAIN_CON" ]
  then
    variable="PRECT"
    varname="RAIN_CON"
    fileinpart=${var}"_daily"
  elif [ "${var}" == "RAIN_GSP" ]
  then
    variable="PRECT"
    varname="RAIN_GSP"
    fileinpart=${var}"_daily"
 else
    varname=${var}
    fileinpart=${var}_"daymean"
  fi
  
  echo "Process "${var}
  ncl 'varname="'${var}'"' \
      'infile="'${filepath}${dataset}'_'${fileinpart}'.grb2"' \
      'plotfile="'${plotpath}'"'   \
      'datainfo="'${dataset}'"'    \
      lats=${nooflats}   \
      lons=${nooflons}   \
      times=${nooftimes} \
      latmin=-10         \
      latmax=10          \
      nDayWin=96         \
      nDaySkip=-30       \
      nSpD=${samperday} \
      mjo_analysis_clivar_wk_spacetime.ncl

done


