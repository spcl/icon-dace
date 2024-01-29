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
# Part 2 in MJO-Analysis Suite
#   * Loops through all variables
#   * Calls NCL-Skript 
#     "mjo_analysis_clivar_anomalies.ncl" 
#     and provides required information on
#     dimensions etc.
#   * NCL-Script:
#     - computes daily anomalies for all
#       relevant variables
#     - writes out anomaly fields as .nc in
#       data directoy
#     - unifies variable names
#     - creates anomaly plots 
#
#-------------------------------------------------------
# DWD, FE 13, Julia Keller, 02/2016
########################################################


for variable in U200 U850 V850 OLR TOT_PREC
do
  if [ "${variable}" == "TOT_PREC" ]
  then
    fileinpart=${variable}"_daily"
  else 
    fileinpart=${variable}"_"${dataext}"mean"
  fi
  fileoutpart=${variable}"_"${dataext}"anom"
  echo "Process "${variable}'"' 
  ncl 'var="'${variable}'"'  \
      'infile="'${filepath}${dataset}'_'${fileinpart}'.grb2"' \
      'outfile="'${filepath}${dataset}'_'${fileoutpart}'.nc"' \
      'plotfile="'${plotpath}'"' \
      'datainfo="'${dataset}'"'  \
      ymdStrt=${StrtDate} \
      ymdLast=${LastDate} \
      lats=${nooflats}    \
      lons=${nooflons}    \
      times=${nooftimes}  \
      mjo_analysis_clivar_anomalies.ncl


done

