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
# Part 3 in MJO-Analysis Suite
#   * Loops through all variables
#   * Calls NCL-Skript 
#     "mjo_analysis_clivar_filtering.ncl"
#     and provides required information on 
#     dimensions etc.
#   * Logical Switches allow to compute/plot 
#     just parts of the procedure
#   * NCL Script:
#     - Filters Daily Anomalies
#     - Determines seasonal variance and 
#       ratio of filtered/unfiltered var
#     - Plots timeseries of filtered area 
#       averaged data
#     - Plots Hovmoeller of filtered data
#     - Plots Lat_Lon Timeseries of filtered 
#       data
#     - Plots Spectra over particual regions
#  
#-------------------------------------------------------
# DWD, FE 13, Julia Keller, 02/2016
########################################################


for variable in OLR TOT_PREC U200 U850 V850  
do
  fileinpart=${variable}"_"${dataext}"anom"
  echo "Process "${variable} 
  ncl 'var="'${variable}'"' \
      'infile="'${filepath}${dataset}'_'${fileinpart}'.nc"' \
      'plotdir="'${plotpath}'"' \
      'datainfo="'${dataset}'"'  \
      ymdStrt=${StrtDate}   \
      ymdLast=${LastDate}   \
      dayint=${timeint}     \
      bpmin=${bandpassmin}  \
      bpmax=${bandpassmax}  \
      bpwgt=${bandpasswgt}  \
      yperhov=5             \
      latmin=-30            \
      latmax=30             \
      lonmin=60             \
      lonmax=300            \
      'vari="'${variplot}'"'      \
      'bptime="'${bptimeplot}'"'  \
      'bparea="'${bpareaplot}'"'  \
      'bphov="'${bphovmplot}'"'   \
      'spectra="'${spectrplot}'"' \
      'wavfreq="'${wavfrqplot}'"' \
      mjo_analysis_clivar_filtering.ncl

done


