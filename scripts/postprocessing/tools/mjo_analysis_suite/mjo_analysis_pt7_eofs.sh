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
# Part 7 in MJO-Analysis Suite
#   * Loops through all variables
#   * Calls NCL-Skript
#     "mjo_analysis_clivar_rmm-index.ncl"
#     and provides required information on
#     dimensions etc.
#   * NCL Script
#     - Computes and plots EOFs
#-------------------------------------------------------
# DWD, FE 13, Julia Keller, 02/2016
########################################################


for variable in  OLR TOT_PREC U200 U850 V850 
do
  fileinpart=${variable}"_"${dataext}"anom"
  echo "Process "${variable} 
  ncl 'var="'${variable}'"' \
      'infile="'${filepath}${dataset}'_'${fileinpart}'.nc"' \
      'plotdir="'${plotpath}'"' \
      'datainfo="'${dataset}'"' \
      ymdStrt=${StrtDate}   \
      ymdLast=${LastDate}   \
      bpmin=${bandpassmin}  \
      bpmax=${bandpassmax}  \
      bpwgt=${bandpasswgt}  \
      spd=${samperday}      \
      latmin=-30             \
      latmax=30             \
      neof=3                \
      mjo_analysis_clivar_eofs.ncl
done


