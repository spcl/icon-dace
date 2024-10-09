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
# Part 6 in MJO-Analysis Suite
#   * Loop through years from Start- to
#     Last Date
#   * For each NH Summer/ NH Winter Season
#     it calls ncl-Skript
#     "mjo_analysis_clivar_rmm-index.ncl"
#     and provides required information on
#     dimensions etc.
#   * NCL Script
#     - Reads PC1/2 as created by  
#       mjo_analysis_clivar_correlations.ncl
#     - Creates RMM-Diagram for each Season
#
#-------------------------------------------------------
# DWD, FE 13, Julia Keller, 02/2016
########################################################

# Extract year from Strt-/LastDate    
yrstrt=${StrtDate:0:4}
yrlast=${LastDate:0:4}
yr=$yrstrt
echo $yrstrt
echo $yrlast

while [ $yr -lt ${yrlast} ]
do
  for seas in summer winter 
  do
    case ${seas} in
      "summer") 
             rmmStrt=${yr}0401
             rmmLast=${yr}0930
       ;;
      "winter") 
             rmmStrt=${yr}1001
             yr=`expr $yr + 1 `
             rmmLast=${yr}0331
       ;;

    esac
    echo $rmmStrt
    echo $rmmLast

    ncl 'infile="'${filepath}${dataset}'_MJO_PC_INDEX.nc"' \
        'plotdir="'${plotpath}'"' \
        ymdStrt=$rmmStrt  \
        ymdLast=$rmmLast  \
        'datainfo="'${dataset}'"' \
        mjo_analysis_clivar_rmm-index.ncl
        
  done

done



