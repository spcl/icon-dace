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
# Part 4 in MJO-Analysis Suite
#   * Calls NCL-Skript
#     "mjo_analysis_clivar_correlations.ncl"
#     and provides required information on
#     dimensions etc.
#   * Logical Switches allow to compute/plot 
#     just parts of the procedure
#   * NCL Script
#     - Computes and plots Lag-Correlation
#     - Computes and plots Cross-Spectra
#     - Computes and plots Multivariate EOFs
#     - Computes and plots MJO Life Cycle 
#  
#-------------------------------------------------------
# DWD, FE 13, Julia Keller, 02/2016
########################################################


ncl 'precname="TOT_PREC"' \
    'olrname="OLR"'       \
    'u850name="U850"'     \
    'u200name="U200"'     \
    'v850name="V850"'     \
    'precfile="'${filepath}${dataset}'_TOT_PREC_'${dataext}'anom.nc"' \
    'olrfile="'${filepath}${dataset}'_OLR_'${dataext}'anom.nc"' \
    'u850file="'${filepath}${dataset}'_U850_'${dataext}'anom.nc"' \
    'u200file="'${filepath}${dataset}'_U200_'${dataext}'anom.nc"' \
    'v850file="'${filepath}${dataset}'_V850_'${dataext}'anom.nc"' \
    'outpath="'${filepath}'"' \
    'plotdir="'${plotpath}'"' \
    'datainfo="'${dataset}'"' \
    ymdStrt=${StrtDate}       \
    ymdLast=${LastDate}       \
    dayint=${timeint}         \
    bpmin=${bandpassmin}      \
    bpmax=${bandpassmax}      \
    bpwgt=${bandpasswgt}      \
    latmin=-30                \
    latmax=30                 \
    latS_IO=-10               \
    latN_IO=5                 \
    lonW_IO=75                \
    lonE_IO=100               \
    latS_bnd=-10              \
    latN_bnd=10               \
    lonW_bnd=80               \
    lonE_bnd=100              \
    neof=2                    \
    'lagcorr="'${lagcorrplot}'"'  \
    'xspectra="'${xspectrplot}'"' \
    'unieofs="'${univareofplot}'"'     \
    'multieofs="'${multivareofplot}'"'     \
    mjo_analysis_clivar_correlations.ncl


