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

###################################
#
# Capsule Shell-Script to analyse
# the Madden-Julian-Oscillation
# in NWP Output or Analysis Data
# Tested for ICON-Clim Experiments
# and NOAA OLR/TRMM/ERA-I Analysis
#
# Predefine parameters for data in 
# upper part of script! 
# 
# Required data:
# - Lat/Lon data
# - 6/12/24h data output
# - Sample name for Input Data:
#   "ICON_80km_10yr_OLR_daymean.grb2"
# 
# 
# DWD, FE 13, Julia Keller, 02/2016
###################################

# Predefine Variables

filepath="/e/uwork/extjkell/icon4mjo/data/10yr_run/processed/"  # Path to input data
plotpath="/e/uwork/extjkell/icon4mjo/data/verification_clim_runs/MJO_EVALUATION/" # Path to directory for plots
dataset="ICON_80km_10yr"  	 # Experiment Identifier
dataext="day"		 # File Extension
nooflats=181		 # Number of Latitudes
nooflons=360		 # Number of Longitudes
nooftimes=3652	   	 # Number of Time Steps
StrtDate=20010101	 # First Date in Data
LastDate=20101231    # Last Date in Data
timeint=1.0			 # Time Interval in Days
samperday=1.0        # Samples per day
bandpassmin=20       # Start of bandpass-window (no of days)
bandpassmax=100      # End of bandpass-window (no of days)
bandpasswgt=201	     # Intervals of bandpass-window

# For most of the sub-scripts it is possible 
# to opt out parts of the analysis procedure
# -> Switch on "true"  = execute plot/analysis 
# -> Switch on "false" = skip plot/analysis 

# Steps in Filtering Routine
variplot="true"
bptimeplot="true"
#bpareaplot="false"
bphovmplot="true" 
spectrplot="true"
wavfrqplot="true"

# Steps in Correlation Routine
lagcorrplot="true"
xspectrplot="true"
univareofplot="true"
multivareofplot="true"

# It now calls the various shell scripts that
# preprocess the data and run the ncl part
# these shell scripts are sourced here and
# have thus access to all variables defined here
echo "----------Mean State------------"
# Compute mean state for summer and winter
. ./mjo_analysis_pt1_mean_state.sh  

echo "-----------Anomalies------------"
# Compute daily anomalies
. ./mjo_analysis_pt2_anomalies.sh 

echo "-----------Filtering------------"
# Apply bandpass filtering and create plots
. ./mjo_analysis_pt3_filtering.sh 

echo "----------Correlations----------"
# Compute various correlations
. ./mjo_analysis_pt4_correlations.sh 

echo "-----Space-Time-Spectra---------"
# Create Wheeler-Kiladis Space Time plots
. ./mjo_analysis_pt5_wk_spacetime.sh  

echo "----------RMM-Index-------------"
# Compute RMM-Index for various phases
. ./mjo_analysis_pt6_rmm-index.sh 

echo "-----------EOFs-----------------"
# Compute and plot EOF fields    
. ./mjo_analysis_pt7_eofs.sh 

