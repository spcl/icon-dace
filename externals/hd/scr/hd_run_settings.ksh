#!/bin/ksh
#
# hd_run_settings.ksh - Sub-script to run_hdmodel.ksh with the main settings for running the HD model
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
#
# ***** HD and forcing Experiment nos. and HD settings ********************************
#
EXPINP=57004                      # Exp. no. of forcing - used for CCLM, HydroPy, Remo
                                   #   ERA5 (55053/54) & JSBACH forcing (25288,25410)
#EXP=70${EXPINP}
EXP=7055204

typeset -Z4 YYYY
YYYY=1979        # First year of simulation
#
# *** Restart or Cold Start
if [ -e ${HDMAIN}/log/${EXP}.year ] ; then
  YYYY=`cat ${HDMAIN}/log/${EXP}.year`
  INEU=0
  cstart='Re-Start'
  if (( ${ICALL} == 1 )) ; then return ; fi
else
  INEU=1
  cstart='Cold-Start'
fi
YEND=1979       # Last year of simulation

IFORCE=1        # Forcing: 1 = HydroPy, 2 = JSBACH-PF, 3 = CCLM, 4 = REMO, 5=WRF, 6=ICON
                # Note that JSBACH forcing must be shifted by cdo remapnn,grid_0_5.txt
HDRES=2         # HD Resolution: 0=0.5 Grad, 1=5 Min, 2= Euro 5 Min with 0.5° or 5 Min. input
                # 3 = Australia, 4=SEA=South East Asia
FORCE_RES=99    # Forcing data resolution (original, without cdo), 99= any non-HD, (for IFORCE=1)
CFORM=nc        # Format of forcing files: 'srv' = Service Format (Default), 'nc' = NetCDCF
ICOUPLE=0       # Coupling type: 0=no, 1=no interpolation, 2=interpolation in HD
DNCOUPLE=${HDFILE}/nemo/hdcouple_hd_vs5_1_to_nemo_imode2.nc

IWORK=1            # Run time: 1=1 year, 2=1 month, 3=year with 30 day months
                   #           4=as 1 but final year with nday_final days
nday_final=212     # Jan-July: 90+91+31
ndate_end=20210731 # end date of run for IWORK=4
MM=01              # Start month
# BIAS CORRECTION
IBC_TYPE=0         # Bias correction type: 0=None, 1 = Mean Bias, 2 = Low, Mid and High Biases   
DN_BCPARA=${HDFILE}/biascor/bias_correction_parameter_7062052_ipol1_2009-2018.nc   # File with Bias correction parameters
IBC_WRITE=0        # Write bias corrected outflow 0/1 = No/Yes

RES_INP='05'       # Inputdata resolution (05, t106, 5min), after interpolation with cdo
case ${IFORCE} in
  1     ) if (( $HDRES != 0 )) && [ $EXPINP -ge 55012 ] ; then RES_INP='5min' ; fi ;;
  3 | 4 | 5 ) if (( $HDRES != 0 )) ; then RES_INP='5min' ; fi ;;
esac
DNREMAP="${EXPINP}_to_${RES_INP}"          # Remap file name
#
# Time steering of run using date_start and date_end - maybe edited for sub-annual simulations
case $IWORK in
  2 ) date_start=${YYYY}${MM}01 ; date_end=${YYYY}${MM}31  ;;  # Editing of last day required 
  3 ) date_start=${YYYY}${MM}01 ; date_end=${YYYY}1230  ;; 
  4 ) date_start=${YYYY}${MM}01 ; date_end=${YYYY}1231     # Editing of last month & day in final year required
      if (( $YYYY == $YEND )) ; then  
        date_end=$ndate_end
      fi ;;  
  * ) date_start=${YYYY}-${MM}-01 ; date_end=${YYYY}-12-31  ;; 
esac
#
# Log output for nhd_diag=99 - Ahr river
xlon1=7.042    # 6.9956   50.5075  --> 2245  474  (-> 217, 258) , i.e. 1 shift east: 7.042
xlat1=50.542
xlon2=7.125    # 7.22     50.55    --> 2246  474   , i.e. 1 shift west: 7.125
xlat2=50.542
cdocom="setmisstoc,0."
#
# User specific settings that will be put as attributes in the output file
# via namelist HDUSER_CTL in file namelist.hduser
HD_USER="Stefan Hagemann"
USER_EMAIL=stefan.hagemann@hereon.de
HD_CONT="${USER_EMAIL}, https://coastmod.hereon.de"
HD_INST="Helmholtz-Zentrum Hereon, Institute of Coastal Systems, Germany"
HD_INSTID="ROR: 03qjp1d79"
#
# *** Settings for HD model resolution
#
IZIP_START=1
case ${HDRES} in
   0  )  IMAP=0              # Remapping type (0=no, 1=vg, 2=05-->5min
         ulimit -s 102400
         IZIP_START=1
         DNPARA="05deg/hdpara_vs1_12.nc"
         nhd_diag=7   # Log output for Elbe river on 0.5°
         HDSTART="05deg/hdstart_05.nc"
#
#        *** Examples if restart files from previous runs are used for initialization.
         case $EXP in
           7025000 ) HDSTART="hdstart/icon/25000/hdrestart_${YYYY}.nc" ;;
           * ) if (( $YYYY > 1901 )) ; then HDSTART="hdstart/7055170/7055170_hdrestart_${YYYY}01.nc" ; fi ;;  # GWSP3 based
         esac

         GRID="grid_0_5.txt"
         DNMAS="05deg/masks_${RES_INP}.nc"
         EXE="hd_05.exe" ;;
   1  )  if [ $RES_INP = '05' ] ; then
           IMAP=2
         else
           IMAP=0
         fi

         ulimit -s 102400

         DNPARA="5min/hdpara_vs5_1.nc"
         HDSTART="5min/hdstart_5min.nc"
         nhd_diag=7   # Log output for Elbe river on 5 Min.
#
#        *** Examples if restart files from previous runs are used for initialization.
         if (( $YYYY > 1979 )) ; then 
            HDSTART="hdstart/7055116/7055116_hdrestart_${YYYY}01.nc" 
            IZIP_START=1
         elif (( $YYYY == 1940 )) ; then 
           HDSTART="5min/7057020_hdrestart_194001.nc"
         fi
         GRID="grid_5min.txt"
         DNMAS="5min/masks_${RES_INP}.nc"
         EXE="hd_5min.exe" ;;
   2  )  echo "HDRES = 2 = euro5min  --> European 5 Min. domain:  -11°W-69°E, 27-72°N"
         if [ $RES_INP = '05' ] ; then
           cdocom="${cdocom} -selindexbox,339,498,37,126"   # = on 5 min: 2029,2988,217,756
           IMAP=2
         else
###           cdocom="setmisstoc,0. -selindexbox,2029,2988,217,756"
           IMAP=0
         fi
         IZIP_START=1
         nhd_diag=7   # 7=Log output for Elbe river on 5 Min., 99 for specified coordinates
         ulimit -s 102400
         DNPARA="euro5min/hdpara_vs5_1_euro5min.nc"
         HDSTART="euro5min/hdstart_euro5min.nc"
#
#        *** Examples if restart files from previous runs are used for initialization.
         case $EXP in
           7062054 ) HDSTART="hdstart/7055192/7055192_hdrestart_${YYYY}01.nc" ;;    # CMOR daily time step
           70620?? | 706220? | 706021[03]) if (( $YYYY <= 1979 )) ; then 
                       HDSTART="hdstart/7062008/7062008_hdrestart_${YYYY}01.nc" 
                     elif (( $YYYY > 2019 )) ; then
                       HDSTART="hdstart/7055145/7055145_hdrestart_${YYYY}01.nc"
          	     else 
                       HDSTART="hdstart/7062029/7062029_hdrestart_${YYYY}01.nc"  # CD3 & HD Vs 5.1
                     fi ;;
           7055037 ) echo 'use default start file' ;;
           * ) if (( $FORCE_RES == 2 )) ; then
                 if (( $YYYY > 1978 )) ; then 
                   HDSTART="hdstart/7055183/7055183_hdrestart_${YYYY}01.nc" ; fi
               else
                 if (( $YYYY > 2019 )) ; then  HDSTART="hdstart/7055145/7055145_hdrestart_${YYYY}01.nc"
                 elif (( $YYYY > 1901 )) ; then 
                   HDSTART="hdstart/7055192/7055192_hdrestart_${YYYY}01.nc"  # GWSP3/WFDE5 based
                 fi  
               fi  ;;
         esac
         GRID="grid_euro5min.txt"
         DNREMAP="${EXPINP}_to_euro5min"
         DNMAS="euro5min/masks_euro${RES_INP}.nc"
         EXE="hd_5min.exe" ;;
   *  )  case ${HDRES} in
           3 ) echo "HDRES = 3 = aus5min   --> Australian 5 Min. domain: 112-154°E, 10-45°S"
##               cdocom="${cdocom} -selindexbox,585,668,201,270"   # = on 5 min: 3505,4008,1201,1620
               cdosel='-sellonlatbox,112,154,-45,-10'
               REG_TAG="aus"    
               xlon1=147.24362     # Log output for gridboxes if nhd_diag=99 
               xlat1=-19.75856
               xlon2=139.6157
               xlat2=-34.3509  
               ## DNPARA="aus/hdpara_vs4d_aus.nc"   # This was used before
               HDSTART_aus="aus/hdstart_aus.nc"
               if (( $YYYY == 1990 )) ; then 
                 HDSTART_aus="hdstart/hdstart_aus_7062010_199001.nc"
               elif (( $YYYY != 1979 )) ; then 
                 HDSTART_aus="hdstart/hdstart_aus_7055051_${YYYY}01.nc" 
               fi   ;;
           4 ) echo "HDRES = 4 = sea5min  --> South Asian 5 Min. domain:  89°E-154°E, 37°N-45°S"
               cdosel='-sellonlatbox,89,154,-45,37'
               REG_TAG="sea5min"    
               xlon1=105.8      # Mekong, Pakse, 105.8      15.1167 
               xlat1=15.1167
               xlon2=105.945    # Mekong, Stung Treng, 105.945    13.533
               xlat2=13.533   ;;
           * ) echo 'Resolution No. ' $HDRES ' not defined --> STOP!' ; exit ;;
         esac

         nhd_diag=99   # 7=Log output for Elbe river on 5 Min., 99 for specified coordinates
         if [ $RES_INP = '05' ] ; then
           cdocom="${cdocom} ${cdosel}"   
           IMAP=2
         else
           IMAP=0
         fi
         IZIP_START=1
         ulimit -s 102400

         DNPARA="5min/${REG_TAG}/hdpara_vs5_1_${REG_TAG}.nc"
         GRID="grid_${REG_TAG}.txt"
         DNREMAP="${EXPINP}_to_${REG_TAG}"

         # Start file
         if [ -e ${HDFILE}/5min/${REG_TAG}/hdstart_${REG_TAG}.nc ] ; then echo "hdstart_${REG_TAG}.nc exists" 
         else
           cdo $cdosel ${HDFILE}/5min/hdstart_5min.nc ${HDFILE}/5min/${REG_TAG}/hdstart_${REG_TAG}.nc
         fi
         HDSTART="5min/${REG_TAG}/hdstart_${REG_TAG}.nc"
         case ${HDRES} in
           4 ) if (( $YYYY > 1950 )) ; then 
                 HDSTART="hdstart/7055191/7055191_hdrestart_${YYYY}01.nc"  # GWSP3/WFDE5 based
               fi  ;;
         esac
  
         if [ -e ${HDFILE}/5min/${REG_TAG}/masks_${REG_TAG}.nc ] ; then echo "masks_${REG_TAG}.nc exists" 
         else
           cdo $cdosel ${HDFILE}/5min/masks_5min.nc ${HDFILE}/5min/${REG_TAG}/masks_${REG_TAG}.nc
         fi
         DNMAS="5min/${REG_TAG}/masks_${REG_TAG}.nc"   # 
         EXE="hd_5min.exe" ;;
esac
#
#
if (( ${IWORK} == 4 )) ; then
  if (( ${YYYY} < ${YEND} )) ; then
    IWORK=1
  fi
fi
echo "End of hd_run_settings.ksh reached"
#
