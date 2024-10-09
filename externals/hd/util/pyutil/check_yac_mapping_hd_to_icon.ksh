#!/bin/ksh
#
# check_yac_mapping_hd_to_icon.ksh - Generates water budgets output to check YAC mapping from HD to ICON-O
# 
# Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#

DNICON=/work/gg0302/g260122/data/icon/trang/icon_grid_ocean.nc

HDPARA=/work/gg0302/g260122/HD/input/05deg/hdpara_vs1_11.nc
IMAP=0    # Mapping: 0: hd_sent mask, 1: dis.nc
#
# Map mask or discharge
case $IMAP in
  0 ) DNHD=hd_sent.nc
      DNIO=hd_weights_on_icon.nc ;;
  * ) DNHD="-mul dis.nc -lec,0 -selvar,FDIR $HDPARA" 
      DNIO=** ;;
esac
#
# HD grid
cdo -s output -fldsum $DNHD > out.log    # m³/s
XVAL1=`cat out.log`
echo "HD volume flux on HD grid:                 " $XVAL1 " m^3/s" 

case $IMAP in
  0 ) cdo -s output -fldsum $DNIO > out.log   # Test with m³/s 
      XVAL2=`cat out.log`
      cdo output -mulc,100. -div -fldsum $DNIO -fldsum $DNHD > out.log
      XVAL3=`cat out.log`
      echo "Mapped volume flux on ICON grid            " $XVAL2 " m^3/s" 
      echo "Mapped percentage: " $XVAL3 "%" ;;

  1 ) echo 'Note that the sink contribution is missing in this calculation'
      cdo -s output -fldsum -mul dis.nc -eqc,5 -selvar,FDIR $HDPARA > out.log
      XVAL4=`cat out.log`
      echo "HD volume flux into sinks:                 " $XVAL4 " m^3/s" 
      cdo -s output -fldsum -mul dis.nc -eqc,-1 -selvar,FDIR $HDPARA > out.log
##      cdo -s output -fldsum -mul dis.nc -lec,-0.5 -selvar,FDIR $HDPARA > out.log
      XVAL5=`cat out.log`
      echo "HD volume flux into open ocean:            " $XVAL5 " m^3/s" 
      # neg. discharges
      cdo -s output -fldsum -mul $DNHD -ltc,0. dis.nc > out.log    # m³/s
      XVAL6=`cat out.log`
      echo "Negative ocean discharges on HD grid:     " $XVAL6 " m^3/s"  
      # neg. inflows
##      cdo -s output -fldsum -mul $DNIO -ltc,0. $DNIO > out.log   # Test with m³/s 
##      XVAL7=`cat out.log`
##      echo "Negative inflows on ICON grid:            " $XVAL7 " m^3/s"  ;;
esac


