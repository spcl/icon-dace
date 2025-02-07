#!/bin/ksh
#
# check_yac_mapping_icon_to_hd.ksh - Generates water budgets output to check YAC mapping from ICON-A to HD
# 
# Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#

#DNICON=/work/gg0302/g260122/data/icon/trang/icon_grid_0030_R02B05_G.nc
DNICON=/work/gg0302/g260122/data/icon/trang/icon_grid_atmos.nc

HDPARA=/work/gg0302/g260122/HD/input/05deg/hdpara_vs1_11.nc
DNL=fr_land_plus_lake.nc   # Fractional Mask
cfmin=0.05                 # Trang's land fraction threshold
DNHD=icon_to_hd_mapped.nc
DNC=cdo_mapped.nc
DNREC=hd_receive.nc
#
cdo gec,$cfmin $DNL lsm_sent_by_trang.nc 

IMAP=0    # Mapping: 0: hd_sent mask, 1: run.nc
case $IMAP in
  0 ) DNI=lsm_sent_by_trang.nc     # 0/1 mask
      cvar='area' ; cunit='km²' ; ufak=1.E-6 ;;
  * ) DNI="-mul run.nc lsm_sent_by_trang.nc"  
      cvar='runoff' ; cunit='m³/s' ; ufak=1.E-3 ;;
esac

#
# ICON grid
cdo -s -w mul $DNI -selvar,cell_area $DNICON t1.nc        # --> m^2 or mm*m²/s
cdo -s output -fldsum -mulc,$ufak t1.nc > out.log         # --> km² or m³/s
XVAL1=`cat out.log`
case $IMAP in
  0 ) cdo -s -w mul $DNL -selvar,cell_area $DNICON t4.nc               # --> m^2 or mm*m²/s
      cdo -s -w mul -gtc,0. $DNL -selvar,cell_area $DNICON t3.nc  ;;
  1 ) cdo -s -w mul -mul run.nc $DNL -selvar,cell_area $DNICON t4.nc   # --> m^2 or mm*m²/s
      cdo -s -w mul run.nc -selvar,cell_area $DNICON t3.nc ;;
esac
cdo -s output -fldsum -mulc,$ufak t4.nc > out.log 
XVAL1b=`cat out.log`
cdo -s output -fldsum -mulc,$ufak t3.nc > out.log 
XVAL1c=`cat out.log`
#
# ********* HD grid ***************************
cdo -s zonsum $DNHD zon.nc                 
cdo -s mul zon.nc -selvar,AREA $HDPARA t2.nc           # --> m^2 or mm*m²/s
cdo -s output -fldsum -mulc,$ufak t2.nc > out.log 
XVAL2=`cat out.log`
#
cdo output -mulc,100. -div -fldsum t2.nc -fldsum t1.nc > out.log
XVAL3=`cat out.log`
#
# Multiply with hd_receive.nc-mask
cdo -s zonsum -mul $DNHD -gtc,0. $DNREC zon.nc
cdo -s mul zon.nc -selvar,AREA $HDPARA t2.nc           # --> m^2 or mm*m²/s
cdo -s output -fldsum -mulc,$ufak t2.nc > out.log 
XVAL6=`cat out.log`
#
# Check contribution over HD ocean points
cdo -s zonsum -mul $DNHD -lec,0. -selvar,FDIR $HDPARA zon.nc
cdo -s mul zon.nc -selvar,AREA $HDPARA t2.nc           # --> m^2 or mm*m²/s
cdo -s output -fldsum -mulc,$ufak t2.nc > out.log 
XVAL7=`cat out.log`
#
# Check contribution over HD open ocean points
cdo -s zonsum -mul $DNHD -eqc,-1 -selvar,FDIR $HDPARA zon.nc
##cdo -s zonsum -mul $DNHD -lec,-0.5 -selvar,FDIR $HDPARA zon.nc
cdo -s mul zon.nc -selvar,AREA $HDPARA t2.nc           # --> m^2 or mm*m²/s
cdo -s output -fldsum -mulc,$ufak t2.nc > out.log 
XVAL8=`cat out.log`
#
# Check contribution over HD land points
cdo -s zonsum -mul $DNHD -gec,1 -selvar,FDIR $HDPARA zon.nc
cdo -s mul zon.nc -selvar,AREA $HDPARA t2.nc           # --> m^2 or mm*m²/s
cdo -s output -fldsum -mulc,$ufak t2.nc > out.log 
XVAL9=`cat out.log`
#
# Check contribution from negative runoffs
cdo -s zonsum -mul $DNHD -ltc,0 $DNHD zon.nc
cdo -s mul zon.nc -selvar,AREA $HDPARA t2.nc           # --> m^2 or mm*m²/s
cdo -s output -fldsum -mulc,$ufak t2.nc > out.log 
XVAL10=`cat out.log`
#
# Check contribution from negative runoffs over the open ocean
cdo -s mul -ltc,0 $DNHD -eqc,-1 -selvar,FDIR $HDPARA t11.nc
cdo -s zonsum -mul $DNHD t11.nc zon.nc
cdo -s mul zon.nc -selvar,AREA $HDPARA t2.nc           # --> m^2 or mm*m²/s
cdo -s output -fldsum -mulc,$ufak t2.nc > out.log 
XVAL11=`cat out.log`
#
# cdo mapping
cdo remapycon,$DNHD $DNI $DNC
cdo -s zonsum $DNC zon.nc                 
cdo -s mul zon.nc -selvar,AREA $HDPARA t3.nc           # --> m^2 or mm*m²/s
cdo -s output -fldsum -mulc,$ufak t3.nc > out.log 
XVAL4=`cat out.log`
#
##cdo remapcon,$DNHD $DNL $DNC
##cdo -s zonsum $DNC zon.nc                 
##cdo -s mul zon.nc -selvar,AREA $HDPARA t5.nc           # --> m^2
##cdo -s output -fldsum -mulc,1.E-6 t5.nc > out.log 
##XVAL4b=`cat out.log`
#
cdo output -mulc,100. -div -fldsum t3.nc -fldsum t1.nc > out.log
XVAL5=`cat out.log`

echo "Gridbox $cvar on ICON Grid with fr_land > 0.05: " $XVAL1 $cunit 
echo "Mapped $cvar on HD Grid:                        " $XVAL2 $cunit 
echo "Mapped percentage: " $XVAL3 "%" 
echo "Cdo-Masked $cvar mapped on HD Grid:             " $XVAL4 $cunit 
echo "Mapped percentage: " $XVAL5 "%" 
echo "Land $cvar on ICON Grid multiplied with fr_land:" $XVAL1b $cunit 
echo "Total box $cvar with land on ICON Grid:         " $XVAL1c $cunit 
echo "Mapped $cvar on hd_receive mask:                " $XVAL6 $cunit 
echo "Mapped $cvar over all HD ocean boxes:           " $XVAL7 $cunit 
echo "Mapped $cvar over HD open ocean boxes:          " $XVAL8 $cunit 
echo "Mapped $cvar over HD land boxes:                " $XVAL9 $cunit 
echo "Mapped $cvar with negative runoff:              " $XVAL10 $cunit 
echo "Mapped $cvar with negative runoff over ocean:   " $XVAL11 $cunit 
##echo "Cdo-Land area mapped on HD Grid:               " $XVAL4b " km^2" 
rm t?.nc zon.nc out.log
