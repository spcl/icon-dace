#!/bin/ksh
#
# fhdbal.com - Water Balance calculation using HD output and input 
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
# ********* Wasserbilanzierung aus HD-Modell Output
# 
##EXP=7055171
EXP=7025000
#Y1=1980 ; Y2=1980
#Y1=1980 ; Y2=2019
Y1=2450 ; Y2=2450
#
DRUN=/scratch/g/g260122/tmp                # Run Directory
HDFILE=/work/gg0302/g260122/HD/input       # HD data Directory
RES_FOR=0               # Forcing resolution: 0=as HD, 1=ICON
#
case $EXP in
  ??55171  ) EXPFOR=57004 ; IRES=0 ; DT=1 ; nrf=4    # HD and riverflow steps per day
             DFORC=/work/gg0302/g260122/HD/forcing/$EXPFOR   # Forcing data directory
             DATA=/work/gg0302/g260122/HD/output/${EXP}      # HD output data directory
             ANNDIS=${DATA}/${EXP}/ann_${EXP}_1979-2019.nc   # Annual discharges
             ANNRUN=${DFORC}/ann_${EXPFOR}_1979-2019.nc ;;   # Annual runoffs
  ??25000  ) EXPFOR=25000 ; IRES=0 ; DT=4 ; nrf=1    # HD and riverflow steps per day
             DFORC=/work/gg0302/g260122/HD/forcing/icon/$EXPFOR   # Forcing data directory
             DATA=/scratch/g/g260122/hd/${EXP}/out     # HD run data directory
             ANNDIS=${DATA}/ann_${EXP}_${Y1}.nc       # Annual discharges
             ANNRUN=${DFORC}/ann_${EXPFOR}_${Y1}.nc          # Annual runoffs
             DNICON=/work/gg0302/g260122/data/icon/trang/icon_grid_0012_R02B04_G.nc
             ANNMAPRUN=${DATA}/ann_mapped_runoff_2450.nc
             RES_FOR=1  ;;
         * ) IRES=2 ; DT=1 ; nrf=48 ;;
esac
#
#
Y2P1=`expr ${Y2} + 1`
#
DNRES1=${EXP}_hdrestart_${Y1}01.nc
DNRES2=${EXP}_hdrestart_${Y2P1}01.nc
echo "restart files: $DNRES1 $DNRES2" 
#
case $IRES in
  0  ) HDPARA=${HDFILE}/05deg/hdpara_vs1_11.nc   ;;
  2  ) HDPARA=${HDFILE}/euro5min/hdpara_vs5_1_euro5min.nc   ;;
esac
#
cd $DRUN
if [ -e ${DATA}/${DNRES1}.gz ];then
  cp -p ${DATA}/${DNRES1}.gz .
  cp -p ${DATA}/${DNRES2}.gz .
  gunzip -f ${DNRES1}.gz
  gunzip -f ${DNRES2}.gz
elif [ -e ${DATA}/${DNRES1} ];then
  cp -p ${DATA}/${DNRES1} .
  cp -p ${DATA}/${DNRES2} .
else
  echo "Neither ${DNRES1} nor ${DNRES1}.gz exist --> Stop!" ; exit
fi
#
YEAR=$Y1
NGES=0       # Number of HD time steps
NDAY_GES=0
while [ $YEAR -le $Y2 ] ; do
  case $YEAR in
    190[48] | 192[048] | 194[048] | 196[048] | 198[048] | 200[048] | 202[048] ) NDAY=366  ;;
    191[26] | 193[26] | 195[26] | 197[26] | 199[26] | 201[26] ) NDAY=366  ;;
    * ) NDAY=365  ;;
  esac
  NDAY_GES=`expr $NDAY_GES + $NDAY`
  let "NSTEP = $DT * $NDAY"
  NGES=`expr $NGES + $NSTEP`
  YEAR=`expr $YEAR + 1`
done
let "NY = $Y2P1 - $Y1"

#
#  Unit OF, BF: " m^3/s*day = " --> XF / NDAY = " m^3/s"
#  Unit RF: m3 d s-1
#     XFvol = XF * 86400. [m³]
#     DXFvol = (XF2-XF1) * 86400 [m³]
#     Flow DXFvol/dt = DXFvol / (NGES*86400) [m³/s] = DXF/NGES 
cdo -s sub ${DNRES2} ${DNRES1} diff_restart.nc
cdo -s -divc,$NGES -varssum -selvar,FLFMEM,FGMEM diff_restart.nc diff_lgf.nc
cdo -s -divc,$NGES -varssum -selcode,711/715 diff_restart.nc diff_rf.nc

echo "No. in time period: Steps = $NGES, Days = $NDAY_GES, Years: $NY"
cdo -s fldsum diff_lgf.nc tglob.nc
cdo -s output tglob.nc > out.log
XVAL1_lgf=`cat out.log`
#
# *** Gesamtdifferenz der Speicher, und Inflow per Gridbox/4. for 1st?? and last time step.
# *** Faktor fuer m^3/s --> km^3/year: 365.25 * 86400 / 10^9
if [ $NY -eq 1 ] ; then 
  let "FAKY = $NDAY * 86400 * 1.E-9"
else
  let "FAKY = 365.25 * 86400 * 1.E-9"
fi
echo "Factor: $FAKY"
cdo output -mulc,$FAKY tglob.nc > out.log
XVAL2_lgf=`cat out.log`
echo "Global LGF storage loss/gain flow: " $XVAL2_lgf " km^3/a = " $XVAL1_lgf " m^3/s"
#cdo -s output -fldsum -mulc,$ufak t1.nc > out.log         # --> km² or m³/s

cdo -s fldsum diff_rf.nc tglob.nc
cdo -s output tglob.nc > out.log
XVAL1_rf=`cat out.log`
cdo -s output -mulc,$FAKY tglob.nc > out.log
XVAL2_rf=`cat out.log`
echo "Global RF storage loss/gain flow:  " $XVAL2_rf " km^3/a = " $XVAL1_rf " m^3/s"
#
# Discharge
cdo -s timmean -selyear,${Y1}/${Y2} $ANNDIS dis.nc
cdo -s lec,0. -selvar,FDIR $HDPARA m0.nc
cdo -s eqc,5. -selvar,FDIR $HDPARA m5.nc
cdo -s -mul dis.nc m0.nc t1.nc
cdo -s -mul dis.nc m5.nc t2.nc
cdo -s add -fldsum t1.nc -fldsum t2.nc dglob.nc
cdo -s output dglob.nc > out.log
XVAL3=`cat out.log`
cdo output -mulc,$FAKY dglob.nc > out.log
XVAL4=`cat out.log`
echo "Global discharge flow into ocean/sink: " $XVAL4 " km^3/a = " $XVAL3 " m^3/s"       
#
cdo -s mul dis.nc -eqc,-1 -selvar,FDIR $HDPARA t2.nc
cdo -s output -fldsum t2.nc > out.log
XVAL3b=`cat out.log`
cdo output -mulc,$FAKY -fldsum t2.nc > out.log
XVAL4b=`cat out.log`
echo "Discharge flow into open ocean:        " $XVAL4b " km^3/a = " $XVAL3b " m^3/s"       
#
# Total Runoff forcing
cdo -s timmean -selyear,${Y1}/${Y2} $ANNRUN run.nc
case $RES_FOR in
  1  ) cdo -s -w mul -mulc,0.001 run.nc -selvar,cell_area $DNICON t1.nc   # -->  m³/s
       if [ -e ${ANNMAPRUN} ];then
         cdo -s mulc,0.001 -zonsum $ANNMAPRUN tz.nc
         ncwa -O -a time,lon tz.nc zon.nc
         cdo -s mul zon.nc -selvar,AREA $HDPARA t2.nc 
         cdo -s mulc,0.001 -zonsum -mul $ANNMAPRUN -lec,0. -selvar,FDIR $HDPARA tz.nc
         ncwa -O -a time,lon tz.nc zon.nc
         cdo -s mul zon.nc -selvar,AREA $HDPARA toc.nc        # all ocean 
         cdo -s mulc,0.001 -zonsum -mul $ANNMAPRUN -eqc,-1 -selvar,FDIR $HDPARA tz.nc
         ncwa -O -a time,lon tz.nc zon.nc
         cdo -s mul zon.nc -selvar,AREA $HDPARA too.nc        # open ocean
       fi ;;
  *  ) cdo -s mulc,0.001 -zonsum run.nc zon.nc                # mm/s --> m/s 
       cdo -s mul zon.nc -selvar,AREA $HDPARA t1.nc ;;        # --> m^3/s
esac

cdo -s output -fldsum t1.nc > out.log         
XVAL5=`cat out.log`
cdo output -mulc,$FAKY -fldsum t1.nc > out.log
XVAL6=`cat out.log`
echo "Global total runoff:        " $XVAL6 " km^3/a = " $XVAL5 " m^3/s"
#
cdo -s output -fldsum t2.nc > out.log         
XVAL5b=`cat out.log`
cdo output -mulc,$FAKY -fldsum t2.nc > out.log
XVAL6b=`cat out.log`
echo "Global mapped total runoff: " $XVAL6b " km^3/a = " $XVAL5b " m^3/s"
# Check contribution over HD ocean points
cdo -s output -fldsum toc.nc > out.log         
XVAL5c=`cat out.log`
cdo output -mulc,$FAKY -fldsum toc.nc > out.log
XVAL6c=`cat out.log`
echo "Mapped runoff over ocean:   " $XVAL6c " km^3/a = " $XVAL5c " m^3/s"
# Check contribution over HD open ocean points
cdo -s output -fldsum too.nc > out.log         
XVAL5d=`cat out.log`
cdo output -mulc,$FAKY -fldsum too.nc > out.log
XVAL6d=`cat out.log`
echo "Mapped runoff - open ocean: " $XVAL6d " km^3/a = " $XVAL5d " m^3/s"
echo " "
#
# Water in restart flow fields
# Water in the FINFL is added to the system in the first river flow time step
# However, the riverfow input is already constrained to the river flow time step per day 
#   (i.e. HD time step) --> No division by nrf necessary
cdo -s -selvar,FINFL $DNRES1 tinf.nc           # m^3/s in 1 day for inflow
cdo -s output -divc,$NGES -fldsum tinf.nc > out.log     # TODO check
XVAL7a=`cat out.log`
cdo -s output -mulc,$FAKY -divc,$NGES -fldsum tinf.nc > out.log
XVAL7b=`cat out.log`
echo "       Riverflow field input: " $XVAL7b " km^3/a = " $XVAL7a " m^3/s"   

cdo -s -selvar,FINFL $DNRES2 tof.nc           # m^3/s in 1 day for outflow
cdo -s -mul tof.nc m0.nc t1.nc
cdo -s -mul tof.nc m5.nc t2.nc
cdo -s add t1.nc t2.nc t3.nc
cdo -s sub tof.nc t3.nc tout.nc    # no extra outflow at ocean/sink boxes as it is regarded in HD output calc. 
cdo -s output -divc,$NGES -fldsum tout.nc > out.log     # TODO check
XVAL9a=`cat out.log`
cdo -s output -mulc,$FAKY -divc,$NGES -fldsum tout.nc > out.log
XVAL9b=`cat out.log`
echo "Riverflow field inland output: " $XVAL9b " km^3/a = " $XVAL9a " m^3/s"   

# For the difference in inflows it must be I-O, not as for the storages (Send -Sanf) --> Factor -1
cdo -s mulc,-1. -selvar,FINFL diff_restart.nc tinf.nc           # m^3/s in 1 day for diff. in inflows
cdo -s output -divc,$NGES -fldsum tinf.nc > out.log    # TODO check
XVAL8a=`cat out.log`
cdo -s output -mulc,$FAKY -divc,$NGES -fldsum tinf.nc > out.log
XVAL8b=`cat out.log`
echo "Riverflow diff. I - O in restart inflows " $XVAL8b " km^3/a = " $XVAL8a " m^3/s"   
echo " "
#
# Summary
echo "No. of days in time period: $NDAY_GES,  No. of years: $NY" > hdbal.log
echo "Global LGF storage loss/gain flow:     " $XVAL2_lgf " km^3/a = $XVAL1_lgf m^3/s" >> hdbal.log
echo "Global RF storage loss/gain flow:      " $XVAL2_rf " km^3/a = $XVAL1_rf m^3/s" >> hdbal.log
echo "Riverflow diff. I - O in rest. inflows " $XVAL8b " km^3/a = $XVAL8a m^3/s" >> hdbal.log
echo "Global discharge flow into ocean/sink: " $XVAL4 " km^3/a = $XVAL3 m^3/s" >> hdbal.log     
echo "Global total runoff:                   " $XVAL6 " km^3/a = $XVAL5 m^3/s"  >> hdbal.log 
echo "-----------------------------------------"  >> hdbal.log
# DW = Inflow - Outflow - Dstorage + DInres = 0?
let "DW = $XVAL6 - $XVAL4 - $XVAL2_lgf - $XVAL2_rf + $XVAL8b"
echo "Inflow - Outflow - Dstorage = 0?:      " $DW " km^3/a"  >> hdbal.log
#
cat hdbal.log 
exit

