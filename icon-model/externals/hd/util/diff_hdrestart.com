#!/bin/ksh
#
# diff_hdrestart.com - Calculate differences and related pseudo flows from two HD restart files 
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
# ********* Calculate difference and related pseudo flows from two HD restart files
#           Test version
# 
Y1=1861 ; Y2=1980
#Y1=1861 ; Y2=1890
#Y1=1951 ; Y2=1980
#Y1=1980 ; Y2=1980

#
nrf=4    # riverflow steps per day
DATA=/work/gg0302/g260122/data/icon/trang/couple        # Data directory with HD restart files
DRUN=/scratch/g/g260122/tmp                # Run Directory
HDFILE=/work/gg0302/g260122/HD/input       # HD data Directory
#
DNRES1=hdrestart_18610101T000000Z.nc
DNRES2=hdrestart_19810101T000000Z.nc
IWORK=1
ANNDIS=${DATA}/hdmeanflow.nc        # Annual discharges
ANNRUN=${DATA}/runoff_total.nc      # Annual total runoff on land
ANNFRESH=${DATA}/ann_FrshFlux_Runoff.nc # Annual freshwater inflow into ocean 

##ANNFRESH=${DATA}/ann_wrong_fresh.nc
IMR=2     # Runoff-Model: 1: HydroPy, 2: ICON
IMO=2     # Ocean Model: 0: no, 1:MPI-OM, 2: ICON-O
#
Y2P1=`expr ${Y2} + 1`
#
echo "'restart files: $DNRES1 $DNRES2" 
IRES=0
case $IMR in
   1 ) echo 'HydroPy runoff'  ;;
   2 ) DNICONA=/work/gg0302/g260122/data/icon/trang/icon_grid_0030_R02B05_G.nc 
       DNLAND=${DATA}/fr_land_plus_lake.nc  
       DNMISSO=${DATA}/hd_missing_in_icono.nc  ;;
esac
case $IMO in
   1 ) echo 'MPI-OM ocean model runoff not adapted up to now'  ;;
   2 ) DNICONO=/work/gg0302/g260122/data/icon/trang/icon_grid_0035_R02B06_O.nc  
       DNOCEAN=/work/gg0302/g260122/data/icon/trang/R2B6_0035/cell_sea_land_mask.nc  ;;
   * ) echo 'No ocean model runoff utilized'  ;;
esac
case $IRES in
  0  ) HDPARA=${HDFILE}/05deg/hdpara_vs1_11.nc   ;;
  2  ) HDPARA=${HDFILE}/euro5min/hdpara_vs5_1_euro5min.nc   ;;
esac
#
# Compute number of days
YEAR=$Y1
NGES=0
while [ $YEAR -le $Y2 ] ; do
  case $YEAR in
    180[48] | 182[048] | 184[048] | 186[048] | 188[048] ) NDAY=366  ;;
    181[26] | 183[26]  | 185[26]  | 187[26]  | 189[26]  ) NDAY=366  ;;
    190[48] | 192[048] | 194[048] | 196[048] | 198[048] | 200[048] | 202[048] ) NDAY=366  ;;
    191[26] | 193[26]  | 195[26]  | 197[26]  | 199[26]  | 201[26] ) NDAY=366  ;;
    * ) NDAY=365  ;;
  esac
  NGES=`expr $NGES + $NDAY`
  YEAR=`expr $YEAR + 1`
done
let "NY = $Y2P1 - $Y1"
#
cd $DRUN

#
#  Unit OF, BF: " m^3/s*day = " --> XF / NDAY = " m^3/s"
#  Unit RF: m3 d s-1
#     XFvol = XF * 86400. [m³]
#     DXFvol = (XF2-XF1) * 86400 [m³]
#     Flow DXFvol/dt = DXFvol / (NGES*86400) [m³/s] = DXF/NGES 
cdo -s sub ${DATA}/${DNRES2} ${DATA}/${DNRES1} diff_restart.nc
cdo -s -divc,$NGES -varssum -selvar,FLFMEM,FGMEM diff_restart.nc diff_lgf.nc
cdo -s -divc,$NGES -varssum -selcode,711/715 diff_restart.nc diff_rf.nc

echo "No. of days in time period: $NGES,  No. of years: $NY"
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
#
cdo -s fldsum diff_rf.nc tglob.nc
cdo -s output tglob.nc > out.log
XVAL1_rf=`cat out.log`
cdo -s output -mulc,$FAKY tglob.nc > out.log
XVAL2_rf=`cat out.log`
echo "Global RF storage loss/gain flow:  " $XVAL2_rf " km^3/a = " $XVAL1_rf " m^3/s"
#
# Discharge
if [ $IWORK -eq 1 ] || [ $IWORK -eq 3 ]; then 
  cdo -s timmean -selyear,${Y1}/${Y2} -selvar,friv $ANNDIS dis.nc
  cdo -s lec,0. -selvar,FDIR $HDPARA m0.nc
  cdo -s eqc,5. -selvar,FDIR $HDPARA m5.nc
  cdo -s -mul dis.nc m0.nc t1.nc
  cdo -s -mul dis.nc m5.nc t2.nc
  cdo -s add t1.nc t2.nc dglob.nc
###     cdo -s mul dglob.nc $DNMISSO t3.nc ; mv t3.nc dglob.nc   # Only regard missing HD in ICON-O boxes

  cdo -s output -fldsum dglob.nc > out.log
  XVAL3=`cat out.log`
  echo "Global discharge flow into ocean/sink: " $XVAL3 " m^3/s"        # m^3/s
  cdo -s output -mulc,$FAKY -fldsum dglob.nc > out.log
  XVAL4=`cat out.log`
  echo "Global discharge flow into ocean/sink: " $XVAL4 " km^3/a"       
fi
#
# Total Runoff forcing
if [ $IWORK -eq 1 ] || [ $IWORK -eq 4 ]; then 
  case $IMR in
  1 ) cdo -s timmean -selyear,${Y1}/${Y2} $ANNRUN run.nc
      cdo -s mulc,0.001 -zonsum run.nc zon.nc                # mm/s --> m/s 
      cdo -s mul zon.nc -selvar,AREA $HDPARA t1.nc           # --> m^3/s
      cdo -s output -fldsum t1.nc > out.log  ;;
  2 ) cdo -s timmean -selyear,${Y1}/${Y2} -shifttime,-1y $ANNRUN run.nc
      cdo -s divc,86400 -divc,365.25 -mulc,0.001 run.nc r1.nc  # mm/a -> m/s
      cdo -s -w mul r1.nc -selvar,cell_area $DNICONA t1.nc        # --> m^3/s
      cdo -s output -fldsum t1.nc > out.log  ;;
  esac

  XVAL5a=`cat out.log`
  cdo output -mulc,$FAKY -fldsum t1.nc > out.log
  XVAL5b=`cat out.log`
  echo "  Global total runoff all: " $XVAL5a " m^3/s    = " $XVAL5b " km^3/a" 
  case $IMR in
    2) cdo -s -w mul t1.nc $DNLAND t2.nc        # Pay regard to land+lake fractions
       cdo -s output -fldsum t2.nc > out.log 
       XVAL5a_cland=`cat out.log`
       cdo output -mulc,$FAKY -fldsum t2.nc > out.log
       XVAL5b_cland=`cat out.log`
       echo "Global total runoff cland:  " $XVAL5a_cland " m^3/s    = " $XVAL5b_cland " km^3/a"

#       cdo -s -w mul t1.nc -gtc,0.5 $DNLAND t3.nc        # Trangs setup v2
       cdo -s -w mul t1.nc -gtc,0.05 $DNLAND t3.nc        # Trangs setup v3
       cdo -s output -fldsum t3.nc > out.log 
       XVAL5a_tr=`cat out.log`
       cdo output -mulc,$FAKY -fldsum t3.nc > out.log
       XVAL5b_tr=`cat out.log`
       echo "Global total runoff Trang: " $XVAL5a_tr " m^3/s    = " $XVAL5b_tr " km^3/a"

       cdo -s -w mul t3.nc $DNLAND t4.nc        # Trangs setup & fr_land
       cdo -s output -fldsum t4.nc > out.log 
       XVAL5a_tr_cl=`cat out.log`
       cdo output -mulc,$FAKY -fldsum t4.nc > out.log
       XVAL5b_tr_cl=`cat out.log`
       echo "Global total runoff Tr_cl: " $XVAL5a_tr_cl " m^3/s    = " $XVAL5b_tr_cl " km^3/a" 

       # Global runoff mapped on HD
#       DNR_HD=/work/gg0302/g260122/data/icon/trang/couple/ann_runoff_1861-1980_mapped_to_hd.nc
       DNR_HD=/work/gg0302/g260122/data/icon/trang/couple/icon_to_hd_mapped.nc
       if [ -e ${DNR_HD} ];then
         cdo -s divc,86400 -divc,365.25 -mulc,0.001 -zonsum ${DNR_HD} zon5.nc   # mm/a --> m/s 
         cdo -s mul zon5.nc -selvar,AREA $HDPARA t5.nc                          # --> m^3/s
         cdo -s output -fldsum t5.nc > out.log
         XVAL5a_onhd=`cat out.log`
         cdo output -mulc,$FAKY -fldsum t5.nc > out.log
         XVAL5b_onhd=`cat out.log`
         echo "Global total runoff on HD: " $XVAL5a_onhd " m^3/s    = " $XVAL5b_onhd " km^3/a"
       fi ;;
  esac
  echo " "
fi
# Water in restart flow fields
# Water in the FINFL is added to the system in the first river flow time step
# However, the riverfow input is already constrained to the river flow time step per day 
#   (i.e. HD time step) --> No division by nrf necessary
if [ $IWORK -eq 1 ] || [ $IWORK -eq 2 ]; then 
  cdo -s -selvar,FINFL ${DATA}/$DNRES1 tinf.nc           # m^3/s in 1 day for inflow
  cdo -s output -divc,$NGES -fldsum tinf.nc > out.log
  XVAL7a=`cat out.log`
  cdo -s output -mulc,$FAKY -divc,$NGES -fldsum tinf.nc > out.log
  XVAL7b=`cat out.log`
  echo "       Riverflow field input: " $XVAL7b " km^3/a = " $XVAL7a " m^3/s"   

  cdo -s -selvar,FINFL ${DATA}/$DNRES2 tof.nc           # m^3/s in 1 day for outflow
  cdo -s -mul tof.nc m0.nc t1.nc
  cdo -s -mul tof.nc m5.nc t2.nc
  cdo -s add t1.nc t2.nc t3.nc
  cdo -s sub tof.nc t3.nc tout.nc    # no extra outflow at ocean/sink boxes as it is regarded in HD output calc. 
  cdo -s output -divc,$NGES -fldsum tout.nc > out.log
  XVAL9a=`cat out.log`
  cdo -s output -mulc,$FAKY -divc,$NGES -fldsum tout.nc > out.log
  XVAL9b=`cat out.log`
  echo "Riverflow field inland output: " $XVAL9b " km^3/a = " $XVAL9a " m^3/s"   

# For the difference in inflows it must be I-O, not as for the storages (Send -Sanf) --> Factor -1
  cdo -s mulc,-1. -selvar,FINFL diff_restart.nc tinf.nc           # m^3/s in 1 day for diff. in inflows
  cdo -s output -divc,$NGES -fldsum tinf.nc > out.log
  XVAL8a=`cat out.log`
  cdo -s output -mulc,$FAKY -divc,$NGES -fldsum tinf.nc > out.log
  XVAL8b=`cat out.log`
  echo "Riverflow diff. I - O in restart inflows " $XVAL8b " km^3/a = " $XVAL8a " m^3/s"   
fi
echo " "
#
# Freshwater inflow into ocean
case $IMO in
  2 ) cdo -s timmean -selyear,${Y1}/${Y2} $ANNFRESH fresh.nc  # m/s  
      cdo -s mul fresh.nc -selvar,cell_area $DNICONO t1.nc        # --> m^3/s
      cdo -s div t1.nc $DNOCEAN t2.nc ; mv t2.nc t1.nc
      cdo -s output -fldsum t1.nc > out.log  ;;
  * ) echo " " ;;
esac
if [ $IMO -gt 1 ] ; then 
  XVAL6a=`cat out.log`
  cdo output -mulc,$FAKY -fldsum t1.nc > out.log
  XVAL6b=`cat out.log`
  echo "Global freswater inflow: " $XVAL6a " m^3/s = " $XVAL6b " km^3/a"        # 
  echo " "
fi
#
# Summary
echo "No. of days in time period: $NGES,  No. of years: $NY" > hdbal.log
echo "Global LGF storage loss/gain flow:     " $XVAL2_lgf " km^3/a = $XVAL1_lgf m^3/s" >> hdbal.log
echo "Global RF storage loss/gain flow:      " $XVAL2_rf " km^3/a = $XVAL1_rf m^3/s" >> hdbal.log
echo "Riverflow diff. I - O in rest. inflows " $XVAL8b " km^3/a = " $XVAL8a " m^3/s" >> hdbal.log

if [ $IWORK -eq 1 ] || [ $IWORK -eq 3 ]; then
  echo "Global discharge flow into ocean/sink: " $XVAL4 " km^3/a = $XVAL3  m^3/s" >> hdbal.log     
fi
if [ $IMO -gt 1 ] ; then 
  echo "Global freswater inflow:               " $XVAL6b " km^3/a = $XVAL6a  m^3/s" >> hdbal.log 
fi
if [ $IWORK -eq 1 ] || [ $IWORK -eq 4 ]; then
  echo "Global total runoff:                   " $XVAL5b " km^3/a = $XVAL5a  m^3/s" >> hdbal.log 
fi
echo "-----------------------------------------"  >> hdbal.log
# DW = Inflow - Outflow - Dstorage + DInres = 0?
case $IWORK in
  1 ) let "DW = $XVAL5b - $XVAL4 - $XVAL2_lgf - $XVAL2_rf + $XVAL8b"  ;;
  2 ) let "DW = $XVAL8b - $XVAL2_lgf - $XVAL2_rf"  ;;
  3 ) let "DW = -1 * $XVAL4"  ;;
  4 ) let "DW = $XVAL5b "  ;;
  * ) echo "DW undefined --> Exit" ; exit  ;;
esac
echo "Inflow - Outflow - Dstorage = 0?:      " $DW " km^3/a"  >> hdbal.log
#
cat hdbal.log 
exit

