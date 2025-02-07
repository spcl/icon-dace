#!/bin/ksh
#
# prepare_hdforcing.ksh - Preparation of the HD model forcing data
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Author: Stefan Hagemann (Hereon, Germany)
# Contact: stefan.hagemann@hereon.de
#_________________________________________
#
#
# ******* Preparation of HD forcing data ***
# 
#   This file is sourced from within the HD run script run_hdmodel.ksh
#
#   Currently, the script is written in a way that it expects the following file names 
#   for the forcing data (Note that this can be modified in the script, too.):
#   NetCDF format: hdforcing.nc (DNRUN) with the variables runoff and drainage representing 
#	               surface runoff and drainage (subsurface runoff), respectively.
#   Service format: DNRUN (e.g. runoff.srv) for surface runoff 
#                   DNDR (e.g. drainage.srv) for drainage (subsurface runoff)   
#
#   Input parameter of script:
#         IFORCE, EXPINP, DNREMAP, GRID, YEAR, iform, DNRUN, DNDR
#         HD directories (already exported in run script using setuserdir.com
#
#   Output parameter using export: ndt_day, UFAK
#         ndt_day = Number of forcing time steps per day => model time step
#            UFAK = Unit factor that needs to be applied to the forcing data so that 
#                   the unit becomes [m/s].
#
DATA=${HDFORCING}/${EXPINP}
case ${IFORCE} in
   1  )  DATA=${HDFORCING}/$EXPINP         # Input data directory 
         UFAK=0.001
         ndt_day=1          # Number of time steps per day
         case $EXPINP in
           *     ) if [ -e $DNREMAP ] ; then 
                     echo 'Remap file exists for forcing ' $EXPINP 
                   else
                     cdo -f nc genycon,${HDMAIN}/grid/$GRID ${DATA}/${EXPINP}_daily_${YEAR}.nc $DNREMAP
                   fi
                   cdocom="$cdocom -remap,${HDMAIN}/grid/$GRID,$DNREMAP" 
                   if [ $iform -eq 0 ] ; then
                     cdo -f srv $cdocom -selvar,qs  ${DATA}/${EXPINP}_daily_${YEAR}.nc $DNRUN
                     cdo -f srv $cdocom -selvar,qsb ${DATA}/${EXPINP}_daily_${YEAR}.nc $DNDR
                     echo 'SRV forcing files created' 
                   elif [ $iform -eq 1 ] ; then
                     cdo $cdocom -selvar,qs,qsb  ${DATA}/${EXPINP}_daily_${YEAR}.nc $DNRUN 
                     ncrename -h -O -v qs,runoff $DNRUN 
                     ncrename -h -O -v qsb,drainage $DNRUN 
                     echo 'NetCDF forcing file created: ' $DNRUN 
                   fi ;;
         esac  ;;
   2  )  UFAK=0.001
         ndt_day=1          # Number of time steps per day
         if [ $iform -eq 0 ] ; then
           cdo -f srv $cdocom -remapnn,${HDMAIN}/grid/grid_0_5.txt ${DATA}/day_${EXPINP}_suru_${YEAR}.sz $DNRUN
           cdo -f srv $cdocom -remapnn,${HDMAIN}/grid/grid_0_5.txt ${DATA}/day_${EXPINP}_drain_${YEAR}.sz $DNDR 
           echo 'SRV forcing files created' 
         elif [ $iform -eq 1 ] ; then
           echo 'NetCDF forcing file operation needs to defined --> STOP' 
           exit
         fi ;;
   3  )  if [ $EXPINP -le 62003 ] ; then DATA=/work/mh0231/m214046/cclm/$EXPINP ; fi
         UFAK=2.777777E-7              # Input is mm/h
         ndt_day=24          # Number of time steps per day
         case $EXPINP in
           62001 ) GRIDCCLM=${HDMAIN}/grid/grid_cclm_cd3_011.txt ;;   # CoastDat3
           62002 ) GRIDCCLM=${HDMAIN}/grid/grid_cclm_cd3_022.txt ;;
           62003 ) GRIDCCLM=${HDMAIN}/grid/grid_cclm_cd3_044.txt ;;
           62007 ) GRIDCCLM=${HDMAIN}/grid/grid_cclm_openfred.txt ;;
           62008 ) GRIDCCLM=${HDMAIN}/grid/grid_cclm_cd2.txt  ;;      # CoastDat2
           62009 ) GRIDCCLM=${HDMAIN}/grid/grid_cclm_openfred17.txt ;; 
           62010 ) GRIDCCLM=${HDMAIN}/grid/grid_cclm_aus.txt 
                   UFAK=4.62963E-8         # Input: mm/6h --> 1/1000/21600
                   ndt_day=4  ;;           # Number of time steps per day
           6201[12345] ) GRIDCCLM=${HDMAIN}/grid/grid_cclm_aus.txt 
                   UFAK=1.157407E-8        # Input: mm/day --> 1/1000/86400
                   ndt_day=1  ;;           # Number of time steps per day
           62018 ) GRIDCCLM=${HDMAIN}/grid/grid_cclm_cd3.txt  ;;      # CoastDat3
           62030 ) GRIDCCLM=${HDMAIN}/grid/grid_cclm_rea6.txt  ;;     # REA6
           62031 ) GRIDCCLM=${HDMAIN}/grid/grid_icon_nukleus.txt ;;   # ICON-CLM Nukleus
           62032 ) ndt_day=4 ; UFAK=4.62963E-8         # Input: mm/6h --> 1/1000/21600
	           GRIDCCLM=${HDMAIN}/grid/grid_cclm_ha.txt  ;;       # Ha: CoastDat3 incl. Boundary zone
	   *     ) echo 'CCLM exp. nr. not defined --> Abbruch!' ; exit ;;
         esac
         echo "Unit Factor: " $UFAK
#
#        *** Interpolation
#        *** remapcon not working without origin and the rotated Pole coordinates in grid description
         if [ -e $DNREMAP ] ; then 
           echo 'Remap file exists for forcing ' $EXPINP 
         else
           cdo -f nc genycon,${HDMAIN}/grid/$GRID -setgrid,$GRIDCCLM ${DATA}/${EXPINP}_suru_${YEAR}.nc $DNREMAP
         fi
         cdocom="$cdocom -remap,${HDMAIN}/grid/$GRID,$DNREMAP" 

         if [ $iform -eq 0 ] ; then
           cdo -f srv $cdocom ${DATA}/${EXPINP}_suru_${YEAR}.nc $DNRUN
           cdo -f srv $cdocom ${DATA}/${EXPINP}_drain_${YEAR}.nc $DNDR 
           echo 'SRV forcing files created' 
         elif [ $iform -eq 1 ] ; then
           cdo -O merge -setvar,runoff ${DATA}/${EXPINP}_suru_${YEAR}.nc \
                        -setvar,drainage ${DATA}/${EXPINP}_drain_${YEAR}.nc trun_dr.nc 
##           cdo $cdocom -setvar,runoff ${DATA}/${EXPINP}_suru_${YEAR}.nc trun.nc
##           cdo $cdocom -setvar,drainage ${DATA}/${EXPINP}_drain_${YEAR}.nc tdr.nc 
           cdo $cdocom trun_dr.nc $DNRUN
           rm trun_dr.nc
           echo 'NetCDF forcing file created: ' $DNRUN 
         fi ;;
   4  )  UFAK=0.001          # Input is mm/s = kg m-2 s-1
         ndt_day=1          # Number of time steps per day
         case $EXPINP in
           6021? ) GRIDREG=${HDMAIN}/grid/grid_remo_e11.txt     ;;    # Remo Cordex Europe 0.11°
           *     ) echo 'REMO Exp. nr. not defined --> Abbruch!' ; exit ;;
         esac
         echo "Unit Factor: " $UFAK
#
#        *** Interpolation
#        *** remapcon not working without origin and the rotated Pole coordinates in grid description
         if [ -e $DNREMAP ] ; then 
           echo 'Remap file exists for forcing ' $EXPINP 
         else
           cdo -f nc genycon,${HDMAIN}/grid/$GRID -setgrid,$GRIDREG ${DATA}/${EXPINP}_suru_${YEAR}.nc $DNREMAP
         fi
         cdocom="$cdocom -remap,${HDMAIN}/grid/$GRID,$DNREMAP" 

         if [ $iform -eq 0 ] ; then
           cdo -f srv $cdocom ${DATA}/${EXPINP}_suru_${YEAR}.nc $DNRUN
           cdo -f srv $cdocom -sub ${DATA}/${EXPINP}_run_${YEAR}.nc ${DATA}/${EXPINP}_suru_${YEAR}.nc $DNDR 
           echo 'SRV forcing files created' 
         elif [ $iform -eq 1 ] ; then
           cdo sub ${DATA}/${EXPINP}_run_${YEAR}.nc ${DATA}/${EXPINP}_suru_${YEAR}.nc ${EXPINP}_drain_${YEAR}.nc
           cdo -O merge -setvar,runoff ${DATA}/${EXPINP}_suru_${YEAR}.nc \
                        -setvar,drainage ${EXPINP}_drain_${YEAR}.nc trun_dr.nc 
           cdo $cdocom trun_dr.nc $DNRUN
           rm trun_dr.nc
           echo 'NetCDF forcing file created: ' $DNRUN 
         fi ;;

   *  ) echo 'Forcing No. ' $IFORCE ' not defined --> STOP!' ; ndt_day=-9999 
        export ndt_day=$ndt_day  ; exit ;;
esac
export ndt_day=$ndt_day
export UFAK=$UFAK

