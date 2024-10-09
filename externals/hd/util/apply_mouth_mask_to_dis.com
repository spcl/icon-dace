#!/bin/ksh
#
# apply_mouth_mask_to_dis.com - Applies the river mouth mask to discharges on the HD grid
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
# *** Script that applies the mouth mask to HD discharge on the HD grid. 
#     Hence, a file is generated that comprises discharge only at the river mouths
#     and (option) at internal sinks.
#
set -ex
#
# *********** Data for Program *************************
# *** This needs to be EDITED  *************************
EXP=7062018   # Exp. number of HD model run     
ISINK=0       # Regard sinks in addition to mouths: 1= yes
# HD parameter file
HDPARA=/work/gg0302/g260122/HD/input/hdpara_vs4d_euro5min.nc
# Directory with HD discharge
DATA=/work/gg0302/g260122/HD/output/$EXP       # Discharge data
YBEG=1979     # Start year of simulation  
YEND=2016     # End year of simulation  
ISEP=1        # HD output in separate files per year (No/yes = 0/1)
DNRIV=${DATA}/${EXP}_meanflow          # ${DNRIV}_YYYY.nc for ISEP=1

DOUT=/work/gg0302/g260122/HD/output/${EXP}/mouth    # Target directory for ocean data

# ******************************************************
set +e
mkdir $DOUT
set -e
cd $TMPDIR
cdo setvar,mask_mouth -eqc,0 -selvar,FDIR $HDPARA mask_mouth.nc
if [ $ISINK -eq 1 ]; then
  cdo -eqc,5 -selvar,FDIR $HDPARA mask_sink.nc
  cdo add mask_mouth.nc mask_sink.nc mask.nc
  cdo setvar,mask_mouth_and_sinks mask.nc mask_mouth.nc
  rm mask.nc
fi

if [ $ISEP -eq 0 ]; then
  cdo setgrid,mask_mouth.nc -selyear,$YBEG/$YEND ${DNRIV} dis.nc
  cdo -b 32 -f nc4 -z zip_2 ifthen mask_mouth.nc dis.nc ${DOUT}/${EXP}_dis_at_mouth_${YEAR}_${YBEG}-${YEND}.nc
  rm dis.nc
else
  YEAR=$YBEG
  while [ $YEAR -le $YEND ] ; do
    DNAM=${DNRIV}_${YEAR}.nc
    cdo setgrid,mask_mouth.nc -ifthen mask_mouth.nc ${DNAM} dis_at_mouth_${YEAR}.nc
    YEAR=`expr $YEAR + 1`
  done
  cdo -b 32 -f nc4 -z zip_2 -O mergetime dis_at_mouth_????.nc ${DOUT}/${EXP}_dis_at_mouth_${YBEG}-${YEND}.nc
fi
#

