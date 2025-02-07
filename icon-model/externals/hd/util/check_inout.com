#!/bin/ksh
#
# check_inout.com - In/Output check for HD model
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
# ********* In/Output check of HD-Modell - Test version

HDPARA=../hdpara.nc
DNDIS=7055203_meanflow_1940.nc
DNINP=../hdforcing.nc
#
# Discharge
cdo -s lec,0. -selvar,FDIR $HDPARA m0.nc
cdo -s eqc,5. -selvar,FDIR $HDPARA m5.nc
cdo -s -mul $DNDIS m0.nc t1.nc
cdo -s -mul $DNDIS m5.nc t2.nc
cdo -s add -fldsum t1.nc -fldsum t2.nc dglob.nc
echo 'Discharge' > out.log
##cdo infov -seltimestep,1/5 dglob.nc  
cdo -s outputf,%15.6f -seltimestep,1/5 dglob.nc >> out.log
#
# Runoff and Drainage
cdo -s mulc,0.001 -zonsum $DNINP zon.nc                # mm/s --> m/s 
cdo -s mul zon.nc -selvar,AREA $HDPARA t1.nc           # --> m^3/s
cdo -s fldsum t1.nc drun.nc
echo 'Surface runoff' >> out.log
cdo -s outputf,%15.6f -seltimestep,1/5 -selvar,runoff drun.nc >> out.log
echo 'Drainage' >> out.log
cdo -s outputf,%15.6f -seltimestep,1/5 -selvar,drainage drun.nc >> out.log
cat out.log
#
exit
