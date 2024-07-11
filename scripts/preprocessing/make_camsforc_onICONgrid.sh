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


#############################################
#  This file prepares the CAMS forecasted aerosols netcdf file for ICON with irad_aero = 8 
#  
#  The user needs to download CAMS data for 11 aerosols aermr01,...,aermr11 i.e. from MARS database
#  The data should be downloaded to one netCDF file that contain all aerosols and timesteps together.
#  This file should be named  $CAMSdir/CAMS_aero_${dathh}.nc where $CAMSdir is the location
#  of the CAMS files and $dathh is the date and hour of the initial forecast time.
#
#  The 3D pressure data of the CAMS forecast is also needed. Unfortunately, this data is not explicitly available in the model output
#  and needs to be calculated. Therefore the 2D data of the natural logarithm of surface pressure - lnsp needs to 
#  be downloaded also. This time the grib format is needed and also separated file for each timestep 
#  $CAMSdir/CAMS_lnsp_${dathh}_${n}.grb where "n" is the forecast hour (3 hours resolution).
#  The 3D pressure is calculated by compute_pressure_on_ml.py which is based on compute_geopotential_on_ml.py given in: 
#  https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height   
#
#  CDO, python3 and ecmwf-toolbox is needed for this script to work
#
#  history:
#  Harel Muskatel    2024-01  first version (E-mail: muskatelh@ims.gov.il)
#############################################

module load cdo
module load python3
module load ecmwf-toolbox

# grid name
GRIDNAME=R02B10_DOM01 # an example, edit here, edit here

# define range of forecast
frange=24 # an example, edit here

# define start date and hour yyyymmddhh
dathh=2024010100 # an example, edit here 

# define the location of CAMS files 
CAMSdir=/scratch/.../.../$dathh/CAMS_files # an example, edit here

# A loop on all forecast range to calculate 3D half-level pressure from 2D lnsp
for n in  `seq 0 3 $frange`; do
	python3 compute_pressure_on_ml.py $CAMSdir/CAMS_lnsp_${dathh}_${n}.grb -o $CAMSdir/CAMS_pres0_${dathh}_${n}.grb
done 

# A loop to change the variable name to "pres" and to create netCDF file
for n in  `seq 0 3 $frange`; do
	grib_set -s shortName=pres $CAMSdir/CAMS_pres0_${dathh}_${n}.grb $CAMSdir/CAMS_pres_${dathh}_${n}.grb
	cdo -f nc copy $CAMSdir/CAMS_pres_${dathh}_${n}.grb $CAMSdir/CAMS_pres_${dathh}_${n}.nc
done

# combine all forecast range to one file for pressure
cdo mergetime $CAMSdir/CAMS_pres_${dathh}_*.nc $CAMSdir/CAMS_pres_${dathh}.nc

# combine aerosols 3D data file with pressure 3D data file (total 12 variables)
cdo merge $CAMSdir/CAMS_pres_${dathh}.nc $CAMSdir/CAMS_aero_${dathh}.nc $CAMSdir/CAMS_forc_${GRIDNAME}.nc

# name of the CAMS forecast file including pressure
sourceFile=$CAMSdir/CAMS_forc_${GRIDNAME}.nc

# name of the target grid file
TARGETGRID=/.../.../CONST/${GRIDNAME}.nc  # an example, edit here

# name of the output file with the interpulated CAMS forecast on ICON grid
outFILE=$CAMSdir/CAMS_aero_${GRIDNAME}.nc

GRIDNUM=`cdo sinfov ${TARGETGRID} | grep nvertex=3 | awk '{print $1}'` 
cdo -s -r -P 4 remapbic,"${TARGETGRID}:${GRIDNUM}" ${sourceFile} t1.nc

cdo mul -gec,0.0 t1.nc t1.nc t2.nc
cdo add -mulc,0.0 -ltc,0.0 t1.nc t2.nc ${outFILE}
rm -rf t1.nc t2.nc

echo 'Done
