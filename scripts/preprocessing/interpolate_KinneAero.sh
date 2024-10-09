#!/bin/ksh

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

#_____________________________________________________________________________
# Interpolate ozone data for ICON runs with transient ozone
# applicable for irad_o3 = 5
# the interpolated aerosol data needs to be linked to the run directory
#_____________________________________________________________________________
set -ex

# ICON atmosphere grid ID:
RES=R2B6_0024
# Begin and end year for interpolation:
BYEAR=2000
EYEAR=2000
# dwd, levante
site='dwd'
#site='levante'

#_____________________________________________________________________________
case $site in
  'dwd')
     DATADIR='/hpc/uwork/icon-sml/Ozone/pool/data/ECHAM6/input/r0006/T127/ozone'   # raw input data
     # Please save the interpolated data under icon-sml account and follow the data structure
     # If you don't have permission, please contact kristina.froehlich@dwd.de:
     OUTDIR='/hpc/uwork/icon-sml/Ozone/${RES}'                                     # output directory
     TARGETGRID='/hpc/rwork0/routfor/routfox/icon/grids/public/edzw/icon_grid_0012_R02B04_G.nc'
   ;;
   'levante')
     DATADIR='/pool/data/ECHAM6/input/r0006/T127/ozone'   # raw input data
     # Please save the interpolated data in ICON pool and follow the data structure
     # If you don't have permission, please contact daniel.reinert@dwd.de:
     RESDIR='/pool/data/ICON/grids/public/edzw/...'       # output directory
     TARGETGRID='/pool/data/ICON/grids/public/edzw/icon_grid_0012_R02B04_G.nc'
     source /sw/etc/profile.levante
     module load cdo/2.0.6-gcc-11.2.0
     module load nco/5.0.6-gcc-11.2.0
   ;;
esac

OUTPUTBASE='bc_ozone_historical'        #base name of the output file
INPUTBASE='T127_ozone_historical'
REMAP_WEIGHTS=wgtdis_T127_to_target.nc
SOURCEGRID=${DATADIR}/T127_ozone_historical_1850.nc

YEAR=$BYEAR
mkdir -p $OUTDIR
cd ${DATADIR}
cdo gencon,$TARGETGRID $SOURCEGRID $REMAP_WEIGHTS
cdo gendis,$TARGETGRID -random,t127grid $REMAP_WEIGHTS
while [ $YEAR -le $EYEAR ]; do
  IFILE=${INPUTBASE}_${YEAR}.nc
  OFILE=${OUTDIR}/${OUTPUTBASE}_${YEAR}.nc
  #cdo -f nc4 -P 8 remap,$TARGETGRID,$REMAP_WEIGHTS ${IFILE} ${OFILE}
  #cdo -f nc4c -P 8 remap,$TARGETGRID,$REMAP_WEIGHTS ${IFILE} ${OFILE}
  cdo -f nc5 -P 8 remap,$TARGETGRID,$REMAP_WEIGHTS ${IFILE} ${OFILE}
  YEAR=$(( YEAR + 1 ))
done
rm $REMAP_WEIGHTS

#case $site in
#  'dwd')
#     chmod -R o+rx /hpc/uwork/icon-sml/Ozone
#     chmod -R g+rx /hpc/uwork/icon-sml/Ozone
#  ;;
#esac

