#!/usr/bin/env sh

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

set -x

#==============================================================================
# Usage:
#   MODEL=icon GRID=R2B04 OMIP_DIR=. FORCE=1 LEV=L20 ./preproc4phc.sh
#   MODEL=mpiom GRID=GR15 OMIP_DIR=. FORCE=1 LEV=L40 ./preproc4phc.sh
#==============================================================================
# CONFIGURATION
#==============================================================================
# input data
         OMIP_POOL_DIR='/pool/data/MPIOM/setup/omip_365_era15'
              OMIP_DIR=${OMIP_DIR:-${OMIP_POOL_DIR}}
            OMIP_FILES='[0-9,a-k,m-z]*.nc'
# input grids
 ICON_GRID_DIR_DEFAULT='/pool/data/ICON/ocean_data/ocean_grid'
         ICON_GRID_DIR=${ICON_GRID_DIR:-${ICON_GRID_DIR_DEFAULT}}

MPIOM_GRID_DIR_DEFAULT="/pool/data/MPIOM/${GRID}"
        MPIOM_GRID_DIR=${MPIOM_GRID_DIR:-${MPIOM_GRID_DIR_DEFAULT}}
#------------------------------------------------------------------------------
# compute $GRID_FILE/$GRID based on given value of GRID_FILE
case "${MODEL}" in
 icon)
   if [ "x$GRID_FILE" = 'x' ]; then
     GRID_FILE=$(ls ${ICON_GRID_DIR}/icon${GRID}*etop*planet.nc | head -1)
   else
     GRID=$(basename $GRID_FILE .nc)
   fi
   ;;
 mpiom)
   if [ "x$GRID_FILE" = 'x' ]; then
     GRID_FILE="${MPIOM_GRID_DIR}/${GRID}s.nc"
   else
     extname=$(echo $GRID_FILE | rev | cut -d '.' -f 1 | rev)
     GRID=$(basename $GRID_FILE .${extname})
   fi
   ;;
 *) echo "Unsupported model! Use 'icon' or 'mpiom'."; exit 1;;
esac
#==============================================================================
# output
              MODEL=${MODEL:-icon}
               GRID=${GRID:-R2B04}
             TARGET=${TARGET:-omip}
                LEV=${LEV:-forcing}
            TIMEAVG=${TIMEAVG:-daily}
         _test_file="${OMIP_DIR}/runoff.nc"
TARGET_MODEL_OUTPUT=${TARGET}_${MODEL}_${GRID}_${LEV}-${TIMEAVG}.nc
              MERGE=${MERGE:-0}
#==============================================================================
# internals
                CDO=${CDO:-cdo}
             CDODEV=${CDODEV:-cdo-dev}
            THREADS=${THREADS:-8}
              FORCE=${FORCE:-0}
#==============================================================================
# remapping setup
remapOperator_general=genbil
remapOperator_special=genbic
targetGrid=./cell_grid-${GRID}-${MODEL}.nc
case "${MODEL}" in
 icon)
   gridSelect=ifs2icon_cell_grid
   [[ ! -f ${targetGrid} ]] && ${CDO} -f nc -selname,${gridSelect} ${GRID_FILE} ${targetGrid}
   ;;
 mpiom)
   [[ ! -f ${targetGrid} ]] && ${CDO} -f nc -sethalo,1,1 -random,${GRID_FILE} ${targetGrid}
   ;;
 *) echo "Unsupported model! Use 'icon' or 'mpiom'."; exit 1;;
esac
# grids/weights for horizontal interpolation
targetWeight_general=./cell_weight-${remapOperator_general}-${TARGET}-to-${MODEL}-${GRID}.nc
targetWeight_special=./cell_weight-${remapOperator_special}-${TARGET}-to-${MODEL}-${GRID}.nc
#==============================================================================
# MAIN
#==============================================================================
# horiz. interpolation for ICON or MPIOM target grid
# some preprocessing:
# * if the targeWeight are not present, they have to be created by with the targeGrid
# * if the targeGrid is not present, it has to be selection from the file
#   provides with the GRID variable
#
# create the target land-sea-mask
lsmFile='lsm4omip.nc'
$CDO gtc,1.0e-4 ${OMIP_DIR}/land_sea_mask.ECMWF.nc ${lsmFile}

if test ! -f "${targetWeight_general}" -o "$FORCE" -eq 1 ;then
  if test ! -f "${targetGrid}" ; then
    echo "GRID variable has to been set correctly!"
    exit 1
  fi
  $CDO -P ${THREADS} ${remapOperator_general},${targetGrid} -fillmiss -ifnotthen ${lsmFile} -seltimestep,1 ${_test_file} ${targetWeight_general}
fi
if test ! -f "${targetWeight_special}" -o "$FORCE" -eq 1 ;then
  if test ! -f "${targetGrid}" ; then
    echo "GRID variable has to been set correctly!"
    exit 1
  fi
  $CDO -P ${THREADS} ${remapOperator_special},${targetGrid} -fillmiss -ifnotthen ${lsmFile} -seltimestep,1 ${_test_file} ${targetWeight_special}
fi
#==============================================================================
# jobs files for use of GNU parallel
[[ -f jobs ]] && rm jobs
# array for output files for later merge
typeset -A  oFiles
typeset -A _oFiles
# cleanup old intermediate stuff
[[ -f $(ls -1 remapped_*.nc) ]] && rm remapped_*nc
[[ -f $(ls -1 fillmiss_*.nc) ]] && rm fillmiss_*nc
# loop over source files
for file in $(ls ${OMIP_DIR}/${OMIP_FILES}); do
  fileBasename=$(basename $file)
  if echo ${fileBasename} | grep -E '(east_west|north_south)_stress' >/dev/null; then
    targetWeight=${targetWeight_special};
  else
    targetWeight=${targetWeight_general};
  fi

  _oFile=fillmiss_$fileBasename
   oFile=remapped_$fileBasename
  echo -n "$CDO settaxis,2001-01-01,12:00:00,1day -fillmiss -ifnotthen ${lsmFile}  $file ${_oFile};" >> jobs
  echo -n "$CDO remap,${targetGrid},${targetWeight} ${_oFile} ${oFile}" >> jobs
  echo "" >> jobs

  oFiles+=" $oFile"
 _oFiles+=" $_oFile"
done
# perform parallel processing of ${THREADS} processes
cat jobs | parallel -j ${THREADS}

# merge together
if [ $MERGE = 0 ]; then
  [[ -f ${TARGET_MODEL_OUTPUT} ]] && rm ${TARGET_MODEL_OUTPUT}
  $CDO merge ${oFiles} _${TARGET_MODEL_OUTPUT}
  #==============================================================================
  # time averaging
  case "${TIMEAVG}" in
    daily)
      mv _${TARGET_MODEL_OUTPUT} ${TARGET_MODEL_OUTPUT}
      ;;
    monthly)
      cdo monmean _${TARGET_MODEL_OUTPUT} ${TARGET_MODEL_OUTPUT}
      oFiles+=" _${TARGET_MODEL_OUTPUT}"
      ;;
    annual)
      cdo yearmean _${TARGET_MODEL_OUTPUT} ${TARGET_MODEL_OUTPUT}
      oFiles+=" _${TARGET_MODEL_OUTPUT}"
      ;;
    *)
      echo "Wrong value for TIMEAVG! Use only: daily, monthly or annual."
      exit 1
  esac
fi
#==============================================================================
# some postprocessing for ICON
if [ "x${MODEL}" = 'xicon' ] ; then
  ncrename -d Time,time -v Time,time ${TARGET_MODEL_OUTPUT}
fi
#==============================================================================
# clean up
for file in ${TEMPFILES} ${_oFiles} ${oFiles}; do
  [[ -f $file ]] && rm $file
done
#==============================================================================
