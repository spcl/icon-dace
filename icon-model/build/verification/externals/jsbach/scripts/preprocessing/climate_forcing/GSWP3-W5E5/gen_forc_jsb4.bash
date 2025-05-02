#!/bin/bash

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2024, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

set -e
# -----------------------------------------------------------------------------------------------------
# This program prepares jsbach standalone data from observations for use in ICON-LAND. It generates a 
# the climate*.nc file at given resolution and time step. Start in target directory.
# -----------------------------------------------------------------------------------------------------
#
cdo="cdo -s"

if [ "$1" == "-h" -o "$1" == "--help" -o ${#@} -eq 0 ] ; then
  echo ""
  echo "Usage: `basename $0` [-h] -d [gswp3-w5e5] input.nc"
  echo "with:"
  echo "  -d, --data   choose dataset: [gswp3-w5e5]"
  echo "optionally provide:"
  echo "  -g, --grid   target grid file, e.g. icon_grid_0043_R02B04_G.nc"
  echo "  -s, --step   target time step [daily]"
  echo ""
  exit 0
fi

CONTACT=$USER # Replace with contact information
INST="?"      # Replace with organisation

PARACHECK=true
DATA="missing"
GRID="/pool/data/ICON/grids/public/mpim/0043/icon_grid_0043_R02B04_G.nc"
STEP="daily"

while $PARACHECK; do
  case $1 in
    -d | --data ) DATA=$2 && shift 2;;
    -g | --grid ) GRID=$(ls $2) && shift 2;;
    -s | --step ) STEP=$2 && shift 2;;
    *           ) PARACHECK=false;;
  esac
done

[[ "$DATA" == "missing"   ]] && echo "ERROR: Provide dataset name with -d option" && exit 1

# Catch not implemented time step request
if   [ $STEP == "daily" ]; then
  true
else
  echo "Timestep $STEP not available for $DATA" && exit 1
fi

# Create target directory
GTAG="$(basename $GRID | cut -d"_" -f4)/$(basename $GRID | cut -d"_" -f3)"
DTRG=$(echo "${GTAG}/$DATA" | tr [a-z] [A-Z])
test -d $DTRG || mkdir -p $DTRG

# -----------------------------------------------------------------------------------------------------
# Set variables depending on dataset
case $DATA in
  gswp3-w5e5 )
    VERSION="20211021; $(date)"
    REFERENCE="The GSWP3-W5E5 dataset is based on GSWP3 v1.09 (Kim 2017) and W5E5 v2.0 (Cucchi et al. 2020, Lange et al. 2021) and provided by ISIMIP (https://doi.org/10.48364/ISIMIP.982724)";;
  *     )
    echo "*** ERROR: Unknown dataset $DATA" && exit 1 ;;
esac

# -----------------------------------------------------------------------------------------------------
# WRITE compute script depending on time resolution
JSBFORC="convert_${DATA}_${STEP}.txt"
if [ ! -f $JSBFORC ]; then
  case $DATA in
    gswp3-w5e5)
      echo "longwave  = rlds;"            >> ${JSBFORC}
      echo "shortwave = rsds;"            >> ${JSBFORC}
      echo "precip    = pr;"              >> ${JSBFORC}
      echo "wspeed    = sfcwind;"         >> ${JSBFORC}
      echo "tmin      = tasmin - 273.15;" >> ${JSBFORC}
      echo "tmax      = tasmax - 273.15;" >> ${JSBFORC}
      echo "qair      = huss;"            >> ${JSBFORC}
      ;;
    *     )
      echo "*** ERROR: Unexpected dataset $DATA" && exit 1 ;;
  esac
fi

PARTAB="partab.txt"
if [ ! -f $PARTAB ]; then
cat > $PARTAB << EOF
&parameter
  name = longwave
  standard_name = surface_downwelling_longwave_flux_in_air
  long_name = "Surface Downwelling Longwave Radiation"
  units = "W m-2"
/
&parameter
  name = shortwave
  standard_name = surface_downwelling_shortwave_flux_in_air
  long_name = "Surface Downwelling Shortwave Radiation"
  units = "W m-2"
/
&parameter
  name = precip
  standard_name = precipitation_flux
  long_name = "Precipitation"
  units = "kg m-2 s-1"
/
&parameter
  name = wspeed
  standard_name = wind_speed
  long_name = "Near-Surface Wind Speed"
  units = "m s-1"
/
&parameter
  name = tmin
  standard_name = air_temperature
  long_name = "Daily Minimum Near-Surface Air Temperature"
  units = "degC"
/
&parameter
  name = tmax
  standard_name = air_temperature
  long_name = "Daily Maximum Near-Surface Air Temperature"
  units = "degC"
/
&parameter
  name = qair
  standard_name = specific_humidity
  long_name = "Near-Surface Specific Humidity"
  units = "kg kg-1"
/
EOF
fi

# -----------------------------------------------------------------------------------------------------
# Loop over forcing data files to get variables and years
unset VARLIST TIMLIST
DATAFILES=$@
for FILE in $DATAFILES ; do
  case $DATA in
    gswp3-w5e5 )
      # Provided in decadal files
      ID=$(basename $FILE)
      VAR=$(echo "print('$ID'.split('_')[2])" | python)
      TIM=$(echo "print('$ID'.split('_')[5])" | python)_$(echo "print('$ID'.split('_')[6].split('.')[0])" | python);;
    * ) echo "*** ERROR: Unexpected dataset $DATA" && exit 1 ;;
  esac
  [[ $VARLIST =~ (^|[[:space:]])"$VAR"($|[[:space:]]) ]] || VARLIST="$VARLIST $VAR"
  [[ $TIMLIST =~ (^|[[:space:]])"$TIM"($|[[:space:]]) ]] || TIMLIST="$TIMLIST $TIM"
done

echo ""
echo "Available variables:    $VARLIST"
echo "Available time periods: ${TIMLIST//_/-}"
echo ""
# -----------------------------------------------------------------------------------------------------
# Generate remap weights
WEIGHTS=weights.nc
${cdo} genycon,$GRID $FILE $WEIGHTS

# -----------------------------------------------------------------------------------------------------
# Loop over time and build forcing files
for TIME in $TIMLIST; do
  OUT=`echo climate_${DATA}_${STEP}_${TIME}.nc | tr [A-Z] [a-z]`
  echo "Processing $OUT"

  # Merge files into one yearly file with all variables (requires dataset specific pre-processing)
  temp1=$(mktemp -u --tmpdir=.)
  case $DATA in
    gswp3-w5e5 )
      # decadal files needs to be merged in splitted in single years
      unset YEARFILES
      for FILE in $DATAFILES ; do
        [[ $FILE =~ $TIME ]] && YEARFILES="$YEARFILES  $FILE"
      done
      $cdo merge $YEARFILES $temp1
      $cdo splityear $temp1 ${temp1}_ && rm $temp1
      SINGLEYEARS=$(ls ${temp1}_*)
      ;;
    * ) echo "*** ERROR: Unexpected dataset $DATA in pre-processing" && exit 1 ;;
  esac

  for SFILE in $SINGLEYEARS ; do
    YEAR=$(echo "print('$SFILE'.split('.')[-2].split('_')[1])" | python)
    echo "Processing $YEAR"
    OUT=`echo climate_${DATA}_${STEP}_${YEAR}.nc | tr [A-Z] [a-z]`
    # Forcing dataset dependent preprocessing for target time steps
    IN=$SFILE
    # Convert variables and remap to target resolution
    $cdo -a -b F64 -f nc4c remap,${GRID},${WEIGHTS} -setpartabn,$PARTAB -exprf,${JSBFORC} $IN $OUT && rm $IN

    # Write File metainfo
    ncatted -h -a title,global,o,c,"ICON-LAND standalone forcing based on $DATA" $OUT
    ncatted -h -a institution,global,o,c,"${INST}"                               $OUT
    ncatted -h -a contact,global,o,c,"${CONTACT}"                                $OUT
    ncatted -h -a version,global,o,c,"${VERSION}"                                $OUT
    ncatted -h -a reference,global,o,c,"${REFERENCE}"                            $OUT
    ncatted -h -a history,global,d,c,""                                          $OUT

    # Move to target directory
    mv $OUT ${DTRG}/
  done

done

test -e $PARTAB && rm $PARTAB
test -e $JSBFORC && rm $JSBFORC
test -e $WEIGHTS && rm $WEIGHTS
