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

set +ex

#module load  nco

##############################################################################################
### This program collects a list of variables e.g. for carbon pools from an ICON restartfile
### and packs them into a new file ###
###
### Note: Tested with nco version 4.6.7 ... older versions might not work.
###       As of March 2018, it doesn't work on MPI workstations because nco is too old.
###       It does work on mistral.
##############################################################################################

### User Interface
INFILE="JS045_restart_atm_19870201T000000Z.nc"
OUTFILE="ic_land_carbon.nc"

PROCESS="carbon"
VARNAME_LIST="c_green    c_woods     c_reserve 
              c_acid_ag1 c_water_ag1 c_ethanol_ag1 c_nonsoluble_ag1
              c_acid_bg1 c_water_bg1 c_ethanol_bg1 c_nonsoluble_bg1 c_humus_1
              c_acid_ag2 c_water_ag2 c_ethanol_ag2 c_nonsoluble_ag2
              c_acid_bg2 c_water_bg2 c_ethanol_bg2 c_nonsoluble_bg2 c_humus_2"
TILE_LIST="box land veg pft01 pft02 pft03 pft04 pft05 pft06 pft07 pft08 pft09 pft10 pft11"

### Delete existing output files ?
if [ -f  ${OUTFILE} ];then
   echo "The output file ${OUTFILE} already exists. Shall I delete it?"
   echo "RETURN =File will be deleted."
   echo "XXX    =Every other input stops this script here."
   read y
   if [ ! ${y} ] ; then 
      rm ${OUTFILE}
      echo "The existing file ${OUTFILE} was deleted."
      echo ""
   else
      exit
   fi
fi

### collect vars from restart file
for VARNAME in ${VARNAME_LIST} ; do
   echo -e "\033[01;31;40m"
   echo "=== ${VARNAME} ==="
   tput sgr0

   for TILE in ${TILE_LIST} ; do
      STREAM_NAME="${PROCESS}_${VARNAME}_${TILE}"
      echo "${STREAM_NAME}"

      cdo selvar,${STREAM_NAME}                      ${INFILE}    ZWERG1
      ncks -A                                        ZWERG1       ZWERG2
   done
done

ncks -C -x -v time            ZWERG2       ZWERG3      # Delete variable time
ncwa -O -a time               ZWERG3       ${OUTFILE}  # Delete dimension time
ncrename -d ncells_2,ncells ${OUTFILE}                 # This needs a recent nco version!

### Clean up
rm ZWERG1  ZWERG2  ZWERG3

### Finish
echo -e "\033[00;31;48m" 
echo ""
echo "Script reached the end without interruption !"
echo ""
tput sgr0
exit

