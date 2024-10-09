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

#_____________________________________________________________________________
# Script to construct a fractional mask
#
# It takes the mask of the (finer) ocean grid and remaps
# it onto the (coarser) atmosphere grid with conservative remapping.

set -e

#oceGridID=0036
#oceRes=R02B04

# atmGridID=0043
#atmRes=R02B04

atmGridFile=${icon_grid}
oceMaskFile=${icon_grid_oce}
outputPath=${output_root_dir}/${grid_label}/fractional_mask

wrkdir=${work_dir}
if ! mkdir -p $wrkdir ; then
  echo "$0: cannot create work dir $wrkdir"
  exit 1
fi
cd $wrkdir

cdo="cdo -s -f nc4 -b F64"

echo '------------------------------------------------------------------'
echo "$0: Working directory: $(pwd)"
echo "$0: Output  directory: ${outputPath}"

ln -s ${oceMaskFile}  ./icon_mask_${oceGridID}_${oceRes}_G.nc
ln -s ${atmGridFile}  ./icon_grid_${atmGridID}_${atmRes}_G.nc

$cdo -b F64 -gtc,0 -selvar,cell_sea_land_mask icon_mask_${oceGridID}_${oceRes}_G.nc cell_sea_land_mask.nc
$cdo setmisstoc,1 -remapycon,icon_grid_${atmGridID}_${atmRes}_G.nc cell_sea_land_mask.nc fractional_lsm_${atmGridID}_${oceGridID}.nc

mkdir -p ${outputPath}
mv fractional_lsm_${atmGridID}_${oceGridID}.nc ${outputPath}/
echo "$0:"
echo "$0: Fractional land sea mask:  ${outputPath}/fractional_lsm_${atmGridID}_${oceGridID}.nc"
echo "$0:  done"
echo "------------------------------------------------------------------"

rm -r $wrkdir
exit 0

