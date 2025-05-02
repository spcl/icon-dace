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

#------------------------------------------------------------------------------
# The extpar file used here is NOT the one just generated in exptpar4jsbach!!
#   That would be too easy ...
#
# We use here the extpar file read by the NWP atmosphere at run time. It
# includes several land variables. For consistency - and to avoid misinter-
# pretations, theses variables are replaced here by the corresponding
# bc_land-file data.
#------------------------------------------------------------------------------
set -e
prog=$(basename $0)
extparfile=${nwp_extpar_file}
extpar_name=${extparfile##*/}
extpar_jsb=${bc_file_dir}/${extpar_name%.nc}_jsb.nc

bc_frac=${bc_file_dir}/bc_land_frac.nc
bc_phys=${bc_file_dir}/bc_land_phys.nc
bc_sso=${bc_file_dir}/bc_land_sso.nc
bc_soil=${bc_file_dir}/bc_land_soil.nc

icon_grid=${icon_grid}

cdo="cdo -s -f nc4"

tmpdir=${work_dir}/extpar_tmp
mkdir -p ${tmpdir}; rm -f ${tmpdir}/*.nc
cd ${tmpdir}

echo "=============================================================="
echo "===  Adaptation of NWP extpar file ${extparfile}"
echo "=============================================================="
if [[ ${extparfile} == "" ]]; then
  echo "ERROR: NWP extpar file ${nwp_extpar_file} does not exist!"
  exit 1
fi

# Variables from bc_land_fract
$cdo selvar,notsea ${bc_frac}  notsea.nc
$cdo setvar,FR_LAND -sub -selvar,notsea ${bc_frac} \
                         -selvar,lake   ${bc_frac}       new_FR_LAND.nc
$cdo setvar,ICE          -selvar,glac   ${bc_frac}       new_ICE.nc
$cdo setvar,FR_LAKE      -selvar,lake   ${bc_frac}       new_FR_LAKE.nc

# Variables from bc_land_phys
# Z0 should be identical, as it was used from extpar in jsbach4_ini_files_from_extpar.sh

# Variables from bc_land_sso
# The orographic parameters are read from "initial_extpar_file" in jsbach4_ini_files_from_gauss.sh.
# If the nwp_extpar file was used as initial_extpar_file (default) a replacement has no effect.
# We however need to adapt these variables to changes of the land sea mask:
# Set orography parameters of complete ocean cells to zero (this is the value of existing
# complete ocean cells)
for var in SSO_GAMMA SSO_OROMAX SSO_OROMIN SSO_SIGMA SSO_STDH SSO_THETA topography_c; do
  if [[ $(cdo showvar ${extparfile} | grep $var ) != "" ]]; then
    $cdo ifthen notsea.nc -selvar,${var} ${extparfile} ${var}.tmp
    $cdo setmisstoc,0 ${var}.tmp                          new_${var}.nc
  fi
done

# Variables from bc_land_soil
$cdo setvar,ROOTDP    -selvar,root_depth   ${bc_soil}     new_ROOTDP.nc

# Adaptation of FAO soil types
if [[ $(cdo showvar ${extparfile} | grep 'SOILTYP' ) != "" ]]; then
  # soiltyp 9 (sea water) and 10 (sea ice): all land points must not have 9 or 10
  $cdo setmisstoc,5 -setmissval,9  -selvar,SOILTYP ${extparfile} SOILTYP.tmp
  $cdo setmisstoc,5 -setmissval,10         SOILTYP.tmp           SOILTYP.tmp2
  # soiltyp 1 (glacier): set 1 according to new glacier mask
  $cdo setmisstoc,5 -setmissval,1          SOILTYP.tmp2          SOILTYP.tmp3
  $cdo setname,SOILTYP -ifthenelse new_ICE.nc new_ICE.nc SOILTYP.tmp3 new_SOILTYP.nc
fi

# Replace all adapted variables
$cdo -O merge new_*.nc new_vars.nc
$cdo setgrid,${icon_grid} -replace ${extparfile} new_vars.nc ${extpar_jsb}

# Remove variables that should not be included to avoid incorrect usage or confusion
# - Remove cell_sea_land_mask in case it exists, because its existence
#   triggers mask corrections in the NWP atmosphere at run time - which we do not
#   want after the adaptations done here.
# - ALB, ALNID, ALUVD not used with 'albedo_type == MODIS'; would need corrections
# - BULK_DENS, SUB_BULK_DENS - currently not used at all; would need adaptation
# - EMISS, EMIS_RAD - not used in current setups (itype_lwemiss = 0); re-check if needed
# - NDVI: not used; would need corrections
# - W_SNOW: currently not used, IFS analysis data
for var in cell_sea_land_mask ALB ALNID ALUVD BULK_DENS NDVI SUB_BULK_DENS W_SNOW; do
  if [[ $(cdo showvar ${extpar_jsb} | grep -w $var) != "" ]]; then
    $cdo delname,$var ${extpar_jsb} ${extpar_jsb}.tmp
    mv ${extpar_jsb}.tmp ${extpar_jsb}
  fi
done

echo "-------------------------------------------------------------------"
echo "$prog: Modified extpar file:"
echo "     ${extpar_jsb}"
echo "----------------------------------------------------------------"

cd ..
rm -rf ${tmpdir}
