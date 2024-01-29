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
extparfile=${coupled_extpar_file}
extpar_file=${extparfile##*/}
extpar_jsb=${bc_file_dir}/${extpar_file%.nc}_jsb.nc

bc_frac=${bc_file_dir}/bc_land_frac_1850.nc
bc_phys=${bc_file_dir}/bc_land_phys_1850.nc
bc_sso=${bc_file_dir}/bc_land_sso_1850.nc

cdo="cdo -s"

tmpdir=${work_dir}/extpar_tmp
mkdir -p ${tmpdir}; rm -f ${tmpdir}/*.nc

$cdo setvar,cell_sea_land_mask -selvar,notsea ${bc_frac} ${tmpdir}/cell_sea_land_mask.nc
$cdo setvar,FR_LAND -sub -selvar,notsea ${bc_frac} \
                         -selvar,lake   ${bc_frac}       ${tmpdir}/FR_LAND.nc
$cdo setvar,ICE          -selvar,glac   ${bc_frac}       ${tmpdir}/ICE.nc
$cdo setvar,FR_LAKE      -selvar,lake   ${bc_frac}       ${tmpdir}/FR_LAKE.nc
$cdo setvar,Z0 -selvar,roughness_length ${bc_phys}       ${tmpdir}/Z0.nc
$cdo setvar,topography_c -selvar,oromea ${bc_sso}        ${tmpdir}/topography_c.nc
$cdo setvar,SSO_STDH     -selvar,orostd ${bc_sso}        ${tmpdir}/SSO_STDH.nc
$cdo setvar,SSO_THETA    -selvar,orothe ${bc_sso}        ${tmpdir}/SSO_THETA.nc
$cdo setvar,SSO_GAMMA    -selvar,orogam ${bc_sso}        ${tmpdir}/SSO_GAMMA.nc
$cdo setvar,SSO_SIGMA    -selvar,orosig ${bc_sso}        ${tmpdir}/SSO_SIGMA.nc
$cdo setmisstoc,5 -setmissval,9  -selvar,SOILTYP ${extparfile} ${tmpdir}/SOILTYP.nc
$cdo setmisstoc,290      -selvar,T_SEA  ${extparfile}    ${tmpdir}/T_SEA.nc

$cdo -O merge ${tmpdir}/*.nc ${tmpdir}/new_vars.nc
$cdo replace ${extparfile} ${tmpdir}/new_vars.nc ${extpar_jsb}

echo "-------------------------------------------------------------------"
echo "$prog: Modified extpar file:"
echo "     ${extpar_jsb}"
echo "-------------------------------------------------------------------"

rm -rf ${tmpdir}
