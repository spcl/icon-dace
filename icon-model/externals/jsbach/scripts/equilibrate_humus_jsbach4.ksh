#!/bin/ksh

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

#----------------------------------------------------------------------------------------------
#
# Script to set the amount of carbon in the YASSO humus pools closer to equilibrium.
#
# Depending on the simulation set-up the script uses different aggregation levels of the in- and output
# to the humus pool: for a simulation without natural land cover changes the values per pft are used,
# for a simulation with natural land cover changes the box value is used.
#
# NOTE: script has only been tested for a limited number of configurations
#       (annual output, jsbach standalone in ECHAM and ICON environment)
#       and it is well possible that file-name structures or the like need to be further adapted
#       for other setups -- or even generalise handling of timestamps in restartfile names etc.
#
# Note:
# - uses ta variables, not per canopy variables to avoid numerical issues when applying natural LCC
#  (+ in addition an average over the pfts is used when applying natural LCC -- using an aggregated
#     variable from an upper layer of the tile hierarchy might catch potential numeric instabilities
#     potentially caused by fractional changes and resulting implicit scaling)
# - works on the pft level - not yet below!
# - The other carbon pools need to be close to equilibrium before ${year_start}.
# - Read and write permission in the restart directory of ${exp} (see end of script) are required.
# - This script needs 4 digit year dates.
#
#
# julia.nabel@mpimet.mpg.de 2020
# Based on the script "equilibriate_humus.ksh" developed by Thomas Raddatz for jsbach3 (Dec 2018)
#----------------------------------------------------------------------------------------------

# parameters that need to be adapted

#----- tested in echam env
exp=run_for_pft_normalMaxLAI_eqhumus_branch # experiment ID
mod=jsb #model short
restartPrefix=restart_${exp}_${mod}_carbon_
restartSuffix=1231
path_restart=./ # the path, where the simulation stores the restart files

#----- tested in icon env
#exp=jsbalone_R2B4_carbon_setting_as_in_echam_runs
#mod=jsbalone
#restartPrefix=${exp}_restart_${mod}_
#restartSuffix=0101T000000Z
#path_restart=./

work=${path_restart}
year_start=2500 #first year of period to determine average fluxes
year_end=2502 #last year of that period
nyear_restart=1 # number of years within one restart cycle

with_nat_LCC=false

pfts=( 01 02 03 04 05 06 07 08 09 10 11 )

#----------------------------------------------------------------------------------------------

cdo="/sw/rhel6-x64/cdo/cdo-1.9.8-magicsxx-gcc64/bin/cdo -s"
rm=/bin/rm
cp=/bin/cp
mv=/bin/mv
expr=/usr/bin/expr

set -e

#copy final restart file
restartToModify=${path_restart}/${restartPrefix}${year_end}${restartSuffix}.nc
${cp} ${restartToModify} ${path_restart}/ORG_${restartPrefix}${year_end}${restartSuffix}.nc
${cp} ${restartToModify} ${work}/tmp_restart.nc

# loop over pfts
for ipft in ${pfts[@]}; do
  echo "pft${ipft}"

  # loop over years in the period to determine average flux to humus pools and average decomposition rate
  if [[ ${with_nat_LCC} == true ]]; then
    # if with natural LCC: collect the information only once and for the veg tile
    fortile=veg
    if [[ ${ipft} == ${pfts[0]} ]]; then
      calcFactors=true
    else
      calcFactors=false
    fi
  else
    # collect it for the specific pft
    fortile=pft${ipft}
    calcFactors=true
  fi

  if [[ ${calcFactors} == true ]]; then
    year=${year_start}
    while [ ${year} -le ${year_end} ] ; do
      # extract variables
      ${cdo} -selvar,carbon_c_decomp_humus_1_sum_ta_${fortile} ${path_restart}/${restartPrefix}${year}${restartSuffix}.nc ${work}/decH1_${year}.nc
      ${cdo} -selvar,carbon_c_decomp_humus_2_sum_ta_${fortile} ${path_restart}/${restartPrefix}${year}${restartSuffix}.nc ${work}/decH2_${year}.nc
      ${cdo} -selvar,carbon_c_into_humus_1_sum_ta_${fortile} ${path_restart}/${restartPrefix}${year}${restartSuffix}.nc ${work}/inH1_${year}.nc
      ${cdo} -selvar,carbon_c_into_humus_2_sum_ta_${fortile} ${path_restart}/${restartPrefix}${year}${restartSuffix}.nc ${work}/inH2_${year}.nc

      (( year = year + 1 ))
    done

    ${cdo} -mergetime ${work}/decH1_*.nc ${work}/decH1_${year_start}-${year_end}.nc
    ${cdo} -mergetime ${work}/decH2_*.nc ${work}/decH2_${year_start}-${year_end}.nc
    ${cdo} -mergetime ${work}/inH1_*.nc ${work}/inH1_${year_start}-${year_end}.nc
    ${cdo} -mergetime ${work}/inH2_*.nc ${work}/inH2_${year_start}-${year_end}.nc

    ${cdo} -timmean ${work}/decH1_${year_start}-${year_end}.nc ${work}/decH1_${year_start}-${year_end}_mean.nc
    ${cdo} -timmean ${work}/decH2_${year_start}-${year_end}.nc ${work}/decH2_${year_start}-${year_end}_mean.nc
    ${cdo} -timmean ${work}/inH1_${year_start}-${year_end}.nc ${work}/inH1_${year_start}-${year_end}_mean.nc
    ${cdo} -timmean ${work}/inH2_${year_start}-${year_end}.nc ${work}/inH2_${year_start}-${year_end}_mean.nc

    rm ${work}/decH?_${year_start}-${year_end}.nc ${work}/inH?_${year_start}-${year_end}.nc

    # calculate factor to equilibrate humus pools
    ${cdo} -div ${work}/inH1_${year_start}-${year_end}_mean.nc ${work}/decH1_${year_start}-${year_end}_mean.nc ${work}/factor_humus_1.nc
    ${cdo} -div ${work}/inH2_${year_start}-${year_end}_mean.nc ${work}/decH2_${year_start}-${year_end}_mean.nc ${work}/factor_humus_2.nc

    # clean
    ${rm} ${work}/decH?_*.nc ${work}/inH?_*.nc
  fi

  # extract humus pools
  ${cdo} -selvar,carbon_c_humus_1_pft${ipft} ${path_restart}/${restartPrefix}${year_end}${restartSuffix}.nc ${work}/c_H1.nc
  ${cdo} -selvar,carbon_c_humus_2_pft${ipft} ${path_restart}/${restartPrefix}${year_end}${restartSuffix}.nc ${work}/c_H2.nc

  # calculate equilibrium humus pools
  ${cdo} -mul ${work}/c_H1.nc ${work}/factor_humus_1.nc ${work}/c_H1_eq.nc
  ${cdo} -mul ${work}/c_H2.nc ${work}/factor_humus_2.nc ${work}/c_H2_eq.nc

  # replace humus pools
  ${cdo} -replace ${work}/tmp_restart.nc ${work}/c_H1_eq.nc ${work}/tmp2_restart.nc
  ${rm} ${work}/tmp_restart.nc
  ${cdo} -replace ${work}/tmp2_restart.nc ${work}/c_H2_eq.nc ${work}/tmp_restart.nc
  ${rm} ${work}/tmp2_restart.nc

  # clean
  ${rm} ${work}/c_H?_eq.nc ${work}/c_H?.nc

  if [[ ${with_nat_LCC} == false ]]; then
    ${rm} ${work}/factor_humus_?.nc
  fi
done

# replace restart file
${mv} ${work}/tmp_restart.nc ${restartToModify}

exit
