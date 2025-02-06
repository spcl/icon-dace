#!/bin/ksh

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# -------------------------------------------------------
# Calculate water budget from output
#
# run as: ./water_budget.s /scratch/ms/de/dei2/temp/dei2_146_2012010100_DOM01_ML_0001.nc
#
# info:
#   * requires native netcdf output (with output_grid=.true.):
#        &output_nml
#        filetype                     =  4                        ! output format: 2=GRIB2, 4=NETCDFv2
#        dom                          = -1                        ! write all domains
#        output_time_unit             =  1                        ! 1: seconds
#        output_bounds                =  0., 86400., 3600.        ! start, end, increment
#        steps_per_file               =  40
#        mode                         =  1  ! 1: forecast mode (relative t-axis), 2: climate mode (absolute t-axis)
#        taxis_tunit                  =  3  ! output time unit, 3: hours, 2:min, 1:sec
#        include_last                 = .TRUE.
#        output_filename              = '${basename}'             ! file name base
#        ml_varlist                   = 'tot_prec', 'acclhfl_s', 'accqhfl_s', 'tqv', 'tqc', 'tqi', 'tqr', 'tqs'
#        output_grid                  = .TRUE.
#        remap                        = 0                         ! 1: latlon,  0: native grid
#       /
#
#   * tip: split GRIB file into typeOfLevel's
#      grib_copy dei2_146_2012010100_DOM01_ML_0001.grb dei2_146_2012010100_DOM01_ML_0001.[typeOfLevel].grb
#      ./water_budget.s /scratch/ms/de/dei2/temp/dei2_146_2012010100_DOM01_ML_0001.surface.grb
#
#   * variable names:
#     variable   netcdf name   grib short name
#     ........................................
#     qfl        accqhfl_s     ??
#     lhf        acclhfl_s     slhf
#     tot_prec   tot_prec      tp     
#     tqv        tqv           tciwv
#     tql        tqc           param69.1.0
#     tqi        tqi           param70.1.0
#     tqr        tqr           tcolr
#     tqs        tqs           tcols 
#                     
# Martin Koehler, May 2014
# -------------------------------------------------------

ml_file=${1}

module load cdo
LANG=''          # for math with 1.0 not 1,0

set -A date     `cdo infov ${ml_file} | grep -w "tqv" | awk '{print $3}'`
set -A time     `cdo infov ${ml_file} | grep -w "tqv" | awk '{print $4}'`
set -A lhf      `cdo output -fldmean -selname,acclhfl_s  ${ml_file}`
set -A qfl      `cdo output -fldmean -selname,accqhfl_s  ${ml_file}`
set -A tot_prec `cdo output -fldmean -selname,tot_prec   ${ml_file}`   
set -A tqv      `cdo output -fldmean -selname,tqv        ${ml_file}`   
set -A tql      `cdo output -fldmean -selname,tqc        ${ml_file}`
set -A tqi      `cdo output -fldmean -selname,tqi        ${ml_file}`
set -A tqr      `cdo output -fldmean -selname,tqr        ${ml_file}`
set -A tqs      `cdo output -fldmean -selname,tqs        ${ml_file}`  

float lhf1
float pme
float dqt
float budget

echo '___________________________________________________________________________________________________________________________________________________________________________'
echo 'date         time        lhf               tot_prec      tqv         tql           tqi           tqr           tqs           P-E               del qt         budget'
echo '                         [mm]              [mm]          [mm]        [mm]          [mm]          [mm]          [mm]          [mm/h]            [mm/h]         [mm/h]'
echo '___________________________________________________________________________________________________________________________________________________________________________'
for ((nn=0; nn<=24; nn+=1)) ; do

 #lhf1=${lhf[$nn]}/2500800.0
  lhf1=${qfl[$nn]}
  if [[ ${nn} -eq 0 ]] ; then
    pme=0.0
    dqt=0.0
    budget=0.0
  else                     # + ${lhf[$nn]}/2500800.0            - ${lhf[$nn-1]}/2500800.0
    pme=$(( ${tot_prec[$nn]} + ${qfl[$nn]} - ${tot_prec[$nn-1]} - ${qfl[$nn-1]} ))
    dqt=$(( ${tqv[$nn]}+${tql[$nn]}+${tqi[$nn]}+${tqr[$nn]}+${tqs[$nn]} - ${tqv[$nn-1]} - ${tql[$nn-1]} - ${tqi[$nn-1]} - ${tqr[$nn-1]} - ${tqs[$nn-1]} ))
    budget=$(( ${pme} + ${dqt} ))
  fi

  printf '%-12s %-11s %-17s %-13s %-11s %-13s %-13s %-13s %-13s %-17s %-14s %-15s \n' \
    ${date[$nn]} ${time[$nn]} ${lhf1}     ${tot_prec[$nn]}                                  \
    ${tqv[$nn]}  ${tql[$nn]}  ${tqi[$nn]} ${tqs[$nn]}      ${tqr[$nn]} ${pme} ${dqt} ${budget}
done
echo '____________________________________________________________________________________________________________________________________________________________________________'


exit

