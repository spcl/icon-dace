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

#
# run as: many.error.s 20120101 00 20120102 00 1 R2B04 25
# -------------------------------------------------------
set -x
	
inidate=${1}
initime=${2}
verdate=${3}
vertime=${4}
ndays=${5}
res=${6}
expnum=${7}
expref=${8}
echo "Arguments: inidate="${inidate}" initime="${initime}" verdate="${verdate}" vertime="${vertime}" ndays="${ndays}" res="${res}" expnum="${expnum}" expref="${expref}

#expnum=34

#              00h      24h      48h      120h
#set -A inidate 20110101 20110101 20110101 20110101
#set -A initime 00       00       00       00      
#set -A verdate 20110101 20110102 20110103 20110106
#set -A vertime 00       00       00       00

#set -A inidate 20110101 
#set -A initime 00       
#set -A verdate 20110102 
#set -A vertime 00       

#mkdir -p "/fe1-daten/"$USER"/plots/icon/nwp.exp"${expnum}
#scriptdir="/fe1-daten/"$USER"/metview/ICON/"
scriptdir="./"
cd ${scriptdir}

metview=metview
#metview=/usr/local/apps/Metview/metview4_dev
module swap metview/4.5.0

met_job=met.job.all.$nstart
\rm -rf $met_job

integer nt
while [[ $nt < ${#inidate[*]} ]]; do

  # -------------------------------------------------------

  set -A vars TQV             TQC            TQI              CAPE               \
              TCC             HCC            MCC              LCC                \
              TQ1             TQ2            TQ3              TQ4       TQ5      \
              PS              PMSL                                               \
              ACCSOB_S        ACCTHB_S       ACCSOB_T         ACCTHB_T           \
              ACCLHFL_S       ACCSHFL_S                                          \
              TOT_PREC        RAIN_GSP       SNOW_GSP         RAIN_CON  SNOW_CON \
              T_G             T_2M           U_10M            V_10M              \
              SP_10M          DZH_SNOW                                           \
              FR_SEAICE       H_SEAICE       T_SEAICE         SST

#             QV_2M           QV_S           T_GT_tile_1      T_S_tile_1         \
#             W_I_tile_1      T_SNOW_tile_1  WTOT_SNOW_tile_1 WLIQ_SNOW_tile_1   \
#             RHO_SNOW_tile_1 H_SNOW         W_SNOW    
##            Z0  T_G 
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done
  set -A vars TQV             TQC            TQI              CAPE               \
              TCC             HCC            MCC              LCC                \
              TQ1             TQ2            TQ3                                 \
              ACCSOB_S        ACCTHB_S       ACCSOB_T         ACCTHB_T           \
              ACCLHFL_S       ACCSHFL_S      TOT_PREC         PS                 \
              T_G             T_2M           U_10M            V_10M              \
              SP_10M 
#             H_SNOW_tile_1   RHO_SNOW_tile_1 
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc diff     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc rms      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc diff_ref ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc ctr      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done

  # -------------------------------------------------------
 
  set -A vars  SP_10M  T_2M  T_G  ACCLHFL_S  ACCSHFL_S  ACCSOB_T  ACCTHB_T     #comparison to observations,     QV_2M  
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc diff_obs ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc ifs_obs  ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done
  
  # -------------------------------------------------------
 
  set -A vars  W_SO  T_SO     #W_SO_ICE_tile_1 
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc ctr      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc diff     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc rms      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var sfc diff_ref ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done
  
  # -------------------------------------------------------
 
  set -A vars  T_SNOW_M  RHO_SNOW_M  DZH_SNOW_M
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var 114 snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var 114 ctr      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var 114 diff     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var 114 rms      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error $expnum $expref $var 114 diff_ref ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done
  
  # -------------------------------------------------------
  	
  set -A vars T  U  V  Q1  Q2  Q3  QV  QC  QI  CC  P  QR  QS QTVAR EXTRA3D1 EXTRA3D2 #QR QS  QTVAR  O3  P               
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var ml diff     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var ml rms      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var ml diff_ref ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var ml snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var ml ctr      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var ml diff     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var ml rms      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var ml diff_ref ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#   echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var ml snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#   echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var ml ctr      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done
    
  # -------------------------------------------------------
  	
  set -A vars ttendcds qtendcds utendcds vtendcds \
              ttendts  qtendt   utendts  vtendts  \
              ttends   utends   vtends            \
              ewgd     nsgd                       \
              ttendsw  ttendlw  O3
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var ml snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var ml diff     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var ml diff_ref ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var ml snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done
    
  # -------------------------------------------------------

  set -A vars T  QV  U  V  FI
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var pl diff     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var pl rms      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var pl diff_ref ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var pl snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var pl ctr      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var pl diff     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var pl rms      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
    echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var pl diff_ref ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#   echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var pl snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#   echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var pl ctr      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done

  set -A vars  Q1  Q2  Q3  QC  QI  CC # QR QS
  for var in ${vars[*]}
  do
    echo ${metview} -b ${scriptdir}zonal.error $expnum $expref $var pl snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#   echo ${metview} -b ${scriptdir}map.error   $expnum $expref $var pl snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
  done

  # -------------------------------------------------------

  set -A vars T  U  V  P  Q1  Q2  Q3  QV  QC  QI  CC  # CRWC CSWC
#  for var in ${vars[*]}
#  do
#    echo ${metview} -b ${scriptdir}map.error  $expnum $expref $var zl diff     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#    echo ${metview} -b ${scriptdir}map.error  $expnum $expref $var zl diff_ref ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#    echo ${metview} -b ${scriptdir}map.error  $expnum $expref $var zl snap     ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#    echo ${metview} -b ${scriptdir}map.error  $expnum $expref $var zl ctr      ${inidate[nt]} ${initime[nt]} ${verdate[nt]} ${vertime[nt]} ${ndays} ${res} >> $met_job
#  done


  nt=nt+1
done  # end loop date

exit
