#!/bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# ----------------------------------------------------------------------------
#
# loop over dates to
# - setup init_SCM.nc and
# - submit SCM run scripst for GOAmazon case on NEC
# - module load cdo nco (do before submitting)
# - # nohup loop_SCM_GOAmazon.s 2008 >2008.stdout 2>&1 &
# - cd workdir/2015
# - nohup /hpc/uhome/mkoehler/icon/icon-nwp-test1/scripts/scm/run/loop_SCM_GOAmazon.s 2015 &
#
# schedule:
# - reformat intput data: get_GASS_DCP_GoAmazon.py
# - loop_scm_GOAMazon.s
#   - loop over days, select input data for that day
#   - run SCM:            run_SCM_GASS_DCP_SGP_nec
#   - reformat output:    GOAmazon_reformat_output.py
#
# quick test of data:
#   - ensemble mean
#     cdo ensavg scm_GASS_flx_ML_2014*T000000Z.nc temp.mean.nc
#     cdo fldavg temp.mean.nc scm_out_ML_2014.mean.nc
#   - deaccumulation
#     file=scm_out_ML_2014.mean.nc
#     nstep=`cdo -s ntime $file`
#     cdo sub -seltimestep,2/$nstep $file -seltimestep,1/`expr $nstep - 1` $file  out.nc
#     cdo mulc,48 out.nc scm_out_ML_2014.mean2.nc
#
#-----------------------------------------------------------------------------

set -ex

# setup

do_init=1     # input processing
do_icon=1     # ICON SCM run
do_post=0     # postprocessing for GASS output

# loop over years

export YYYY=$1       # 2008
YYYYlast=$1          # 2008

while (( ${YYYY} <= ${YYYYlast} )) ; do

EXPNAME=SCM_GASS_DCP_GoAmazon
#EXPNAME=SCM_GASS_DCP_SGP
EXPDIR=${WORK}/run-icon/scm/${EXPNAME}/${YYYY}
mkdir -p ${EXPDIR}/scm-out ${EXPDIR}/scm-GASS-format
ICONDIR=${HOME}/icon/icon-nwp-test1

# SCM data directory (grids, init data, extpar)
SCMDATA=/hpc/uwork/mkoehler/scm/data       # at DWD on NEC


cd ${EXPDIR}

# loop over dates

# GoAmazon
date=${YYYY}022000            #2014010100
datelast=${YYYY}022000        #2014123100

# SGP
#date=${YYYY}042800       
#datelast=${YYYY}083000

# test
#date=${YYYY}072400
#datelast=${YYYY}072400

while (( ${date} <= ${datelast} )) ; do

  start_date=`date_calc.py --date=${date} -a printfmt`
  end_date=`date_calc.py   --date=${date} -a printfmt -s -5`
# end_date="2014-03-26T21:00:00Z"

  echo $start_date > ${EXPDIR}/start_date   # 5day case study
  echo $end_date   > ${EXPDIR}/end_date     # ...
# echo $YYYY       > ${WORK}/run-icon/scm/${EXPNAME}/YYYY
  echo $YYYY       > ${EXPDIR}/YYYY

  YYYYMMDD=`echo ${date} | cut -c 1-8`
  YYYY_MM_DD=`echo ${start_date} | cut -c 1-10`


  # --- cut period from input (on lce or rcl)

  if [ $do_init = 1 ] ; then
    cdo seldate,`echo ${start_date} | cut -c 1-19`,`echo ${end_date} | cut -c 1-19`             \
      ${SCMDATA}/init_data/init_${EXPNAME}_${YYYY}.nc init_SCM.nc

    # fix dtime to 20min:
    cdo settaxis,${YYYY_MM_DD},00:00:00,20minutes init_SCM.nc out.nc
    mv -f out.nc init_SCM.nc

    ncks -O -v longitude,latitude,FR_LAND,PLCOV_MX,LAI_MX,ROOTDP,RSMIN,SOILTYP,Z0,EMIS_RAD,TOPO \
      ${SCMDATA}/init_data/init_${EXPNAME}_${YYYY}.nc init_SCM_2.nc
    ncks -A init_SCM_2.nc init_SCM.nc
    ncks -O -d time,0,360 init_SCM.nc init_SCM.nc 
#   ncks -O -d time,0,351 init_SCM.nc init_SCM.nc  # last forecast only to 21UTC (2014032200)
    ncks -O -3 init_SCM.nc out.nc
    ncrename -d time,nt out.nc
    ncks -O -4 out.nc init_SCM.nc
    ncatted -O -a units,time,a,c,'seconds since '`echo $start_date | cut -c 1-10`' 0:00:00 0:00' init_SCM.nc

#   \mv -f init_SCM.nc ${EXPDIR}
 
  # clean
    rm -rf init_SCM_2.nc out.nc

  fi


  # --- submit SCM

  if [ $do_icon = 1 ] ; then
    qsubw ${ICONDIR}/scripts/scm/run/run_SCM_GASS_DCP_GoAmazon_nec       # qsubw waits till completion
#   qsubw ${ICONDIR}/scripts/scm/run/run_SCM_GASS_DCP_SGP_nec            # qsubw waits till completion

  # move SCM output data to scm-out

    mv ${EXPDIR}/scm_GASS_out_PL_${YYYYMMDD}T000000Z.nc \
       ${EXPDIR}/scm_GASS_out_ML_${YYYYMMDD}T000000Z.nc \
       ${EXPDIR}/scm_GASS_flx_ML_${YYYYMMDD}T000000Z.nc \
       ${EXPDIR}/scm_out_ML_${YYYYMMDD}T000000Z.nc      \
       ${EXPDIR}/scm-out
  fi


  # --- prepare data to GASS format

  if [ $do_post = 1 ] ; then
    python3 ${ICONDIR}/scripts/scm/run/GOAmazon_reformat_output.py $YYYYMMDD
  fi


  # end loop
  
  date=`date_calc.py --date=${date} -s -1`

done

  YYYY=`expr $YYYY + 1`
done


