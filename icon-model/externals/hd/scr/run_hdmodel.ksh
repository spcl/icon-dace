#!/bin/ksh
#
# run_hdmodel.ksh - HD model runscript on Levante
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
#
#################################################
# batch settings for HD model run on Levante at DKRZ
#################################################
#
### Batch Queuing System is SLURM
#SBATCH --output=HD5.run.o%j
#SBATCH --error=HD5.run.o%j
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stefan.hagemann@hereon.de
#SBATCH --account=gg0302
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=1             # Specify max. number of tasks to be invoked
# #SBATCH --mem=5G
# #SBATCH --mem-per-cpu=500      # Specify real memory required per CPU in MegaBytes
# #SBATCH --mem-per-cpu=2000      # Specify real memory for global 5 min run or levante
#SBATCH --mem-per-cpu=4000      # Specify real memory for global 5 min run or levante
#SBATCH --time=01:00:00         # Set a limit on the total run time (0:30 for reg., 1:20 for G)
###############################################################################
#
# Levante numbers - Global 5 Min. mem_per_cpu: 3500, time: 0:40
# Levante numbers - Euro 5 Min. mem_per_cpu: 2000, time: 0:10
#
CMPI=2   # 1 = MPI, 2= OpenMPI
if (( ${CMPI} == 1 )) ; then
  export MEMORY_AFFINITY=MCM
  export MP_PRINTENV=YES
  export MP_LABELIO=YES
  export MP_INFOLEVEL=0
  export MP_EAGER_LIMIT=64k
  export MP_BUFFER_MEM=64M,256M
  export MP_USE_BULK_XFER=NO
  export MP_BULK_MIN_MSG_SIZE=128k
  export MP_RFIFO_SIZE=4M
  export MP_SHM_ATTACH_THRESH=500000
  export LAPI_DEBUG_STRIPE_SEND_FLIP=8

  # Important on Levante
  export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
  export I_MPI_PMI=pmi2

  export ELG_BUFFER_SIZE=$(echo '100*1000*1000' | bc)
  export SCT_EVENTCOUNTERS=0
elif (( ${CMPI} == 2 )) ; then
  export OMPI_MCA_pml="ucx"
  export OMPI_MCA_btl=self
  export OMPI_MCA_osc="pt2pt"
  export UCX_IB_ADDR_TYPE=ib_global
  # for most runs one may or may not want to disable HCOLL
  export OMPI_MCA_coll="^ml,hcoll"
  export OMPI_MCA_coll_hcoll_enable="0"
  export HCOLL_ENABLE_MCAST_ALL="0"
  export HCOLL_MAIN_IB=mlx5_0:1
  export UCX_NET_DEVICES=mlx5_0:1
  export UCX_TLS=mm,knem,cma,dc_mlx5,dc_x,self
  export UCX_UNIFIED_MODE=y
  export HDF5_USE_FILE_LOCKING=FALSE
  export OMPI_MCA_io="romio321"
  export UCX_HANDLE_ERRORS=bt
fi
#export SCT_NESTEDTIMERS=0
set -e
date
HOME=$(pwd)
#
# *********** Include user settings from sub scripts ***
#
# Read paths definitions from file setuserdir.com
#     HDMAIN              # HD main Directory that includes script, grid and log directories.
#     HDDIR               # Run Directory
#     HDFILE              # HD input data Directory, e.g. HD parameter and start files.
#     HDFORCING           # HD directory with forcing data in subdirectories
if [ -e setuserdir.com ] ; then        # Script was started in $HDMAIN/scr
  source ./setuserdir.com
elif [ -e ../scr/setuserdir.com ] ; then  # Script was started in $HDMAIN/exp. directory
  source ../scr/setuserdir.com
else
  echo 'Script setuserdir.com neither exist in present directory nor in the HDMAIN subdirectory scr!'
  echo 'Hence, run script was likely neither run in the HD main directory HDMAIN/exp. nor its HDMAIN/scr.'
  exit
fi
if (( ${ierrset} == 1 )) ; then
  echo "This is the end! Some basic directories are missing." 
  exit
fi
#
# *** Settings for HD run
EXP=EXPERIMENT
case $EXP in
  ''|*[!0-9]*) echo "run hd_run_settings.ksh in ${HDMAIN}/scr"
               if [ -e ${HDMAIN}/scr/hd_run_settings.ksh ] ; then
                 # Export variable ICALL used by script hd_run_settings.ksh
                 export ICALL=1
                 source ${HDMAIN}/scr/hd_run_settings.ksh
               else
                 echo "Script hd_run_settings.ksh does not exist in ${HDMAIN}/scr!"
                 exit
               fi  ;;
           *)  echo "Experiment = $EXP" 
               if [ -e ${HDMAIN}/${EXP}/hd_run_settings_${EXP}.ksh ] ; then
                 export ICALL=1
                 source ${HDMAIN}/${EXP}/hd_run_settings_${EXP}.ksh
               else
                 echo "Script hd_run_settings_${EXP}.ksh does not exist in ${HDMAIN}/${EXP}!"
                 exit
               fi  ;;
esac
echo $cstart ' for Exp. ' $EXP ' and Year ' $YYYY

if (( ${INEU} == 1 )) ; then
  set +e 
  mkdir ${HDMAIN}/$EXP
  set -e 
  cd ${HDMAIN}/$EXP
  cp -p ${HDMAIN}/scr/hd_run_settings.ksh ./hd_run_settings_${EXP}.ksh
  cp -p ${HDMAIN}/scr/hd_post.ksh .
#
# copy standard run script in experiment directory and assign experiment no. to file name 
  cp -p ${HDMAIN}/scr/run_hdmodel.ksh ${HDMAIN}/${EXP}/run_hdmodel_${EXP}.ksh
  sed -i "s/EXPERIMENT/${EXP}/"  ${HDMAIN}/${EXP}/run_hdmodel_${EXP}.ksh
else  
  # Export variable ICALL used by script hd_run_settings_${EXP}.ksh
  export ICALL=2
  source ${HDMAIN}/${EXP}/hd_run_settings_${EXP}.ksh
  if [ -e ${HDMAIN}/${EXP}/run_hdmodel_${EXP}.ksh ];then 
    echo "${HDMAIN}/${EXP}/run_hdmodel_${EXP}.ksh exists"
  else
    cp -p ${HDMAIN}/scr/run_hdmodel.ksh ${HDMAIN}/${EXP}/run_hdmodel_${EXP}.ksh
    sed -i "s/EXPERIMENT/${EXP}/"  ${HDMAIN}/${EXP}/run_hdmodel_${EXP}.ksh
  fi
fi
#
# Settings Summary
echo $cstart ' for Exp. ' $EXP ' and Year ' $YYYY
echo 'Resolution of forcing: ' $RES_INP ' with Exp. ID' $EXPINP
echo 'Start file: '$HDSTART ' in year' $YYYY
#
# ******* Paths setting ***
#
HDOUT=${HDDIR}/${EXP}/out                  # HD Output dir.
#
# *** Make output directory
cd $HDDIR
rm -rf coredir*
#rm -f HD.err* HD.out*
set +e 
mkdir $EXP
cd $EXP
mkdir out
###rm -f runoff.srv drainage.srv
set -e
#
#
# Input format of forcing files with surface runoff and drainage (subsurface runoff) data
case $CFORM in
  srv ) DNRUN='runoff.srv'
        DNDR='drainage.srv'
        iform=0 ;;
  nc  ) DNRUN='hdforcing.nc'
        DNDR=$DNRUN
        iform=1 ;;
    * ) echo 'File format ' $CFORM ' not defined --> STOP!' ; exit ;;
esac
echo 'File format ' $CFORM ' --> iform = ' $iform
#
# *** Monthly runs --> read month info
if (( ${IWORK} == 2 )) ; then
  if (( ${INEU} == 0 )) ; then
    MM=`cat ${HDMAIN}/log/${EXP}.month`
  fi
  echo HD RESTART ' for Exp. ' $EXP ' and Month ' $MM
  cdocom="${cdocom} -selmon,$MM"
fi
#
# ***** Prepare HD forcing data ***************************************************************
#
# Copy forcing preparation scripts to run directory to avoid that changes in these scripts affect running simulations
if (( ${INEU} == 1 )) ; then
  cp -p ${HDMAIN}/scr/prepare_hdforcing.ksh ${HDDIR}/${EXP}/
  cp -p ${HDMAIN}/${EXE} ${HDDIR}/${EXP}/${EXP}_${EXE}
fi
#
# Export variables used by script prepare_hdforcing.ksh
export HDRES=$HDRES
export FORCE_RES=$FORCE_RES
export IFORCE=$IFORCE
export EXPINP=$EXPINP
export DNREMAP=$DNREMAP
export GRID=$GRID
export YEAR=$YYYY 
export iform=$iform
export DNRUN=$DNRUN
export DNDR=$DNDR

source ${HDDIR}/${EXP}/prepare_hdforcing.ksh
echo "Number of forcing time steps per day: $ndt_day"
echo "Unit Factor: $UFAK"
if (( ${ndt_day} <= 0 )) ; then
  echo "Number of forcing time steps per day < 0 --> Error" 
  exit
fi
#
# ******* Linking of necessary HD model files ******************************
if (( ${INEU} == 1 )) ; then
  set +e ; rm grid_hd.txt hdpara.nc  masks.nc hdstart.nc rmp_hd.nc ; set -e
fi
if [ -e ${HDOUT}/${EXP}_hdrestart_${YYYY}${MM}.nc ];then
  set +e ; rm hdstart.nc ; set -e
  ln -s ${HDOUT}/${EXP}_hdrestart_${YYYY}${MM}.nc hdstart.nc
fi
if [ -e hdpara.nc ];then echo 'hdpara.nc exists' ; else
  ln -sf ${HDFILE}/${DNPARA} hdpara.nc
fi
#
if [ -e hdstart.nc ];then echo 'hdstart.nc exists' ; else
  if (( ${IZIP_START} == 1 )) ; then
    if [ -e ${HDFILE}/${HDSTART}.gz ];then
      cp -p ${HDFILE}/${HDSTART}.gz ./hdstart.nc.gz
      gunzip hdstart.nc.gz
    elif [ -e ${HDFILE}/${HDSTART} ];then
      ln -sf ${HDFILE}/${HDSTART} hdstart.nc
    fi
  else
    ln -sf ${HDFILE}/${HDSTART} hdstart.nc
  fi
fi
if [ -e grid_hd.txt ];then echo 'grid_hd.txt exists' ; else
  ln -sf ${HDMAIN}/grid/${GRID} grid_hd.txt
fi
if [ -e masks.nc ] ; then echo 'masks.nc exists' ; else
  ln -sf ${HDFILE}/${DNMAS} masks.nc
fi
#
# *** if remap via Veronikas mapping
if (( ${IMAP} == 1 )) ; then
  set +e ; rm rmp_hd.nc ; set -e
  ln -sf ${HDFILE}/rmp_hd_${RES_INP}.nc rmp_hd.nc
fi
#
# ******* Handling of time info *****************************************

if (( ${IWORK} == 1 )) ; then
  ndays=365  
elif (( ${IWORK} == 2 )) ; then
  case ${MM} in
    "01" ) ndays=31;;
    "02" ) ndays=28;;
    "03" ) ndays=31;; 
    "04" ) ndays=30;; 
    "05" ) ndays=31;;
    "06" ) ndays=30;;
    "07" ) ndays=31;; 
    "08" ) ndays=31;; 
    "09" ) ndays=30;;
    "10" ) ndays=31;;
    "11" ) ndays=30;;
    "12" ) ndays=31;;
  esac
elif (( ${IWORK} == 3 )) ; then
  ndays=360  
elif (( ${IWORK} == 4 )) ; then
  ndays=$nday_final
fi

if (( ${IWORK} < 3 )) ; then
  #---- leap year  
  year_4=$((${YYYY}/4))
  year_1=$((${YYYY}-$year_4*4))
  if (( $year_1 == 0 )) ; then
    if (( ${MM} == 02 )) ; then ndays=29 ; fi
    if (( ${IWORK} == 1 )) ; then ndays=366 ; fi
  fi

  year_100=$((${YYYY}/100))
  year_1=$((${YYYY}-$year_100*100))
  if (( $year_1 == 0 )); then
    if (( ${MM} == 02 )) ; then ndays=28 ; fi
    if (( ${IWORK} == 1 )) ; then ndays=365 ; fi
  fi

  year_400=$((${YYYY}/400))
  year_1=$((${YYYY}-$year_400*400))
  if (( $year_1 == 0 )) ; then
    if (( ${MM} == 02 )) ; then ndays=29 ; fi
    if (( ${IWORK} == 1 )) ; then ndays=366 ; fi
  fi
fi
#
# ******* Calculation based on the number of time steps per day
nstep=$((${ndays}*${ndt_day}))
dt=$((86400/${ndt_day}))
echo $ndays ' --> Number of timesteps dt =' $dt ' : ' $nstep
#
# Time steering of run using date_start and date_end - supercedes month and year settings
if [ "$date_start" != "" ] && [ "$date_end" != "" ] ; then
  YYYY=0
  MM=0
fi
#
# ******* namelist.hd *************************************************
if (( $ICOUPLE != 2 )) ; then
   lcoup_out=.FALSE.
else
   lcoup_out=.TRUE.
fi
if (( $IBC_WRITE == 0 )) ; then
   lbc_write=.FALSE.
else
   lbc_write=.TRUE.
fi
#
rm -f namelist.hd
cat > namelist.hd << end_hdalone_ctl
&HD_CTL
  YEAR1 = ${YYYY}
  MONTH1 = ${MM}
  date_start = "$date_start"
  date_end = "$date_end"
                       !!   1h     2h     4h       6h      8h     12h     24h
  NSTEP = ${nstep}     !! 8760 ! 4380 !  2190 !  1460 !  1095 !   730 !   365
  delta_time = ${dt}   !! 3600 ! 7200 ! 14400 ! 21600 ! 28800 ! 43200 ! 86400
  runoff_file =   $DNRUN
  drainage_file = $DNDR
  forcing_freq = 0   !! 0: stepwise, 1: daily
  IOUT = 6           !! 5: monthly, 6: daily
  OUT_DATAPATH = "${HDOUT}"
  UFAKRU = ${UFAK}
  coupling_type = ${ICOUPLE}        ! 0=no, 1=no interpolation, 2=interpolation in HD
  coupling_file = "${DNCOUPLE}"
  lcoupling_out = ${lcoup_out}      ! Write discharge on ocean grid dep. on IOUT (couling_type 2 only)
  iform_input = ${iform}            ! Format Input files: 0 = SRV (Default), 1 = NetCDCF
  ibc_type = ${IBC_TYPE}            ! Bias correction type: 0=None, 1 = Mean Bias, 2 = Low, Mid and High Biases 
  dn_bcpara = "${DN_BCPARA}"
  lbc_write = ${lbc_write}
/
end_hdalone_ctl
#.................... NOTE for IOUT of HD ..................................
#!     ***   IOUT = Mittelungsartvariable, d.h. ueber wieviel Zeitschritte
#!     ***          1   30-Day Averages   --> NT = 30 * hd_steps_per_day
#!     ***          2   Decadal Averages  --> NT = 10 * hd_steps_per_day
#!     ***          3   Weekly Averages   --> NT = 7  * hd_steps_per_day
#!     ***          4   Monthly Averages ohne Schaltjahre
#!     ***          5   Monthly Averages inklusive Schaltjahre
#!     ***          6   Daily Output
#!     ***          -1  Mittelung ueber restliche Zeitschritte und Return
#
# *********** namelist.hdset *************************************************
rm -f namelist.hdset
cat > namelist.hdset << end_hydrology_ctl
&HYDROLOGY_CTL
  ldebughd =  .FALSE.
!  ldebughd =  .TRUE.
  diag_water_budget = .FALSE.
!!! locean = .FALSE. : closure of water budget for ocean coupling
  locean = .FALSE.
  nhd_diag = ${nhd_diag}
!!! lhd_highres = .TRUE.: CALL hd_highres_write --> /output/hd/hd_YYYY_MM_02_hd_highres.nc
  lhd_highres = .FALSE.
  fllog1 = ${xlon1} 
  fblog1 = ${xlat1}
  fllog2 = ${xlon2}
  fblog2 = ${xlat2}  
  nremap = ${IMAP}
  lhd_rout = .FALSE.
! Factors to allow sensitivity studies on 5 Min.
  fk_rfk = 1.
  fk_lfk = 1.
  fk_gfk = 1.
! Discharge dependent riverflow velocity (irf_vel <> 0)
  irf_vel = 1
  qrf_ref = 1000.
/
end_hydrology_ctl
#
#! nremap = Type of Interpolation from input (atmospheric) grid to HD grid
#!          0   Input = Output
#!          1   using HDMAP routine by Veronika (default)
#!          2   0.5 degree to 5 Min./
#
# *********** namelist.hdset *************************************************

rm -f namelist.hduser
cat > namelist.hduser << end_hduser_ctl
&HDUSER_CTL
  hd_user = "${HD_USER}"
  hd_cont = "${HD_CONT}"
  hd_inst = "${HD_INST}"
  hd_instid = "${HD_INSTID}"
/
end_hduser_ctl

#################################################
# run the program
#################################################

if [ -e ${HDDIR}/${EXP}/${EXP}_${EXE} ] ; then 
  ls -al ${HDDIR}/${EXP}/${EXP}_${EXE}
else
  ls -al ${HDMAIN}/${EXE}
  cp -p ${HDMAIN}/${EXE} ${HDDIR}/${EXP}/${EXP}_${EXE}
fi

echo "----- start HD with ${EXP}_${EXE}"
if (( ${CMPI} == 1 )) ; then
  srun --propagate=STACK --mpi=pmi2 ${HDDIR}/${EXP}/${EXP}_${EXE}
  status=$?
  if [ $status -ne 0 ] ; then
    echo "HD model crashed for ${EXP}"
    ls -lat
    exit $status
  fi
else
  srun --propagate=STACK ${HDDIR}/${EXP}/${EXP}_${EXE}
  status=$?
  if [ $status -ne 0 ] ; then
    echo "HD model crashed for ${EXP}"
    ls -lat
    exit $status
  fi
fi
echo "----- HD finished"
echo "--- Post-processing for HD is called now ---"
cd ${HDMAIN}/$EXP
./hd_post.ksh $EXP
#
exit
#
