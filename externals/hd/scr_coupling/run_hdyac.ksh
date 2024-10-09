#!/bin/ksh
#
# run_hdyac.ksh - HD model runscript in a coupled system with YAC
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Author: Ha Ho-Hagemann (Hereon, Germany)
# Contact: ha.hagemann@hereon.de
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
#SBATCH --mail-user=ha.hagemann@hereon.de
#SBATCH --account=gg0302
#SBATCH --partition=shared     # Specify partition name
#SBATCH --ntasks=3             # Specify max. number of tasks to be invoked
####SBATCH --ntasks=1             # Specify max. number of tasks to be invoked
# #SBATCH --mem=5G
# #SBATCH --mem-per-cpu=500      # Specify real memory required per CPU in MegaBytes
#SBATCH --mem-per-cpu=2000      # Specify real memory for global 5 min run or levante
#SBATCH --time=01:30:00        # Set a limit on the total run time (1:20 for G)
###############################################################################
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

WRKDIR=$(pwd)
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
CURRENT_DATE="0"

while [[ -n $1 ]] ; do
 CURRENT_DATE=$1
 shift
done

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** Include user settings ***
#
if [ -e ./job_settings ] ; then 
  source ./job_settings
else
  echo 'Script job_settings does not exist in present directory scr!'
  exit
fi

ierrset=0
if [ -e $HDMAIN ];then echo "$HDMAIN exists" ; else
  echo "$HDMAIN does NOT exists!"
  echo "   --> You need to create it"
  ierrset=1
fi
if [ -e $HDFILE ];then echo "$HDFILE exists" ; else
  echo "$HDFILE does NOT exists!"
  echo "   --> You need to create it"
  ierrset=1
fi
if [ -e $HDDIR ];then echo "$HDDIR exists" ; else
  echo "$HDDIR does NOT exists!"
  echo "   --> You need to create it"
  ierrset=1
fi
if [ -e $HDFORCING ];then echo "$HDFORCING exists" ; else
  echo "$HDFORCING does NOT exists!"
  echo "   --> You need to create it"
  ierrset=1
fi
if [ ${ierrset} == 1 ]; then
  echo "Some HD directories does not exist --> TERMINATION"
  export ierrset=1
else
  export ierrset=0
fi

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** Running mode defined ***
#
set -e
#
echo "YYYY:" ${YYYY}

if [ ${CURRENT_DATE} == "0" ];then
 runstate="start"
else
 if [ ${CURRENT_DATE} == ${YYYY} ];then
  runstate="start"
 else
  runstate="restart"
 fi
fi

echo "runstate:" ${runstate}

if [ ${runstate} == "start" ];then
 if [ -e ${WRKDIR}/log/${EXP}.year ] ; then
  rm -f ${WRKDIR}/log/${EXP}.year
 fi
fi

if [ -e ${WRKDIR}/log/${EXP}.year ] ; then
  YYYY=`cat ${WRKDIR}/log/${EXP}.year`
  INEU=0
else
  INEU=1
fi
echo 'Exp. ' $EXP ' and Year ' $YYYY

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** Paths setting ***
#
HDOUT=${HDDIR}/${EXP}/out          # HD Output dir

if [ $EXP != "hd_gl05deg" ] && [ $EXP != "hd_eu5min" ] ;then
 echo "please choose the supported option of scale and resolution"
 echo "see job_settings: EXPINP | IFORCE | HDRES | RES_INP"
 exit
fi

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** Prepare coupling.yaml ***
#
if [ $runmod == "CPL" ];then
 ./coupling_yaml_generate.sh ${runstate}	# ==> coupling_${EXP}.yaml
fi

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** Make output directory and set simulation year(s) ***
#
cd $HDDIR
rm -rf coredir*
#rm -f HD5.run.o*
set +e 
mkdir $EXP
cd $EXP
mkdir out

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** Linking executable directory of HD ***
#
if [ -e ${exe_dir} ];then
 unlink ${exe_dir}
fi

ln -sf $WRKDIR/BUILD/${exe_dir} .

if [ $runmod == "CPL" ];then
 if [ -e toy_ocean.exe ];then
  unlink toy_ocean.exe
 fi
 if [ -e toy_land.exe ];then
  unlink toy_land.exe
 fi
 ln -sf $WRKDIR/BUILD/toy_ocean.exe .
 ln -sf $WRKDIR/BUILD/toy_land.exe .
fi

if [ $runmod == "CPL" ];then
 if [ -e coupling.yaml ];then
  unlink coupling.yaml
 fi

 cp ${WRKDIR}/coupling_${EXP}.yaml coupling_${YYYY}.yaml
 ln -sf ${WRKDIR}/coupling_${EXP}.yaml coupling.yaml
fi

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** Settings for forcing data ***
#
case ${IFORCE} in
  1     ) if (( $HDRES != 0 )) && [ $EXPINP -ge 55012 ] ; then RES_INP='5min' ; fi ;;
  3 | 4 ) if (( $HDRES != 0 )) ; then RES_INP='5min' ; fi ;;
esac
echo 'Resolution of forcing: ' $RES_INP
DNREMAP="${EXPINP}_to_${RES_INP}"
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

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** Settings for HD model resolution ***
#
IZIP_START=1
case ${HDRES} in
   0  )  IMAP=0              # Remapping type (0=no, 1=vg, 2=05-->5min
         ulimit -s 102400
         IZIP_START=0
         DNPARA="05deg/hdpara_vs1_11.nc"
         HDSTART="05deg/hdstart_05.nc"
         nhd_diag=7   # Log output for Elbe river on 0.5°
         GRID="grid_0_5.txt"
         DNMAS="05deg/masks_${RES_INP}.nc"
         EXE="${exe_dir}/hd_05.exe" ;;
   1  )  if [ $RES_INP = '05' ] ; then
           IMAP=2
         else
           IMAP=0
         fi
         DNPARA="5min/hdpara_vs5_1.nc"
         HDSTART="5min/hdstart_5min.nc"
         nhd_diag=7   # Log output for Elbe river on 5 Min.
#
#        *** Examples if restart files from previous runs are used for initialization.
         if (( $YYYY > 1979 )) ; then 
            HDSTART="hdstart/7055116/7055116_hdrestart_${YYYY}01.nc" 
            IZIP_START=1
         fi
         GRID="grid_5min.txt"
         DNMAS="5min/masks_${RES_INP}.nc"
         EXE="${exe_dir}/hd_5min.exe" ;;
   2  )  if [ $RES_INP = '05' ] ; then
           cdocom="${cdocom} -selindexbox,339,498,37,126"   # = on 5 min: 2029,2988,217,756
           IMAP=2
         else
###           cdocom="setmisstoc,0. -selindexbox,2029,2988,217,756"
           IMAP=0
         fi
         IZIP_START=1
         nhd_diag=7   # Log output for Elbe river on 5 Min.
         ulimit -s 102400
         DNPARA="euro5min/hdpara_vs5_1_euro5min.nc"
         HDSTART="euro5min/hdstart_euro5min.nc"
#
#        *** Examples if restart files from previous runs are used for initialization.
         case $EXP in
           70620?? | 706220? | 706021[03]) if (( $YYYY <= 1979 )) ; then 
                       HDSTART="hdstart/hdstart_euro5min_7062008_${YYYY}01.nc" 
                       IZIP_START=0
          	     else 
                       HDSTART="hdstart/hdstart_euro5min_7062029_${YYYY}01.nc"  # CD3 & HD Vs 5.1
                     fi ;;
           7055037 ) echo 'use default start file' ;;
           * ) if (( $YYYY > 1901 )) ; then HDSTART="hdstart/hdstart_euro5min_7062029_${YYYY}01.nc" ; fi ;;  # GWSP3 based
         esac
         GRID="grid_euro5min.txt"
         DNREMAP="${EXPINP}_to_euro5min"
         DNMAS="euro5min/masks_euro${RES_INP}.nc"
         EXE="${exe_dir}/hd_euro.exe" ;;
   3  )  if [ $RES_INP = '05' ] ; then
           IMAP=2
           cdocom="${cdocom} -selindexbox,585,668,201,270"   # = on 5 min: 3505,4008,1201,1620
         else
           IMAP=0
         fi
         ulimit -s 102400
         DNPARA="aus/hdpara_vs4d_aus.nc"
         HDSTART="aus/hdstart_aus.nc"
         if (( $YYYY == 1990 )) ; then 
           HDSTART="hdstart/hdstart_aus_7062010_199001.nc"
         elif (( $YYYY != 1979 )) ; then 
           HDSTART="hdstart/hdstart_aus_7055051_${YYYY}01.nc" 
         fi
         GRID="grid_aus.txt"
         DNMAS="aus/masks_aus${RES_INP}.nc"       # Comprises Input data dimensions!
         nhd_diag=100
         xlon1=147.24362     # Log output for gridboxes if nhd_diag=100 
         xlat1=-19.75856
         xlon2=139.6157
         xlat2=-34.3509  
         EXE="${exe_dir}/hd_aus.exe" ;;
    * ) echo 'Resolution No. ' $HDRES ' not defined --> STOP!' ; exit ;;
esac
#
  echo 'Start file: '$HDSTART ' in year' $YYYY
#
# *** Monthly runs --> read month info
if (( ${IWORK} == 2 )) ; then
  if (( ${INEU} == 0 )) ; then
    MM=`cat ${WRKDIR}/log/${EXP}.month`
  fi
  echo HD RESTART ' for Exp. ' $EXP ' and Month ' $MM
  cdocom="${cdocom} -selmon,$MM"
fi

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** Prepare HD forcing data ***
#
# Copy forcing preparation scripts to run directory to avoid that changes in these scripts affect running simulations
if (( ${INEU} == 1 )) ; then
  cp -p ${HDMAIN}/scr/prepare_hdforcing.ksh ${HDDIR}/${EXP}/.
fi
#
# Export variables used by script prepare_hdforcing.ksh
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
# *** Linking of necessary HD model files ***
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

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** time info ***
#
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
fi

if (( ${IWORK} != 3 )) ; then
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
# *** Calculation based on the number of time steps per day ***
nstep=$((${ndays}*${ndt_day}))
dt=$((86400/${ndt_day}))
echo $ndays ' --> Number of timesteps dt =' $dt ' : ' $nstep

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** namelist.hd ***
#
if (( $ICOUPLE != 2 )) ; then
   lcoup_out=.FALSE
else
   lcoup_out=.TRUE
fi
rm -f namelist.hd
cat > namelist.hd << end_hdalone_ctl
&HDALONE_CTL
  YEAR1 = ${YYYY}
  MONTH1 = ${MM}
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

echo ----- start HD with ${EXE}
if (( ${CMPI} == 1 )) ; then
  srun --propagate=STACK --mpi=pmi2 ${HDMAIN}/${EXE}
else
 if [ $runmod == "UNCPL" ];then
# Stand-alone HD:
  srun --propagate=STACK ${EXE}
 else
# Coupled HD + toy_ocean + toy_land
  echo "0 ${EXE}" > mpmd.lst
  echo "1 toy_ocean.exe" >> mpmd.lst
  echo "2 toy_land.exe" >> mpmd.lst

  srun -l \
    -l --hint=nomultithread --distribution=block:cyclic \
    --multi-prog mpmd.lst  || {
  exit 1
  } 
 fi

fi
echo "----- HD finished"

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# *** Postprocessing for HD ***
echo "--- Post-processing for HD ---"
cd ${WRKDIR}
./post_HD.sh ${runstate}

