#!/bin/ksh -l
#   
# convert_inflow.com - Script to convert HD discharge to an ocean grid with predefined options
# 
# Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
# *** Script/Program to generate a transformation matrix that associates each source (e.g. HD)
#     mouth point with a corresponding mouth point on the target ocean grid. 
#     Here, the source coordinates are directly linked to the target mouth points, 
#     i.e. no remapping required.
#
# Parse command line parameters
#
set -e
#
while [[ -n $1 ]] ; do
  case $1 in
    --id        | -i    ) EXP=$2                 ; shift    ;;
    --bgc       | -b    ) IBGC=$2                ; shift    ;;
    --mode      | -m    ) IMODE=$2               ; shift    ;;
    --ybeg      | -y    ) YBEG=$2                ; shift    ;;
    --yend      | -z    ) YEND=$2                ; shift    ;;
    --source    | -s    ) ISRC=$2                ; shift    ;;
    --umask     | -u    ) dn_src_mouth=$2        ; shift    ;;
    --ocean     | -o    ) IOCEAN=$2                          ;;
    -*                  ) err_strg="\nERROR: invalid option $1 !\n"
                          [[ "$1" != $( echo "$1" | tr -d = ) ]] && \
                          err_strg=$err_strg"  Please use blank insted of '=' to separate option's argument.\n"
                          usage ;;
    *                   ) export CPL_MOD=$1 ;;
  esac
  shift
done
#
# *** Set Experiment ID / List 
if [ "$IOCEAN" = "" ] && [ "$ISRC" = "" ]; then
   echo "Neither source model nor Ocean model are set. Please use the following options"
   echo "     -m <Mode No.>,     1=prescr. mouth mask, 2=LS mask (Def), 3=1&2, 4=pre-1&2"
   echo "     -i <7-digit ID>,   ID = Source exp. no., Default: Last Exp. no. "
   echo "     -b <BGC mode>,     0=No data conversion, 1=Discharge only (Def.), 2=1 or more tracers"
   echo "                        3=Bias corrected discharge, 4=Bias corrected tracer flow (planned)"
   echo "     -s <Source>,       1=HD4, 2=HD5 (Def), 3=HD1.10, 4=mHm, 5=Ute, 6=MPIOM-BGC, 7=MPIOM"
   echo "     -u <mouth file>    optional: mouth mask on source grid"
   echo "     -o <Ocean target>, 1=NEMO3.6, 2=ECOSMO3, 3=SCHISM, 4=ECOSMO2, 5=ICON, 6=Nils, 7=ICON-C "
   echo "                        8=NEMO4.0, 9=TRIM,10=NEMO-med, 11=HD4, 12=HD5, 13=HD1.10, 14=MOM, 15=NEMO-NSS"
   echo "     -y  <First Year>   e.g. YYYY, Default: 1979, if 0 use preset values"
   echo "     -z  <Last Year>    e.g. YYYY, Default: <First Year>"
   exit
fi

#
# *** Set Options for possible ocean model mouth points 
set +x
if [ "$IMODE" = "" ]; then
   echo "ocean model mouth option is not set. Please use Option -n <IMODE>"
   echo "  --> Default setting IMODE = 2 is used"
   echo "IMODE = 1    Use existing mask with potential ocean model mouth points on ocean model grid"
   echo "              e.g. coast_ocean_NEMO.nc"
   echo "IMODE = 2    Generate mask of coastal ocean points from ocean model sea mask (Default)"
   echo "IMODE = 3    Use both an existing mask and the coastal ocean points derived from ocean model sea mask"
   echo "IMODE = 4    Use both an existing primary mask and secondary mask of coastal ocean points"
   IMODE=2
else
   case $IMODE in
     1  ) echo "Existing mouth mask is used"  ;;
     2  ) echo "Mouths are generated from ocean model sea mask"  ;;
     3  ) echo "Combine existing mouth mask & mouths generated from ocean model sea mask"  ;;
     4  ) echo "Combine existing primary mouth mask & secondary coastal ocean mask"  ;;
     *  ) echo "ERROR - Option does not exist: " $IMODE ; exit
   esac
fi
#
# *** Set Options for discharge only or for one/several bgc tracers 
CFLOW='discharge'
if [ "$IBGC" = "" ]; then
   echo "BGC tracer option is not set. Please use Option -b <IBGC>"
   echo "  --> Default setting IBGC = 1 is used: Discharge only "
   echo "                      IBGC = 0:         No data conversion"
   echo "                      IBGC = 2:         One or more tracers"
   echo "                      IBGC = 3:         Bias Corrected Discharge"
##   echo "                      IBGC = 4:         Bias Corrected Tracer Flow for one or more tracers"
   IBGC=1
else
   case $IBGC in
     0  ) echo "No data conversion"  ;;
     1  ) echo "Conversion for discharge only"  ;;
     2  ) echo "Conversion for one or more bgc tracers"
          CFLOW='tracer_flow'  ;;
     3  ) echo "Conversion for bias coorected discharge"
          CFLOW='discharge_bc' ;;
##     4  ) echo "Conversion for Bias corrected tracer flow"  ;;
     *  ) echo "ERROR - Option -b IBGC does not exist: " $IBGC ; exit
   esac
fi
#
# *** Set Begin and End year 
if [ "$YBEG" = "" ] ; then echo "First year is not set --> use default: 1979" ; Y1=1979 ; YBEG=0 ; fi
if [ "$YEND" = "" ] ; then 
  echo "Last year is not set --> use first year"  
  if [ $YBEG -eq 0 ] ; then Y2=$Y1 ; else ; Y2=$YBEG ; fi
fi
#
#  *** Set Options for source grid with river mouths.
if [ "$ISRC" = "" ]; then
   echo "Source grid option is not set. Please use Option -s <ISRC>"
   echo "  --> Default setting ISRC = 2 is used"
   echo "ISRC = 1    HD Model Vs. 4"
   echo "ISRC = 2    HD Model Vs. 5.1"
   echo "ISRC = 3    HD Model Vs. 1.11"
   echo "ISRC = 4    mHm Vs. 2"
   echo "ISRC = 5    Utes standard input data available from 1940-2016"
   echo "ISRC = 6    MPIOM BGC data, e.g. Fabrice's BGC data"
   echo "ISRC = 7    MPIOM runoff flux data, e.g. REMO-MPIOM"
   ISRC=2  ; HDVS=vs5_1 ; CHD=hd5_1
else
   case $ISRC in
     1  ) echo "Source data: HD Model Vs. 4"      ; CHD=hd4 ; HDVS=vs4    ;;
     2  ) echo "Source data: HD Model Vs. 5"      ; CHD=hd5_1 ; HDVS=vs5_1  ;;
     3  ) echo "Source data: HD Model Vs. 1.11"   ; CHD=hd1_11 ; HDVS=vs1_11 ;;
     4  ) echo "Source data: mHm Vs. 2"           ; CHD=mhm ; HDVS=mhm    ;;
     5  ) echo "Source data: Utes standard input" ; CHD=usid ; HDVS=ute    ;;
     6  ) echo "Source data: MPIOM TP04 data"     ; CHD=mpiom ; HDVS=TP04    ;;
     7  ) echo "Source data: MPIOM data"          ; CHD=mpiom_uhh ; HDVS=UHH    ;;
     *  ) echo "ERROR - Source grid option does not exist: " $ISRC ; exit  ;;
   esac
fi
#
#  *** Is mask file for the mouths on the source grid provided.
if [ "$dn_src_mouth" = "" ]; then
   echo " No mouth mask on source grid provided - defintion in script expected "
else
   if [ -s $dn_ocean ] ; then 
     echo " Mouth mask on source grid: " $dn_src_mouth
   else
     echo " Mouth mask on source grid does not exist: " $dn_src_mouth
     exit
   fi
fi
#
#  *** Set Options for ocean model.
if [ "$IOCEAN" = "" ]; then
   echo "Ocean model option is not set. Please use Option -o <IOCEAN>"
   echo "  --> Default setting IOCEAN = 0 is used"
   echo "IOCEAN = 1    NEMO Vs. 3.6"
   echo "IOCEAN = 2    ECOSMO Vs. 3"
   echo "IOCEAN = 3    SCHISM"
   echo "IOCEAN = 4    ECOSMO Vs. II"
   echo "IOCEAN = 5    ICON"
   echo "IOCEAN = 6    Nils Nordsee-Modell"
   echo "IOCEAN = 7    ICON-Coast"
   echo "IOCEAN = 8    NEMO Vs. 4.0"
   echo "IOCEAN = 9    TRIM"
   echo "IOCEAN = 10   NEMO-med7KM"
   echo "IOCEAN = 11    HD Model Vs. 4"
   echo "IOCEAN = 12    HD Model Vs. 5.1"
   echo "IOCEAN = 13    HD Model Vs. 1.11"
   echo "IOCEAN = 14   MOM from IOW"
   echo "IOCEAN = 15   NEMO-NSS"
   IOCEAN=0
else
   case $IOCEAN in
     1  ) echo "NEMO Vs. 3.6 is used as ocean model"  ;;
     2  ) echo "ECOSMO Vs. 3 is used as ocean model"  ;;
     3  ) echo "SCHISM is used as ocean model"  ;;
     4  ) echo "ECOSMO Vs. II is used as ocean model"  ;;
     5  ) echo "ICON is used as ocean model"  ;;
     6  ) echo "Nils Nordee-Modell is used as ocean model"  ;;
     7  ) echo "ICON-Coast is used as ocean model"  ;;
     8  ) echo "NEMO Vs. 4.0 is used as ocean model"  ;;
     9  ) echo "TRIM is used as ocean model"  ;;
    10  ) echo "NEMO-med7km is used as ocean model"  ;;
    11  ) echo "Mouths of HD Model Vs. 4"        ;;
    12  ) echo "Mouths of HD Model Vs. 5"        ;;
    13  ) echo "Mouths of HD Model Vs. 1.10"     ;;
    14  ) echo "MOM is used as ocean model"  ;;
    15  ) echo "NEMO-NSS is used as ocean model"  ;;
     *  ) echo "ERROR - Ocean model Option does not exist: " $IOCEAN ; exit  ;;
   esac
fi
if [ "$EXP" = "" ]; then
   echo "Exp. No. ID is not set. Please use Option -i <ID>"
     # Exp. number of HD model run, e.g. 7055056, 7055086, 7055045 (HD1.10), 7056101 Utes data    
     # 7056100 E-Hype 7056002 mhm, bgc002
   case $ISRC in
     1  ) EXP=7055056  ;;  # HD vs4 Europe
     2  ) EXP=7055115  ;;  # HD vs5_0 Europe
     3  ) EXP=7055045  ;;  # HD vs1_10
     4  ) EXP=7056002  ;;  # mhm
     5  ) EXP=7056101  ;;  # Utes std. input data
     6  ) echo "No default for MPIOM TP04 data --> EXIT" ; exit  ;;
     7  ) echo "No default for MPIOM data --> EXIT" ; exit  ;;
   esac
fi
echo "    --> Exp. No. ' $EXP ' is used"
set -x
#
# *********** Data for Program *************************
# *** This needs to be EDITED  *************************
DSRC=/home/g/g260122/hdmodel/util      # Dir. with src files
DIN=/work/gg0302/g260122/HD/input
DRUN=/work/gg0302/g260122/HD/data
#
# *** Conversion of discharge or bgc flow in one go? 
ICONV=1  # 0=NO conversion, 1=daily data, 2=monthly climatology, 3=sequence of 5 bgc fields
         # 4 = annual sequence of several bgc variables and time steps, e.g. monthly
ISEPIN=1      # Inflow in separate files per year (No/yes = 0/1)
ISEPOUT=1     # Flow on ocean in separate files per year (No/yes = 0/1)
ICAL=0        # Calendar: 0=normal, 1=360 day years
IMACH=2       # Machine: 1 Mistral, 2 Levante
#
DATA=/work/gg0302/g260122/HD/output/$EXP   # Discharge data directory

case $IBGC in
  0 ) ICONV=0  ;;
  1 ) DNRIV=${DATA}/${EXP}_meanflow_YYYY.nc  ;;    # Def. HD name ${EXP}_meanflow_YYYY.nc for ISEPIN=1
  2 ) DNRIV=${DATA}/${EXP}_tracer_flow_YYYY.nc  ;;    # Def. HD name ${EXP}_meanflow_YYYY.nc for ISEPIN=1
  3 ) DNRIV=${DATA}/bc/${EXP}_bcflow_YYYY.nc  ;;    # Def. HD name ${EXP}_meanflow_YYYY.nc for ISEPIN=1
esac
#
case $EXP in 
  7055045 ) CHD=hd1_10 ; CDIS=jsbach_$CHD ; Y1=1979 ; Y2=2009  ;;
  7055055 ) CHD=hd4 ; CDIS=jsbach_$CHD ; Y1=1979 ; Y2=1979 ;;
  7055056 ) CHD=hd4 ; CDIS=gswp3_mpihm_$CHD ; Y1=1901 ; Y2=2014 ;;
  7055086 ) CHD=hd5b ; CDIS=era5_hydropy_$CHD ; HDVS=vs5b  ;;
  7055104 ) CHD=hd5b ; CDIS=era5_hydropy_$CHD ; HDVS=vs5b  ;;
  7055112 ) CHD=hd5_0 ; CDIS=eobs20_hydropy_$CHD ; HDVS=vs5_0 ; Y1=1950 ; Y2=2019  ;;
  7055115 ) CHD=hd5_0 ; CDIS=era5_hydropy_$CHD ; HDVS=vs5_0 ; Y1=1979 ; Y2=2018  ;;
  7055119 ) CHD=hd1_10 ; CDIS=gswp3_hydropy_$CHD ; Y1=1979 ; Y2=2014  ;;
  7055136 ) CHD=hd5_1 ; CDIS=era5_hydropy_$CHD ; Y1=1979 ; Y2=2018  ;;
  7055144 ) CHD=hd5_1 ; CDIS=gswp3_hydropy_$CHD ; Y1=1979 ; Y2=2014  ;;
  7055145 ) CHD=hd5_1 ; CDIS=wfde5_hydropy_$CHD ; Y1=1979 ; Y2=2014  ;;
  7055156 ) CHD=hd5_1 ; CDIS=gswp3_hydropy_$CHD ; Y1=1901 ; Y2=2014 ;;
  7055157 ) CHD=hd1_11 ; CDIS=gswp3_hydropy_$CHD ; Y1=1901 ; Y2=2014 ;;
  7056002 ) CHD=mhm ; CDIS=${CHD}_vs2 
            DRUN=/work/gg0302/g260122/mhm/data  
            DATA=/work/gg0302/g260122/mhm/output/$EXP  
            case $IBGC in
              1 ) Y1=1950 ; Y2=2019
                  DNRIV=${DATA}/${EXP}_meanflow_YYYY.nc  ;;   
              2 ) Y1=1960 ; Y2=2016
                  DNRIV=${DATA}/${EXP}_tracer_flow_YYYY.nc  ;;
              * ) echo "ERROR - Option IBGC not valid for mhm: " $IBGC ; exit
            esac  ;;
  7056101 ) CHD=ute ; CDIS=${CHD}_stdinput ; ISEPIN=0 ; CFLOW='inflow' ;; 
  7060210 ) Y1=1950 ; Y2=2005 ; CDIS=remo_hadgem2_${CHD} ; ICAL=1    ;; 
  706021[12] ) Y1=2006 ; Y2=2099 ; CDIS=remo_hadgem2_${CHD} ; ICAL=1 ;;
  7060213 ) Y1=1950 ; Y2=2005 ; CDIS=remo_mpiesm_${CHD}     ;; 
  706021[456] ) Y1=2006 ; Y2=2100 ; CDIS=remo_mpiesm_${CHD} ;;
  7062008 ) CHD=hd4 ; CDIS=coastdat2_$CHD  ;;
  7062018 ) CHD=hd4 , CDIS=coastdat3_$CHD  ;;
  ??56100 ) CHD=ehype ; CDIS=${CHD}_on_hd5 
            if [ $EXP -eq 7156100 ] ; then CFLOW='nflow' ; fi
            if [ $EXP -eq 7256100 ] ; then CFLOW='pflow' ; fi
            DATA=/work/gg0302/g260122/bgc/ehype/${EXP}
            DNRIV=${DATA}/${EXP}_discharge_on_hd_YYYY.nc ;;
   bgc00? ) ISEPIN=0
            CFLOW="${EXP}_flow" 
            DATA=/work/gg0302/g260122/HD/bgc
            case $ISRC in
              3 ) CHD=hd1 ; CDIS=on_${CHD}
                  case $EXP in
                    bgc001 ) DNBGC=${DATA}/hd/bgc001_flow_mpiom_tp04_on_hdvs1_10_1850.nc ;;
                    bgc002 ) DNBGC=${DATA}/hd/bgc002_flow_mpiom_tp04_on_hdvs1_10_1900-2010.nc  ;;
                  esac ;;
###                  DNBGC=${DATA}/${EXP}_fabrice_flow_on_hd_vs1.nc ;; # 
              6 ) CHD=mpiom ; CDIS=${CHD}_tp04
                  case $EXP in
                    bgc001 ) DNBGC=${DATA}/river_input_TP04_preindustrial.nc ;;
                    bgc002 ) DNBGC=${DATA}/river_input_TP04_DIN_DIP_1900-2010.nc  ;;
                  esac ;;
              * ) echo 'Source format incompatible with Exp. ' $EXP ' --> Stop!' ; exit ;;  
            esac
            DNRIV=${TMPDIR}/bgc_inflow.nc
            cdo setmisstoc,0. $DNBGC $DNRIV
            case $EXP in
              bgc001 ) ICONV=3 ; Y1=1850 ; Y2=$Y1 ;;
              bgc002 ) ICONV=4  
                       Y1=1900 ; Y2=$Y1 ;;    # All years are in single file -->set only 1 year
            esac ;;
    * ) CDIS=${CHD}_${EXP} ;;
esac
#
# Use preset values for YBEG and YEND
if [ $YBEG -eq 0 ] ; then
  YBEG=$Y1
  YEND=$Y2
fi
#
if [ $ICONV -eq 2 ] ; then
  ISEPIN=0
  DNRIV=${DATA}/mean_${EXP}_${YBEG}-${YEND}.nc  # 
fi

#
cd $DRUN

case $ISRC in
   1   )   DNHDPARA=${DIN}/euro5min/hdpara_vs4d_euro5min.nc ; cgr='5-Min.'  ;;    # cgr is currently a dummy
   2   )   DNHDPARA=${DIN}/euro5min/hdpara_${HDVS}_euro5min.nc ; cgr='5-Min.'  ;;   # e.g. HDVS=vs5_0 
#   3   )   DNHDPARA=${DIN}/05deg/hdpara_vs1_10_ext.nc ; cgr='0.5' ;; 
   3   )   DNHDPARA=${DIN}/05deg/hdpara_vs1_11.nc ; cgr='0.5' ;; 
   4   )   DNHDPARA=/work/gg0302/g260122/mhm/input/mask_catchment_full_${CDIS}.nc ; cgr='0.0625' ;; 
   5   )   DATA=/work/gg0302/g260122/bgc/ute
           DNHDPARA=${DATA}/river_loads_and_discharge_NSBS_ECOSMO1.nc
           DNRIV=$DNHDPARA
           ISEPIN=0 ;;
   6   )   DNHDPARA=${DATA}/mouth_${EXP}_${CDIS}.nc 
           cdo timsum -nec,0. $DNBGC ts.nc
           cdo gtc,0. -varssum ts.nc $DNHDPARA  ;;
   7   )   if [ "$dn_src_mouth" = "" ]; then
              echo 'Mouth mask of source grid not specified (Option -u) --> Script terminates!'
              exit
           fi 
           DNHDPARA=${dn_src_mouth}   ;; 
esac
#
# ******* End of source related definitions
OVS=''
case $IOCEAN in
  1 | 8 ) OM='nemo'
        if [ $IOCEAN -eq 8 ] ; then OVS='_vs4' ; fi
        DNAM=NEMO${OVS}_grid_GCOAST.nc                    # Mask file if mouths are generated (0 for mouth)
        DN_OCEANMOUTH=${DIN}/${OM}/gcoast_mouth.nc  ;;    # Vs. 3.3: coast_ocean_NEMO.nc
  10 )  OM='nemo' ; OVS='-med7km'
        DNAM=lsm_med7km.nc                          ;;    # NEMO Sea-Land Mask (1/0) 
  15 )  OM='nemo' ; OVS='-nss'
        DNAM=nemo_nss_mask.nc                       ;;    # NEMO Sea-Land Mask (1/0) 
  2  )  OM='ecosmo'  ; OVS='_hr3' ; DNAM='oceanmask_ecosmo_hr3.nc' ;;  # Mask file if mouths are generated
  3  )  OM='schism'  ; DNAM='schism_river_points_db.nc'      # Mask file if mouths are generated
        IMODE=1 ; echo "!!! For SCHISM, only IMODE = 1 is used !!!" 
        DN_OCEANMOUTH=${DIN}/${OM}/$DNAM ;;
##        DN_OCEANMOUTH=${DIN}/${OM}/schism_mouth.nc ;;
  4  )  OM='ecosmo'  ; OVS='_ii' ; DNAM='oceanmask_ecosmo_ii.nc'  # Mask file if mouths are generated
        if [ $IMODE -eq 1 ] ; then OVS=${OVS}_mouth_ute ; fi
        DN_OCEANMOUTH=${DIN}/${OM}/mask_ute_e10.nc ;;
  5  )  OM='icon'   ; OVS='omip_60km' ; DNAM='icon_river_points.nc'   
        IMODE=1 ; echo "!!! For ICON, only IMODE = 1 is used !!!" 
##        DN_OCEANMOUTH=${DIN}/${OM}/omip_60km_surface_cell_sea_land_mask.nc ;; # (-1 = coastal ocean)
        DN_OCEANMOUTH=${DIN}/${OM}/omip_60km_lsm_c_lev6.nc ;; # (-1 = coastal ocean)
  6  )  OM='nsea'  ; OVS='3' ; DNAM=''        # Vs. 1 was nsea
        IMODE=1 ; echo "!!! For NSEA model, only IMODE = 1 is used !!!" 
        DN_OCEANMOUTH=${DIN}/${OM}/${OM}${OVS}_coastline_mask.nc ;;
  7  )  OM='icon' ; OVS='_omip_244_500' ; DNAM='icon_river_points.nc'   
        DN_OCEANMOUTH=${DIN}/${OM}/omip_244_500_lsm_c_lev8.nc  # (-1 = coastal ocean)
        case $IMODE in
          1 ) echo 'IMODE = 1 uses only ' $DN_OCEANMOUTH ;;
          4 ) DN_OMOUTH_PRESET=${DIN}/${OM}/bgc002_on_hd1_inflow_mask_on_iconomip_244_500.nc ;;
          * ) echo "!!! For ICON, only IMODE = 1 & 4 are allowed -> STOP !!!" ; exit ;;
        esac ;;
  9  )  OM='trim'  ; DNAM='trim_seamask_rea6_lidia.nc' ;;  # Mask file if mouths are generated
 14  )  OM='mom'   ; DNAM='iow_mom_slm.nc' ;;              # MOM Sea-Land Mask (1/0) 

 1[123] ) case $IOCEAN in
            11 )  dnhd=${DIN}/hdpara_vs4d_euro5min.nc ; OVS=vs4    ; area_min=10000.  ;; 
            12 )  dnhd=${DIN}/hdpara_vs5_1_euro5min.nc ; OVS=vs5_1    ; area_min=10000.  ;;
            13 )  dnhd=${DIN}/hdpara_vs1_11.nc    ; OVS=vs1_11 ; area_min=50000. ;; 
          esac
        OM='hd'
        DN_OCEANMOUTH=${DIN}/mouth_${OM}_${OVS}.nc  
        cdo eqc,0 -selvar,FDIR $dnhd $DN_OCEANMOUTH 
        case $IMODE in
          1 ) echo 'IMODE = 1 uses only ' $dnhd ;;
          4 ) DN_OMOUTH_PRESET=${DIN}/mouth_${OM}_${OVS}_large_rivers.nc
              cdo mul -eqc,0 -selvar,FDIR $dnhd -selvar,CAT_AREA $dnhd ts.nc
              cdo gec,$area_min ts.nc $DN_OMOUTH_PRESET
              echo 'What about discharge crit.' ;;
          * ) echo "!!! For HD mouths, only IMODE = 1 & 4 are allowed -> STOP !!!" ; exit ;;
        esac ;; 
  *  ) echo "ERROR - Ocean model Option does not exist: " $iocean ; exit  ;;
esac
DN_OMASK=${DIN}/${OM}/$DNAM         # Ocean mask: 1= ocean, 0=land and missing values
#
# Target directory for ocean data: sub directory of input data dir.
case $IBGC in
  3  ) DATA_OUT=${DATA}/bc/$OM  
       if [ $IOCEAN -eq 7 ] ; then DATA_OUT=${DATA}/bc/$OM-coast ; fi  ;; 
  *  ) DATA_OUT=${DATA}/$OM    
       if [ $IOCEAN -eq 7 ] ; then DATA_OUT=${DATA}/$OM-coast ; fi   ;;  
esac

#DATA_OUT=$TMPDIR  # for testing


if [ $ISEPOUT -eq 1 ] ; then 
  DNOUTFLOW=${DATA_OUT}/${CFLOW}_${CDIS}_on_${OM}${OVS}_YYYY.nc
else  
  if [ $YBEG -eq $YEND ] ; then
    DNOUTFLOW=${DATA_OUT}/${CFLOW}_${CDIS}_on_${OM}${OVS}_${YBEG}.nc
  else
    DNOUTFLOW=${DATA_OUT}/${CFLOW}_${CDIS}_on_${OM}${OVS}_${YBEG}-${YEND}.nc
  fi
fi
#
# ******************************************************
#
# *** Compiling
case $IMACH in
  1 ) source /sw/rhel6-x64/etc/profile.mistral     # Necessary for module command
      [[ -n `whence ifort` ]] && module unload intel
      module load intel/18.0.4
      LPATH="/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/lib"
      LIBS="-lnetcdf"
      NC_INCLUDE="-I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/include"
      NC_LIB="`/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/bin/nf-config --flibs`"
      NC_INC1="`/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-gcc48/bin/nc-config --cflags`" 
      NC_INC2="`/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/bin/nf-config --fflags`"   ;;

# Levante: requires module load netcdf-fortran/4.5.3-intel-oneapi-mpi-2021.5.0-intel-2021.5.0   for netcdf   
#          module load netcdf-c/4.8.1-intel-oneapi-mpi-2021.5.0-intel-2021.5.0                  for netcdf-c
  2 ) module load netcdf-fortran/4.5.3-intel-oneapi-mpi-2021.5.0-intel-2021.5.0
      module load netcdf-c/4.8.1-intel-oneapi-mpi-2021.5.0-intel-2021.5.0
      LPATH="/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/lib" 
      LIBS="-lnetcdf"
        NC_INCLUDE="-I/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/include"
      NC_LIB1="`/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/bin/nf-config --flibs`"
      NC_LIB="$NC_LIB1 -Wl,-rpath,$LPATH"
      NC_INC1="`/sw/spack-levante/netcdf-c-4.8.1-7dq6g2/bin/nc-config --cflags`" 
      NC_INC2="`/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/bin/nf-config --fflags`"  ;;
esac 
NC_INCLUDE="$NC_INC2 $NC_INC1"
#
F90=ifort
#
# ******************************************************
#
# *** Write namelist
if [ $ICONV -gt 0 ] ; then
cat > convert_inflow_ctl.nml << EOF
&convert_inflow_ctl
  dn_inflow = "$DNRIV"
  dn_outflow = "$DNOUTFLOW"
  iconv = $ICONV              ! Method of conversion
  isep_in = $ISEPIN           ! Input data in separate files (No/Yes=0/1)
  isep_out = $ISEPOUT         ! Output data in separate files (No/Yes=0/1)
  ybeg = $YBEG                ! Start year of data
  yend = $YEND                ! End year of data
/
EOF
  if [ -e $DATA_OUT ] ; then echo 'Output directory exists: ' $DATA_OUT ; else
  mkdir $DATA_OUT
  fi
else
  set +e
  rm convert_inflow_ctl.nml
  set -e
fi
#
set +e
rm -f convert.exe gen_mouth.exe
##rm ${OM}_mask.nc ${OM}_grid.nc mask_hd.nc
##rm hdmouth_mask.nc rivmouth_source.nc hdmouth_on_${OM}.nc
rm coast_oceangrid.nc
set -e

pwd
#
# *** Select/Generate mouths on ocean model grid.
case $IMODE in 
 1 | 4 ) case $IOCEAN in
        [18] ) cdo setvar,mask_coast_ocean -eqc,0 ${DN_OCEANMOUTH} coast_oceangrid.nc  ;;
#          10 ) cdo setvar,mask_coast_ocean -eqc,0 ${DN_OCEANMOUTH} coast_oceangrid.nc  ;;
        [34] ) ln -s ${DN_OCEANMOUTH} coast_oceangrid.nc  ;;
           5 ) cdo setvar,mask_coast_ocean -eqc,-1 ${DN_OCEANMOUTH} coast_oceangrid.nc  ;;
           6 ) cdo chname,mask,mask_coast_ocean ${DN_OCEANMOUTH} coast_oceangrid.nc  ;;
           7 ) cdo setvar,mask_coast_ocean -eqc,-1 ${DN_OCEANMOUTH} coast_oceangrid.nc  ;;
      1[123] ) cdo setvar,mask_coast_ocean ${DN_OCEANMOUTH} coast_oceangrid.nc  ;;
           * ) echo 'IMODE=[1,4] not defined for IOCEAN = ' $IOCEAN ' --> ABBRUCH!!!' ; exit ;;
         esac  ;;
     3 ) DN_OMOUTH_PRESET=$DN_OCEANMOUTH  ;;
esac
# *** At this point, coast_oceangrid.nc comprises prescribed (1,4: from $DN_OCEANMOUTH) 
# *** coastal ocean points. For IMODE=[2,3], these will be calculated below.
# *** If IMODE=[3,4], this becomes the secondary mask and the primary mask is taken either 
# *** from $DN_OCEANMOUTH or prescribed, respectively.
case $IMODE in 
 2 | 3 ) ${F90} ${DSRC}/mo_grid.f90 ${DSRC}/generate_mouth.f90 -o gen_mouth.exe $NC_INCLUDE $NC_LIB
         ln -sf ${DN_OMASK} ocean_mask.nc
#
#        *** Run program that generates potential mouths as coastal ocean points from ocean model sea mask 
#        *** that neighbor land points or that neighbor missing values at the NSEW side.
#        ***  --> Output: coast_oceangrid.nc
         ./gen_mouth.exe
         if [ -s coast_oceangrid.nc ] ; then 
           echo "potential river mouths on ocean grid created from ocean mask !"
           cp -p coast_oceangrid.nc coast_ocean_${OM}.nc
         else
           echo "coast_oceangrid.nc does not exists = Failure in gen_mouth.exe --> TERMINATION of SCRIPT !"
           exit
         fi 
         if [ $IMODE -eq 3 ] ; then 
           cdo setvar,nemo_mask_preset -selvar,nmhd_msk ${DN_OMOUTH_PRESET} coast_ocean_${OM}_preset.nc
           ncks -a -v nemo_mask_preset -A coast_ocean_${OM}_preset.nc coast_ocean_${OM}.nc
           ln -sf coast_ocean_${OM}.nc coast_oceangrid.nc
         fi
         rm gen_mouth.exe  ;;   
  4 )    cp -p coast_oceangrid.nc coast_ocean_${OM}.nc
         cdo setvar,nemo_mask_preset ${DN_OMOUTH_PRESET} coast_ocean_${OM}_preset.nc
         ncks -a -v nemo_mask_preset -A coast_ocean_${OM}_preset.nc coast_ocean_${OM}.nc
         ln -sf coast_ocean_${OM}.nc coast_oceangrid.nc  ;;
esac
#
# *** Allocate HD mouths
${F90} ${DSRC}/mo_grid.f90 ${DSRC}/mo_time.f90 ${DSRC}/mo_flow_inout.f90 ${DSRC}/mo_interpol.f90 ${DSRC}/mo_convert.f90 ${DSRC}/convert_discharge.f90 -o convert.exe $NC_INCLUDE $NC_LIB
#
case $ISRC in
  [123] ) cdo setctomiss,0. -eqc,0. -selvar,FDIR $DNHDPARA hdmouth_mask.nc
          cdo setvar,FMOUTH hdmouth_mask.nc rivmouth_source.nc
          rm hdmouth_mask.nc ;;
    [4] ) cdo eqc,0. -selvar,FDIR $DNHDPARA mrm_mouth.nc
          cdo eqc,5. -selvar,FDIR $DNHDPARA mrm_sinks.nc
          cdo add mrm_mouth.nc mrm_sinks.nc mrmmouth_mask.nc
          cdo setvar,FMOUTH -setctomiss,0. -setmissval,-9999. mrmmouth_mask.nc rivmouth_source.nc
          rm mrm_mouth.nc mrm_sinks.nc ;;
    5   ) cdo setvar,FMOUTH -gtc,-90. -selvar,lon $DNHDPARA t3.nc
          cdo selvar,lon $DNHDPARA t1.nc
          cdo ifthenelse t3.nc t1.nc -setctomiss,0. t3.nc t4.nc
          cdo selvar,lat $DNHDPARA t2.nc
          ncrename -h -O -d lat_dim,lon_dim t2.nc
          cdo ifthenelse t3.nc t2.nc -setctomiss,0. t3.nc t5.nc
          cdo -O merge t[345].nc rivmouth_source.nc
          rm t[12345].nc ;;
   [67] ) cdo setvar,FMOUTH $DNHDPARA rivmouth_source.nc  ;;
    *   ) echo 'Treatment of mouth mask needs to be specified -> STOP' ; exit ;;
esac
#
# *** Run the Program to create data that will be used for the remapping of mouth points
# ***  --> Creates files hd_to_ocean_mouth.nc, hdmouth_on_oceangrid.nc
./convert.exe $IMODE $ISRC $IOCEAN 
set +x
DNOUT=`ls hdcouple*to*imode?.nc`

if [ $IOCEAN -eq 1 ] || [ $IOCEAN -eq 8 ] ; then 
  cdo gtc,0.5 -setmisstoc,0. -selvar,FMOU_HD_TO_NEMO $DNOUT hd_used_mouths_for_${OM}${OVS}.nc
  cdo lec,0 hd_used_mouths_for_${OM}${OVS}.nc ${DIN}/hd_${HDVS}_used_mouths_${OM}${OVS}_inverted_for_oasis.nc
  rm hd_used_mouths_for_${OM}${OVS}.nc
fi
#
if [ $ISRC -eq 6 ] ; then 
  mv $DNOUT ${DIN}/${OM}/${EXP}_$DNOUT
else
  mv $DNOUT ${DIN}/${OM}/$DNOUT
fi
echo 'Coupling file: ' ${DIN}/${OM}/$DNOUT ' generated'
rm coast_oceangrid.nc
rm convert.exe
#
# *** 360-day calendar?
if [ $ICAL -eq 1 ] ; then
  YEAR=$YBEG
  while [ $YEAR -le $YEND ] ; do
    DNOUTFLOW=${DATA_OUT}/${CFLOW}_${CDIS}_on_${OM}${OVS}_${YEAR}.nc
    cdo --no_history -f nc4 -z zip_2 settaxis,${YEAR}-01-01,12:00:00,1d -setcalendar,360_day $DNOUTFLOW tmp.nc
    mv tmp.nc $DNOUTFLOW
    echo " Calendar corrected for " $DNOUTFLOW 
    YEAR=`expr $YEAR + 1`
  done
fi
#
# Inflow conversion?
if [ $ICONV -gt 0 ] ; then
  if [ $IOCEAN -eq 5 ] || [ $IOCEAN -eq 7 ]; then
    if [ $ICONV -gt 1 ] ; then
      echo 'For ICON model, ICON grid setting must be written to converted discharge files.'
      echo "  e.g,  cdo -z zip_2 setgrid,$DN_OCEANMOUTH"
    fi
  fi
fi
#
#

