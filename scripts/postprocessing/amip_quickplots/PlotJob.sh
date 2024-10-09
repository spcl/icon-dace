#!/bin/sh

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

#
set -e
#####################################################
#
# Please adjust the following variables in the script
# explanation see below
#
# The following modules are loaded :
# - cdo/1.7.2-gcc48 (mistral and other maschine)
# - ncl/6.2.1-gccsys (mistral)
# - ncl/6.2.0-precompiled (other maschine)
# maybe you must change it; check it with:
# "module avail cdo" and "module avail ncl"
#

TAB=0

ATM_2d=1
ATM_3d=0

SINGLE=0
PAGE=1

TYP=ANN


EXP=mag0153_amip
YY1=1979
YY2=2008
DATDIR=/work/mh0287/users/marco/experiments/mag0153_amip_post
NAME=timmean_monmean_${EXP}


# atm_RES= ICON-grid resolution r2b4 or r2b6
atm_RES=r2b4

COMMENT='amip 160 km '      #   only for Single-Plots

WORKDIR=/mnt/lustre01/work/mh0081/m214091/

MODELDIR=~/icon-aes/
#####MODELDIR=/pool/data/ICON/post/
#####MODELDIR=/home/zmaw/m214091/contrib_quick/
#######################################################
#
# cell=filled triangles cont=filled contour 
# default= 1: model=cell, ERAinterim=cont, model-ERAinterim=cont
default=1
# otherwise default = 0
# if cell=1 all plots are filled triangles
# if cell=0 all plots are filled contour
  cell=1
if [ "$default" = "1" ];then
   cell=-1
fi

#
# ERAinterim (time frame: 1979-2008)
ERAystrt=$YY1
ERAylast=$YY2
export ERAystrt
echo ERAystrt path $ERAystrt
export ERAylast
echo ERAylast path $ERAylast
# era_RES= atmospheric grid resolution 63 
#          (used for ERAinterim-data and table program)
era_RES=63
# oce_RES= ocean grid resolution GR15 GR30 TP04 TP10
oce_RES=GR15
# Lev 47 possible
LEV=47
# ERAinterim path for trigangle data
   ERAinDir=/pool/data/ICON/post/${atm_RES}_amip/ERA-Interim/
export ERAinDir




QUELLE=${MODELDIR}/scripts/postprocessing/amip_quickplots/

export QUELLE
echo QUELLE path $QUELLE
echo $QUELLE

GrdDir=/pool/data/ICON/post/${atm_RES}_amip
GrdInfoFile=${GrdDir}/${atm_RES}_amip.nc
export GrdInfoFile
echo GrdInfoFile path $GrdInfoFile

  PLTDIR=${WORKDIR}/${EXP}_${TYP}
  export PLTDIR
  if [ ! -d ${PLTDIR} ] ; then
    mkdir ${PLTDIR}
    echo ${PLTDIR}
  fi

cd ${PLTDIR}
pwd


# Load modules 
MODULES=

    case `hostname` in
    mlogin*|mistral*)
        CDO_MODULE=cdo/1.7.2-gcc48;;
    *)  CDO_MODULE=cdo/1.7.2-gccsys;;
    esac
    MODULES="$MODULES $CDO_MODULE"

    case `hostname` in
    mlogin*|mistral*)
        NCL_MODULE=ncl/6.2.1-gccsys;;
    *)  NCL_MODULE=ncl/6.2.0-precompiled;;
    esac
    MODULES="$MODULES $NCL_MODULE"

    . $MODULESHOME/init/ksh
    module unload cdo
    module unload ncl
    module load $MODULES

which cdo 
which ncl


#
#--------------- TABLE --------------------------
if [ "$TAB" = "1" ]
then

echo $TYP $DATDIR $NAME $EXP $YY1 $YY2

  ${QUELLE}/TABLE.job $TYP $NAME $EXP $YY1 $YY2 $DATDIR 


fi
echo '####################################################'
echo  you find your table on
echo ${PLTDIR}
echo '#####################################################'
#
#--------------- QUICKPLOTS ---------------------
#
if [ "$ATM_2d" = "1" ]
then

MEANTIME="(${YY1}-${YY2})"
cat >var.txt << eof00
$NAME
$TYP
$EXP
$MEANTIME
$COMMENT
$PLTDIR
eof00
#
#

#---prepare seasonal amd timaverage from 2d-ERA-iterim
if [ ! -s "ERAin_${atm_RES}_atm_2d_${ERAystrt}_${ERAylast}_${TYP}.nc" ]
then
 ${QUELLE}/PREPAREera $ERAystrt $ERAylast $TYP $atm_RES 
fi


${QUELLE}/PREPAREatm_2d $TYP $NAME $atm_RES $DATDIR

if [ "$PAGE" = "1" ]
then
  nclsh  ${QUELLE}/atm_2d_page.ncl -default=${default} -cell=${cell}
fi
if [ "$SINGLE" = "1" ]
then

  nclsh ${QUELLE}/atm_2d_single.ncl -default=${default} -cell=${cell}

fi

set +e
rm Ubusy_*.nc   var.txt  sea.nc


set -e

echo '####################################################'
echo  you find your plots in
echo ${PLTDIR}
echo '#####################################################'

fi

############################################################
if [ "$ATM_3d" = "1" ]
then
echo ATM_dyn

MEANTIME="(${YY1}-${YY2})"
cat >var.txt << eof00
$NAME
$TYP
$EXP
$MEANTIME
$COMMENT
$PLTDIR
eof00

cp ${QUELLE}/partab .

#---prepare seasonal amd timaverage from 3d-ERA-iterim
if [ ! -s "ERAin_T63_atm_3d_zon_${ERAystrt}_${ERAylast}_${TYP}.nc" -o \
     ! -s "ERAin_T63L47_atm_3d_zon_${ERAystrt}_${ERAylast}_${TYP}.nc" ]
then
 ${QUELLE}/PREPAREera_3d $ERAystrt $ERAylast $TYP $atm_RES 
fi

#---prepare seasonal amd timaverage from 2d-ERA-iterim
if [ ! -s "ERAin_${atm_RES}_atm_2d_${ERAystrt}_${ERAylast}_${TYP}.nc" ]
then
 ${QUELLE}/PREPAREera $ERAystrt $ERAylast $TYP $atm_RES $WORKDIR
fi

#---
${QUELLE}/PREPAREatm_3d_logp $TYP $NAME $atm_RES $era_RES $ERAylast $LEV $DATDIR 

if [ "$PAGE" = "1" ]
then
  ncl   ${QUELLE}/atm_3d_logp_page.ncl
fi
if [ "$SINGLE" = "1" ]
then
  ncl   ${QUELLE}/atm_3d_logp_single.ncl
fi
rm -f Ubusy_*.nc  Uatm_dyn_pl_log Uatm_dyn_pl


#---
${QUELLE}/PREPAREatm_3d $TYP $NAME $atm_RES $era_RES $ERAylast $DATDIR

if [ "$PAGE" = "1" ]
then
  nclsh  ${QUELLE}/atm_3d_linp_page.ncl 
  nclsh  ${QUELLE}/atm_3d_map_page.ncl -default=${default} -cell=${cell}
fi
if [ "$SINGLE" = "1" ]
then
  nclsh  ${QUELLE}/atm_3d_linp_single.ncl 
  nclsh  ${QUELLE}/atm_3d_map_single.ncl -default=${default} -cell=${cell}
fi

rm -f Ubusy_*.nc   var.txt var1.txt Uatm_dyn_pl_log Uatm_dyn_pl


echo '####################################################'
echo  you find your plots in
echo ${PLTDIR}
echo '#####################################################'

fi
exit

#Please adjust the following variables in the script:
#
# EXP= experiment number, appears in the caption of the plots
#
# COMMENT= the comment appears in the subtitle of the plots
#          maximum length 20 characters (it works only for single plots)
# TYP= average to compare with ERAinterim-data(1979-1999)or (1979-2008)
#      ANN(annual), DJF(Dec-Feb), MAM(mar-may)  JJA(jul-aug), SON(sep-nov),
#      JAN ... DEC
#
# YY1= start date, appears in the caption of the plots
# YY2= end date, appears in the caption of the plots
#                                
#      
# NAME= XXX name of data files (XXX_atm_2d_ml.nc and XXX_atm_3d_ml.nc)
# WORKDIR= working directory 
# DATDIR = directory for input data XXX_atm_2d_ml.nc and XXX_atm_3d_ml.nc
# MODELDIR= model directory
#
#
#
# ATM_3d= 1 plot atmosphere data
#          0 no plot of atmospheric data
# ATM_2d= 1 plot surface data
#          0 no plot of surface data
#
#       the plot program expects the following two files:
#               XXX_atm_2d_ml.nc (surface data, containing at least:
#                           variable: 
#                                 clwvi Liquid water + ice content
#                                 clt   total cloud cover     
#                                 psl   sea level pressure    
#                                 tas   2 m temperature       
#                                 ts    surface temperature   
#                                 tauu  zonal wind stress     
#                                 prw   column water vapor    
#                                       vertical integral of cloud liquid water
#                                 pr    total precipitation   
#                                                          
#       the interpolation from model level to pressure level computes this programm automatically 
# 
#              XXX_atm_3d_ml.nc (atmosphere data, pressure levels 
#                         in hPa:  1000,925,850,775,700,600,500,400,300,250,
#                                   200,150,100,70,50,30,10
#                         containing at least:
#                             variable: ta  temperature           
#                                       ua  zonal wind            
#                                       va  meridional wind       
#                                       hus specific humidity     
#                                           velocity potential
#                                       clw cloud liquid water    
#                                       cli cloud ice             
#                                       zg  geopotential height
#                                       hur relative humidity
#                                       cl  cloud cover
#
#               XXX_atm_3d_ml.nc (atmosphere data, pressure levels (47) in hPa:
#                         100900,99500,97100,93900,90200,86100,81700,77200,
#                         72500,67900,63300,58800,54300,49900,45700,41600,
#                         37700,33900,30402,27015,23833,20867,18116,15578,
#                         13239,11066,9102,7406,5964,4752,3743,2914,2235,1685,
#                         1245,901,637,440,296,193,122,74,43,23,11,4,1
#
#                         
#                         containing at least:
#                             variable: ta  temperature
#                                       ua  zonal wind
#                                       va  meridional wind
#        

