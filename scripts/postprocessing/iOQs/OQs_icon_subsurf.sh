#! /usr/bin/env bash

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

set -e

set -x # DEBUG

PROGRAM=$(basename $0)

TMP_DIR=${PROGRAM}_$Y1-$Y2
[ -e $TMP_DIR ] && rm -rf $TMP_DIR
mkdir -p $TMP_DIR
cd $TMP_DIR

$BINDIR/OQs_icon_input_mask.sh $GRID $LEV

ref_SAL=$POOL/$GRID/initial_state.nc 
ref_TEM=$POOL/$GRID/initial_state.nc 

CATFILE=subsurfdata.nc

CODE='t_acc,s_acc'

for YEAR in $(seq $Y1 $CHUNK $Y2) 
do

    INFILE=$(printf $DATADIR/${EXPID}/${EXPID}_oce_def_%04.0f0101T000000Z.nc $YEAR)
    #INFILE=$(printf $DATADIR/${EXPID}/${EXPID}_%04.0f0101T000000Z.nc $YEAR)

    if [ -e $INFILE ]
    then
        $CDO cat -timmean -selvar,$CODE $INFILE $CATFILE
    else
        echo Warning ! $INFILE is missing, skipping to next file ...
    fi
done

$CDO sinfo $CATFILE # DEBUG
inbase=$(printf ${EXPID}_icon_subsurf_%04.0f0101_%04.0f1231 $Y1 $Y2)
infile=$inbase.nc
outfile=$inbase.pdf
$CDO timmean $CATFILE $infile

for basin in glo atl_arc indopacific 
do
    $CDO -f nc -zonmean -remapnn,r720x360 -ifthen mask_$basin.nc -selvar,t_acc  $infile  thetao_${basin}.$infile
    $BINDIR/plotsec_ncl --title="Sea_water_potential_temperature_${basin}" --cstring="thetao_${basin}" --lstring="[C]" --rstring="${Y1}-${Y2}" --var=t_acc --min=-2 --max=30 --inc=2 --pal=GMT_haxby  thetao_${basin}.$infile

    $CDO -f nc -zonmean -remapnn,r720x360 -ifthen mask_$basin.nc -sub -selvar,t_acc  $infile -selvar,T $ref_TEM thetao-ano_${basin}.$infile
    $BINDIR/plotsec_ncl --title="Sea_water_potential_temperature_bias_${basin}" --cstring="thetao_${basin}" --lstring="[C]" --rstring="${Y1}-${Y2}" --var=t_acc --min=-5 --max=5 --inc=0.5 --pal=BlueWhiteOrangeRed thetao-ano_${basin}.$infile

    $CDO -f nc -zonmean -remapnn,r720x360 -ifthen mask_$basin.nc -selvar,s_acc $infile  so_${basin}.$infile
    $BINDIR/plotsec_ncl --title="Sea_water_salinity_${basin}" --cstring="so_${basin}" --lstring="[1e-3]" --rstring="${Y1}-${Y2}" --var=s_acc --min=32 --max=37 --inc=0.25 --pal=GMT_haxby so_${basin}.$infile

    $CDO -f nc -zonmean -remapnn,r720x360 -ifthen mask_$basin.nc -sub -selvar,s_acc $infile  -selvar,S $ref_SAL so-ano_${basin}.$infile
    $BINDIR/plotsec_ncl --title="Sea_water_salinity_bias_${basin}" --cstring="so_${basin}" --lstring="[1e-3]" --rstring="${Y1}-${Y2}" --var=s_acc --min=-1 --max=1 --inc=0.1 --pal=BlueWhiteOrangeRed so-ano_${basin}.$infile
done

cd -

mkdir -p ${EXPID}_${Y1}-${Y2}

mv $TMP_DIR/*.$outfile ${EXPID}_${Y1}-${Y2}
mv $TMP_DIR/*.lst_ ${EXPID}_${Y1}-${Y2}


rm -rf $TMP_DIR
