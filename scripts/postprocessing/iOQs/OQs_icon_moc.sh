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

TMP_DIR=${EXPID}_${PROGRAM}_$Y1-$Y2
[ -e $TMP_DIR ] && rm -rf $TMP_DIR
mkdir -p $TMP_DIR
cd $TMP_DIR

$BINDIR/OQs_icon_input_mask.sh $GRID $LEV

CATFILE=moc.nc

CODE='777,778,779'

for YEAR in $(seq $Y1 $CHUNK $Y2) 
do

    INFILE=$(printf $DATADIR/$EXPID/MOC.%04.0f-01-01T00:00:00.000 $YEAR)

    if [ -e $INFILE ]
    then
#        $CDO cat -timmean -selcode,$CODE $INFILE $CATFILE
        $CDO cat -selcode,$CODE $INFILE $CATFILE
    else
        echo Warning ! $INFILE is missing, skipping to next file ...
    fi
done

$CDO sinfo $CATFILE # DEBUG
inbase=$(printf ${EXPID}_icon_moc_%04.0f0101_%04.0f1231 $Y1 $Y2)
infile=$inbase.nc
outfile=$inbase.pdf

$CDO -r -f nc -setgrid,r1x180 -timmean $CATFILE $infile

$CDO -r -f nc -setname,"global_moc" -setunit,"Sv" -mulc,1.025e-9 -ifthen mask_glo_zon.nc -selcode,777 -setgrid,r1x180 $infile  moc_glo.$infile
ncrename -h -O -d lev,depth -v lev,depth  moc_glo.$infile
#$BINDIR/plotsec_ncl --title="Global_ocean_meridional_overturning_streamfunction" --unit="Sv" --lstring="[Sv]" --cstring="global_moc" --rstring="${Y1}-${Y2}" --code=777 --min=-24 --max=24 --inc=3 --pal=BlueWhiteOrangeRed moc_glo.$infile
$BINDIR/plotsec_ncl --title="Global_ocean_meridional_overturning_streamfunction" --cstring="global_moc" --rstring="${Y1}-${Y2}" --code=777 --min=-24 --max=24 --inc=3 --pal=BlueWhiteOrangeRed moc_glo.$infile

$CDO -r -f nc -setname,"atlantic_moc" -setunit,"Sv" -mulc,1.025e-9 -ifthen mask_atl_arc_zon.nc -selcode,778  -setgrid,r1x180 $infile  moc_atl.$infile
ncrename -h -O -d lev,depth -v lev,depth moc_atl.$infile
#$BINDIR/plotsec_ncl --title="Atlantic_ocean_meridional_overturning_streamfunction" --unit="Sv" --lstring="[Sv]" --cstring="atlantic_moc" --rstring="${Y1}-${Y2}" --code=778 --min=-24 --max=24 --inc=3 --pal=BlueWhiteOrangeRed moc_atl.$infile
$BINDIR/plotsec_ncl --title="Atlantic_ocean_meridional_overturning_streamfunction" --cstring="atlantic_moc" --rstring="${Y1}-${Y2}" --code=778 --min=-24 --max=24 --inc=3 --pal=BlueWhiteOrangeRed moc_atl.$infile

$CDO -r -f nc -setname,"indopacific_moc" -setunit,"Sv" -mulc,1.025e-9 -ifthen mask_indopacific_zon.nc -selcode,779 -setgrid,r1x180 $infile  moc_indopacific.$infile
ncrename -h -O -d lev,depth -v lev,depth moc_indopacific.$infile
#$BINDIR/plotsec_ncl --title="IndoPacific_ocean_meridional_overturning_streamfunction" --unit="Sv" --lstring="[Sv]" --cstring="indopacific_moc" --rstring="${Y1}-${Y2}" --code=779 --min=-24 --max=24 --inc=3 --pal=BlueWhiteOrangeRed moc_indopacific.$infile
$BINDIR/plotsec_ncl --title="IndoPacific_ocean_meridional_overturning_streamfunction" --cstring="indopacific_moc" --rstring="${Y1}-${Y2}" --code=779 --min=-24 --max=24 --inc=3 --pal=BlueWhiteOrangeRed moc_indopacific.$infile

cd -

mkdir -p ${EXPID}_${Y1}-${Y2}
mv $TMP_DIR/*.$outfile ${EXPID}_${Y1}-${Y2}
mv $TMP_DIR/*.lst_ ${EXPID}_${Y1}-${Y2}

rm -rf $TMP_DIR
