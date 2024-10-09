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

CATFILE=surfdata.nc
CATFILE2=surfdata39.nc

CODE='h_acc,t_acc,s_acc,h_sqr'
CODE2='hi_acc,hs_acc,conc_acc,mld'

for YEAR in $(seq $Y1 $CHUNK $Y2) 
do

    INFILE=$(printf $DATADIR/${EXPID}/${EXPID}_oce_min_%04.0f0101T000000Z.nc $YEAR)
    #INFILE=$(printf $DATADIR/${EXPID}/${EXPID}_%04.0f0101T000000Z.nc $YEAR)

    if [ -e $INFILE ]
    then
        $CDO cat -timmean -selvar,$CODE -sellevidx,1 $INFILE $CATFILE &
        $CDO cat -selmon,3,9 -selvar,$CODE2 $INFILE $CATFILE2 &
        wait
    else
        echo Warning ! $INFILE is missing, skipping to next file ...
    fi
done

$CDO sinfo $CATFILE # DEBUG
$CDO sinfo $CATFILE2 # DEBUG
inbase=$(printf ${EXPID}_icon_surf_%04.0f0101_%04.0f1231 $Y1 $Y2)
infile=$inbase.nc
outfile=$inbase.pdf
inbase2=$(printf ${EXPID}_icon_surf39_%04.0f0101_%04.0f1231 $Y1 $Y2)
infile2=$inbase2.nc
outfile2=$inbase2.pdf
$CDO -timmean $CATFILE $infile &
$CDO -ymonmean $CATFILE2 $infile2 &
wait
rm -f $CATFILE $CATFILE2

#$BINDIR/plot2d_ncl  --cstring=barotropic_streamfunction --lstring=[Sv] --rstring="${Y1}-${Y2}" --code=27 --min=-100 --max=250 --inc=10 --scal=1.025e-9 --type=0 $infile
#mv $outfile psi.$outfile

$CDO -remapnn,${PLOTGRID_2D} -ifthen mask_glo_surf.nc -selvar,h_acc  $infile zos.$infile
$BINDIR/plot2d_ncl --title="Sea_surface_height" --lev=-1 --lstring="[m]" --cstring=zos --rstring="${Y1}-${Y2}" --var=h_acc --min=-2 --max=2 --inc=0.2 --pal=BlueWhiteOrangeRed zos.$infile

#optional
X=$($CDO showvar -selvar,h_sqr $infile | xargs echo -n )
if [ "$X" == "h_sqr" ] ; then
$CDO -remapnn,${PLOTGRID_2D} -ifthen mask_glo_surf.nc -selvar,h_sqr  $infile zos_sqr.$infile
$CDO -log10 -setname,h_var -sub -timavg -selvar,h_sqr zos_sqr.$infile -sqr -timavg -selvar,h_acc zos.$infile variance_zos.$infile
$BINDIR/plot2d_ncl --title="Log_of_sea_surface_height_variance" --lev=-1 --lstring="[log10]" --cstring=zos --rstring="${Y1}-${Y2}" --var=h_var --min=-3 --max=0 --inc=0.3 --pal=WhiteBlueGreenYellowRed --proj=SPS variance_zos.$infile
fi

#$BINDIR/plot2d_ncl  --lstring=t20d --code=109 --min=0 --max=300 --inc=30 --pal=WhiteBlueGreenYellowRed $infile
#mv $outfile t20d.$outfile

#$BINDIR/plot2d_ncl  --lstring=hfsd --code=170 --min=-150 --max=150 --inc=15 --pal=BlueWhiteOrangeRed $infile
#mv $outfile hfsd.$outfile

#$BINDIR/plot2d_ncl  --lstring=wfonocor --code=165 --min=-1e-4 --max=1e-4 --inc=1e-5 --pal=BlueWhiteOrangeRed $infile
#mv $outfile wfonocor.$outfile

$CDO -remapnn,${PLOTGRID_2D} -ifthen mask_glo_surf.nc -sellevidx,1 -selvar,t_acc  $infile tos.$infile
$BINDIR/plot2d_ncl  --title="Sea_surface_temperature" --lstring="[C]" --cstring=tos --rstring="${Y1}-${Y2}" --var=t_acc --min=-2 --max=30 --inc=2 --pal=GMT_haxby tos.$infile

#### compute anomalies to Levitus ####
$CDO -f nc -sub tos.$infile -remapnn,${PLOTGRID_2D} -selvar,T -sellevidx,1 $ref_TEM tos-ano.$infile
$BINDIR/plot2d_ncl  --title="Sea_surface_temperature_bias"  --lstring="[C]" --cstring=tos --rstring="${Y1}-${Y2}" --var=t_acc --min=-5 --max=5 --inc=0.5 --pal=BlueWhiteOrangeRed tos-ano.$infile


$CDO -remapnn,${PLOTGRID_2D} -ifthen mask_glo_surf.nc -sellevidx,1 -selvar,s_acc  $infile sos.$infile
$BINDIR/plot2d_ncl  --title="Sea_surface_salinity" --lstring="[1e-3]" --cstring=sos --rstring="${Y1}-${Y2}" --var=s_acc --min=32 --max=37 --inc=0.25 --pal=WhiteBlueGreenYellowRed sos.$infile

#### compute anomalies to Levitus ####
$CDO -f nc -sub sos.$infile -remapnn,${PLOTGRID_2D} -selvar,S -sellevidx,1 $ref_SAL sos-ano.$infile
$BINDIR/plot2d_ncl  --title="Sea_surface_salinity_bias" --lstring="[1e-3]" --cstring=sos --rstring="${Y1}-${Y2}" --var=s_acc --min=-2 --max=2 --inc=0.2 --pal=BlueWhiteOrangeRed sos-ano.$infile


#### Seasonal plots ####

$CDO -remapnn,${PLOTGRID_2D} -ifthen mask_glo_surf.nc -sellevidx,1 -mul -selvar,hi_acc  $infile2 -selvar,conc_acc  $infile2 sivol.$infile2

$CDO -O -selmon,3 sivol.$infile2 sivol-NH3.$infile2
$BINDIR/plot2d_ncl --title="Sea_ice_equivalence_thickness_NH_March" --lstring="[m]" --cstring=hi_acc --rstring="${Y1}-${Y2}" --lev=1 --var=hi_acc --min=0 --max=4 --inc=0.25 --pal=WhiteBlueGreenYellowRed --proj=NPS sivol-NH3.$infile2
$CDO -O -selmon,3 sivol.$infile2 sivol-SH3.$infile2
$BINDIR/plot2d_ncl --title="Sea_ice_equivalence_thickness_SH_March" --lstring="[m]" --cstring=hi_acc --rstring="${Y1}-${Y2}" --lev=1 --var=hi_acc --min=0 --max=4 --inc=0.25 --pal=WhiteBlueGreenYellowRed --proj=SPS sivol-SH3.$infile2

$CDO -O -selmon,9 sivol.$infile2 sivol-NH9.$infile2
$BINDIR/plot2d_ncl --title="Sea_ice_equivalence_thickness_NH_September" --lstring="[m]" --cstring=hi_acc --rstring="${Y1}-${Y2}" --lev=1 --var=hi_acc --min=0 --max=4 --inc=0.25 --pal=WhiteBlueGreenYellowRed --proj=NPS sivol-NH9.$infile2
$CDO -O -selmon,9 sivol.$infile2 sivol-SH9.$infile2
$BINDIR/plot2d_ncl --title="Sea_ice_equivalence_thickness_SH_September" --lstring="[m]" --cstring=hi_acc --rstring="${Y1}-${Y2}" --lev=1 --var=hi_acc --min=0 --max=4 --inc=0.25 --pal=WhiteBlueGreenYellowRed --proj=SPS sivol-SH9.$infile2



$CDO -remapnn,${PLOTGRID_2D} -ifthen mask_glo_surf.nc -sellevidx,1 -selvar,conc_acc  $infile2 siconc.$infile2

cdo -O -selmon,3 siconc.$infile2 siconc-NH3.$infile2
$BINDIR/plot2d_ncl --title="Sea_ice_concentration_NH_March" --lstring="[frac]" --cstring=siconc --rstring="${Y1}-${Y2}" --lev=1 --var=conc_acc --min=0 --max=1 --inc=0.05 --pal=WhiteBlueGreenYellowRed --proj=NPS siconc-NH3.$infile2
cdo -O -selmon,3 siconc.$infile2 siconc-SH3.$infile2
$BINDIR/plot2d_ncl --title="Sea_ice_concentration_SH_March" --lstring="[frac]" --cstring=siconc --rstring="${Y1}-${Y2}" --lev=1 --var=conc_acc --min=0 --max=1 --inc=0.05 --pal=WhiteBlueGreenYellowRed --proj=SPS siconc-SH3.$infile2

cdo -O -selmon,9 siconc.$infile2 siconc-NH9.$infile2
$BINDIR/plot2d_ncl --title="Sea_ice_concentration_NH_September" --lstring="[frac]" --cstring=siconc --rstring="${Y1}-${Y2}" --lev=1 --var=conc_acc --min=0 --max=1 --inc=0.05 --pal=WhiteBlueGreenYellowRed --proj=NPS siconc-NH9.$infile2
cdo -O -selmon,9 siconc.$infile2 siconc-SH9.$infile2
$BINDIR/plot2d_ncl --title="Sea_ice_concentration_SH_September" --lstring="[frac]" --cstring=siconc --rstring="${Y1}-${Y2}" --lev=1 --var=conc_acc --min=0 --max=1 --inc=0.05 --pal=WhiteBlueGreenYellowRed --proj=SPS siconc-SH9.$infile2



$CDO -remapnn,${PLOTGRID_2D} -ifthen mask_glo_surf.nc -sellevidx,1 -mul -selvar,hs_acc  $infile2 -selvar,conc_acc  $infile2 sicsno.$infile2

cdo -O -selmon,3 sicsno.$infile2 sicsno-NH3.$infile2
$BINDIR/plot2d_ncl --title="Snow_equivalence_thickness_NH_March" --lstring="[m]" --cstring=sicsno --rstring="${Y1}-${Y2}" --lev=1 --var=hs_acc --min=0 --max=1 --inc=0.05 --pal=WhiteBlueGreenYellowRed --proj=NPS sicsno-NH3.$infile2
cdo -O -selmon,3 sicsno.$infile2 sicsno-SH3.$infile2
$BINDIR/plot2d_ncl --title="Snow_equivalence_thickness_SH_March" --lstring="[m]" --cstring=sicsno --rstring="${Y1}-${Y2}" --lev=1 --var=hs_acc --min=0 --max=1 --inc=0.05 --pal=WhiteBlueGreenYellowRed --proj=SPS sicsno-SH3.$infile2

cdo -O -selmon,9 sicsno.$infile2 sicsno-NH9.$infile2
$BINDIR/plot2d_ncl --title="Snow_equivalence_thickness_NH_September" --lstring="[m]" --cstring=sicsno --rstring="${Y1}-${Y2}" --lev=1 --var=hs_acc --min=0 --max=1 --inc=0.05 --pal=WhiteBlueGreenYellowRed --proj=NPS sicsno-NH9.$infile2
cdo -O -selmon,9 sicsno.$infile2 sicsno-SH9.$infile2
$BINDIR/plot2d_ncl --title="Snow_equivalence_thickness_SH_March" --lstring="[m]" --cstring=sicsno --rstring="${Y1}-${Y2}" --lev=1 --var=hs_acc --min=0 --max=1 --inc=0.05 --pal=WhiteBlueGreenYellowRed --proj=SPS sicsno-SH9.$infile2


$CDO -remapnn,${PLOTGRID_2D} -ifthen mask_glo_surf.nc -sellevidx,1 -selvar,mld  $infile2 mld.$infile2

cdo -O -selmon,3 mld.$infile2 mld3.$infile2
$BINDIR/plot2d_ncl --title="Mixed_layer_depth_March" --lstring="[m]" --cstring=zmld --rstring="${Y1}-${Y2}" --lev=-1 --var=mld --min=0 --max=5000 --inc=250 --pal=WhiteBlueGreenYellowRed mld3.$infile2
cdo -O -selmon,9 mld.$infile2 mld9.$infile2
$BINDIR/plot2d_ncl --title="Mixed_layer_depth_September" --lstring="[m]" --cstring=zmld --rstring="${Y1}-${Y2}" --lev=-1 --var=mld --min=0 --max=5000 --inc=250 --pal=WhiteBlueGreenYellowRed mld9.$infile2

cd -

mkdir -p ${EXPID}_${Y1}-${Y2}

mv $TMP_DIR/*.$outfile ${EXPID}_${Y1}-${Y2}
mv $TMP_DIR/*.$outfile2 ${EXPID}_${Y1}-${Y2}
mv $TMP_DIR/*.lst_ ${EXPID}_${Y1}-${Y2}

rm -rf $TMP_DIR
