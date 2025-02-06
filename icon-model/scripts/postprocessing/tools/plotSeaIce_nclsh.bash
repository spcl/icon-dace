#!/usr/bin/env bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

set -ex

# Calculates sea ice thickness and snow times ice concentration for NH and SH
# Plots thickness and concentration for March/September using nclsh
# Arranges 4 plots for SH and NH together using convert
#  - SLO 29.04.2014
#  - update path etc.
#  - needs to select a year with 12 monthly timesteps
#  - now global datfile with S, T to bottom, NH and SH datafiles with upper levels only
#    calculation of nh/sh data is (normally) done in postproc_omip.bash

# Erweiterungen

# USAGE: plotSeaice_nclsh.bash datafilename {nh|sh} year [briefname expfilename]

#expfile=$1    # long filename
#selyear=$2    # year to select from output file - 12 data sets, Jan-Dec
expname=$(basename $(pwd))
datfile=$1
suffix=$2     # nh or sh
selyear=$3
brfname=${4:-$expname}
datglob=dat.$brfname.icets.glb.${selyear}y.nc
expfile=${5:-$datglob}

#echo $selyear
#brfname=${expname##.}
#brfname=$expname
#brfname=sice.wGMR
#datfile=dat.r16207.def.icets.$suffix.${selyear}y.nc
#datfile=dat.$expname.icets.$suffix.${selyear}y.nc
#datfile=dat.$brfname.icets.$suffix.${selyear}y.nc

calc=n     # calculate nh or sh ice + ts variables
ncsl=y     # plot using nclsh
conv=y     # convert plots 4 in 1
mvpl=y     # move plots of sea ice to xplot.seaice

vi=hi      # thickness only
vi=hixc    # thickness times concentration


SHOWGRID=" "
SHOWGRID=" -showGrid "

pldens=90
pldens=120

nclsh="nclsh /pool/data/ICON/tools/icon_plot.ncl -altLibDir=/pool/data/ICON/tools"

if [[ $calc = y ]]; then

  # calculate variables ice/snow time concentration, hiXconc and hsXconc
  cdo -selvar,wet_c,t_acc,s_acc,hi_acc,hs_acc,conc_acc,ice_u_acc,ice_v_acc \
      -selyear,$selyear $expfile xscr.icets.nc
  cdo chname,hi_acc,hixc_acc -mul -selvar,hi_acc xscr.icets.nc -selvar,conc_acc xscr.icets.nc xscr.hiXconc.nc
  cdo chname,hs_acc,hsxc_acc -mul -selvar,hs_acc xscr.icets.nc -selvar,conc_acc xscr.icets.nc xscr.hsXconc.nc
  rm -f xscr.hihsXconc.nc
  cdo merge xscr.hiXconc.nc xscr.hsXconc.nc xscr.hihsXconc.nc

  # calculate magnitude of velocity
  cdo chname,ice_u_acc,ice_magv_acc -sqrt \
     -add -sqr -selvar,ice_u_acc xscr.icets.nc -sqr -selvar,ice_v_acc xscr.icets.nc xscr.magv.nc

  # merge
  rm -f $datglob
  cdo merge xscr.icets.nc xscr.hihsXconc.nc xscr.magv.nc $datglob
  rm -f xscr.*nc

  # select upper 5 levels of S and T and collect to NH and SH
  datnhem=dat.$brfname.icets.nh.${selyear}y.nc
  datshem=dat.$brfname.icets.sh.${selyear}y.nc
  cdo -sellonlatbox,0,360,50,90 -sellevidx,1/5 $datglob $datnhem
  cdo -sellonlatbox,0,360,-90,-50 -sellevidx,1/5 $datglob $datshem

  if [[ "$suffix" = "nh" ]]; then datfile=$datnhem; fi
  if [[ "$suffix" = "sh" ]]; then datfile=$datshem; fi
  
fi  # calculate
  
# plot using icon_plot.ncl:
# timestep=2/8: assumes to be March/Sept - monthly data in the file

if [[ $ncsl = y ]]; then

  ts=2       # timestep
  vari=${vi}_acc
  iceplot=pl.$brfname.$suffix.${selyear}y.$vi.$ts

  TIT="MPIOM - ice thickness"
  TIT="ICON - ice thickness"

  if [[ "$suffix" = "nh" ]]; then mapt=NHps; fi
  if [[ "$suffix" = "sh" ]]; then mapt=SHps; fi

  PLOTLVth="0.01,0.1,0.2,0.3,0.5,0.7,1,1.5,2,3,4,6,9"
# PLOTLVco="0.01,0.05,0.1,0.3,0.5,0.7,0.8,0.9,0.99"
# now in %
# PLOTLVco="1,5,10,20,30,50,70,80,90,99"
  PLOTLVco="1,5,10,20,30,50,70,80,95,98"
  
  $nclsh \
     -iFile=$datfile -oFile=$iceplot -timeStep=$ts -mapType=$mapt -varName=$vari \
     -plotLevs=$PLOTLVth \
     -colormap=WhiteBlueGreenYellowRed \
     $SHOWGRID \
     -tStrg="$TIT" \
  
  convert -density $pldens $iceplot.eps $iceplot.png
  
  ts=8
  iceplot=pl.$brfname.$suffix.${selyear}y.$vi.$ts

  $nclsh \
     -iFile=$datfile -oFile=$iceplot -timeStep=$ts -mapType=$mapt -varName=$vari \
     -plotLevs=$PLOTLVth \
     -colormap=WhiteBlueGreenYellowRed \
     $SHOWGRID \
     -tStrg="$TIT" \

  convert -density $pldens $iceplot.eps $iceplot.png
  
  ts=2
  vc=conc
  varc=${vc}_acc
  iceplot=pl.$brfname.$suffix.${selyear}y.$vc.$ts
  TIT="MPIOM - ice concentration"
  TIT="ICON - ice concentration"
  
  $nclsh \
     -iFile=$datfile -oFile=$iceplot -timeStep=$ts -mapType=$mapt -varName=$varc \
     -plotLevs=$PLOTLVco -scaleFactor=100 \
     -colormap=WhiteBlueGreenYellowRed \
     $SHOWGRID \
     -tStrg="$TIT" \
  
  convert -density $pldens $iceplot.eps $iceplot.png
  
  ts=8
  iceplot=pl.$brfname.$suffix.${selyear}y.$vc.$ts
  
  $nclsh \
     -iFile=$datfile -oFile=$iceplot -timeStep=$ts -mapType=$mapt -varName=$varc \
     -plotLevs=$PLOTLVco -scaleFactor=100 \
     -colormap=WhiteBlueGreenYellowRed \
     $SHOWGRID \
     -tStrg="$TIT" \
  
  convert -density $pldens $iceplot.eps $iceplot.png

fi  # ncsl

if [[ $conv = y ]]; then

  x1dat=xscr1.png
  x2dat=xscr2.png
  splnam=pl.$brfname.$suffix.${selyear}y
  quadnhigh=quadpl.high.$brfname.$suffix.${selyear}y.$vi.png
  quadnlow=quadpl.low.$brfname.$suffix.${selyear}y.$vi.png
  
  convert +append $splnam.${vc}.2.png $splnam.${vc}.8.png $x1dat
  convert +append $splnam.${vi}.2.png $splnam.${vi}.8.png $x2dat
  convert -geometry 1200x1200 -append $x1dat $x2dat $quadnhigh
# convert -geometry  700x700  -append $x1dat $x2dat $quadnlow
  rm $x1dat $x2dat

fi  # conv

psst=n     # plot sst and convert
# SST
if [[ $psst = y ]]; then

  ts=2       # timestep
  vit=t      # temperaure plots in upper row
  # there could be other variable than thickness plotted below
  vari=${vit}_acc
  iceplot=pl.$brfname.$suffix.${selyear}y.$vit.$ts
  if [[ "$suffix" = "nh" ]]; then mapt=NHps; fi
  if [[ "$suffix" = "sh" ]]; then mapt=SHps; fi

  PLOTLVtm="-2,-1.7,-1.5,-1,0,1,2,3,5,7,10"
  
  $nclsh \
     -iFile=$datfile -oFile=$iceplot -timeStep=$ts -mapType=$mapt -varName=$vari \
     -plotLevs=$PLOTLVtm \
     -colormap=ViBlGrWhYeOrRe \
     $SHOWGRID \
  
  convert -density $pldens $iceplot.eps $iceplot.png
  
  ts=8
  iceplot=pl.$brfname.$suffix.${selyear}y.$vit.$ts

  $nclsh \
     -iFile=$datfile -oFile=$iceplot -timeStep=$ts -mapType=$mapt -varName=$vari \
     -plotLevs=$PLOTLVtm \
     -colormap=ViBlGrWhYeOrRe \
     $SHOWGRID \

  convert -density $pldens $iceplot.eps $iceplot.png

  x1dat=xscr1.png
  x2dat=xscr2.png
  splnam=pl.$brfname.$suffix.${selyear}y
  quadnhigh=quadpl.high.$brfname.$suffix.${selyear}y.t.$vi.png
  convert +append $splnam.$vit.2.png $splnam.$vit.8.png $x1dat
  convert +append $splnam.$vi.2.png $splnam.$vi.8.png $x2dat
  convert -geometry 1200x1200 -append $x1dat $x2dat $quadnhigh
  rm $x1dat $x2dat

fi  # psst

if [[ $mvpl = y ]]; then  # move to xplot.seaice
  mkdir -p xplot.seaice
  mv pl.$brfname.* xplot.seaice
fi  # mvpl
  
exit
  

  
  
  
  
  
  
  
