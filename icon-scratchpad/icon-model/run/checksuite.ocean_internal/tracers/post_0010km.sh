#!/bin/bash

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

filepath=../../../experiments/oce_tracer_bubble_0010km;
cd $filepath; \
mkdir plots; \
cd plots; \

for file in ../oce_*.nc; do \
 pwd
 fname=$(basename "$file"); \
 cdo -splitsel,1 $file $fname; \
 if [ ! -f tst.nc ] 
  then
      f=$fname'000000.nc';
      cdo -f nc -sellonlatbox,-55,-20,-45,45 -remapbil,r3600x1800 -topo tst.nc; \
      cdo -P 16 gencon,tst.nc $f AtlBoxTOregularBox.nc; \
      cdo -remap,tst.nc,AtlBoxTOregularBox.nc $f remapnn_r3600x1800_$f; \
 fi
done


for file in ../oce_*+; do \
 fname=$(basename "$file")
 cdo -splitsel,1 $file $fname; \
done

for file in oce_*.nc; do \
 cdo remap,tst.nc,AtlBoxTOregularBox.nc $file remapnn_r3600x1800_$file; \
 nclsh /pool/data/ICON/tools/icon_plot.ncl -iFile=$file -varName=t_acc  -oFile=t_$file -oType=png -secLC=-40,10 -secRC=-40,30 -rStrg='-' -resolution=r3600x1800; \
 nclsh /pool/data/ICON/tools/icon_plot.ncl -iFile=$file -varName=u_acc  -oFile=u_$file -oType=png -secLC=-40,10 -secRC=-40,30 -rStrg='-' -resolution=r3600x1800; \
 nclsh /pool/data/ICON/tools/icon_plot.ncl -iFile=$file -varName=v_acc  -oFile=v_$file -oType=png -secLC=-40,10 -secRC=-40,30 -rStrg='-' -resolution=r3600x1800; \
 nclsh /pool/data/ICON/tools/icon_plot.ncl -iFile=$file -varName=w_acc  -oFile=w_$file -oType=png -secLC=-40,10 -secRC=-40,30 -rStrg='-' -resolution=r3600x1800; \
done

