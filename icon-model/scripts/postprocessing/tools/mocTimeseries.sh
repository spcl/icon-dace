#!/usr/bin/env bash

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

set -x

# file patterns:
#  - 'MOC*Z' : if more than one file, files are catted to one annual mean file 'moc*'
#  - moc*    : one annual mean file to be plotted
mocPattern='MOC*Z'
if [[ ! -z "$1" ]];then
  mocPattern=$1
fi

mocVar=var778 #AMOC

# calculate the moc timeseries in 1000 depth
# cat, yearmean, netcdf of MOC to moc.ym.nc - grid is 1x180 values, i.e. 1 deg lat-resolution
files=$(echo $mocPattern | wc -w)
if [[ $files -eq 0 ]]; then
  echo "no file[s] available"
  exit
fi
if [[ $files -ge 1 ]]; then
  cdo -f nc -r -setgrid,r1x180 -yearmean -cat "${mocPattern}"  moc.ym.nc
  mocPattern=moc.ym.nc
fi

# "moc.gnuplot.dat" - plot only:
if [[ "$mocPattern" != "moc.gnuplot.dat" ]]; then

  #  select AMOC at 1000m depth and write ascii data for gnuplot
  levs=$(cdo nlevel -selvar,$mocVar $(ls -1 $mocPattern | head -n 1))
  if [[ "$levs" = "20" ]]; then
    # old: 40-60N mean value
 #  cdo outputkey,date,value -mulc,1.e-9 -fldmean -selname,$mocVar -sellonlatbox,0,1,40,60 \
 #       -sellevel,1000 $mocPattern > moc.gnuplot.dat
    # new: at 26N, i.e. index 90+26 
    cdo outputkey,date,value -mulc,1.e-9 -selname,$mocVar -selindexbox,1,1,116,116 \
         -sellevel,1000 $mocPattern > moc.gnuplot.dat
  fi
  if [[ "$levs" = "40" ]]; then
    cdo outputkey,date,value -mulc,1.e-9 -selname,$mocVar -selindexbox,1,1,116,116 \
         -sellevel,1020 $mocPattern > moc.gnuplot.dat
  fi
# else
#   cdo outputkey,date,value -mulc,1.e-9 -fldmean -selname,$mocVar -sellonlatbox,0,1,40,60 \
#        -sellevel,1020 $mocPattern > moc.gnuplot.dat
#   cdo outputkey,date,value -mulc,1.e-9 -fldmean -selname,$mocVar -sellonlatbox,0,1,40,60 \
#        -intlevel,1000 -sellevel,900/1100 -setgrid,r1x180 $mocPattern > moc.gnuplot.dat
# fi
fi

gnuplot -p <<EOF
#set terminal postscript color
#set output 'moc-tims.ps'
set terminal png size 1200,600
set output 'moc-tims.png'

set size 0.98,0.95
set origin 0.02,0.02

set style data lines
set grid
set title font "Times-Bold,16"
set title  "Atlantic MOC (1000m depth at 26N)"
#set title  "Atlantic MOC (1000m depth, 40-60N)"
#set title  "Atlantic MOC"
#set title  "Atlantic MOC ($mocPattern)"
#set title  "N-Atlantic Meridional Overturning Circulation (MOC)"
#set title  "N-Atlantic Meridional Overturning Circulation (MOC)" font "Helvetica-Bold,17"

set timefmt x "%Y-%m-%d"
set xdata time
set format x "%Y"
#set xtics 20
#set mxtics 5
set xlabel "simulation years" font "Times-Bold,14

#set xrange ["2005":"2010"]
#set yrange [8.1:27.8]
#set ytics 4
#set mytics 4   # minor tickmarks
set ylabel "[SV]" font "Times-Bold,14

#set key left
#set key left font "Helvetica-Bold,10"
set key left font "Helvetica,11"
plot 'moc.gnuplot.dat' using 1:2 w l t "${mocPattern}" lw 3 lt -1
#plot 'moc.gnuplot.dat' using 1:2 w l t "${mocPattern}, 1000m/40N-60N" lw 3 lt -1
#plot 'moc.gnuplot.dat' using 1:2 w l t
#plot 'moc.gnuplot.dat' using 1:2 w l t 'AMOC, 1000m, 40N-60N' lw 3 lt -1
EOF

#gwenview moc-tims.png
