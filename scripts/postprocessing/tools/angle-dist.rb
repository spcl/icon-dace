#!/usr/bin/env ruby

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#
# @author: ralf mueller, ralf.mueller@dkrz.de
#
# Plot angle distribution of an icon grid file
#
# Usage:
# ./angle-dist.rb <input-File> [<numberOfBins>]
#
# Requirements: ruby,ruby-netcdf,ruby-gsl

require 'numru/netcdf'
include NumRu
require 'gsl'

d2r       = 45.0/Math.atan(1.0)
maxDIV    = 1000
plotSetup = lambda {|title| "-L 'Angle distribution, #{title}' --bitmap-size 800x600 --page-size a4,xsize=20cm,ysize=15cm -f 0.03 -g 3"}
plot      = lambda {|hist,title,flag|
  ['ps','svg','png'].each {|fmt|
    hist.graph("-C -T #{fmt}  #{title} > #{flag}-angle-dist.#{fmt}")
  }
}
preview   = false

iFile = ARGV[0].nil? ? "iconR2B04-grid_spr0.90.nc" : ARGV[0]
nBINS = ARGV[1].nil? ? 300                         : ARGV[1].to_i

g4    = NetCDF.open(iFile)
x     = g4.var('meridional_normal_primal_edge').get
y     = g4.var('zonal_normal_primal_edge').get
yDIVx = y/x

# clean up to small/large values: set them to 0
yDIVx = (yDIVx*(yDIVx.abs < maxDIV))

# compute degrees
angles = yDIVx.collect {|v| Math.atan(v)*d2r}

# create histogram
histogram = GSL::Histogram.alloc(nBINS,-100.0,100.0)
histogram.increment(angles.to_gv)
histogram.graph("-C -T X -g 3") if preview
plot.call(histogram,plotSetup['ICON grid R2B04'],'icon-grid')

#==========================================================
# create plot for a pure lonlat grid with the same number of gridpoints
gridSize     = x.size
anglesLONLAT = GSL::Vector.alloc(gridSize)
(0...gridSize/2).each {|i| anglesLONLAT[i] = 0.0}
(gridSize/2...gridSize).each {|i| anglesLONLAT[i] = 90.0}

histogramLONLAT = GSL::Histogram.alloc(nBINS,-100.0,100.0)
histogramLONLAT.increment(anglesLONLAT)
histogramLONLAT.graph("-C -T X -g 3") if preview
plot.call(histogramLONLAT,plotSetup['Regular grid'],'regular-grid')
