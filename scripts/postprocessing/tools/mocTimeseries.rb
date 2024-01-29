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

require 'jobqueue'
require 'cdo'

#==============================================================================
# USAGE:
# <TITLE='title'> <PLOTFILE='image.png'> ./mocTimeseries.rb '../oce_mpiom_r11698_3600_1.0E-17_30/MOC*Z'
#
# @author: ralf.mueller , ralf.mueller@dkrz.de
#==============================================================================

mocPattern      = ARGV[0].nil?      ? 'MOC*Z'    : ARGV[0]
mocVar          = ARGV[1].nil?      ? 'var778'   : ARGV[1]
title           = ENV['TITLE'].nil? ? mocPattern : ENV['TITLE']
plotFile        = ENV['PLOTFILE']
q               = JobQueue.new
lock            = Mutex.new
Cdo.debug       = true
Cdo.forceOutput = !ENV['FORCE'].nil?

# postprocessing/expanding
mocPattern = [File.expand_path(File.dirname(mocPattern)),File.basename(mocPattern)].join(File::SEPARATOR)
# calculate the moc timeseries in 1000 depth

mocFiles = Dir.glob(mocPattern).sort

levels = Cdo.nlevel(:input => "-seltimestep,1 #{mocFiles[0]}")[0].to_i

# process each moc file in parallel and join them at the end
mocOutputFiles = []
mocFiles.each_with_index {|file,i|
  next if i >= 84
  q.push {
    levelSelection = 20 == levels ? "-sellevel,1000" : "-intlevel,1000 -sellevel,900/1100"
    ofile = "_moc_processing_#{i}"
    Cdo.fldmean(:input   => " -selname,#{mocVar} -mulc,1.e-9 -sellonlatbox,0,1,40,60 #{levelSelection} -yearmean -setgrid,r1x180 #{file}",
                :options => '-f nc -r',
                :output  => ofile)
    lock.synchronize { mocOutputFiles << [i,ofile]}
  }
}
q.run
pp mocOutputFiles
# merge together
mocComplete = Cdo.outputkey('date,value',:input => " -cat #{mocOutputFiles.sort.transpose[1].join(' ')}")
File.open("moc.dat","w") {|f| f << mocComplete.join("\n")}
IO.popen("LD_LIBRARY_PATH=/usr/lib gnuplot -p <<EOF
#{plotFile.nil? ? "" : "set terminal png;set output '#{plotFile}'"}
set timefmt x '%Y-%m-%d'
set xdata time
set grid
set format x '%Y'
set title '#{title}' font 'arial,18pt'
plot 'moc.dat' using 1:2 w l t 'AMOC, 1000m, 50N (interpol. from 40N-60N)'
EOF").read
