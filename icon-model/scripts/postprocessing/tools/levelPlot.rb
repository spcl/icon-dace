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

require 'extcsv'
require 'cdo'
require 'shellwords'
require 'unifiedPlot'

# Usage
#
# ./levelPlot.rb <inputICONfile> <variableName> <operation:optional> <plotfile:optional>
#
#
# operations can be operators file fldmin (default), fldmax, fldavg, ... (everything that has only one input and one output stream)

ifile      = ARGV[0]
varname    = ARGV[1]
operation  = ARGV[2].nil? ? 'fldmin' : ARGV[2]
plotfile   = ARGV[3]
createPlot = (not ARGV[3].nil?)

#_plotfile  = "plot_#{varname}-#{Time.new.strftime("%Y%m%d-%H%M%S")}"

#Cdo.debug = true
Cdo.setCdo('/home/ram/src/cdo/trunk/cdo/build/bin/cdo') if 'thingol' == `hostname`.chomp
#Cdo.checkCdo

# Temporal file for text output
dataFile = MyTempfile.path

# read the date
IO.popen("echo 'date|time|depth|#{varname}' > #{dataFile}")
Cdo.outputkey('date,time,level,value', 
              :input => "-#{operation} -selname,#{varname} #{ifile} >>#{dataFile}")
unit = Cdo.showunit(:input => "-selname,#{varname} #{ifile}").first

# postprocessing for correct time values
data = []
File.open(dataFile).each_with_index {|line,lineIndex|
  next if line.chomp.empty?
  _t = line.chomp.strip.gsub(/ +/,'|').split('|')
  if 0 == lineIndex then
    data << _t
    next
  end
  data << _t
}
icon = ExtCsv.new("array","plain",data.transpose)

# Create datetime column for timeseries plot
icon.datetime = []
icon.date.each_with_index{|date,i| icon.datetime << [date,icon.time[i]].join(' ') }
icon.datacolumns << "datetime"

# Plot data with automatic splitting by depth
unless icon.datacolumns.include?(varname)
  warn "Variable cannot be found!"
  exit -1
end
data = []
icon.depth.uniq.each {|level|
  obj = icon.selectBy(:depth => level)
  data << {:x => (0...obj.size).to_a, :y => obj.send(varname),:title => "#{level}m"}
}
UnifiedPlot.linePlot(data,plotConf: {:title => 'salinity fldmax for all depths over the last timesteps',
                    :label_position => 'inside left top',:xlabel => "timestep",:ylabel => 'salinity [psu]'},
:oType => 'png',:oName => 'sal_fldmax_2029' )
