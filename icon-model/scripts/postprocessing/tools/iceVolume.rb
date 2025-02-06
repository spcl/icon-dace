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

require 'cdo'
require 'cdp'
require 'fileutils'
require 'jobqueue'
require 'gsl'

#=============================================================================== 
# check input
if ARGV[0].nil?
  warn "no input files given"
  exit(-1)
else
  puts 'Please remember, that the last file in the list is expected to be an ICON grid file'
  sleep 1
end
files     = ( ARGV.size > 1 ) ? ARGV : Dir.glob(ARGV[0])
# check files
files.each {|file|
  warn "Cannot read file '#{file}'" unless File.exist?(file)
}

q         = JobQueue.new([JobQueue.maxnumber_of_processors,20].min)
p         = JobQueue.new([JobQueue.maxnumber_of_processors,20].min)
lock      = Mutex.new
#=============================================================================== 
# setup of CDO on different machines
Cdp.setCDO
Cdp.setDebug
Cdo.forceOutput = ! ENV['FORCE'].nil?
Cdo.debug       = ! ENV['DEBUG'].nil?
#=============================================================================== 
#=============================================================================== 
# helper method for plotting ice volume and extent for NH and SH
def plot(nhIceVolume,shIceVolume,nhIceExtent,shIceExtent,oType=nil,oTag=nil)
  volumeOutput, extentOutput = '',''
  volumeOutput = "-T #{oType} > iceVolume_#{oTag}.#{oType}" unless (oType.nil? and oTag.nil?)
  extentOutput = "-T #{oType} > iceExtent_#{oTag}.#{oType}" unless (oType.nil? and oTag.nil?)

  size = nhIceVolume.size
  GSL::Vector.graph([GSL::Vector.linspace(0,size-1,size),nhIceVolume],
                    [GSL::Vector.linspace(0,size-1,size),shIceVolume],
                    "-C -g 3 -X 'timesteps' -Y 'Ice Volume [km^3]' -L 'IceVolume (red:NH, green:SH)' #{volumeOutput}")
  GSL::Vector.graph([GSL::Vector.linspace(0,size-1,size),nhIceExtent],
                    [GSL::Vector.linspace(0,size-1,size),shIceExtent],
                    "-C -g 3 -X 'timesteps' -Y 'Ice Extent [km^2]' -L 'IceExtent (red:NH, green:SH)' #{extentOutput}")
end
#=============================================================================== 
# compute the experiments from the data directories and link the corresponding files
#=============================================================================== 
# MAIN
gridfile, experimentFiles, experimentAnalyzedData = Cdp.splitFilesIntoExperiments(files)
iceHeight                                         = 'p_ice_hi'
iceConcentration                                  = 'p_ice_concSum'
#   process the experiments results
experimentFiles.each {|experiment, files|
  files.each {|file|
    q.push {
      iceDiagFile = "iceDiag_#{File.basename(file)}"
      volumeFile  = "iceVolume_#{File.basename(file)}"
      extentFile  = "iceExtent_#{File.basename(file)}"
      nhFile      = "nh_#{File.basename(file)}"
      shFile      = "sh_#{File.basename(file)}"

      Cdo.setname('ice_volume',
                  :input => " -divc,1e9 -mul -mul -selname,#{iceHeight} #{file}  -selname,cell_area #{gridfile}  -selname,#{iceConcentration} #{file}",
                  :output => volumeFile)

      Cdo.setname('ice_extent',
                  :input => " -divc,1e6 -mul -selname,cell_area #{gridfile}  -gtc,0.15 -selname,#{iceConcentration} #{file}",
                  :output => extentFile)

      Cdo.merge(:input => [volumeFile,extentFile].join(' '),
                :output => iceDiagFile)

      # split the file in norther and sourtern hemisphere
      Cdo.fldsum(:input => "-sellonlatbox,-180,180,50,90 #{iceDiagFile}",  :output => nhFile) 
      Cdo.fldsum(:input => "-sellonlatbox,-180,180,-50,-90 #{iceDiagFile}",:output => shFile) 

      lock.synchronize {experimentAnalyzedData[experiment] << [nhFile,shFile] }
    }
  }
}
q.run
# merge data together for NH and SH
q.clear
iceOutputsfiles = []
experimentAnalyzedData.each {|experiment,files|
  files.transpose.each_with_index {|hemisphereFiles,i|
    q.push {
      tag   = ['nh','sh'][i]
      ofile = [experiment,'iceDiagnostics',tag].join('_') + '.nc'
      FileUtils.rm(ofile) if File.exist?(ofile)
      Cdo.cat(:input => hemisphereFiles.sort.join(' '), :output => ofile, :options => '-r')
      lock.synchronize{ iceOutputsfiles << ofile }
    }
  }
}
q.run

nhFile, shFile = iceOutputsfiles

nhIceVolume = Cdo.readArray(nhFile,'ice_volume').flatten.to_gv
shIceVolume = Cdo.readArray(shFile,'ice_volume').flatten.to_gv
nhIceExtent = Cdo.readArray(nhFile,'ice_extent').flatten.to_gv
shIceExtent = Cdo.readArray(shFile,'ice_extent').flatten.to_gv

plot(nhIceVolume,shIceVolume,nhIceExtent,shIceExtent)
#plot(nhIceVolume,shIceVolume,nhIceExtent,shIceExtent,'png','iceTest')
