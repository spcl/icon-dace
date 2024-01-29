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
require 'jobqueue'

# USAGE =======================================================================
# ./omip2icon.rb <ifile0,ifile1,...,ifileN> <resolution> <targetGrid> <targetWeight> <timeIntegration> <nProcs>
#
# <ifile0..N> should be a comma separated list of the original omip data files
# <resolution> can be R2B02, R2B04 or any thring you have a valid grid file
#              for. This is only used to 
#   * tag the outout file
#   * look for the right defaults for grid and weights
# <targetGrid> is a singe ICON grid (cell-only)
# <targetWeight> the weights for the above grid, this optional, weights will be
#                precomputed, if they are not given
# <timeIntegration> is an optional operator for temporal averaging, if left
#                   out, output is daily like in the original
# <nProcs> is the number of parallel processing threads (default:8)
#============================================================================== 

# CONFIG PRESETS ==============================================================
#   target horizontal icon resolution; will influece the default values for
#   grid and weight file
RESOLUTION            = 'R2B02'
#   target averaging method; omip has daily values, other possible values are
#   monmean, seasmean or what ever you think is usefull
TIMEINTERVAL_OPERATOR = ''
# COMMANDLINE OPTIONS ========================================================= 
#   list of omip intput files; in the default case additional land-sea-mask
#   files are ignored a comma separated list can be given as first argument
#   instead
iFiles                = ARGV[0].nil? ? Dir.glob("./orig/*nc").delete_if {|f| 
                                              f =~ /land_sea_mask_larger_continents/ or 
                                              f =~ /land_sea_mask.ECMWF/
                                                            }          : ARGV[0].split(',')
#
#   spacial resolution
resolution            = ARGV[1].nil? ? RESOLUTION                      : ARGV[1]
#
#   target grid file
#   use 'cdo selname,ifs2icon_cell_grid <icon-grid file> <ofile>' for creating another on
targetGridFile        = ARGV[2].nil? ? "cell_grid-#{resolution}.nc"    : ARGV[2]
#
#   weights for the target grid; this is optional and will be computed
#   automatically is it's not there
targetGridweightsFile = ARGV[3].nil? ? "cell_weight-#{resolution}.nc"  : ARGV[3]
#
#   operator for the temporal resolution of the output
timeIntervalOperator  = ARGV[4].nil? ? TIMEINTERVAL_OPERATOR           : ARGV[4]
#
#   number of parallel processes to run; depends on the hardware
nWorkers              = ARGV[5].nil? ? 8                               : ARGV[5]
#============================================================================== 

# lets work in debug mode if given by the user, i.e. DEBUG=1
Cdo.debug = ENV['DEBUG'].nil? ? false : true

# FILE CHECKING ===============================================================
# targetGridFile is mandatory, whereas the weights can be computed if there are
# not given
[iFiles,targetGridFile].flatten.each {|f|
  unless File.exist?(f)
    warn "Input file '#{f}' cannot be found!"
    exit -1
  end
}
unless File.exist?(targetGridweightsFile) then
  puts "WeightFile '#{targetGridweightsFile}' is cannot be found. It will be computed ..."

  # precompute the weight
  Cdo.gencon(targetGridFile,:in => iFiles[0],:out => targetGridweightsFile,:options => "-P #{nWorkers}")
end
# get the filename with the land-sea-mask; warn if missing
#   lsm is used to limit the horizontal interpolation to the ocean area, this
#   is esp. important for wind stress, but here it;s done for all variables
lsmFile = iFiles.find {|v| v =~ /land_sea_mask\.nc/}
if lsmFile.nil?
  warn "#================================================================================" 
  warn "Land-Sea-Mask (land_sea_mask.nc) file is MISSING! Going on without respecting continents ..."
  sleep 1
  fillFromWater_Operator = ''
else
  fillFromWater_Operator = "-fillmiss -ifthen #{lsmFile}"
  # remove the lsm from the input files
  iFiles.delete(lsmFile)
end
#============================================================================== 

# PROCESSING ================================================================== 
# proccess variables in parallel with the following steps:
# * correct the time axes
# * mask out the land points with original land sea mask
# * fill the missing values
# * remap to the icon grid

# create a queue with a predifined number of workers and put the jobs for each
# file into it
jq     = JobQueue.new(nWorkers)
# catch the output files, fill the queue and run it
oFiles = []
iFiles.each {|file|
  oFile = "remapped_#{File.basename(file)}"
  oFiles << oFile
  jq.push {
    Cdo.remap(targetGridFile,targetGridweightsFile,
              :in => "#{fillFromWater_Operator} -settaxis,2001-01-01,12:00:00,1day #{file}",
              :out => oFile)
  }
}
jq.run

# Merge all the results together
timeIntervalTag = timeIntervalOperator == '' ? 'daily' : timeIntervalOperator
Cdo.merge(:in => oFiles.sort.join(" "),:out => "omip4icon-#{resolution}-#{timeIntervalTag}.nc")
#============================================================================== 
