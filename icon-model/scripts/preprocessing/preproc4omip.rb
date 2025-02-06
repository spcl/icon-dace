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
require 'pp'

# USAGE =======================================================================
#   MODEL=icon  GRID=../ocean_grid/iconR2B04-ocean_etopo40_planet.nc RESOLUTION=R2B04 OMIP_DIR=. FORCE=1 ./preproc4omip.rb
#   MODEL=mpiom GRID=../../mpiom/GR30L20_fx.nc  RESOLUTION=GR15 OMIP_DIR=. FORCE=1 ./preproc4omip.rb
#============================================================================== 

#==============================================================================
# CONFIGURATION
#==============================================================================
# input
OMIP_POOL_DIR = '/pool/data/MPIOM/setup/omip_365_era15'
     OMIP_DIR = ENV.has_key?('OMIP_DIR') ? ENV['OMIP_DIR'] : OMIP_POOL_DIR
  OMIP_FILES  = Dir.glob("#{OMIP_DIR}/*.nc").delete_if {|f|
    /cell_/.match(f) or /land_sea_mask.nc/.match(f)
  }.sort
#==============================================================================
# output
             TARGET = 'omip'
              MODEL = ENV.has_key?('MODEL')      ? ENV['MODEL']      : 'icon'
         RESOLUTION = ENV.has_key?('RESOLUTION') ? ENV['RESOLUTION'] : 'R2B04'
TARGET_MODEL_NOMISS = "#{TARGET}_#{MODEL}_#{RESOLUTION}-nomiss.nc"
TARGET_MODEL_OUTPUT = "#{TARGET}_#{MODEL}_#{RESOLUTION}.nc"
              FORCE = ENV.has_key?('FORCE')
#==============================================================================
# internals
                   CDO = 'cdo-dev'
               THREADS = ENV.has_key?('THREADS') ? ENV['THREADS'] : 8.to_s

                 GRID  = ENV['GRID']
 ICON_GRID_DIR_DEFAULT = '/pool/data/ICON/ocean_data/ocean_grid'
         ICON_GRID_DIR = ENV.has_key?('ICON_GRID_DIR') ? ENV['ICON_GRID_DIR'] : ICON_GRID_DIR_DEFAULT

MPIOM_GRID_DIR_DEFAULT = '/pool/data/MPIOM/input'
        MPIOM_GRID_DIR = ENV.has_key?('MPIOM_GRID_DIR') ? ENV['MPIOM_GRID_DIR']  : MPIOM_GRID_DIR_DEFAULT
# grids/weights for horizontal interpolation
remapConfig = {
  'icon' => {
    'remapOperator'=> 'genbic',
    'gridSelect'   => 'ifs2icon_cell_grid',
  },
  'mpiom' => {
    'remapOperator'=> 'genbic',
    'gridSelect'   => 'area',
  },
  'targetGrid'   => "#{MODEL}_cell_grid-#{RESOLUTION}.nc",
  'targetWeight' => "#{MODEL}_#{TARGET}-cell_weight-#{RESOLUTION}.nc",
}
#============================================================================== 
#   operator for the temporal resolution of the output
#timeIntervalOperator  = ARGV[4].nil? ? TIMEINTERVAL_OPERATOR           : ARGV[4]
#
#   number of parallel processes to run; depends on the hardware
nWorkers              = 4

targetGrid   = remapConfig['targetGrid']
targetWeight = remapConfig['targetWeight']
#============================================================================== 
# get the filename with the land-sea-mask; warn if missing
#   lsm is used to limit the horizontal interpolation to the ocean area, this
#   is esp. important for wind stress, but here it;s done for all variables
#
# TODO:
# special handling of real lsm from /pool/data/MPIOM/setup/omip_365_era15/land_sea_mask.ECMWF.nc
#   land(i+1,j+1) = MERGE(1._sp, dummy(i, j), dummy(i, j) .GT. 1.e-4_sp)
# ==>> cdo gtc,1,0e-4 real_lsm lsm
#
filterLSM = lambda {|iFiles|
  lsmFile = iFiles.find {|v| v =~ /land_sea_mask.ECMWF.nc/}
  if lsmFile.nil?
    warn "#================================================================================" 
    warn "Land-Sea-Mask (land_sea_mask.ECMWF.nc) file is MISSING! Going on without respecting continents ..."
    sleep 1
  else
    iFiles.delete(lsmFile)
  end
  return lsmFile
}
fillWaterOperator = lambda {|file,lsmFile|

  operator = ! /runoff/.match(file).nil? ? '' : "-fillmiss -ifthen #{lsmFile}"

  pp ['file:',file,'|operator:',operator].join

  return operator
}

# lets work in debug mode if given by the user, i.e. DEBUG=1
Cdo.debug = ENV['DEBUG'].nil? ? false : true
Cdo.forceOutput = false

# FILE CHECKING ===============================================================
# targetGridFile is mandatory, whereas the weights can be computed if there are
# not given
[OMIP_FILES].flatten.each {|f|
  unless File.exist?(f)
    warn "Input file '#{f}' cannot be found!"
    exit -1
  end
}
# preproc grids/weight for omip data
# compute integer lsm with MPIOMs threshold
sourceLsm = Cdo.gtc(1.0e-4,:input => filterLSM[OMIP_FILES], :output => '_lsm.nc')
if not File.exist?(targetWeight ) or FORCE then 
  if not File.exist?(targetGrid) or FORCE then
    if not File.exist?(GRID) then
      puts "GRID variablen has to be set correctly!"
      exit 1
    else
      Cdo.selname(remapConfig[MODEL]['gridSelect'],:input => GRID, :output => targetGrid)
    end
  end
  Cdo.send(remapConfig[MODEL]['remapOperator'],targetGrid,
           :options => " -P #{THREADS}",
           :input => "#{fillWaterOperator[OMIP_FILES[0],sourceLsm]} #{OMIP_FILES[0]}",
           :output => targetWeight)
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
jq     = JobQueue.new(8)
fmq    = JobQueue.new(8)
# catch the output files, fill the queue and run it
oFiles = []
OMIP_FILES.each {|file|
  oFile = "remapped_#{File.basename(file)}"
  _oFile = "fillmiss_#{File.basename(file)}"
  fmq.push {
    Cdo.fillmiss(:input => "-ifthen #{sourceLsm} #{file}", :output => _oFile)
  }
  oFiles << oFile
  puts oFile
  jq.push {
    Cdo.remap(targetGrid,targetWeight,
              :input => "#{_oFile}",
              :output => oFile)
  }
}
fmq.run
jq.run

# Merge all the results together
#timeIntervalTag = timeIntervalOperator == '' ? 'daily' : timeIntervalOperator
Cdo.settaxis('2001-01-01,12:00:00,1day',
             :input => Cdo.merge(:input => oFiles.sort.join(" ")),
             :output => "#{TARGET_MODEL_OUTPUT}")

#============================================================================== 
