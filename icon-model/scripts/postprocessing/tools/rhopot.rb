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
require 'socket'
require 'iconPlot'

#==============================================================================
# some helper methods =========================================================
def dbg(msg); pp msg unless ENV['DEBUG'].nil? ; end
def myPlotter
  plotFile = 'thingol' == Socket.gethostname \
           ? '/home/ram/src/git/icon/scripts/postprocessing/tools/icon_plot.ncl' \
           : '/pool/data/ICON/tools/icon_plot.ncl'
  plotter  = 'thingol' == Socket.gethostname \
           ? IconPlot.new(ENV['HOME']+'/local/bin/nclsh', plotFile, File.dirname(plotFile),'png','qiv',true,true) \
           : IconPlot.new("/home/zmaw/m300064/local/bin/nclsh", plotFile, File.dirname(plotFile), 'png','display',true,true)
  [plotter,plotFile]
end
#------------------------------------------------------------------------------
def initFilename(experiment)
  ENV['INIT'].nil? ? "initial_#{experiment}.nc" : ENV['INIT']
end
#------------------------------------------------------------------------------
def secPlot(ofile,experiment,secPlots,lock,temp='t_acc',sal='s_acc',rho='rhopoto',plotDir=".")
  plotDir << '/' unless '/' == plotDir[-1]

  title = (true) ? experiment : '"ICON Ocean, Mimetic-Miura, L40"'

  plotter, plotFile = myPlotter

  im = plotter.scalarPlot(ofile,plotDir+'T_'+     File.basename(ofile,'.nc'),temp,
                          :tStrg => title, :bStrg => '" "',
                          :hov => true,
                          :minVar => -3.0,:maxVar => 3.0,:withLines => false,:lStrg => 'T',
                          :numLevs => 20,:rStrg => 'Temperature', :colormap => "BlWhRe")
  lock.synchronize {(secPlots[experiment] ||= []) << im }
  im =  plotter.scalarPlot(ofile,plotDir +'S_'+     File.basename(ofile,'.nc'),sal,
                           :tStrg => title, :bStrg => '" "',
                           :hov => true,
                           :minVar => -0.2,:maxVar => 0.2,:withLines => false,:lStrg => 'S',
                           :numLevs => 16,:rStrg => 'Salinity', :colormap => "BlWhRe")
  lock.synchronize {(secPlots[experiment] ||= []) << im }
  im = plotter.scalarPlot(ofile,plotDir+'rhopot_'+File.basename(ofile,'.nc'),rho,
                          :tStrg => title, :bStrg => '"  "',
                          :hov => true,
                          :minVar => -0.6,:maxVar => 0.6,:withLines => false,:lStrg => 'rhopot',
                          :numLevs => 24,:rStrg => 'Pot.Density', :colormap => "BlWhRe")
  lock.synchronize {(secPlots[experiment] ||= []) << im }
end
#------------------------------------------------------------------------------
# plot the first 3 levels if possbile
def horizPlot(ifile,timesteps,varnames,experiment,plots,lock,plotDir=".")
  plotter, plotFile = myPlotter

  # preselect the timestep
  timestepData = Cdo.seltimestep(timesteps.join(','), input: ifile,output: [timesteps.join('-'),File.basename(ifile)].join('_'))
  # precompute levels
  levelsOfVars = {}
  varnames.each {|varname| (1..3).each {|levidx| (levelsOfVars[varname] ||= [] ) << Cdo.showlevel(input: " -sellevidx,#{levidx} -selname,#{varname} #{timestepData}").first }}

  timesteps.each {|ts| 
    levelsOfVars.each {|varname,levels|
      levels.each_with_index {|level,levidx|
        im = plotter.scalarPlot(timestepData,"#{varname}_#{level}m_TS#{ts}_#{File.basename(ifile,'.nc')}",varname,
                          :tStrg => experiment, :bStrg => '" "',:maskName => "wet_c",:maskFile => ifile,
                          :levIndex => levidx+1, :timeStep => ts-1, :tStrg => "#{varname} ")
        lock.synchronize {(plots[experiment] ||= []) << im }
      }
    }
  }
end
#------------------------------------------------------------------------------
def cropSecPlots(secPlots,plotDir='.')
  plotDir << '/' unless '/' == plotDir[-1]
  q = JobQueue.new
  secPlots.each {|exp,files|
    q.push {
      cropfiles = []
      files.each {|sp|
        cropfile = plotDir+"crop_#{File.basename(sp)}"
        cmd      = "convert -resize 60%  -crop 650x560+50+50 #{sp} #{cropfile}"
        puts cmd
        system(cmd)
        cropfiles << cropfile
      }
      images = cropfiles[-3,3]
      system("convert -append #{(images.grep(/crop_T/) + images.grep(/crop_S/) + images.grep(/crop_rho/)).join(' ')} #{plotDir}#{exp}.png")
    }
  }
  q.run
  #system("convert +append #{(cropfiles.grep(/crop_T/) + cropfiles.grep(/crop_S/) + cropfiles.grep(/crop_rho/)).join(' ')} exp.png")
  #system("display #{cropfiles.join(' ')}") #if 'thingol' == Socket.gethostname
end
#------------------------------------------------------------------------------
def cropMapPlots(plots,plotDir='.')
  plotDir << '/' unless '/' == plotDir[-1]
  q = JobQueue.new
  plots.each {|exp,files|
    q.push {
      cropfiles = []
      files.each {|sp|
        cropfile = plotDir+"crop_#{File.basename(sp)}"
        cmd      = "convert -crop 820x450+80+140 -resize 80% #{sp} #{cropfile}"
        puts cmd
        system(cmd)
        cropfiles << cropfile
      }
#      images = cropfiles[-3,3]
#      system("convert +append #{(images.grep(/crop_T/) + images.grep(/crop_S/) + images.grep(/crop_rho/)).join(' ')} #{plotDir}#{exp}.png")
    }
  }
  q.run
  #system("convert +append #{(cropfiles.grep(/crop_T/) + cropfiles.grep(/crop_S/) + cropfiles.grep(/crop_rho/)).join(' ')} exp.png")
  #system("display #{cropfiles.join(' ')}") #if 'thingol' == Socket.gethostname
end
#------------------------------------------------------------------------------
def computeRhopot(ifile,ofile,tempName,salName)
  changeNames = ''
  changeNames << '-chname' if 't' != tempName.downcase or 's' != salName.downcase
  changeNames << ',' << [tempName,'t'].join(',') unless tempName.downcase == 't'
  changeNames << ',' << [salName,'s'].join(',')  unless salName.downcase  == 's'

  ifileName = ifile
  ifile     = "#{changeNames} #{ifile}" unless changeNames.empty?
  if Cdo.version < "1.6.0"
    Cdo.rhopot(0,:input => ifile,:output => ofile)
  else
    # remove all codes so that adisit can use names
    #  use ncatted (fast in-place edit) + mv (new filename to mark, that code
    #  attribute is removed)
    if Cdo.showcode(:input => " -seltimestep,1 #{ifileName}")[0].split.map(&:to_i).reduce(0,:+) >= 0
      cmd = "ncatted -O -a code,,d,, #{ifileName}"
      dbg(cmd)
      puts IO.popen(cmd).read
    end
    Cdo.rhopot(0,:input => "-adisit #{ifile}",:output => ofile)
  end
end
#==============================================================================
######### M A I N   S C R I P T ###############################################
#==============================================================================
# USAGE:
#   MASK=mask.nc GRID=grid.nc <INIT=initFile.nc> rhopot.rb <list-of-icon-ocean-result-files>
#------------------------------------------------------------------------------
# check input
if ARGV[0].nil?
  warn "no input files given"
  exit(-1)
end

files            = ( ARGV.size > 1 ) ? ARGV : Dir.glob(ARGV[0])
maskFile         = ENV["MASK"].nil? ? "mask.nc" : ENV["MASK"]
maskVar          = ENV["MASKVAR"].nil? ? "wet_c" : ENV["MASKVAR"]
gridFile         = ENV["GRID"].nil? ? "grid.nc" : ENV["GRID"]
expName          = ENV["EXP"]
initFile         = ENV['INIT']
tempName,salName = 't_acc','s_acc'
# check files
if files.empty?
  warn "no files given"
  exit 1
end
files.each {|file|
  warn "Cannot read file '#{file}'" unless File.exist?(file)
}
unless File.exist?(maskFile) 
  warn "Cannot open maskfile '#{maskFile}' - try to compute it"
  Cdo.selname(maskVar,input: " -seltimestep,1 "+files[0], output: maskFile)
  exit -1 unless File.exist?(maskFile)
end
pp files unless ENV['DEBUG'].nil?
#------------------------------------------------------------------------------
#q         = JobQueue.new([JobQueue.maxnumber_of_processors,16].min)
q         = JobQueue.new
lock      = Mutex.new

Cdp.setCDO
Cdp.setDebug
Cdo.forceOutput = !ENV['FORCE'].nil?
Cdo.debug       = true
diff2init       = true
def plot?
  not /(thingol|thunder)/.match(Socket.gethostname).nil?
end
#------------------------------------------------------------------------------
# compute the experiments from the data directories and link the corresponding
# files
experimentFiles, experimentAnalyzedData = Cdp.splitFilesIntoExperiments(files,expName)
pp experimentFiles if Cdo.debug

#------------------------------------------------------------------------------
# create a mask file if it is not given via ENV['MASK']

# process the files
#   start with selectiong the initial values from the first timestep
experimentFiles.each {|experiment, files|
  pp files unless ENV['DEBUG'].nil?
  q.push {
    initFile    = initFilename(experiment)
    # create a separate File with the initial values
    if not File.exist?(initFile) or not Cdo.showname(:input => initFile).flatten.first.split(' ').include?("rhopot")
      initTS     = Cdo.selname([tempName,salName].join(','),
                               :input => "-seltimestep,1 #{files[0]}",
                               :options => '-r -f nc',
                               :output => "initTS_#{experiment}")
      initRhopot = computeRhopot(initTS,"initRhopot_#{experiment}",tempName,salName)
      Cdo.merge(:input => [initTS,initRhopot].join(' '),:output => initFile)
    end
  }
}
q.run
# compute meaked weight
maskedAreaWeights = Cdp.maskedAreaWeights("cell_area",
                                          gridFile,
                                          "wet_c",
                                          maskFile,
                                          "maskedAeraWeightsFrom#{File.basename(maskFile)}")
#   process the experiments results
experimentFiles.each {|experiment, files|
  files.each {|file|
    q.push {
      maskedYMeanFile = "masked_#{File.basename(file)}"
      fldmeanFile     = "fldmean_#{File.basename(file)}"
      rhopotFile      = "rhopot_#{File.basename(file)}"
      mergedFile      = "T-S-rhopot_#{File.basename(file)}"
      diffFile        = "T-S-rhopot_diff2init_#{File.basename(file)}"
      initFile        = initFilename(experiment)

      Cdo.div(:input  => " -selname,#{[tempName,salName].join(',')} #{file} #{maskFile}",
              :output => maskedYMeanFile)
      # compute rhopot
      computeRhopot(maskedYMeanFile,rhopotFile,tempName,salName)

      Cdo.merge(:input => [maskedYMeanFile,rhopotFile].join(' '), :output => mergedFile)
#      Cdo.sub(:input => [mergedFile,initFile].join(' '),:output => diffFile)
      Cdo.fldsum(:input => "-mul -sub #{[mergedFile,initFile].join(' ')} #{maskedAreaWeights}", :output => fldmeanFile,:options => '-r -f nc')
      lock.synchronize {experimentAnalyzedData[experiment] << fldmeanFile }
    }
  }
}
q.run

# merge all yearmean data (T,S,rhopot) into one file per experiment
q.clear
secPlots, mapPlots = {}, {}
experimentAnalyzedData.each {|experiment,files|
  tag      = diff2init ? 'diff2init' : ''
  ofile    = [experiment,'T-S-rhopot',tag].join('_') + '.nc'
  FileUtils.rm(ofile)  if File.exist?(ofile)
  Cdo.cat(:input => files.sort.join(' '), :output => ofile, :force => plot?)

  yearmean = "yearmean"
  ymfile   = [yearmean,ofile].join("_")
  FileUtils.rm(ymfile) if File.exist?(ymfile)
  Cdo.settunits('years',:input => "-yearmean #{ofile}", :output => ymfile,:force => plot?)

  q.push { secPlot(ymfile,experiment,secPlots,lock) if plot? }
  q.push { horizPlot(experimentFiles[experiment][0],[1,2,3],['t_acc','s_acc'],experiment,mapPlots,lock) if plot? }
}
q.run
cropSecPlots(secPlots) if plot?
cropMapPlots(mapPlots) if plot?

