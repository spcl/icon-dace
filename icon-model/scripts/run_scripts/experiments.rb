# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

require 'pp'

# Methods for generating Namelist files
module NmlLib
  INDENT = '    '
  def NmlLib.create(name,values)
    {name => values}
  end
  def NmlLib.write(nml,filename)
    File.open(filename,"a") {|f|
      nml.each {|k,v|
        f << "&#{nml.keys.first}\n"
          v.each {|kk,vv|
            f << INDENT
            f << if vv.kind_of?(Array)
              "#{[kk,vv.join(',')].join(' = ')}\n"
            else 
              "#{[kk,vv].join(' = ')}\n"
            end
          }
        f << '/' << "\n"
      }
    }
  end
end

# Basic structure of what a testcase is
class TestCase
  attr_accessor :initialization, :forcing, :physics, :algorithms
  def initialize(name,initialization,forcing,algorithms)
    @name           = name
    @initialization = initialization
    @forcing        = forcing
    @algorithms     = algorithms
    @namelists      = []
  end
end
module TestCases
  initNml = NmlLib.create( :ocean_forcing_and_init_nml,{
    :iforc_oce          =>10,
    :iforc_omip         =>10,
    :forcing_timescale          =>10,
    :iforc_stat_oce     => 3,
    :init_oce_prog      =>10,
    :itestcase_oce      =>  50,
    :idiag_oce          =>   1,
    :temperature_relaxation  => 2,
    :relaxation_param   =>   1.0,
    :irelax_3d_T        =>   0,
    :relax_3d_mon_T     =>   1.0,
    :irelax_3d_S        =>   0,
    :relax_3d_mon_S     =>   0.0
  })
  dynNml = NmlLib.create('ocean_dynamics_nml',{
  :n_zlev             =>  10,
  :dzlev_m      =>  [50.0,  100.0,  200.0,  300.0,  450.0,  600.0,  800.0, 1000.0, 1000.0, 1000.0],
  :iswm_oce           =>   0  ,
  :idisc_scheme       =>   1  ,
  :l_inverse_flip_flop=>'.FALSE.',
  :l_rigid_lid        =>'.FALSE.',
  :ab_beta            =>   0.7,
  :ab_gam             =>   0.7,
  :solver_tolerance   =>1.0E-6,
  :i_bc_veloc_lateral =>   0  ,
  :i_bc_veloc_top     =>   1  ,
  :i_bc_veloc_bot     =>   1  ,
  :basin_center_lat   =>  30.0,
  :basin_center_lon   =>   0.0,
  :basin_width_deg    =>  60.0,
  :basin_height_deg   =>  60.0,
  :coriolis_type      =>   1  ,
  :expl_vertical_velocity_diff => 1,
  :expl_vertical_tracer_diff   => 1
  })

  TC33 = TestCase.new('33',initNml,nil,dynNml)
end

# Implements the main functionality for running an experiment:
# * 
class Experiment
  @@defaulCompiler = 'gcc'
  @@defaultMpi     = 'openmpi'
  def initialize(name,testcase,resolution,verticalLevels,nProcs=8,nThreads=4)
    @name       = name
    @testcase   = testcase
    @resolution = resolution
    @nProcs     = nProcs
    @nThreads   = nThreads
  end
  def run(host,compiler=@@defaulCompiler,mpi=@@defaultMpi)
    puts "# create exp dir"
    puts "# setup and checkinitial data"
    puts "# setup and check forcing date"
    puts "# write name lists"
    pp @testcase
    puts "# load modules"
    puts "# run the model"
  end

  def setExperimentDir
  end
  def setupInputData
  end
  def setupForcingData
  end
  def checkFile(filename)
    checkHorizontalResolution
    checkVerticalLevels
  end
  def checkHorizontalResolution
  end
  def checkVerticalLevels
  end
  def createNamelists
  end
  private :setupForcingData, :setupInputData, :setExperimentDir, :checkFile, :checkHorizontalResolution, :checkVerticalLevels
end
