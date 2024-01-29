# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# This is a prototype of a statistics module. Target language is fortran 95
require 'narray'

class DataStatistics

  # There are multiple cases of initialization
  # * use mode (Fixnum) to define the shape of the array, where the data is
  #   stored - values, which should be added must have the same shape. If mode
  #   is missing, then 
  # * the target argument has to be used. It defines an external NArray, where
  #   the data is accumulated. In addition to these two options
  # * the source paramters can be used to define an external NArray, from which
  # the data is read during the addtition to the stats object
  def initialize(target: [],source: [])
    @accumulation = target
    @source       = source
    @max, @min, @mean, @median = [],[],[],[]
    @numberOfAccumulations = 0
  end
  attr_accessor :accumulation, :max, :min ,:mean, :median, :numberOfAccumulations, :source

  # add values (singe, list or field) to the accumulation
  def add(values=nil)
    if @source.empty? and values.nil?
      warn "Nothing to add!"
      raise ArgumentError,
            "Both @source attribute and values parameters are not given!"
    end
    if (not @source.empty? and not values.nil?) then
      warn "To much to add!"
      raise ArgumentError,
            "Both @source attribute and values parameters are given!"
    end

    toBeAdded = values.nil? ? @source[0] : values

    @accumulation << toBeAdded
    nextValue = @accumulation.inject(0,&:+)
    @accumulation.clear
    @accumulation << nextValue

    @max          << toBeAdded.max  ; @min     << toBeAdded.min
    @mean         << toBeAdded.mean ; @median  << toBeAdded.median

    @numberOfAccumulations += 1
  end
  def variance
  end
  def meanHorizontal
  end
end
class DataStatisticsCollection
  def initialize
    @statistics =[]
  end
  def add(statistic)
    @statistics << statistic
  end
  def update
    @statistics.each {|statistic|
      statistic.add
    }
  end
end

