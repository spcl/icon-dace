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

require 'pp'
require "numru/netcdf"
require 'numru/netcdf_miss'
include NumRu

# 
# This script splits the glocal cell grid into northern and southern hemisphere.
# @author: ralf mueller , ralf.mueller@dkrz.de
#
# Requirements: ruby, ruby modules: ruby-netcdf
def splitHemisphere(iFilename,varname,lon,lat)
  iFile     = NetCDF.open(iFilename)
  lats      = iFile.var(lat).get
  iceValues = iFile.var(varname).get_with_miss

  # compute location indices and corresponding values
  nhIndeces = (lats>0.0).where
  shIndeces = (lats<0.0).where
  iceValuesNH = iceValues[nhIndeces,true,1..-1]
  iceValuesSH = iceValues[shIndeces,true,1..-1]


  # create output
  nhFile,shFile = "nh_#{iFilename}","sh_#{iFilename}"
  [nhFile,shFile].each_with_index {|file,i|
    puts "Creating '#{file}' ...."
    indeces = [nhIndeces,shIndeces][i]
    f = NetCDF.create(file)
    iFile.each_dim {|dim| 
      next if ['clon','clat','ncells'].include?(dim.name) or 
      f.def_dim(dim.name,dim.length)
    }
    ["clon","clat","ncells"].each {|hdim| 
      f.def_dim(hdim,indeces.size)
    }

    iFile.each_var{|var|
      a      = {var.name => var.dim_names}
      newvar = f.def_var( var.name, var.ntype, var.dim_names )
      var.each_att{|att| newvar.put_att( att.name, att.get )}
    }
    f.enddef
    iFile.each_var{|var| 
      #puts var.name
      case var.name 
      when varname
	f.var(var.name).put(var.get[indeces,true,true])
      when 'p_ice_concSum'
	f.var(var.name).put(var.get[indeces,true])
      when 'cell_area'
	f.var(var.name).put(var.get[indeces])
      when lon,lat
	f.var(var.name).put(var.get[indeces])
      when lon+'_vertices',lat+'_vertices'
	f.var(var.name).put(var.get[true,indeces])
      else
	f.var(var.name).put(var.get)
      end
    }
    f.close
  }

  [nhFile,shFile]
end


if $0 == __FILE__; then

  iFilename = ARGV[0]

  varname    = ARGV[1]
  lon,lat   = 'clon','clat'

  splitHemisphere(iFilename,varname,lon,lat)
end

