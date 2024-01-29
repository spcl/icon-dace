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

require 'facets/array/nonuniq'
require 'numru/netcdf'
include NumRu

# = MAIN =======================================================================
ifile = ARGV[0]
if ifile.nil? or not File.exist?(ifile) 
  warn "Cound not read input file '#{ifile}'!"
  warn "Usage:\n\t./findEdges.rb <icon-grid-file>"
  exit(1)
end

ifileHandle  = NetCDF.open(ifile,"r")
cellEdges    = ifileHandle.var('edge_of_cell').get
cellVertices = ifileHandle.var('vertex_of_cell').get
edgeVertices = ifileHandle.var('edge_vertices').get

# R2B04: list of cellsof which their common edge define the path
cellPairsLists= {
  :gibraltar      => [[4543,4485]],
  :denmarkStrait  => [[4690,4695]],
  :drakePassage   => [
    [18481,18495],
    [18483,18493],
    [18642,18649],
    [18651,18643],
    [19186,19191]],
}

puts '#==========================================================================='
cellPairsLists.each {|location,cellPairs|
  puts '============================================================================='
  puts location
  puts ["cells: ", cellPairs.join(' ')].join

  verts = []
  cellPairs.each {|cellpair|
    cellpair.each {|cell| verts << cellVertices[cell,0..-1].to_a }
  }
  puts ["common verts: ", verts.flatten.nonuniq.join(',')].join

  # get the corresponding edges
  edges = []
  cellPairs .each {|cellpair|
    cellpair.each {|cell| edges << cellEdges[cell,0..-1].to_a }
  }
  puts ["common edges: ", edges.flatten.nonuniq.join(',')].join
}
  puts '============================================================================='
