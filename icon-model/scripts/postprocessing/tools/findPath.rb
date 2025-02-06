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
require 'facets/array/nonuniq'
require 'numru/netcdf'
include NumRu

module IconPathsAlongCells
  # R2B04: list of cellsof which their common edge define the path
  CellPairsLists= {
    :gibraltar      => [[4543,4485]],
    :denmarkStrait  => [[4690,4695]],
    :drakePassage   => [
      [18481,18495],
      [18483,18493],
      [18642,18649],
      [18651,18643],
      [19186,19191]],

    :indonesianThroughflow_00 => [
      [6864,6866],
      [6868,6869],
      [6872,6875],
      [12064,12067],
      [12072,12073],
      [12152,12154]],

    :indonesianThroughflow_01 => [
      [6855,6857],
      [6858,6869],
      [6878,6879],
      [11265,11269],
      [11461,11462],
      [11470,11471],
      [11509,11510],
      [11518,11519],
      [12106,12107],
      [12109,12110],
      [12133,12134],
      [12142,12143]],

    :scotlandIceland => [
      [4682,4702],
      [4238,4270],
      [4268,4271],
      [4256,4259]],

    :mozambique => [
      [5809,5813],
      [5810,5815],
      [11155,11163]],

    :framStrait => [ #spitsbergen - greenland
      [690,697],
      [4755,4765]],

    :beringStrait => [
      [2642,2647]],

    :barentsOpening => [ #spitsgerben - norway
      [562,567],
      [560,561],
      [572,575]],

    :agulhas => [ #eastOfsouth africa - coutner current to acc
      [20384,20385]],
    :agulhas_long => [ #longer
      [20384,20385],
      [20388,20391]],
    :agulhas_longer => [ #longer
      [20384,20385],
      [20388,20391],
      [15788,15790]],
  }

  def IconPathsAlongCells.getEdgesAndVerts(ifile)
    ifileHandle  = NetCDF.open(ifile,"r")
    cellEdges    = ifileHandle.var('edge_of_cell').get
    cellVertices = ifileHandle.var('vertex_of_cell').get
    edgeVertices = ifileHandle.var('edge_vertices').get

    paths = {}

    puts '#==========================================================================='
    CellPairsLists.each {|location,cellPairs|

      verts, edges = [], []
      cellPairs.each {|cellpair|
        pair_verts, pair_edges = [], []
        cellpair.each {|cell|
          pair_edges << cellEdges[cell,0..-1].to_a
          pair_verts << cellVertices[cell,0..-1].to_a
        }
       # pp pair_edges
       # puts
        pair_verts.flatten.nonuniq.each {|v| verts << v}
        pair_edges.flatten.nonuniq.each {|e| edges << e}
      }

      paths[location] = {verts: verts.uniq, edges: edges.uniq, cells: cellPairs.flatten}

      puts '============================================================================='
      puts location
      puts ["cells:        ", cellPairs.join(' ')].join
      puts ["common verts: ", paths[location][:verts].join(',')].join
      puts ["common edges: ", paths[location][:edges].join(',')].join
    }
    puts '============================================================================='

    return paths
  end
end

if $0 == __FILE__ then
  # = MAIN =======================================================================
  ifile = ARGV[0]

  if ifile.nil? or not File.exist?(ifile) 
    warn "Cound not read input file '#{ifile}'!"
    warn "Usage:\n\t./findPath.rb <icon-grid-file>"
    exit(1)
  end

  a = IconPathsAlongCells.getEdgesAndVerts(ifile)
end
