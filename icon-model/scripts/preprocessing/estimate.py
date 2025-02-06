#! /usr/bin/env python

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

#------------------------------------------------------------------------------
import getopt, sys, math

#"R<grid_refine_root>B<grid_refine_level>

#------------------------------------------------------------------------------

cells=0
verts=0
edges=0

#------------------------------------------------------------------------------
vertical_levels=96
#------------------------------------------------------------------------------
# set machine parameters
#nodes_list=[ 32.0 ]
#opnemp_threads_list=[1.0, 4.0, 8.0, 16.0, 32.0]
#opnemp_threads_list=[ 8.0 ]
nproma_list=[  8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 18.0 ]

#mpi_procs_list=[ 128.0, 256.0, 512.0, 541.0, 640.0, 768.0,  1082.0 ]
#mpi_nodes_list=[ 2.0, 4.0, 8.0, 16.0, 24.0, 34.0 ]
omp_threads=2.0
mpi_procs_pernode=32
# 64=compute decompositions SMT mode
# 32=compute decompositions ST mode
total_threads_per_node=64.0


def display_thread_load(message, entities, threads, block_size):
    #blocks_own_cell   = round(owned_cells/nproma)
    #blocks_own_cell   = (owned_cells/nproma)
    #blocks_total_cell = round(domain_cells  / nproma)
    #blocks_own_edges  = round(owned_edges   / nproma)
    #blocks_own_verts  = round(owned_verts   / nproma)
  
  blocks = (entities / block_size )
  blocks_per_thread = blocks / threads
  blocks_i = round(blocks + 0.25)
  blocks_per_thread_i = round(blocks_per_thread)
  #iterations = int(entities / blocks_per_thread_i)
  remain_block = blocks - (blocks_per_thread_i * threads)
  #------------------------------------------------------------------------------
  print message, " blocks=",blocks, " blocks_per_thread=", blocks_per_thread,\
  blocks_per_thread_i, " remain_block=", remain_block, remain_block*block_size

def compute_decomposition(mpi_procs, opnemp_threads):
  owned_cells = round((cells / mpi_procs) + 0.499999)
  halo_cells = 4.4 * math.sqrt(owned_cells)
  domain_cells = owned_cells + halo_cells
  owned_edges = edges / mpi_procs
  owned_verts = verts / mpi_procs


  #------------------------------------------------------------------------------
  print "==============================================================="
  print "= decomposition info"
  print "  nodes=",mpi_procs/32.0,  "  mpi_procs=",mpi_procs, " openmp threads=", opnemp_threads
  print "  Per domain:  own cells=",owned_cells, " halo cells=", halo_cells, "total_cells=", domain_cells
  print "  Per domain:  own edges=", owned_edges, " own vertices=", owned_verts
  print "=========================="
  #------------------------------------------------------------------------------
  # compute nproma parameters
  for nproma in nproma_list:
    #------------------------------------------------------------------------------
    print "= threads info, nproma=", nproma
    display_thread_load("own cells",   owned_cells, opnemp_threads, nproma)
    display_thread_load("total cells", domain_cells, opnemp_threads, nproma)
    display_thread_load("own edges",   owned_edges, opnemp_threads, nproma)
    display_thread_load("own verts",   owned_verts, opnemp_threads, nproma)
    print "=========================="

#------------------------------------------------------------------------------

def check_nproma(grid_name):

  print "==============================================================="
  print "= "
  print "= grid_name:",grid_name
  print "  cells=",cells, " edges=", edges, " vertices=", verts
  print "  vertical_levels=", vertical_levels
  print "= "
  print "==============================================================="
  #------------------------------------------------------------------------------
  #for mpi_procs in mpi_procs_list:
      #compute_decomposition(mpi_procs, omp_threads)
  for mpi_nodes in mpi_nodes_list:
      compute_decomposition(mpi_nodes * mpi_procs_pernode, omp_threads)


#------------------------------------------------------------------------------
# compute size of grid
def compute_icon_standard_size(grid_refine_root, grid_refine_level):
  cells = (20*grid_refine_root*grid_refine_root)*(4**grid_refine_level)
  verts = (10*grid_refine_root*grid_refine_root)*4**grid_refine_level+2
  edges = (30*grid_refine_root*grid_refine_root)*4**grid_refine_level
  return (cells, verts, edges)

#iconR2B04-grid
#mpi_nodes_list=[ 1.0, 4.0, 8.0, 16.0, 20.0, 27.0 ]
#nproma_list=[   16   16   20   10     8    12   ]
# 1 days = 1.5min on 1 node, 96 levels
#(cells, verts, edges)=compute_icon_standard_size(2, 4)
#check_nproma("iconR2B04-grid")

# iconR2B05-grid_dec-1082.nc
#mpi_nodes_list=[ 1.0,  5.0, 12.0, 23.0 ]
#nproma_list=  [   9    9    9     9   ]
# 1 days = 2.0min on 1 node, 96 levels
#cells = 25920
#verts = 12962
#edges = 38880
#check_nproma("iconR2B05-grid_dec-1082")

# iconR2B05-grid_dec-362.nc
##mpi_nodes_list=[ 1.0, 6.0, 27.0, 45.0 ]
##nproma_list=  [  10   9    10    12 ]
#mpi_nodes_list=[ 2.0, 6.0, 27.0, 45.0 ]
#nproma_list=  [  9   9    10    12 ]
# 1 days = 2.5min on 1 node, 96 levels
#cells = 34560
#verts = 17282
#edges = 51840
#check_nproma("iconR2B05-grid_dec-362")

# iconR2B05-grid_dec-1922
#mpi_nodes_list=[ 2.0, 6.0, 12.0,  24.0 ]
#nproma_list=  [ 12   12   10     10   ]
# 1 days = 3.5min on 1 node, 96 levels
#cells = 46080
#verts = 23042
#edges = 69120
#check_nproma("iconR2B05-grid_dec-1922")


## iconR2B06-grid_dec-2432
#mpi_nodes_list=[  2.0, 4.0, 12.0, 21.0, 42.0  ]
##nproma_list=  [  12   12   11    11    11   ]
## 1 days = 6min on 1 node, 96 levels
#cells = 58320
#verts = 29162
#edges = 87480
#check_nproma("iconR2B06-grid_dec-2432")

## iconR2B05-grid_dec-642
#mpi_nodes_list=[  4.0, 12.0,  40.0,  80.0  ]
##nproma_list=  [  12    10    12     12  ]
#mpi_nodes_list=[  80.0 ]
## 10 days = 45min on 1 node, 96 levels
#cells = 61440
#verts = 30722
#edges = 92160
#check_nproma("iconR2B05-grid_dec-642")

##iconR2B05-grid
#mpi_nodes_list=[ 4.0, 10.0, 20.0, 40.0,  80.0 107 ]
##nproma_list=[   16   16    16    16      8    12 ]
#mpi_nodes_list=[ 107.0 ]
## 1 days = 7.2min on 1 node, 96 levels
#(cells, verts, edges)=compute_icon_standard_size(2, 5)
#check_nproma("iconR2B05-grid")

## iconR2B06-grid_dec-1442.nc
##mpi_nodes_list=[ 8  20   90   180 ]
##nproma_list=[    9  12   12    12]
## 1 days = 18.5 min on 1 node, 96 levels
#mpi_nodes_list=[ 180.0 ]
#cells = 138240
#verts = 69122
#edges = 207360
#check_nproma("iconR2B06-grid_dec-1442")

# iconR2B07-grid_dec-2432.nc
#mpi_nodes_list=[ 9  19  38  152  228]
#nproma_list=[    9  12  12   12   16]
# 1 days = 37.75min on 1 node, 96 levels
#mpi_nodes_list=[ 32.0 ]
#cells = 233280
#verts = 116642
#edges = 349920
#check_nproma("iconR2B07-grid_dec-2432")

# iconR2B06-grid.nc
#mpi_nodes_list=[ 32 64 128 214 244 ]
#nproma_list=[    16 10  10  12  11 ]
#1 days = 58.25min on 1 node, 96 levels
#(cells, verts, edges)=compute_icon_standard_size(2, 6)
#mpi_nodes_list=[  244.0 ]
#check_nproma("iconR2B06-grid")

# iconR2B07-grid_dec-4322.nc
#mpi_nodes_list=[  34  70   231 240 ]
#nproma =       [  12  12    14   9 ]
# 1 days = 88.5 on 1 node, 96 levels
cells = 414720
verts = 207362
edges = 622080
mpi_nodes_list=[  35.0 ]
check_nproma("iconR2B07-grid_dec-4322.nc")

sys.exit(0)


## iconR2B06-grid_dec-1082.nc
#mpi_nodes_list=[  8.0 ]
#cells = 103680
#verts = 51842
#edges = 155520
#check_nproma("iconR2B06-grid_dec-1082")
##nproma = [ 16 ]


## iconR2B06-grid_dec-1922.nc
#cells = 184320
#verts = 92162
#edges = 276480
#check_nproma("iconR2B06-grid_dec-1922")
##nproma = [ 16 ]




sys.exit(0)

gridpoints=0.0
(cells, verts, edges)=compute_icon_standard_size(2, 6)
gridpoints=round(cells)

levels = 96.0
fields = 1.0

bytesize = 8.0

kilo = 1024.0
mega = kilo * kilo
giga = kilo * mega
tera = kilo * giga
peta = kilo * tera


size_per_3dvar = (gridpoints*levels*bytesize)/giga
size_per_step  = (gridpoints*fields*levels*bytesize)/giga
size_per_day = size_per_step * 4.0 # assuming output every six hours
size_per_year = size_per_day * 365.0 # 

print "Model output estimate:"
print " - grid points per variable and level:     %10d" % gridpoints
print " - size per variable and horizontal slice: %10.3f gb" % size_per_3dvar
print " - size per output time:                   %10.3f gb" % size_per_step  
print " - size per output day:                    %10.3f gb" % size_per_day  
print " - size per output year:                   %10.3f gb" % size_per_year
