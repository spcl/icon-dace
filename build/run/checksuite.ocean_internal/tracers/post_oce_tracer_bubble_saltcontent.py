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

from __future__ import print_function
import unittest,os,tempfile,sys,glob
from stat import *
from cdo import *
cdo = Cdo()
import numpy as np
from matplotlib import pylab as pl
import scipy.io
import datetime
import os.path
import dateutil.parser

def parseOptions():
  
  options = {
              'GRID'        : 'global',                               # shortcut for the used gridtype(global,box,channel)
              'PLOTDIR'     : '/plots',
              'EXP'         : '../../../experiments/oce_tracer_bubble',                    # default experiment name
              'FILEPATTERN' : '../../../experiments/oce_tracer_bubble/oce_tracer_bubble_', # default output file pattern
              'DATE'        : '2001-01-01T00:00:00Z',
              'DATUM'       : '20010101T000000Z',
              'RB'          : 'R2B04',
              'PROCS'       : 8                                       # number of threads/procs to be used for parallel operations
             }
             
  return options

def saltcontent(options):

  gridfile = options['EXP'] + '/icon' + options['RB'] + '-ocean_etopo40_planet.nc'  
  f = options['FILEPATTERN'] + options['RB'] + '_oce_' + options['DATUM'] + '.nc'
#  fid = scipy.io.netcdf_file(f,'r')
#  var = fid.variables
#  s_acc = var['s_acc'][:,:,:].copy() 	# [time, depth, ncells]
  
  salt_content = cdo.vertsum(input = " -fldsum -mulc,1025.022 -mul -selname,s_acc " + str(f) + " -mul -selname,prism_thick_c " + str(f) + " -selname,cell_area " + str(gridfile), returnCdf=True).variables['s_acc'][:] 

  pl.figure(1)
  pl.subplot(211)
  pl.title('cumulative sum of global saltcontent variation')
  pl.plot((salt_content-salt_content.mean()).cumsum())

  pl.subplot(212)
  t = salt_content.flatten()
  pl.title('global saltcontent')
  pl.plot(t)
  pl.show()
  
  return

options = parseOptions()
saltcontent(options)