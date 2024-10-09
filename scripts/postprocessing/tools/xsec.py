#!/usr/bin/env python

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

import pylab as pyl
from cdo import *
from string import rjust
import sys,os
from jobqueue import mpJobQueue

#=============================================================================== 
# input file (expected is ICON ocean model result)
ifile     = sys.argv[1]
# variable name to plot (default T)
if (2 == sys.argv.__len__()):
  varname = 'T'
else:
  varname = sys.argv[2]
# location of the x-section
lon=[0,1]
lat=[-20,20]
# image scaling factor
scale = 5
# display debug stuff the script is started with DEBUG=<some value>
debug = 'DEBUG' in os.environ
#=============================================================================== 
def makeJob(args):
  t,var,shape = args
  ofile       = 'xsect_'+rjust(str(t),4,'0')+'.png'
  v           = var[t,0:shape[1],0,:]
  pyl.subplots_adjust(left=.05,right=.95,bottom=.1,top=.85)
  pyl.imshow(v,interpolation="nearest")
  pyl.savefig(ofile,bbox_inches='tight',dpi=200)
  print('save image for timestep %i'%(t))
#=============================================================================== 

cdo       = Cdo()
cdo.debug = debug
q = mpJobQueue(8,True)

_ifile = cdo.sellonlatbox('%i,%i,%i,%i'%(lat[0],lat[1],lon[0],lon[1]),
                          input=' -remapnn,r360x180  -sellonlatbox,%i,%i,%i,%i  -selname,%s %s'%(
                            lat[0]-1,
                            lat[1]+1,
                            lon[0]-1,
                            lon[1]+1,
                            varname,ifile),
                          options='-P 8')

var = cdo.readCdf(_ifile).variables[varname]

shape = var.shape

pyl.rcParams['figure.figsize']=[shape[3]/scale,shape[1]/scale]

for t in range(0,shape[0]):
  q.push(makeJob,(t,var,shape))

q.run()

#os.system("convert xsect_[0-9][0-9][0-9][0-9].png xsect.gif")
