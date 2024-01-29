#!/usr/bin/env python

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

from cdo import *
import os,sys,math
import matplotlib
matplotlib.use('AGG')
import numpy as np
import mpl_toolkits.basemap as bm
import matplotlib.pyplot as plt
cdo        = Cdo()
cdo.debug  = 'DEBUG' in os.environ

# USAGE {{{ =============================================================================
def usage():
    return """
# USAGE =================================================================================
#   ./calc_psi.py <ifile> VAR=<varname> PLOT=<plotfile> CMAP=<colormap> REMAP=<True> LEVELS=<levels> AREA=<area> ASPECT=<aspect>
#
# default values are:
#   varname  = u_vint_acc
#   plotfile = psi.png'    (other output types:  png, pdf, ps, eps and svg)
#   colormap = jet         (see http://matplotlib.org/examples/color/colormaps_reference.html for more)
#   levels   = -500,-200,-150,-100,-75,-50,-30,-20,-10,-5,0,5,10,20,30,50,75,100,150,200,500
#   remap    = True        (expect icon input, so that remapping to r360x180 is done internally;
#                           can be set to False, false or 0 to disable)
#   area     = global      (if set to another value, the continents will not be drawn)
#   aspect   = auto        (other values: equal or a number)
# =======================================================================================
"""
# }}} ===================================================================================
# INPUT HANDLING {{{ ====================================================================
if len(sys.argv) < 2:
    print("Provide an input file")
    print(usage())
    exit(1)

inputfile = sys.argv[1]
# stop if file cannot be read in
if not os.path.isfile(inputfile):
    print("Cannot read input: "+inputfile)
    print(usage())
    exit(1)

# }}} ===================================================================================
# OPTION HANDLING {{{ ===================================================================
options = {'VAR'     : 'u_vint_acc',
           'REMAP'   : True,
           'PLOT'    : 'psi.png',
           'CMAP'    : 'BrBG',
           'LEVELS'  : [-250,-150,-100,-75,-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,30,40,50,75,100,150,250],
           'AREA'    : 'global',
           'SHOWMAP' : True,
           'TITLE'   : '',
           'DEBUG'   : False,
           'BOX'     : '',
           'WRITEPSI': False,
           'ASPECT'  : 'auto'}

optsGiven = sys.argv[2:]
for optVal in optsGiven:
    key,value    = optVal.split('=')
    if 'LEVELS' == key:
        value = value.split(',')
        if 1 == len(value):
            value = int(value[0])
        else:
            value = map(lambda x: float(x), value)
    if key in ['REMAP','DEBUG']:
        if value in ['false','False','FALSE','0']:
            value = False

    options[key] = value

varName    = options['VAR']
remapInput = options['REMAP']
plotfile   = options['PLOT']
colormap   = options['CMAP']
levels     = options['LEVELS']
area       = options['AREA']
showMap    = options['SHOWMAP']
aspect     = options['ASPECT']
# }}} ===================================================================================
# DATA PREPARATION {{{ ==================================================================

# if the input is NOT global, the remapping has to be limitted to the area
if 'global' != area and True == remapInput:
    # read min/max of input lon/lat
    file_h      = cdo.readCdf(inputfile)
    v           = file_h.variables[varName]
    x,y         = getattr(v,'coordinates').split(' ')
    xArr, yArr  = file_h.variables[x][:], file_h.variables[y][:]
    x_min,x_max = xArr.min() ,xArr.max()
    y_min,y_max = yArr.min() ,yArr.max()

    if 'radian' == getattr(file_h.variables[x],'units'):
        rad2deg = 45.0/math.atan(1.0)
        x_min,x_max = x_min*rad2deg,x_max*rad2deg
        y_min,y_max = y_min*rad2deg,y_max*rad2deg

    x_min,x_max = math.floor(x_min), math.ceil(x_max)
    y_min,y_max = math.floor(y_min), math.ceil(y_max)

# remapnn  to regular 1deg grid - violated divergence-free velocity at ocean boundaries are visible
# remapcon to regular 1deg grid - smoothing on boundary problems
# replace missing value with zero for later summation
if remapInput:
    ifile = cdo.setmisstoc(0.0,input = '-remapcon,r360x180 '+inputfile,options='-P 8')
#   ifile = cdo.setmisstoc(0.0,input = '-remapnn,r360x180 '+inputfile,options='-P 8')
    if 'global' != area:
        ifile = cdo.sellonlatbox(x_min,x_max,y_min,y_max,input = ifile)
else:
    ifile = cdo.setmisstoc(0.0,input = inputfile)

file_h    = cdo.readCdf(ifile)
var       = file_h.variables[varName]
varData   = var[:]
varDims   = var.dimensions
# read in dimensions: expectes is 2d with time axis (time,lat,lon)
a         = map(lambda x: file_h.variables[x][:], varDims)
#times, depth, lats, lons = a[0], a[1], a[2], a[3] # MPIOM psi input
times, lats, lons = a[0], a[1], a[2]

if options['DEBUG']:
    print("# DEBUG ===================================================================")
    print(inputfile)
    print(varName)
    print(plotfile)
    print(colormap)
    print("#==========================================================================")
    print(varData)
    print(varData.shape)
    print(varDims)
    print(times)
    print(lons)
    print(lats)
    print("#==========================================================================")
    print(options)
    print("# DEBUG ===================================================================")

# use first timestep only
if times.size > 1: 
    print('Will only use the first timestep!!!!!!!!!')
#varData = varData[-1,0,:,:] # MPIOM psi input
varData = varData[-1,:,:]

# }}} ===================================================================================
# CALC PSI {{{ ==========================================================================
psi      = np.array(varData)
psi[:,:] = varData[:,:]
# parial sum from south to north
for lat in range(lats.size-2,-1,-1):
    psi[lat,:] = psi[lat,:] + psi[lat+1,:]

erad = 6.371229e6
pi   = 3.141592653
dist = (pi/lats.size)*erad
psi  = -psi * dist * 1.0e-6
#psi  = psi * 1.0e-6 / 1025.0 # MPIOM psi input
# }}} ===================================================================================
# WRITE PSI {{{ ==========================================================================
if options['WRITEPSI']:
  # create a copy of the input data and rename the variable
  psiFileName = cdo.chname('%s,psi'%(varName),input = ifile, output = os.path.dirname(inputfile)+'/psi_remapped.nc',force=True,options='-r')

  try:
    from netCDF4 import Dataset
  except:
    print("Could not load netCDF4 module!")

  psiFile                  = Dataset(psiFileName,'r+')
  psiFileVar               = psiFile.variables['psi']
  psiFileVar.standard_name = "barotropic stream function"
  psiFileVar.long_name     = "psi"
  psiFileVar.units         = "Sv"
  psiFileVar[1,:,:]        = psi[:,:]
  psiFile.history          = psiFile.history + ' changed by calc_psi.py'
  psiFile.close()
# }}} ===================================================================================
# PLOTTING {{{ ==========================================================================
#fig = plt.figure(figsize=(10,5))
fig = plt.figure()


# limit area if BOX is given
if '' != options['BOX']:
  lonMin,lonMax,latMin,latMax = options['BOX'].split(',')
  latMin = math.floor(float(latMin)) 
  lonMin = math.floor(float(lonMin))
  latMax = math.ceil( float(latMax))
  lonMax = math.ceil( float(lonMax))
else:
  latMin = math.floor(lats.min()) 
  lonMin = math.floor(lons.min())
  latMax = math.ceil(lats.max())
  lonMax = math.ceil(lons.max())

matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
lon2d, lat2d = np.meshgrid(lons, lats)
mapproj = bm.Basemap(projection='cyl', 
                     llcrnrlat=latMin, 
                     llcrnrlon=lonMin,
                     urcrnrlat=latMax,
                     urcrnrlon=lonMax)

if 'global' == area or True == showMap:
    mapproj.drawcoastlines(linewidth=.2)
    mapproj.fillcontinents(color='grey',lake_color='k')

mapproj.drawmapboundary(fill_color='0.99')
mapproj.drawparallels(np.array([-80,-60,-40,-20, 0, 20,40,60,80]), labels=[1,1,0,0],fontsize=6,linewidth=0.2)
mapproj.drawmeridians(range(0,360,30), labels=[0,0,0,1],fontsize=6,linewidth=0.2)
lonsPlot, latsPlot = mapproj(lon2d, lat2d)

plt.gca().set_aspect(aspect)

# filled contours
CS = plt.contourf(lonsPlot, latsPlot, psi,
                    levels,
                    extend='both',
                    cmap=colormap)
# contour lines
CSBar = plt.contour(lonsPlot,
                    latsPlot,
                    psi,
                    CS.levels[::1],
                    colors = ('k',),
                    linewidths = (0.5,),
                    origin = 'lower')

# contour line labels
plt.clabel(CSBar, fmt = '%2.1f', colors = 'k', fontsize=6)

# colorbar
if 'global' == area:
    orientation = 'horizontal'
else:
    orientation = 'vertical'

cbar = plt.colorbar(CS,orientation=orientation)
cbar.set_label("Sv")

if ('' == options['TITLE']):
  title = "Bar. Streamfunction for\n"+inputfile
else:
  title = options['TITLE']
plt.suptitle(title,fontsize=9)
plt.title("psi",fontsize=8,loc='left')

fig.savefig(plotfile)
# }}}
# =======================================================================================
# vim:fdm=marker
