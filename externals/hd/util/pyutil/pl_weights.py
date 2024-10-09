#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 
# pl_weights.py - Plot weights and mappings for ICON couplings via YAC with HD 
# 
# Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#

"""
###############################################################################
## Read and Convert YaC weight files from 1D coordinate fields as scatter plots
#    1) Plot mapping & plot weights and write NC file of sending HD boxes for HD --> ICON-O 
#    2) Plot weights and write NC file of receiving HD boxes for ICON-A --> HD 
#
###############################################################################
## Author: Stefan Hagemann (Hereon, Institute of Coastal Systems, GER)
##
"""
import sys,os,subprocess,datetime
import fnmatch
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from util_plot import *
from subprocess import call
from util_var import Variable
import util_inout as uio
from math import degrees
#
import psyplot.project as psy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
# 
import yaml

if os.path.exists("./config_pl_weights.yml"):
    cfg = yaml.safe_load(open("./config_pl_weights.yml"))   # Use path to config file
else:
    print('ERROR: You need to copy a YAML file config_pl_weights.yml to the current directory')
    print('   This file is necessary to steer the program. Example files are provided in ./example/')
    print('   config_pl_weights_hdtoicon.yml    Plot mapping and plot/generate mask of HD boxes sending output to ICON-O')
    print('   config_pl_weights_icontohd.yml    Plot/generate mask of HD boxes receiving input from ICON-A')
    exit()


# Directories     
DATA=cfg['DATA']
DPLOT=cfg['DPLOT']
DNICON=cfg['DNICON']
#Input file
ifile=cfg['IFILE']   # 0: use DNINP, 1: generate file name with cvar
dninp=cfg['DNINP']
iwork=cfg['IWORK']
ireg=cfg['IREGION']  # Region no.: 0 = globe, 10=euro5min
if 'IDOMAIN' in cfg:
  idomain=cfg['IDOMAIN']
else:
  idomain=0
#output files
iout=cfg['IOUT']   # 0: use DNOUT, 1: generate file name with PLTAG and cvar
dnout=cfg['DNOUT']
#
PLTAG=cfg['PLTAG']
LINTER=cfg['LINTER']       # Interactive run: no/yes
# Special plots
# Initial Plot options
xmin=cfg['XMIN']
xmax=cfg['XMAX']
xlow=cfg['XLOW']
xmiss=cfg['XMISS']
ncol=cfg['NCOL']
icoltyp=cfg['ICOLTYP']
if 'IVARBAR' in cfg:
  ivarbar=cfg['IVARBAR']
else:
  ivarbar=0
if 'IBARLOC' in cfg:
  ibarloc=cfg['IBARLOC']        # Location/orientation of colorbar. 0=horiz., 1=vert.
else:
  ibarloc=0
if 'IDATE' in cfg:
  idate=cfg['IDATE']
else:
  idate=1
pathname = os.path.dirname(sys.argv[0])        
os.chdir(os.path.abspath(pathname))


def main():
  print('\n>>>>>>>> Plot Mask is running! <<<<<<<<<<<<')
  print('DATA directory:', DATA)
  print('Note that the following files should be located in/linked from the DATA directory:')
#   
  xmiss=-9999.     # Missing value
#
# ********* Basic settings **************** 
# *** TODO put basic settings in subroutine of util_plot
  opt=PlotOptions()
  opt.pdate = idate 
  opt.ufak = 1.
  opt.xmin = xmin ; opt.xmax = xmax ; opt.xlow = xlow ; opt.ncol = ncol
  opt.icoltyp = icoltyp
  opt.cmap = get_colortype(icoltyp, False)
  opt.ibarloc = ibarloc
  opt.ptitle = PLTAG
  opt.ifile = ifile ; opt.iout = iout ; opt.dninp = dninp ; opt.dnout = dnout

  reg = Region(ireg)
  status = reg.get_region()
  loop = True
  lmap = True
  lmask=[]
#
# Which variable should be read?
# iwork = 1 --> read src_address and plot ICON --> HD mapping file
# iwork = 2 --> read dst_address and plot mask of HD boxes that receive data from ICON-A
  csrc = 'src_address'   # Name of variable in input data file
  cdst = 'dst_address'
  if iwork == 1:
      cvar = csrc   # Name of variable in input data file
      print(' Weight file - HD to ICON: ', dninp)
      print('     ICON ocean grid file: ', DNICON)
  elif iwork == 2:
      cvar = cdst
      print(' Weight file - ICON to HD: ', dninp)
      print('ICON atmosphere grid file: ', DNICON)
  else:
      print('IWORK = ', iwork, ' is NOT defined --> ABBRUCH!')
      exit() 
  print(' ')
  var=Variable(100)
  var.longname = cvar
#
# Grid definition
  if idomain == 0:
      nlon = 720  ; nlat = 360
#    flat = df.variables['lat'][:]
#    flon = df.variables['lon'][:]
      dlat = 0.5
      flat_n = 90.   ; flat_s = -90.
      flon_w = -180. ; flon_e = 180.

  print('Target Grid: ', nlon, ' * ',nlat)

  while (loop):
    if (LINTER):
      plot_menu(opt, reg, lmap)
    else:
      loop = False
#
#   Directory path to observations 
##    DDIR = DATA + os.sep + str(EXP)
    DDIR = DATA
    print ('DIR. DATA:', DDIR)

    if ifile == 0:
      DNAM = DATA + os.sep + opt.dninp
    else:
#
#     *** using fnmatch if only part of the filename is known
      file=fnmatch.filter(os.listdir(DDIR), cvar + '_' + '*' + '.nc')
      print ('Data file: ', file[0])
      DNAM = DATA + os.sep + file[0]

#
#   *** Read weight file
    df = Dataset(DNAM, 'r')
    # Source grid
    vardims=len(df.variables[csrc].dimensions)
    print('SRC vardims: ', vardims,': ', df.variables[csrc].dimensions)

    if vardims == 1:
        fsrc = df.variables[csrc][:]
        fweight = df.variables['remap_matrix'][:]
        fdst = df.variables[cdst][:]
    elif vardims == 2:
        fsrc = df.variables[csrc][:,:] 
        fweight = df.variables['remap_matrix'][:,:]
        fdst = df.variables[cdst][:,:]
    df.close()
    
    print (csrc + ' index ranges from ', fsrc.min(), ' to ', fsrc.max())
    print (cdst + ' index ranges from ', fdst.min(), ' to ', fdst.max())
    print ('Weights range from ', fweight.min(), ' to ', fweight.max())
    if iwork == 1:
        fdata = fsrc
    elif iwork == 2:
        fdata = fdst
    
      
    valmin=fdata.min()
    valmax=fdata.max()
    print (cvar + ' index ranges from ', valmin, ' to ', valmax)
#
    fmask = generate_mask(nlon, nlat, fdata, fweight)

    filetag = dnout

    print('fmask len:', len(fmask), ' ranges from ', fmask.min(), ' to ', fmask.max()) 

    plot_map(nlon, nlat, fmask, var, reg, lmask,
             flat_n, flat_s, flon_w, flon_e, 
             opt, xmiss, filetag, DPLOT)

    #
    # Write NC file with mask
    out_data = uio.Mapdata()
    out_data.nlon = nlon ;  out_data.nlat = nlat
    out_data.vars = ['slm'] ; out_data.units = ['-'] 
    # 
    # Coordinates
    out_data.xlat = []
    for jb in range(nlat):
        out_data.xlat.append(flat_n - ((jb+0.5) * dlat) )
    out_data.xlon = []
    for jl in range(nlon):
        out_data.xlon.append(flon_w + ((jl+0.5) * dlat) )

    dnmas = DPLOT + os.sep + dnout + '.nc' 

    ofile = out_data.openw(dnmas)
    out_data.write(fmask)
    out_data.close(ofile)

    nval = np.count_nonzero(fmask)
    print('NC file written: ', dnmas, ' with ', nval, ' masked boxes') 
    
    jb60s = nlat - int(nlat / 6) 
    nval = np.count_nonzero(fmask[0:jb60s, 0:nlon])
    print('  Masked boxes without Antarctica, north of 60°S: ', nval) 
#
#   *** Plot Mapping for HD --> ICON-O
    if iwork == 1:
        plot_mapping_hd_to_icon(nlon, fdata, fdst, out_data, DNICON, var, reg, opt, DPLOT, PLTAG)
    elif iwork == 2:
        plot_weights_on_icon(fsrc, fweight, DNICON, var, reg, opt, DPLOT)
    else:
        exit()


# ***************************************************
def plot_mapping_hd_to_icon(nlon, fdata, fdst, out_data, DNICON, var, reg, opt, DPLOT, PLTAG):
# ***************************************************
#
#   *** Colormap and Projection
    icoltyp2 = 15
    cstrmap = get_colortype(icoltyp2, False)
    colmap = plt.get_cmap(cstrmap, 5)
    # colmap.set_over('gray')
    cproj = ccrs.PlateCarree()
 
    fig2, ax2 = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
#
#   *** Open ICON dataset
    ds = psy.open_dataset(DNICON)

    #from IPython import embed
    #embed()
    #ds2 = ds.where(ds['cell_sea_land_mask']  != -2, drop=True)
    #ds['cell_sea_land_mask'] = ds.cell_sea_land_mask.load().astype(float)
    #ds['cell_sea_land_mask'].values[ds['cell_sea_land_mask'].values == -2] = np.nan
    indices = np.where(ds['cell_sea_land_mask'].values != -2)[0]
    ds = ds.load()
    ds2 = ds.isel(cell=indices)

#   *** Plot dataset
    plot_icon(ds2.psy.cell_sea_land_mask, ax=ax2, cmap=colmap)
#
#  
    feld = ds['cell_sea_land_mask']
    print('nicon: ', len(feld))

    xlon_icon = feld.clon.values    # 1D np.ndarray of longitudes at the desired locations
    xlat_icon = feld.clat.values    # 1D np.ndarray of latitudes at the desired locations

    xlon_src, xlat_src, xlon_dst, xlat_dst = generate_mapping(fdata, fdst, nlon,   \
            out_data.xlon, out_data.xlat, xlon_icon, xlat_icon)

    plot_mapping_arrows(xlon_src, xlat_src, xlon_dst, xlat_dst, reg, cproj)
#
#   *** Set plot title, axis range & legends and background map settings
    var.longname = 'Mapping'
    set_plot_char(ax2, cproj, opt, reg, var)

    plotfile = DPLOT + os.sep + 'mapping_' + PLTAG + '.pdf'
    plt.savefig(plotfile, dpi=200, format='pdf')

    plt.show()
    plt.close()
    print('Plot --> ' + plotfile)

# ***************************************************
def plot_weights_on_icon(fsrc, fweight, DNICON, var, reg, opt, DPLOT):
# ***************************************************
#
#   *** Colormap and Projection
    icoltyp2 = 0
    cstrmap = get_colortype(icoltyp2, False)
    colmap = plt.get_cmap(cstrmap, 5)
    # colmap.set_over('gray')
    cproj = ccrs.PlateCarree()
 
    fig2, ax2 = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
#
#   *** Open ICON dataset
    print('Open ICON file: ', DNICON)
    ds = psy.open_dataset(DNICON)

##    fmas = ds['cell_area'].copy(deep=True)
    fmas = xr.zeros_like(ds['cell_area'])
    nbox = fmas.size    
    print('ncells = ', nbox)

    fmas.psy.base = ds
    fmas.encoding["coordinates"] = ds.cell_area.encoding["coordinates"]  

##    for i in range(nbox):
##        fmas[i] = 0.
##    fmas * 0.
    nlink = len(fsrc)
    print('nlink = ', nlink, ' = ', len(fweight))
    for i in range(nlink):
##        isrc = fsrc[i]-1   # This is how it is supposed to be.
##        isrc = fsrc[i]     #This looks worse
        isrc = fsrc[i] - 2
#        print(isrc, '. weight: ', fweight[i])
        fmas[isrc] = fmas[isrc] + float(fweight[i])
    print('fmas ranges from ', float(fmas.min()), ' to ', float(fmas.max()) )
    #
    # Write NC data
    icon_data = uio.Mapdata()
    icon_data.nlon = nbox ;  icon_data.nlat = 0
    icon_data.vars = ['weights'] ; icon_data.units = ['-'] 
    # 
    # Coordinates
    icon_data.xlat = ds.variables['clat'][:]
    icon_data.xlon = ds.variables['clon'][:]
 
    dntmp = DPLOT + os.sep + 'weights_plain' + '.nc' 
    dnw = DPLOT + os.sep + 'weights_on_icon' + '.nc' 

    print('fmas= ',fmas.values)

    ofile = icon_data.openw(dntmp)
    icon_data.write(fmas.values)
    icon_data.close(ofile)
    call('cdo setgrid,' + DNICON + ' ' + dntmp + ' ' + dnw , shell=True)
    print('Weights written to ' + dnw)
    call('rm ' + dntmp, shell=True)
#
#   *** Plot dataset 
    plot_icon(fmas, ax=ax2, cmap=colmap)
##    plot_icon(fmas, ax=ax2, cmap=colmap, vmin=opt.xmin, vmax=opt.xmax)    #TODO define formatoptions!
#
#   *** Set plot title, axis range & legends and background map settings
    var.longname = 'ICON-Src'
    set_plot_char(ax2, cproj, opt, reg, var)

    plotfile = DPLOT + os.sep + 'weights_on_icon' + '.pdf'
    plt.savefig(plotfile, dpi=200, format='pdf')

    plt.show()
    plt.close()
    print('Plot --> ' + plotfile)



# ***************************************************
def generate_mask(nlon, nlat, fdata, fweight):
# ***************************************************
#
#   Generate Mask/data file on 2D grid [0:nlon[ [0:nlat[
#   fdata contains indices on 1D grid [0:nlon*nlat[
#   In mo_couling_hd, global cell_index already runs from 0 to nlon*nlat -1. However, the 
#      addresses in the weight file are global_id + 1 --> index_py = fdata - 1.
#
    nbox = len(fdata)
    fmask = np.zeros((nlat, nlon))   # [y, x] Coordinates
   
    for i in range(nbox):
        index_py = fdata[i] - 1   
        jb = int(float(index_py) / float(nlon))
        jl = index_py - jb * nlon 
        if fdata[i] == nlon * nlat:
            print('Gen.Mask confirms to substract 1 from 1Dindex = ', fdata[i], ' --> Py jl,jb ', jl, jb)
        elif jl == nlon or jb == nlat:
            print('Gen.Mask error: jl or jb too large: ', fdata[i], ' --> Py jl,jb ', jl, jb)
            exit()
        elif jl < 0 or jb < 0:
            print('Gen.Mask error: jl or jb negative: ', fdata[i], ' --> Py jl,jb ', jl, jb)
            exit()
        if i == 3175:
            print(i, 'Gen.Mask: 1Dindex = ', fdata[i], ' --> Fort jl,jb ', jl+1, jb+1)
        fmask[jb,jl] = fmask[jb,jl] + fweight[i]

    return fmask

# ***************************************************
def generate_mapping(fsrc, fdst, nlon, xlon_hd, xlat_hd, xlon_icon, xlat_icon):
# ***************************************************
#
#   Generate Mapping coordinates from HD to ICON:
#   HD coordinates: xlon in [0:nlon[, xlat in [0:nlat[
#   In mo_couling_hd, global cell_index already runs from 0 to nlon*nlat -1. However, the 
#   addresses in the weight file are global_id + 1. This may change in future YAC versions.
#
    nbox = len(fsrc)
    nicon = len(xlat_icon)
    xlon_src = []
    xlat_src = []   
    xlon_dst = []
    xlat_dst = []   

    for i in range(nbox):
        index_py = fsrc[i] - 1    # -1 due to YAC definition as global_id+1 (see above)
        jb = int(float(index_py) / float(nlon))
        jl = index_py - jb * nlon 
        xlon_src.append(xlon_hd[jl])
        xlat_src.append(xlat_hd[jb])

        if fdst[i] >= nicon:
            print('i=', i, ' ERROR: fdst[i] ', fdst[i], ' >= nicon=', nicon)
            print('You may have read the wrong weight file?!')

        xlon_dst.append(degrees(xlon_icon[fdst[i]-2]))   # -1 wegen Python starting at 0   before bugfix: is there an impact??
        xlat_dst.append(degrees(xlat_icon[fdst[i]-2]))
#        xlon_dst.append(degrees(xlon_icon[fdst[i]-1]))   # -1 wegen Python starting at 0
#        xlat_dst.append(degrees(xlat_icon[fdst[i]-1]))
        if 41 < xlat_src[i] < 43 and 10 < xlon_src[i] < 12:
            print(i, jl+1, jb+1, 'Mapping: ', xlon_src[i], xlat_src[i], ' -> ', fdst[i], xlon_dst[i], xlat_dst[i]) 

    return xlon_src, xlat_src, xlon_dst, xlat_dst

#***************************************************
if __name__ == "__main__":
    main()
#***************************************************

