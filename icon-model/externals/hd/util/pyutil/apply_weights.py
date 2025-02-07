#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# apply_weights.py - Reading and applying YAC weight files ICON-A and ICON-O coupings with HD
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
## Program reads YaC weight files with 1D coordinate fields and applies the 
#    mapping to agiven NetCDF input file
#    1) HD --> ICON-O 
#    2) ICON-A --> HD 
#
# Moritz: src/tgt_address im AC weight file sind (global_id + 1). Wenn diese vom Benutzer angegeben wurden, 
#         dann hängt es davon ab ob die minimale id bei 0 oder 1 liegt.
#         Falls die ID nicht angegeben wurden, hat YAC die selbst erstellt...
#         dann lässt sich die Datei schwer interpretieren...
#
###############################################################################
## Author: Stefan Hagemann (Hereon, Institute of Coastal Systems, GER)
##
"""
import sys,os,subprocess,datetime
import fnmatch
import numpy as np
from netCDF4 import Dataset
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

if os.path.exists("./config_apply_weights.yml"):
    cfg = yaml.safe_load(open("./config_apply_weights.yml"))   # Use path to config file
else:
    print('ERROR: You need to copy a YAML file config_apply_weights.yml to the current directory')
    print('   This file is necessary to steer the program. Example files are provided in ./example/')
#    print('   config_pl_weights_hdtoicon.yml    Plot mapping and plot/generate mask of HD boxes sending output to ICON-O')
    print('   config_apply_weights_icontohd.yml   Map ICON-A data to the HD grid')
    exit()


# Directories     
DATA=cfg['DATA']
DNICON=cfg['DNICON']
#Input files
dninp=cfg['DNINP']     # File that should be remapped
dnweight=cfg['DNWEI']  # weight file
iwork=cfg['IWORK']     # 1: 1D->2D mapping, 2: 
if 'IGRID' in cfg:
    igrid=cfg['IGRID']
else:
    igrid=0
#output files
iout=cfg['IOUT']   # 0: use DNOUT, 1: generate file name with PLTAG and cvar
dnout=cfg['DNOUT']
#
pathname = os.path.dirname(sys.argv[0])        
os.chdir(os.path.abspath(pathname))


def main():
    print('>>>>>>>> Apply YAC weigths is running! <<<<<<<<<<<<')
#   
    xmiss=-9999.     # Missing value
#
# ********* Basic settings **************** 
    csrc = 'src_address'   # Source grid indices
    cdst = 'dst_address'   # Target grid indices
#
##    var=Variable(100)
##    var.longname = cvar
#
# Target grid definition
    if igrid == 0:
        nlon = 720  ; nlat = 360
#      flat = df.variables['lat'][:]
#      flon = df.variables['lon'][:]
        dlat = 0.5
        flat_n = 90.   ; flat_s = -90.
        flon_w = -180. ; flon_e = 180.

    print('Target Grid: ', nlon, ' * ',nlat)

#
#   Directory path to data 
##    DDIR = DATA + os.sep + str(EXP)
    DDIR = DATA
    print ('DIR. DATA:', DDIR)

    DNAM = DATA + os.sep + dninp
    DNW = DATA + os.sep + dnweight
    print('Weights: ', DNW)
    print('Data: ', DNAM)

#
#   *** Read weight file
    df = Dataset(DNW, 'r')
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
    print (' Sum of weights: ', fweight.sum())
#
    if iwork == 1:
##        idelta = 2        # Python index = dst-index - idelta
        fmap = mapping_1d_to_2d(nlon, nlat, fdst) 
    elif iwork == 2:
        fmap = mapping_1d_to_2d(nlon, nlat, fsrc) 
    print('fmap ranges from: y ', fmap[0,:].min(), ' to ', fmap[0,:].max()) 
    print('                  x ', fmap[1,:].min(), ' to ', fmap[1,:].max()) 
   
#
#   *** Read data file
    in_data = uio.Mapdata()
    ifile = in_data.openr(DNAM)
    varlist = in_data.varinfo(ifile)

    print('varist: ', varlist)

    if len(varlist) == 1:
        cvar = varlist[0]
    else:
        print('More than one variable in NC file --> extend program by selecting variable name via YAML file')
        exit()

    if iwork == 1:
##        fdata = ifile[cvar][:]  # Data in NetCDF structure
        fdata = ifile.variables[cvar][:]
    elif iwork == 2:
        fdata = ifile.variables[cvar][:][:]
        print(' fd(23,23) ->', fdata[23][23] )
    print(fdata)

    nmes = fdata.size
    print('fdata(', nmes,') ranges from: ', fdata.min(), ' to ', fdata.max()) 
    print (' Sum of fdata: ', fdata.sum())

# 
#   *** Mapping
    if iwork == 1:
        fout = remap_1d_to_2d(fsrc, fweight, fmap, fdata[:], nlon, nlat)
        print('fout ranges from: ', fout.min(), ' to ', fout.max()) 
#
#       *** Write remapped file into NC file
        fileout = DATA + os.sep + dnout + '.nc' 
        out_data = uio.Mapdata()
        out_data.nlon = nlon ;  out_data.nlat = nlat
        out_data.vars = [cvar] ; 
        if hasattr(ifile.variables[cvar], 'units'):
            out_data.units = [ifile.variables[cvar].getncattr('units')]
        else:
            out_data.units = ['-']
        # 
        # Coordinates
        if igrid == 0: 
            out_data.xlat = []
            for jb in range(nlat):
                out_data.xlat.append(flat_n - ((jb+0.5) * dlat) )
            out_data.xlon = []
            for jl in range(nlon):
                out_data.xlon.append(flon_w + ((jl+0.5) * dlat) )
        ofile = out_data.openw(fileout)
        out_data.write(fout)
        out_data.close(ofile)
        print('NC file written: ', fileout) 
#
    elif iwork == 2:
#
#       *** Open ICON dataset
        print('Open ICON file: ', DNICON)
        ds = psy.open_dataset(DNICON)
        fmas = ds.variables['cell_sea_land_mask'][:]
##        fmas = xr.zeros_like(ds['cell_sea_land_mask'])
        nicon = fmas.size    
        print('ncells = ', nicon)
#
        fout = remap_2d_to_1d(fmap, fweight, fdst, fdata, nicon)
        print('fout ranges from: ', fout.min(), ' to ', fout.max()) 
        print (' Sum of fout: ', fout.sum())
        write_remapped_data(ds, DATA, fout)
    

# ***************************************************
def mapping_1d_to_2d(nlon, nlat, fdst):
# ***************************************************
#
#   Rewrite the YAC 1D->1D Mapping into a 1D->2D mapping, i.e. on a 2D grid [0:nlon[ [0:nlat[
#   fdata contains indices on 1D grid [0:nlon*nlat[
#
#   In mo_couling_hd, global cell_index already runs from 0 to nlon*nlat -1. However, the 
#      addresses in the weight file are global_id + 1 --> index_py = fdst - 1.
#
    nbox = len(fdst)
    fmap = np.zeros((2, nbox))   # (0,:) -> y, (1::) -> x
   
    for i in range(nbox):
        index_py = fdst[i] - 1   ### see above
        jb = int(float(index_py) / float(nlon))
        jl = index_py - jb * nlon 
        if fdst[i] == nlon * nlat:
            print('mapping_1d_to_2d confirms to substract 1 from 1Dindex = ', fdst[i], ' --> Py jl,jb ', jl, jb)
        elif jl >= nlon or jb >= nlat:
            print('mapping_1d_to_2d error: jl or jb too large: ', fdst[i], ' --> Py jl,jb ', jl, jb)
            exit()
        elif jl < 0 or jb < 0:
            print('mapping_1d_to_2d error: jl or jb negative: ', fdst[i], ' --> Py jl,jb ', jl, jb)
            exit()
        if i == 3175:
            print(i, 'mapping_1d_to_2d: 1Dindex = ', fdst[i], ' --> Fort jl,jb ', jl+1, jb+1)
        fmap[0,i] = jb
        fmap[1,i] = jl

    return fmap

# ***************************************************
def remap_1d_to_2d(fsrc, fweight, fmap, fdata, nlon, nlat):
# ***************************************************
#   Remap the data fdata 1D->2D mapping using the weights fweight and the mapping fmap
#
    nmap = len(fweight)
    fout = np.zeros((nlat,nlon))
#
    for i in range(nmap):
       icell = fsrc[i] - 2        # Test yielded that -2 is necessary for a python index 
##       icell = fsrc[i] - 1 

###       if abs(fdata[icell]) > 0:
       if abs(fweight[i]) > 0:
           jb = int(fmap[0,i])     
           jl = int(fmap[1,i])
##           if icell == 76: 
##               print(icell, ' --> jb = ', jb, ' jl = ', jl)
##               print('        --> fdata = ', fdata[icell], ' fweight = ', fweight[i])

           fout[jb,jl] = fout[jb,jl] + fdata[icell] * fweight[i]
           if jl == 35 and jb == 68:   # Python indices
               print(icell, ' --> jb = ', jb, ' jl = ', jl, ' --> fout= ', fout[jb,jl])
               print('        --> fdata = ', fdata[icell], ' fweight = ', fweight[i])
##           elif jl == 344 and jb == 56:   # Python indices
           elif icell == 4823:   # Python indices
               print(icell, ' --> jb = ', jb, ' jl = ', jl, ' --> fout= ', fout[jb,jl])
               print('        --> fdata = ', fdata[icell], ' fweight = ', fweight[i])

    return fout

# ***************************************************
def remap_2d_to_1d(fmap, fweight, fdst, fdata, nicon):
# ***************************************************
#
#   Generate Mapping coordinates from HD to ICON:
#   HD coordinates: xlon in [0:nlon[, xlat in [0:nlat[
#   In mo_couling_hd, global cell_index already runs from 0 to nlon*nlat -1. However, the 
#   addresses in the weight file are global_id + 1. This may change in future YAC versions.
#
    nbox = len(fmap[0,:])
    fout = np.zeros((nicon))
    print('nbox= ', nbox, ' nicon=', nicon)

    for i in range(nbox):
        jb = int(fmap[0, i])
        jl = int(fmap[1, i])

##        index_py = fsrc[i] - 1    # -1 wegen Python starting at 0
##        jb = int(float(index_py) / float(nlon))
##        jl = index_py - jb * nlon 

        icell = fdst[i] - 2
        fout[icell] = fout[icell] + fdata[jb][jl] * fweight[i]

        if fdata[jb][jl] <= 0:   # Python indices
            print(i, ' --> jb = ', jb, ' jl = ', jl, ' --> fout= ', fout[icell])
            print(icell,' --> fdata = ', fdata[jb][jl], ' fweight = ', fweight[i])

        if jl == 347 and jb == 55:   # Python indices
            print(i, ' --> jb = ', jb, ' jl = ', jl, ' --> fout= ', fout[icell])
            print(icell,' --> fdata = ', fdata[jb][jl], ' fweight = ', fweight[i])

        if fdst[i] >= nicon or fdst[i] < 0:
            print('i=', i, ' ERROR: fdst[i] ', fdst[i], ' outside [0;', nicon,'[')
            print('You may have read the wrong weight file?!')

##        xlon_dst.append(degrees(xlon_icon[fdst[i]-2]))   # -1 wegen Python starting at 0
##        xlat_dst.append(degrees(xlat_icon[fdst[i]-2]))

    return fout

# ***************************************************
def write_remapped_data(ds, DOUT, fout):

##    fmas = ds.variables['cell_sea_land_mask'][:]
##    fmas = xr.zeros_like(ds['cell_sea_land_mask'])
    nbox = fout.size    
    print('ncells = ', nbox)

    # Write NC data
    icon_data = uio.Mapdata()
    icon_data.nlon = nbox ;  icon_data.nlat = 0
    icon_data.vars = ['weights'] ; icon_data.units = ['-'] 
    # 
    # Coordinates
    icon_data.xlat = ds.variables['clat'][:]
    icon_data.xlon = ds.variables['clon'][:]
 
    dntmp = DOUT + os.sep + 'hd_weights_plain' + '.nc' 
    dnw = DOUT + os.sep + 'hd_weights_on_icon' + '.nc' 

##    print('fmas= ',fmas.values)
    print('fout= ',fout)

    ofile = icon_data.openw(dntmp)
##    icon_data.write(fmas.values)
    icon_data.write(fout)
    icon_data.close(ofile)
    call('cdo setgrid,' + DNICON + ' ' + dntmp + ' ' + dnw , shell=True)
    print('Weights written to ' + dnw)
    call('rm ' + dntmp, shell=True)


#***************************************************
if __name__ == "__main__":
    main()
#***************************************************

