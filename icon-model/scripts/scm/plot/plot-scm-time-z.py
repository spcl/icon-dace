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

'''
Purpose
    Read and plot SCM from NetCDF file with Python3
Details
    Time-Vertical plot
Author
    Martin Koehler, DWD
Revision history
    20200701 -- Initial version
'''

from netCDF4 import Dataset, date2num, num2date
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits as mp

# setup

varname  = 'temp'      # temp, u, v, w, clc, rh, tot_qv_dia, tot_qc_dia, tot_qi_dia
#npoint  = 0           # one point (0 is 1st point)
npoint   = -1          # mean of 32 points

#nc_file = '/hpc/uwork/mkoehler/run-icon/scm/SCM_GoAmazon_OP1/scm_GASS_out_PL_20140220T000000Z.nc'
nc_file = '/hpc/uwork/mkoehler/run-icon/scm/SCM_ICON/scm_out_ML_20200101T000000Z.nc'
nc_fid  =  Dataset(nc_file, 'r')

# read file

dates   = num2date(nc_fid.variables['time'][:], nc_fid.variables['time'].units)
vardata = nc_fid.variables[varname][:]
geopot  = nc_fid.variables['geopot'][:]

# time conversion

num_dates = [int(d.strftime('%Y%m%d%H')) for d in dates]    # YYYYMMDDHH e.g. 2019081516
hours     = [(d.day-dates[0].day)*24+(d.hour-dates[0].hour)+(d.minute/60)+(d.second/3600) for d in dates]

# plot

title  = "ICON SCM forecast    %s" % (nc_fid.variables[varname].standard_name)
xtitle = "hours since %s" % (num_dates[0])
ytitle = "height [m]"

fig, ax = plt.subplots(figsize=(8, 5))
fig.subplots_adjust(bottom=0.15, right=1.0)

ax.set_title (title)
ax.set_xlabel(xtitle)
ax.set_ylabel(ytitle)

grav  = 9.80665
if npoint >= 0:
  dataz = vardata[:,:,npoint]
  datay = geopot [0,:,npoint] / grav
else:
  dataz = vardata[:,:,:].mean(axis=2)
  datay = geopot [0,:,:].mean(axis=1) / grav

[X, Y] = np.meshgrid(hours, datay)        # creating 2-D grid
Z      = dataz
Z      = np.swapaxes(Z,0,1)

#cmap = 'jet'                             # color map options
cmap1 = plt.cm.get_cmap('Spectral')
cmap  = cmap1.reversed()

CS = ax.contourf(X, Y, Z, 30, cmap=cmap, extend='both')   # plot filled contour plot

# color bar

cbar = fig.colorbar(CS)
cbar.ax.set_ylabel('['+nc_fid.variables[varname].units+']')
 
# right axis with model level locations

ax2 = ax.twinx()                          # second axes
ax2.set_yticks(datay)
ax2.tick_params(axis='y', labelright=False)

# screen, png, pdf

fig.savefig('time_z_'+varname+'.png')     # png
#fig.savefig('var_z'+varname+'.pdf')      # pdf
plt.show()                                # to screen
