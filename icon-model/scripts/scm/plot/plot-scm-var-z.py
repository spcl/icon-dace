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
    Variable-Vertical plot
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

varname  = 'temp'                # temp, u, tot_qv_dia
#steps   = (0,)                  # time step (0 is 1st step)
steps    = (0,3,6,9,12,15,18,21,24)   # time step (0 is 1st step)

#npoints = (0,)                  # one point (0 is 1st point)
#npoints = range(32)             # plot all 32 points
npoints  = (-1,)                 # mean of 32 points

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

title  = "ICON SCM forecast"
xtitle = "%s [%s]" % (nc_fid.variables[varname].standard_name,\
                      nc_fid.variables[varname].units)
ytitle = "height [m]"

fig, ax = plt.subplots(figsize=(5, 7))
fig.subplots_adjust(bottom=0.15, left=0.2)

ax.set_title (title)
ax.set_xlabel(xtitle)
ax.set_ylabel(ytitle)

grav  = 9.80665

for ns in steps:
  for np in npoints:
    if np >= 0:
      datax = vardata[ns,:,np-1]
      datay = geopot [ns,:,np-1] / grav
    else:
      datax = vardata[ns,:,:].mean(axis=1)
      datay = geopot [ns,:,:].mean(axis=1) / grav
    ax.plot(datax,datay, label=str(hours[ns])+" h")   # (" + str(num_dates[ns]) + ")")
    ax.margins(y=0)

# right axis with model level locations

ax2 = ax.twinx()                          # second axes
ax2.set_yticks(datay)
ax2.tick_params(axis='y', labelright=False)  

# label for each line

ax.legend()

fig.savefig('var_z_'+varname+'.png')      # png
#fig.savefig('var_z'+varname+'.pdf')      # pdf
plt.show()                                # to screen
