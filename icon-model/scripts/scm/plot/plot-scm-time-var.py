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
    Time-Variable plot
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

varname  = 'temp'
level    = 90          # 1000hPa
#npoints = (0,)        # one point (0 is 1st point)
#npoints = range(32)   # plot all 32 points
npoints  = (-1,)       # mean of 32 points

#nc_file = '/hpc/uwork/mkoehler/run-icon/scm/SCM_GoAmazon_OP1/scm_GASS_out_PL_20140220T000000Z.nc'
nc_file = '/hpc/uwork/mkoehler/run-icon/scm/SCM_ICON/scm_out_ML_20200712T000000Z.nc'
nc_fid  =  Dataset(nc_file, 'r')

# read file

dates   = num2date(nc_fid.variables['time'][:], nc_fid.variables['time'].units)
vardata = nc_fid.variables[varname][:]

# time conversion

num_dates = [int(d.strftime('%Y%m%d%H')) for d in dates]    # YYYYMMDDHH e.g. 2019081516
hours     = [(d.day-dates[0].day)*24+(d.hour-dates[0].hour)+(d.minute/60)+(d.second/3600) for d in dates]

# plot

title  = "ICON SCM forecast"
xtitle = "hours since %s" % (num_dates[0])
ytitle = "%s (%s)" % (nc_fid.variables[varname].standard_name,\
                      nc_fid.variables[varname].units)
ticks  = np.arange(0.,max(hours)+0.1,24.)

fig, ax = plt.subplots(figsize=(7, 5))
fig.subplots_adjust(bottom=0.15)

ax.set_title (title)
ax.set_xlabel(xtitle)
ax.set_ylabel(ytitle)

for np in npoints:
  if np >= 0:
    data = vardata[:,level-1,np-1]
  else :
    data = vardata[:,level-1,:].mean(axis=1)
  ax.plot(hours, data)
  ax.set_xticks(ticks)

fig.savefig('time-var_'+varname+'.png')   # png
#fig.savefig('time-var_'+varname+'.pdf')  # pdf
plt.show()                                # to screen
