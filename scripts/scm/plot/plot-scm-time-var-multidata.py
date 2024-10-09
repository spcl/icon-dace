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
from scipy.signal import savgol_filter
import math

# setup

varname      = 'lhfl_s'      # temp, t_2m, qv_2m, u_10m, v_10m, shfl_s, lhfl_s, clct
level        = 1           # 1000hPa, 90 (model level), 1 (e.g. t_2m)
#npoints     = (1,)        # one point (0 is 1st point)
#npoints     = range(32)   # plot all 32 points
npoints      = (-1,)       # mean of 32 points
varname_axis = 'Latent heat flux'      # TCC, T2m, $\mathit{\\theta_v}$
                           # 'Latent heat flux',

metgrm_nvar  = 21          # 20:shfl_s, 21:lhfl_s, 24:t_2m, 26:u_10m, 27:v_10m, 5:q_s, 42:clct


# SCM output:
#nc_file = '/hpc/uwork/mkoehler/run-icon/scm/SCM_GoAmazon_OP1/scm_GASS_out_PL_20140220T000000Z.nc'
nc_file = '/hpc/uwork/mkoehler/run-icon/scm/SCM_ICON_fixcpcv/scm_out_ML_20200712T000000Z.nc'
nc_fid  =  Dataset(nc_file, 'r')

# LES output:
nc4_file = '/hpc/uwork/mkoehler/run-icon/scm/LES_ICON/scm_out_sfc_ML_LES_ICON_2020071200_mean.nc'
nc4_fid  =  Dataset(nc4_file, 'r')

# CRM output:
nc5_file = '/hpc/uwork/mkoehler/run-icon/scm/CRM_ICON_fixcpcv/scm_out_sfc_ML_CRM_ICON_fixcpcv_2020071200_mean.nc'
nc5_fid  =  Dataset(nc5_file, 'r')

# ICON global meteogram output:
nc2_file = '/hpc/uwork/mkoehler/run-icon/experiments/exp_002_2020071200/METEOGRAM_patch001.nc'
nc2_fid  =  Dataset(nc2_file, 'r')
metgrm_2d3d = 'sfcvalues'     # 'sfcvalues' or 'values' for 3D variables
metgrm_np   = 2               # 2:Lindenberg
metgrm_np6  = 1               # 1,2: two points near Lindeneberg
metgrm_dtime= 360             # ICON model time step
metgrm_dtime6= 1.5            # LES LAM time step

# ICON global point output: (cdo remapnn,lon=14.11/lat=52.21 in out)
nc3_file = '/hpc/uwork/mkoehler/run-icon/experiments/exp_002_2020071200/out-Lindenberg.nc'
nc3_fid  =  Dataset(nc3_file, 'r')

# LES LAM meteogram output:
nc6_file = '/hpc/uwork/mkoehler/run-icon/scm/LES_LAM_ICON/METEOGRAM_patch003.nc'
nc6_fid  =  Dataset(nc6_file, 'r')


# read file

dates    = num2date(nc_fid.variables['time'][:] , nc_fid.variables['time'].units)
dates3   = num2date(nc3_fid.variables['time'][:], nc3_fid.variables['time'].units)
dates4   = num2date(nc4_fid.variables['time'][:], nc4_fid.variables['time'].units)

vardata  = nc_fid.variables[varname][:]
#vardata3 = nc3_fid.variables[varname][:]
vardata4 = nc4_fid.variables[varname][:]
vardata5 = nc5_fid.variables[varname][:]
vardata2 = nc2_fid.variables[metgrm_2d3d][:,metgrm_nvar-1,metgrm_np-1]
vardata2 = vardata2.reshape(vardata2.shape[0],1,1)
#vardata6 = nc6_fid.variables[metgrm_2d3d][:,metgrm_nvar-1,metgrm_np6-1]
vardata6 = nc6_fid.variables[metgrm_2d3d][:,metgrm_nvar-1,:].mean(axis=1)
vardata6 = vardata6.reshape(vardata6.shape[0],1,1)

if varname == 'lhfl_s' or varname == 'shfl_s' or varname == 'clct':
  vardata  = vardata.reshape(vardata.shape[0],1,vardata.shape[1])
if varname == 'lhfl_s' or varname == 'shfl_s':
  vardata [0,0,:] = vardata [1,0,:]   # first flux point is diagnostic and wrong
  vardata4[0,0,0] = vardata4[1,0,0]
  vardata5[0,0,0] = vardata5[1,0,0]
vardata4 = vardata4.reshape(vardata4.shape[0],1)
vardata5 = vardata5.reshape(vardata5.shape[0],1)

if varname == 'clct':
  vardata2 = vardata2 * 100.0
  vardata6 = vardata6 * 100.0

print('variable shapes', vardata.shape, vardata2.shape, vardata4.shape, vardata5.shape, vardata6.shape)

#-- time conversion

num_dates  = [int(d.strftime('%Y%m%d%H')) for d in dates]    # YYYYMMDDHH e.g. 2019081516
hours      = [(d.day-dates[0].day)*24 +(d.hour-dates[0].hour) +(d.minute/60)+(d.second/3600) for d in dates]
#hours3     = [(d.day-dates3[0].day)*24+(d.hour-dates3[0].hour)+(d.minute/60)+(d.second/3600) for d in dates3]
hours4     = [(d.day-dates4[0].day)*24 +(d.hour-dates4[0].hour) +(d.minute/60)+(d.second/3600) for d in dates4]
hours2     = nc2_fid.variables['time_step'][:] * metgrm_dtime  / 3600.0
hours6     = nc6_fid.variables['time_step'][:] * metgrm_dtime6 / 3600.0

# smoothing:
#   savitztky/golay filter: least squares to regress a small window of your data onto a polynomial
#   savgol_filter(x, window_length, polyorder, ...)

tsmooth = 1.5         # smoothing in [h]
print( 'smoothing with window width:', math.ceil(tsmooth/(hours2[1]-hours2[0])/2)*2+1 )

vardata   = savgol_filter(vardata , math.ceil(tsmooth/(hours [1]-hours [0])/2)*2+1, 1, axis=0)
vardata2  = savgol_filter(vardata2, math.ceil(tsmooth/(hours2[1]-hours2[0])/2)*2+1, 1, axis=0)
#vardata3  = savgol_filter(vardata3, math.ceil(tsmooth/(hours3[1]-hours3[0])/2)*2+1, 1, axis=0)
vardata4  = savgol_filter(vardata4, math.ceil(tsmooth/(hours4[1]-hours4[0])/2)*2+1, 1, axis=0)
vardata5  = savgol_filter(vardata5, math.ceil(tsmooth/(hours [1]-hours [0])/2)*2+1, 1, axis=0)
vardata6  = savgol_filter(vardata6, math.ceil(tsmooth/(hours6[1]-hours6[0])/2)*2+1, 1, axis=0)

#-- unit formatting:  convert [W m-2] to [$W m^{-2}$] etc
unit = nc_fid.variables[varname].units
for r in (("-", "^{-"), ("1", "1}"), ("2", "2}")):
  unit = unit.replace(*r)
unit = '[$' + unit + '$]'
print('old/new unit format: ', unit)

#-- plot

title  = "ICON SCM forecast"
xtitle = "Simulation time [$h$]"    # since %s [h]" % (num_dates[0])
#ytitle = "%s (%s)" % (nc_fid.variables[varname].standard_name,\
#                      nc_fid.variables[varname].units)
#ytitle = "%s [%s]" % (varname_axis,nc_fid.variables[varname].units)
ytitle = varname_axis + '  ' + unit

dt_ticks = 3.0
ticks  = np.arange(0.,max(hours)+0.1,dt_ticks)

fig, ax = plt.subplots(figsize=(7, 5))
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(left=0.15)

#ax.set_title (title)
ax.set_title('(d)',fontsize=16, fontweight='bold',loc='right')
ax.set_xlabel(xtitle, fontsize=18)
ax.set_ylabel(ytitle, fontsize=18)
ax.set_xticks(ticks)
ax.set_xlim(0,24)
#ax.set_ylim(281,297)   # t_2m
#ax.set_ylim(0.004,0.008)   # qv_2m
#ax.set_ylim(-4,1)   # v_10m
#ax.set_ylim(-300,10)   # lhfl_s
ax.set_ylim(-400,50)   # shfl_s

plt.setp(ax.get_xticklabels(), fontsize=14) # xtick labels fontsize
plt.setp(ax.get_yticklabels(), fontsize=14) # ytick labels fontsize
plt.rc('legend', fontsize=13)               # legend fontsize

# plot SCM line:

for np in npoints:
  if np >= 0:
    data = vardata[:,level-1,np-1]
  else :
    data = vardata[:,level-1,:].mean(axis=1)
  ax.plot(hours, data, linewidth=2.0, label='SCM')

# plot LES line:

ax.plot(hours4, vardata4[:,level-1], linewidth=2.0, label='LES-PER')

# plot CRM line:

ax.plot(hours,  vardata5[:,level-1], linewidth=2.0, label='CRM-PER')

# plot ICON global meteogram line:

ax.plot(hours2, vardata2[:,level-1], linewidth=2.0, label='global')

# plot LES LAM meteogram line:

ax.plot(hours6, vardata6[:,level-1], linewidth=2.0, label='LES-LAM')

# plot ICON global line:

#ax.plot(hours3, vardata3[:,level-1,0,0], linewidth=2.0, label='ICON global')

# label for each line

ax.legend()


fig.savefig('time-var_'+varname+'.png')   # png
fig.savefig('time-var_'+varname+'.pdf')   # pdf
plt.show()                                # to screen
