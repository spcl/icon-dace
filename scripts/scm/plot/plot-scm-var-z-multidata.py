# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

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

varname      = 'u'        # theta_v, temp, u, qv, tot_qv_dia
metgrm_nvar  = 6                  # 5:theta_v, 6:u, 7:v, 12:qv, 18:tot_qv_dia ICON_global
metgrm_nvar2 = 6                  # 5:theta_v, 6:u, 7:v, 10:qv                LEM-LAM
varname_axis = 'zonal wind speed' # '$\mathit{wind speed}$'  q_v, wind speed, $\mathit{\\theta_v}$

#steps   = (0,)                  # time step (0 is 1st step)
#steps    = (0,3,6,9,12,15,18,21,24)   # time step (0 is 1st step)
#steps    = (0,24,48,72,96)   # time step (0 is 1st step)
#steps_mtg= (0,12,24,36,48) 
#steps    = (0,4,8,12)   # time step (0 is 1st step)
#steps_mtg= (0,2,4,6)

# 15UTC for 3 days
# hours_plot = np.array([15,15+24,15+48])
# steps    = hours_plot * 4   # time step (0 is 1st step)
# steps_mtg= hours_plot * 2
# steps4   = hours_plot / 3

# 15UTC for first day
hours_plot = np.array([15])
steps    = hours_plot * 4   # time step (0 is 1st step)
steps_mtg= hours_plot * 2
steps4   = hours_plot / 3
steps6   = np.array([2880])

print(steps)
print(steps_mtg)
print(steps4)


#npoints = (0,)                  # one point (0 is 1st point)
#npoints = range(32)             # plot all 32 points
npoints  = (-1,)                 # mean of 32 points

# SCM output:
#nc_file = '/hpc/uwork/mkoehler/run-icon/scm/SCM_GoAmazon_OP1/scm_GASS_out_PL_20140220T000000Z.nc'
nc_file = '/hpc/uwork/mkoehler/run-icon/scm/SCM_ICON_fixcpcv/scm_out_ML_20200712T000000Z.nc'
nc_fid  =  Dataset(nc_file, 'r')

# LES output:
nc4_file = '/hpc/uwork/mkoehler/run-icon/scm/LES_ICON/scm_out_ml_ML_LES_ICON_2020071200_mean.nc'
nc4_fid  =  Dataset(nc4_file, 'r')

# CRM output:
nc5_file = '/hpc/uwork/mkoehler/run-icon/scm/CRM_ICON_fixcpcv/scm_out_ml_ML_CRM_ICON_fixcpcv_2020071200_mean.nc'
nc5_fid  =  Dataset(nc5_file, 'r')

# ICON global meteogram output:
nc2_file = '/hpc/uwork/mkoehler/run-icon/experiments/exp_002_2020071200/METEOGRAM_patch001.nc'
nc2_fid  =  Dataset(nc2_file, 'r')
metgrm_2d3d = 'values'        # 'sfcvalues' or 'values' for 3D variables
metgrm_np   = 2               # 2:Lindenberg
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
vardata2 = nc2_fid.variables[metgrm_2d3d][:,:90,metgrm_nvar-1,metgrm_np-1]
vardata6 = nc6_fid.variables[metgrm_2d3d][:,:90,metgrm_nvar2-1,:].mean(axis=2)
#vardata3 = nc3_fid.variables[varname][:]
vardata4 = nc4_fid.variables[varname][:]
vardata5 = nc5_fid.variables[varname][:]

# wind speed
# vardata__v  = nc_fid.variables['v'][:]
# vardata2_v = nc2_fid.variables[metgrm_2d3d][:,:90,7-1,metgrm_np-1]
# vardata6_v = nc6_fid.variables[metgrm_2d3d][:,:90,7-1,:].mean(axis=2)
# vardata4_v = nc4_fid.variables['v'][:]
# vardata5_v = nc5_fid.variables['v'][:]
# vardata  = ( vardata  ** 2 + vardata__v ** 2 ) ** 0.5
# vardata2 = ( vardata2 ** 2 + vardata2_v ** 2 ) ** 0.5
# vardata6 = ( vardata6 ** 2 + vardata6_v ** 2 ) ** 0.5
# vardata4 = ( vardata4 ** 2 + vardata4_v ** 2 ) ** 0.5
# vardata5 = ( vardata5 ** 2 + vardata5_v ** 2 ) ** 0.5


geopot   = nc_fid.variables ['geopot'][:]
heights6 = nc6_fid.variables['heights'][:,7,:].mean(axis=1)
heights6 = heights6 - min(heights6) + 10.0
#print('LES LAM heights', heights6.shape, heights6)

print('variable shapes', vardata.shape, vardata2.shape, vardata4.shape, vardata5.shape, vardata6.shape) #vardata3.shape, 

#-- time conversion

num_dates  = [int(d.strftime('%Y%m%d%H')) for d in dates]    # YYYYMMDDHH e.g. 2019081516
hours      = [(d.day-dates[0].day)*24 +(d.hour-dates[0].hour) +(d.minute/60)+(d.second/3600) for d in dates]
hours3     = [(d.day-dates3[0].day)*24+(d.hour-dates3[0].hour)+(d.minute/60)+(d.second/3600) for d in dates3]
hours4     = [(d.day-dates4[0].day)*24 +(d.hour-dates4[0].hour) +(d.minute/60)+(d.second/3600) for d in dates4]
hours2     = nc2_fid.variables['time_step'][:] * metgrm_dtime  / 3600.0
hours6     = nc6_fid.variables['time_step'][:] * metgrm_dtime6 / 3600.0

#-- unit formatting:  convert [W m-2] to [$W m^{-2}$] etc

unit = nc_fid.variables[varname].units
for r in (("-", "^{-"), ("1", "1}"), ("2", "2}")):
  unit = unit.replace(*r)
unit = '[$' + unit + '$]'

#-- plot

title  = "ICON SCM forecast (dot: ICON global, dash: LES)"
#xtitle = "%s [%s]" % (nc_fid.variables[varname].standard_name,\
#                      nc_fid.variables[varname].units)
#xtitle = "%s [%s]" % (varname_axis,nc_fid.variables[varname].units)
xtitle = varname_axis + '  ' + unit

ytitle = "height [$m$"
ytitle_right = "model levels"

fig, ax = plt.subplots(figsize=(6, 6))
fig.subplots_adjust(bottom=0.15, left=0.2)

#ax.set_title (title).set_fontsize(14)
ax.set_title('(c)',fontsize=16, fontweight='bold',loc='right')
ax.set_xlabel(xtitle).set_fontsize(18)
ax.set_ylabel(ytitle).set_fontsize(18)
#ax.set_xlim([285,300])  # theta_v
#ax.set_xlim([290,304])  # theta_v
#ax.set_xlim([0.0,0.007])  # qv
#ax.set_xlim([0,10])     # u
#ax.set_xlim([-10,2])   # v
ax.set_xlim([0,10])    # speed
ax.set_ylim([0,4000])

plt.setp(ax.get_xticklabels(), fontsize=14) # xtick labels fontsize
plt.setp(ax.get_yticklabels(), fontsize=14) # ytick labels fontsize
plt.rc('legend', fontsize=13)               # legend fontsize

grav  = 9.80665

colors = ["mediumblue", "forestgreen", "red", "gold", "purple"]

nc=0
for ns in steps:
  for np in npoints:
    if np >= 0:
      datax = vardata[ns,:,np-1]
      datay = geopot [ns,:,np-1] / grav
    else:
      datax = vardata[ns,:,:].mean(axis=1)
      datay = geopot [ns,:,:].mean(axis=1) / grav
    ax.plot(datax,datay, linewidth=1.5, 
     #color=colors[nc], label=str(hours[ns])+" h")   # (" + str(num_dates[ns]) + ")")
      label= 'SCM')   # (' + str(num_dates[ns]) + ')')
    ax.margins(y=0)
    nc=nc+1

# right axis with model level locations

ax2 = ax.twinx()                          # second axes
ax2.set_yticks(datay)
ax2.tick_params(axis='y', labelright=False)  
#ax2.tick_params(axis='y')  
ax2.set_ylabel(ytitle_right).set_fontsize(14)

# plot LES lines

nc=0
for ns in steps4:
  ns = int(ns)
  datax = vardata4[ns,:,0,0]
  datay = geopot  [ns,:,0] / grav
  ax.plot(datax, datay, linewidth=2.0, 
   #color=colors[nc]), linestyle="dashed"
    label= 'LES-PER')   # (' + str(dates4[ns]) + ')')
  nc=nc+1

# plot CRM lines

nc=0
for ns in steps4:
  ns = int(ns)
  datax = vardata5[ns,:,0,0]
  datay = geopot  [ns,:,0] / grav
  ax.plot(datax, datay, linewidth=2.0, 
   #color=colors[nc]),, linestyle="dashed" linestyle="dashed"
    label= 'CRM-PER')   # (' + str(dates4[ns]) + ')')
  nc=nc+1

# plot ICON global meteogram lines

nc=0
for ns in steps_mtg:
  print(vardata2.shape,geopot.shape)
  datax = vardata2[ns,:]
  datay = geopot  [ns,:,0] / grav
  ax.plot(datax, datay, linewidth=2.0, 
   #color=colors[nc]), linestyle="dotted"
    label= 'global')  # (' + str(hours2[ns]) + ')')
  nc=nc+1

# plot LES LAM meteogram lines

nc=0
for ns in steps6:
  print(vardata6.shape,geopot.shape)
  datax = vardata6[ns,:]
  datay = heights6
  ax.plot(datax[:datax.shape[0]-1], datay[:datax.shape[0]-1], linewidth=2.0, 
   #color=colors[nc]), linestyle="dotted"
    label= 'LES-LAM')   # (' + str(hours6[ns]) + ')')
  nc=nc+1

# label for each line

ax.legend()

fig.savefig('var_z_'+varname+'.png')      # png
fig.savefig('var_z_'+varname+'.pdf')      # pdf
plt.show()                                # to screen
