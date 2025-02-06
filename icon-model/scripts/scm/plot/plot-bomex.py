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

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


infiles=[
         '../../../data/scm/SCM_BOMEX/scm_out_ML_19690626T000000Z.nc',
         '../../../data/scm/CRM_BOMEX/scm_out_ML_19690626T000000Z.nc',
         '../../../data/scm/LEM_BOMEX/scm_out_ML_19690626T000000Z.nc'
         ]

labels=[
        'SCM',
        'CRM',
        'LEM'
]

vars=['qv','temp','u','v']
vars_surf=['shfl_s','lhfl_s','umfl_s','vmfl_s']
times=[4*3600,5*3600] #seconds
ds = [ [] for _ in range(len(infiles)) ]

for ifl in range(0,len(infiles)):
 ds[ifl]=xr.open_dataset(infiles[ifl])
 ds[ifl]['time_s']=np.round(ds[ifl].time%1*3600.*24)


for it in range(0,len(times)):
  fig, ax = plt.subplots(nrows=len(vars), ncols=1,figsize=(6, 4*len(vars)))
  for ifl in range(0,len(infiles)):
    indt=np.argmin(np.abs(times[it]-ds[ifl]['time_s'].values))
    for ivar in range(0,len(vars)):
     ax[ivar].plot(ds[ifl][vars[ivar]][indt,:,:].mean(axis=1),ds[ifl].z_mc.mean(axis=1),label=labels[ifl])
     ax[ivar].set_xlabel(vars[ivar])
     ax[ivar].set_ylabel('z[m]')
     if(ivar==0):
      ax[ivar].legend()
  fig.savefig('bomex-prof-t'+str(times[it])+'.png')   
  plt.close()

fig, ax = plt.subplots(nrows=len(vars), ncols=1,figsize=(6, 4*len(vars_surf)))
for ifl in range(0,len(infiles)):
  for ivar in range(0,len(vars_surf)):
   ax[ivar].plot(ds[ifl]['time_s'].values,ds[ifl][vars_surf[ivar]][:,:].mean(axis=1).values,label=labels[ifl])
   ax[ivar].set_ylabel(vars_surf[ivar])
   ax[ivar].set_xlabel('time[s]')
   if(ivar==0):
    ax[ivar].legend()
fig.savefig('bomex-evol.png')   
plt.close()

for ifl in range(0,len(infiles)):
 ds[ifl].close()
