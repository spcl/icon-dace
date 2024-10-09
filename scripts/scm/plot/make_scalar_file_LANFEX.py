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

from sys import argv
import numpy as np
from matplotlib import pyplot as plt
import netCDF4
import os
import diagnostics as diag

# function to deaccumulate fluxes
def deaccumulate(aX):
    X = np.zeros_like(aX)
    X[1:] =  aX[1:] - aX[:-1]
    return X/(60.0*time[1])
## plotting parameters
nc = 30 # cell number

## read in file
expname = argv[1]
sname = '../'+expname+'/scalars_'+expname+'_ICON_SCM_DWD_DOM01_ML_0001.nc'
print("scalar nc file : ",sname)
sfile = netCDF4.Dataset(sname)
N_t = sfile.variables['time'][:].shape[0]


pname = '../'+expname+'/profs_'+expname+'_ICON_SCM_DWD_DOM01_ML_0001.nc'
pfile = netCDF4.Dataset(pname)
z_mc = sfile.variables['z_mc'][:,nc]


#oufile
oname = 'scalars_'+expname+'_ICON_SCL_DWD_v1.nc'
os.system('mkdir '+expname)

ofile = netCDF4.Dataset(expname+'/'+oname,'w')
ofile.createDimension('time',size=0)
# VARS #######

#time (seconds)
time = sfile.variables['time']
time_o = ofile.createVariable('time',np.float64,('time'))
n_hours = 19.0
dt_in_hours = (n_hours)/(len(time[:])-1)
time_o[:] = [3600.0 *(dt_in_hours) * float(i) for i in range(len(time[:]))] # 60.0  *time[:]
time_o.standard_name = 'time'
time_o.units = 'seconds since 2014-11-24 17:00:00'

# create array of datetimes
import datetime
init_date = datetime.datetime(year=2014,month=11,day=24,hour=17)
date_list = []
for t in range(len(time_o)):
    timedelta = datetime.timedelta(seconds = time_o[t] )
    date_list.append( init_date + timedelta )
    print (date_list[t])

#lwdn (W/m^2)
lwdn_o = ofile.createVariable('lwdn',np.float64,('time'))
lwdn_o[:] = deaccumulate( sfile.variables['accthd_s'][:,nc] ) # Surface down thermal radiationacc. since model start
lwdn_o.standard_name = 'downwelling longwave radiation'
lwdn_o.units= 'W/m**2'

#lwup (W/m**2)
lwup_o = ofile.createVariable('lwup',np.float64,('time'))
lwup_o[:] = deaccumulate( sfile.variables['accthu_s'][:,nc] )
lwup_o.standard_name = 'upwelling longwave radiation'
lwup_o.units= 'W/m**2'

#swdn
swdn_o = ofile.createVariable('swdn',np.float64,('time'))
swdndir = deaccumulate( sfile.variables['accsodird_s'][:,nc] ) # Surface down solar direct rad.acc. since model start
swdndif = deaccumulate( sfile.variables['accsodifd_s'][:,nc] ) # Surface down solar diff. rad. acc. since model start
swdn_o[:] = swdndir[:] + swdndif[:]
swdn_o.standard_name = 'downwelling shortwave radiation'
swdn_o.units= 'W/m**2'

#swup
swup_o = ofile.createVariable('swup',np.float64,('time'))
swup_o[:] = sfile.variables['sou_s'][:,nc]
swup_o.standard_name = 'upwelling shortwave radiation'
swup_o.units= 'W/m**2'

#tstar
tstar_o = ofile.createVariable('tstar',np.float64,('time'))
tstar_o[:] = sfile.variables['t_g'][:,nc]
tstar_o.standard_name = 'surface temperature'
tstar_o.units = 'K'

#shf (W/m**2)
shf_o = ofile.createVariable('shf',np.float64,('time'))
shf_o[:] =  -deaccumulate( sfile.variables['accshfl_s'][:,nc] ) # surface sensible heat flux acc. since model start
shf_o.standard_name = 'Surface sensible heat flux'
shf_o.units = 'W/m**2'

#lhf (W/m**2)
lhf_o = ofile.createVariable('lhf',np.float64,('time'))
lhf_o[:] =  -deaccumulate( sfile.variables['acclhfl_s'][:,nc] ) # surface sensible heat flux acc. since model start
lhf_o.standard_name = 'Surface latend heat flux'
lhf_o.units = 'W/m**2'

# u10m (m/s)
u10m_o = ofile.createVariable('u10m',np.float64,('time'))
u_10m = sfile.variables['u_10m'][:,0,nc]
u10m_o[:] = u_10m[:]
u10m_o.standard_name = '10m zonal wind velocity'
u10m_o.unts = 'm/s'

# v10m (m/s)
v10m_o = ofile.createVariable('v10m',np.float64,('time'))
v_10m = sfile.variables['v_10m'][:,0,nc]
v10m_o[:] = v_10m[:]
v10m_o.standard_name = '10m meridional wind velocity'
v10m_o.unts = 'm/s'

#ustar (m/s)
ustar_o = ofile.createVariable('ustar',np.float64,('time'))
tvm = sfile.variables['tvm'][:,nc] # turbulent transfer velocity for momentum tvm = C_M U_A
ustar_o[:] = np.sqrt( tvm[:] * np.sqrt( u_10m[:]**2 + v_10m[:]**2 )  )
ustar_o.standard_name = 'Surface friction velocity'
ustar_o.units = 'm/s'

#blh (m)
blh_o = ofile.createVariable('blh',np.float64,('time'))
theta_v = sfile.variables['theta_v'][:,:,nc]
for t in range(time[:].shape[0]):
    blh_o[t] = diag.get_zi(theta_v[t,:],z_mc[:])
blh_o.standard_name = 'boundary layer height'
blh_o.units = 'm'

# zct (m) cloud top height
zct_o = ofile.createVariable('zct',np.float64,('time'))
clc = sfile.variables['clc'][:,:,nc]
for t in range(time[:].shape[0]):
    zct_o[t] = diag.cloud_top(clc[t,:],z_mc[:])
zct_o.standard_name = 'Cloud top height'
zct_o.units='m'

#lwp (kg/m**2)
lwp_o = ofile.createVariable('lwp',np.float64,('time'))
tqc_dia = sfile.variables['tqc_dia'][:,nc]
tqr = sfile.variables['tqr'][:,nc]
lwp_o[:] = tqc_dia[:] + tqr[:]
lwp_o.standard_name = 'liquid water path'
lwp_o.units = 'kg/m**2' #sfile.variables['tqc_dia'].units

#tca clould cover
tca_o = ofile.createVariable('tca',np.float64,('time'))
tca_o[:] = sfile.variables['clct'][:,nc]/100.0
tca_o.standard_name = 'total cloud cover'
tca_o.units = '0-1'

#cldsed (kg/m^2s)
cldsed_o = ofile.createVariable('cldsed',np.float64,('time'))
cldsed_o[:] = sfile.variables['cloud_gsp_rate'][:,nc]
cldsed_o.standard_name = 'Cloud droplet sedimentation rate onto the surface'
cldsed_o.units = 'kg/m**2s' #sfile.variables['cloud_gsp_rate'].units

#rainsed (kg/m^2s)
rainsed_o = ofile.createVariable('rainsed',np.float64,('time'))
rainsed_o[:] = sfile.variables['rain_gsp_rate'][:,nc]
rainsed_o.standard_name = 'Rain droplet sedimentation rate onto the surface'
rainsed_o.units = 'kg/m**2s' #sfile.variables['rain_gsp_rate'].units


# optional diagnistics
##qcfl_s (kg/m^2s)
#qcfl_s_o = ofile.createVariable('qcfl_s',np.float64,('time'))
#qcfl_s_o[:] = sfile.variables['qcfl_s'][:,nc]
#qcfl_s_o.standard_name =sfile.variables['qcfl_s'].standard_name
#qcfl_s_o.units = sfile.variables['qcfl_s'].units
##qhfl_s
#qhfl_s_o = ofile.createVariable('qhfl_s',np.float64,('time'))
#qhfl_s_o[:] = sfile.variables['qhfl_s'][:,nc]
#qhfl_s_o.standard_name = sfile.variables['qhfl_s'].standard_name
#qhfl_s_o.units = sfile.variables['qhfl_s'].units
##umfl_s
#umfl_s_o = ofile.createVariable('umfl_s',np.float64,('time'))
#umfl_s_o[:] = sfile.variables['umfl_s'][:,nc]
#umfl_s_o.standard_name = sfile.variables['umfl_s'].standard_name
#umfl_s_o.units = sfile.variables['umfl_s'].units
##qhfl_s
#vmfl_s_o = ofile.createVariable('vmfl_s',np.float64,('time'))
#vmfl_s_o[:] = sfile.variables['vmfl_s'][:,nc]
#vmfl_s_o.standard_name = sfile.variables['vmfl_s'].standard_name
#vmfl_s_o.units = sfile.variables['vmfl_s'].units

##tscrn (K)
#tscrn_o = ofile.createVariable('tscrn',np.float64,('time'))
#tscrn = sfile.variables['t_2m'][:,nc]
#tscrn_o.standard_name = 'Screen level (2.0m) temperature'
#tscrn_o.units = 'm'

#t_2m (K)
t_2m_o = ofile.createVariable('t_2m',np.float64,('time'))
t_2m_o[:] = sfile.variables['t_2m'][:,0,nc]
t_2m_o.standard_name = 'Temperature at 2 meter'
t_2m_o.units = 'K'

#rh_2m (%)
rh_2m_o = ofile.createVariable('rh_2m',np.float64,('time'))
rh_2m_o[:] = sfile.variables['rh_2m'][:,0,nc]
rh_2m_o.standard_name = 'rh_2m' #sfile.variables['rh_2m'].standard_name
rh_2m_o.units = '%' # sfile.variables['rh_2m'].units

#vis (m)
vis_o = ofile.createVariable('vis',np.float64,('time'))
vis_o.standard_name = 'surface visibility'
vis_o.units = 'm'

lat = 52.1/180.0*np.pi
lon = -0.42/180.0*np.pi

import zenith
from vis import calculate_visibility_wrf
tot_qv_dia = sfile.variables['tot_qv_dia'][:,:,nc]
tot_qc_dia = sfile.variables['tot_qc_dia'][:,:,nc]
tot_qi_dia = sfile.variables['tot_qi_dia'][:,:,nc]
qr = sfile.variables['qr'][:,:,nc]
qs = sfile.variables['qs'][:,:,nc]
temp = sfile.variables['temp'][:,:,nc]
pres = sfile.variables['pres'][:,:,nc]
rho = sfile.variables['rho'][:,:,nc]
u = sfile.variables['u'][:,:,nc]
v = sfile.variables['v'][:,:,nc]
vis = np.zeros_like(u_10m)

for t,date in enumerate(date_list):
    year = date.year
    month = date.month
    day = date.day
    hour = date.hour
    day_in_year = date.timetuple().tm_year
    cos_mu0 = zenith.pre_radiation_nwp(hour,year,day_in_year,lat,lon)
    vis[t] = calculate_visibility_wrf(tot_qv_dia[t,:],tot_qc_dia[t,:],tot_qi_dia[t,:],qr[t,:],qs[t,:],
            u[t,:],v[t,:],temp[t,:],pres[t,:],rho[t,:],cos_mu0,tot_qv_dia.shape[1]                      )
vis_o[:] = vis[:]

 #= ofile.createVariable('',np.float64,('time'))
 #= ofile.createVariable('',np.float64,('time'))

### close file
ofile.close()



