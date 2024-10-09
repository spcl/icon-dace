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



num_lev = 10

def set_var(oname,Iname,levtyp='full',name='',unit=''):
    if levtyp=='full':
        ovar = ofile.createVariable(oname,np.float64,('time','zf'))
    if levtyp=='half':
        ovar = ofile.createVariable(oname,np.float64,('time','zh'))
    if isinstance(Iname, str): #basestring):  
        Ivar  =  pfile.variables[Iname]
    else:
        Ivar = Iname
    if name == '':
        ovar.standard_name = Ivar.standard_name
    else:
        ovar.standard_name = name
    if unit =='':   
       ovar.units = Ivar.units
    else:
        ovar.units = unit   
    if isinstance(Iname,str): #basestring): 
        if levtyp=='full':  ovar[:,:] = Ivar[:,-num_lev:,nc]
        if levtyp=='half':  ovar[:,:] = Ivar[:,-(num_lev+1):,nc]
    else:
        #ovar[:,:] = Ivar[:,-num_lev:]
        if levtyp=='full':  ovar[:,:] = Ivar[:,-num_lev:]
        if levtyp=='half':  ovar[:,:] = Ivar[:,-(num_lev+1):]
    return ovar


# plotting parameters
nc = 20 # cell number
cp = 1005.0 # J/K kg
L = 2.5e6 # J/kg
expname = argv[1]
gscp_scheme = "2mom" # expname[0:4] #"1mom"
print (gscp_scheme)
## read in file
sname = '../'+expname+'/scalars_'+expname+'_ICON_SCM_DWD_DOM01_ML_0001.nc'
sfile = netCDF4.Dataset(sname)

pname = '../'+expname+'/profs_'+expname+'_ICON_SCM_DWD_DOM01_ML_0001.nc'
pfile = netCDF4.Dataset(pname)

N_t = pfile.variables['time'][:].shape[0]
z_mc = pfile.variables['z_mc'][:,nc]
z_ifc = pfile.variables['z_ifc'][:,nc]
rho = pfile.variables['rho'][:,:,nc] # kg/m^3

#oufile
oname = 'profs_'+expname+'_ICON_SCM_DWD_v1.nc'
os.system('mkdir '+expname)

ofile = netCDF4.Dataset(expname+'/'+oname,'w')
ofile.createDimension('time',size=0)
ofile.createDimension('zf',size=z_mc[-num_lev:].shape[0])
ofile.createDimension('zh',size=z_ifc[-(num_lev+1):].shape[0])

### VARS
#time
time_o = ofile.createVariable('time',np.float64,('time'))
time = pfile.variables['time']
n_hours = 19.0
dt_in_hours = (n_hours)/(len(time[:])-1)
time_o[:] = [3600.0 *(dt_in_hours) * float(i) for i in range(len(time[:]))] # 60.0  *time[:]
time_o.standard_name = 'time'
time_o.units = 'seconds since 2014-11-24 17:00:00'

# height at cell centers
height_mc_o = ofile.createVariable('zf',np.float64,('zf'))
height_mc_o[:] =  z_mc[-num_lev:]
height_mc_o.standard_name = 'Height at full model levels'
height_mc_o.units = 'm'

# height at half levels
height_ifc_o = ofile.createVariable('zh',np.float64,('zh'))
height_ifc_o[:] =  z_ifc[-(num_lev+1):]
height_ifc_o.standard_name = 'Height at half model levels'
height_ifc_o.units = 'm'

#pressure [Pa]
p_o = set_var('p','pres',name='Pressure',unit='Pa')


# theta [K]
temp = pfile.variables['temp'][:,:,nc]
exner = pfile.variables['exner'][:,:,nc]
theta_v = pfile.variables['theta_v'][:,:,nc]
#theta_o = ofile.createVariable('theta',np.float64,('time','zf'))
#theta_o.standard_name = 'Potential temperature'
#theta_o.units = 'K'
#theta_o[:,:] = temp[:,:] / exner[:,:]

# relhum [%]
rh_o = set_var('rh','rh',name='Relative humidity',unit='%')

#qv [kg/kg]
qv_o = set_var('qv','tot_qv_dia',name='Vapour specific humidity',unit='kg/kg')

#qc [kg/kg]
qc_o = set_var('qc','tot_qc_dia',name='Condensed water specific humidity',unit='kg/kg')

#qr [kg/kg]
set_var('qr','qr',name='Rain water specific humidity',unit='kg/kg')

## special output for 2 moment microphysics
if gscp_scheme == "2mom":
    #qnc [1/m^3]
    qnc = pfile.variables['qnc'][:,:,nc] # 1/kg
    var = rho[:] * qnc[:]
    qnc_o = set_var('nc',var,name='Cloud droplet number concentration',unit='1/m**3')
    
    #qnr [1/m^3]
    qnr = pfile.variables['qnr'][:,:,nc] # 1/kg
    var = rho[:] * qnr[:]
    set_var('nr',var,name='Rain droplet number concentration',unit='1/m**3')
    
    # Rc [m]
    rho_w = 1000.0 # kg/m^3
    #Rc = np.ma.masked_where( qnc < 1.0e-12, (3.0/4.0*qc_o[:]/(rho_w*qnc[:]*np.pi))**(1.0/3.0) )
    #Rc = np.ma.masked_where( qnc < 1.0e-12, (3.0/4.0*qc_o[:]/(rho_w*qnc[:]*np.pi))**(1.0/3.0) )
    Rc = np.ma.masked_where( qnc_o[:] < 1.0e-12, (3.0/4.0*qc_o[:]/(rho_w*qnc_o[:]*np.pi))**(1.0/3.0) )
    set_var('rc_diag',Rc,name='diagnostic radius of cloud droplet',unit='m')
### 

##Rc_rad [m] cloud droplet radius used in radiation
#rc_rad = np.ma.masked_where(qc_o[:,:] < 1e-12, pfile.variables['re_diag'][:,:,nc]*1.0e-6) # mum
#rc_rad = np.ma.masked_where(qc_o[:,:] < 1e-12, pfile.variables['reff_qc'][:,:,nc]) # mum
#rc_rad = pfile.variables['re_diag'][:,:,nc]*1.0e-6
#set_var('rc_rad',rc_rad,name='Mean radius of cloud droplets used in radiation',unit='m')
#var = np.where(qc_o[:,:] > 1e-5, pfile.variables['re_diag'][:,:,nc]*1.0e-6,0.0) # mum
#rc_rad = np.insert(rc_rad,-1,0.0,axis=1)
#set_var('rc_rad',rc_rad,'half',name='Mean radius of cloud droplets used in radiation',unit='m')
#rc_diag_o = ofile.createVariable('rc_rad',np.float64,('time'))

#effective radius
set_var('reff_qc','reff_qc',name='effective radius of cloud droplets',unit='m')
set_var('reff_qr','reff_qr',name='effective radius of rain droplets',unit='m')

# cf
clc = pfile.variables['clc'][:,:,nc] * 0.01
set_var('cf',clc,name='Cloud fraction',unit='0-1')

#ddt_temp_gscp [K/s]
set_var('dt_mic','ddt_temp_gscp',name='Temperature increment due to cloud microphysics',unit='K/s')

# u
set_var('u','u',name='Zonal wind komponent',unit='m/s')

# v
set_var('v','v',name='Meridional wind komponent',unit='m/s')

#T flux 
T_g = pfile.variables['t_g'][:,nc]
p_sfc = pfile.variables['pres_sfc'][:,nc]
ex_ifc = diag.exner_ifc(exner[:,:],z_mc[:],z_ifc[:],p_sfc[:])
rho_ifc = diag.rho_ifc(rho[:,:],z_mc[:],z_ifc[:],p_sfc[:],T_g[:])
#var = - pfile.variables['tetfl_turb'][:,:,nc] /(rho_ifc[:,:] * cp)*ex_ifc[:,:]
var = - pfile.variables['tetfl_turb'][:,:,nc] /(rho_ifc[:,:])*ex_ifc[:,:]
set_var('wt',var,'half',name='Vertical turbulent  temperature flux',unit='Km/s')
#
##qv flux
#var = - pfile.variables['vapfl_turb'][:,:,nc]/(rho_ifc[:,:] * L)
var = - pfile.variables['vapfl_turb'][:,:,nc]/(rho_ifc[:,:] )
#set_var('wq',var,'half',name='Vertical turbulent moisture flux',unit='kg/kg m/s')
qv_o = set_var('wq',var,'half',name='Vertical turbulent moisture flux',unit='kg/kg m/s')

# turbulent liquid flux
var = -pfile.variables['liqfl_turb'][:,:,nc]/rho_ifc[:,:]
set_var('wl',var,'half',name='Vertical turbulent liquid flux',unit='kg/kg m/s')


# TKE [m^2/s^2]
set_var('TKE','tke','half',name='Turbulent kinetic energy',unit='m**2/s**2')

# buouyancy source of TKe [m^2/s^3]
set_var('buoy','ddt_tke_therm','half',name='Buoyancy production of TKE',unit='m**2/s**3')

set_var('shear','ddt_tke_mech','half',name='Shear production of TKE',unit='m**2/s**3')

#set_var('mech','ddt_tke_mech','half',name='Total mechanical production of TKE',unit='m**2/s**3')

# edr [m^2/s^3]
set_var('diss','edr','half',name='Dissipation of TKE',unit='m**2/s**3')

#tkvm [m^2/s]
set_var('km','tkvm','half',name='Eddy diffusivity for momentum', unit='m**2/s')
# tkvh
set_var('kh','tkvh','half',name='Eddy diffusivity for heat', unit='m**2/s')

#dt_turb [K/s]
set_var('dt_turb','ddt_temp_turb',name='Temperature increment due to turbulence',unit='K/s')

#lwup [W/m^2]
set_var('lwup','lwflx_up','half',name='Upwelling long-wave radiation',unit='W/m**2')
#lwdn [W/m^2]
set_var('lwdn','lwflx_dn','half',name='Downwelling long-wave radiation',unit='W/m**2')
#swup [W/m^2]
set_var('swup','swflx_up','half',name='Upwelling short-wave radiation',unit='W/m**2')
#swdn [W/m^2]
set_var('swdn','swflx_dn','half',name='Downwelling short-wave radiation',unit='W/m**2')

#dt_lw [K/s]
set_var('dt_lw','ddt_temp_radlw',name='Temperature increment due to long-wave radiaion',unit='K/s')
#dt_sw [K/s]
set_var('dt_sw','ddt_temp_radsw',name='Temperature increment due to short-wave radiaion',unit='K/s')
#ofile.createVariable('',np.floa,64,('time','zf'))
#ofile.createVariable('',np.float64,('time','zf'))
#ofile.createVariable('',np.float64,('time','zf'))



ofile.close()
