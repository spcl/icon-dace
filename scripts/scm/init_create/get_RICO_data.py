#!/usr/bin/env python

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# prepare SCM/LES input and forcing data for the RICO
# (marine shallow convection) case in netCDF file
# based on microHH setup of RICO case: https://github.com/microhh/microhh
#
# Ivan Bastak Duran, 2020
#-------------------------------------------------------------------

import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import sys, getopt, os       # file handling and operators
from datetime import datetime

#-------------------------------------------------------------------
case = 'gcss' # Original RICO
#case = 'ss08' # Moist RICO from Stevens/Seifert & Seifert/Heus
#case = 'test' # More moist mixed-layer for testing
print ("Beginning of main")

icon_dir = '../../../'
data_dir = 'data/scm/init_data/'
file_out = icon_dir + data_dir + 'init_SCM_RICO.nc'

# specify model levels
hlevels = np.array([0.0,740,2260,2980,3260,4000.,6000.,20000.])
hlevels = hlevels[hlevels.size::-1] # inverse order like model levels ?

# constants
R_d        = 287.06            # dry air gas constant
g          = 9.81              # grav. acceleration
dtdz_str   = 0.04              # stratospheric gradient of potential temperature [K/m]
dtdz_fratm = 0.005             # free atmosphere gradient of potential temperature [K/m]

# variable initialisation
z     = hlevels
thl   = np.zeros(np.size(z))
qt    = np.zeros(np.size(z))
u     = np.zeros(np.size(z))
v     = np.zeros(np.size(z))
tke   = np.zeros(np.size(z))
ug    = np.zeros(np.size(z))
vg    = np.zeros(np.size(z))
wls   = np.zeros(np.size(z))
thlls = np.zeros(np.size(z))
qtls  = np.zeros(np.size(z))

for k in range(z.size):

 # Liquid water potential temperature: same in GCSS and SS08
 if(z[k] < 740.):
     thl[k] = 297.9
 elif(z[k] <= 4000):
     thl[k] = 297.9 + (317.0 - 297.9)/(4000. - 740.) * (z[k] - 740.) 
 else:
    thl[k] = 317+(z[k]-4000.)*dtdz_fratm

 if(case == 'gcss'):
     if(z[k] < 740.):
         qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
     elif(z[k] < 3260.):
         qt[k] = 13.8 + (2.4 - 13.8) / (3260. - 740.) * (z[k] - 740.) 
     elif(z[k] <= 4000):
         qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.) 
     else:
         qt[k] = 0.0

 elif(case == 'ss08'):
     if(z[k] < 740.):
         qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
     elif(z[k] < 3260.):
         qt[k] = 13.8 + (4.4 - 13.8) / (3260. - 740.) * (z[k] - 740.) 
     elif(z[k] <= 4000):
         qt[k] = 4.4 + (3.6 - 4.4)/(4000. - 3260.) * (z[k] - 3260.) 
     else:
         qt[k] = 0.0

 elif(case == 'test'):
     q0 = 18.
     q1 = 15.8
     if(z[k] < 740.):
         qt[k] = q0 + (q1 - q0) / 740. * z[k]
     elif(z[k] < 3260.):
         qt[k] = q1 + (2.4 - q1) / (3260. - 740.) * (z[k] - 740.) 
     elif(z[k] <= 4000):
         qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.) 
     else:
         qt[k] = 0.0

 # Subsidence
 if(z[k] < 2260):
     wls[k] = -0.005 * (z[k] / 2260.)
 else:
     wls[k] = -0.005

 # U and V component wind
 u[k]  = -9.9 + 2.0e-3 * z[k]
 ug[k] = u[k]
 v[k]  = -3.8
 vg[k] = v[k]

 # Advective and radiative tendency thl
 thlls[k] = -2.5 / 86400.

 # Advective tendency qt
 if(z[k] < 2980):
     qtls[k] = -1.0 / 86400. + (1.3456/ 86400.) * z[k] / 2980.
 else:
     qtls[k] = 4e-6

# normalize profiles to SI
qt  /= 1000.
qtls/= 1000.

ep = 287.04 / 461.5 

#surface settings
def esat(T):
 c0 = 0.6105851e+03; c1 = 0.4440316e+02; c2 = 0.1430341e+01; c3 = 0.2641412e-01 
 c4 = 0.2995057e-03; c5 = 0.2031998e-05; c6 = 0.6936113e-08; c7 = 0.2564861e-11 
 c8 = -.3704404e-13 
 x  = max(-80.,T-273.15)
 return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

def qsat(p, T):
 return ep*esat(T)/(p-(1.-ep)*esat(T))

# surface parameters
p_surf       = 101540.  # surface pressure [Pa]
SST = 299.8 
ths = SST / (p_surf/1.e5)**(287.04/1005.)
qs  = qsat(p_surf, SST) 
z0           = 0.0002   #m
CM           = 0.001229 
CQ           = 0.001133 
CH           = 0.001094 
dt_relax     = 7200.0   # s

longitude = 6.417
latitude  = -61.77


# write data into netcdf file

nlev=hlevels.size
nlevTsoil=9
nlevWsoil=8
nnclass_lu=23
print ("write to netCDF file ", file_out)
do=Dataset(file_out,mode='w',format='NETCDF4_CLASSIC')
         
do.description='Initial profiles and forcings for RICO case'
do.createDimension('lev',size=nlev)
do.createDimension('nt', 1)
do.createDimension('levTsoil',size=nlevTsoil)
do.createDimension('levWsoil',size=nlevWsoil)
do.createDimension('nclass_lu',size=nnclass_lu)

#--- coordinates

#time
nc_time=do.createVariable('time',np.float64,'nt')
nc_time.units='s since 2004-12-16T00:00:00Z'

#level
lev=do.createVariable('lev',np.float64,'lev')
lev.units='m'
lev.standard_name='height'
lev[:]=np.arange(nlev)

#soil levels
levTsoil=do.createVariable('levTsoil',np.float64,'levTsoil')
levTsoil.units='level'
levTsoil[:]=np.arange(nlevTsoil)

levWsoil=do.createVariable('levWsoil',np.float64,'levWsoil')
levWsoil.units='level'
levWsoil[:]=np.arange(nlevWsoil)

#atmospheric coordinates
nc_height=do.createVariable('height',np.float64,('nt','lev'))
nc_height.units='meter above surface'

#soil coordinates
nc_depth_T=do.createVariable('depth_T',np.float64,('levTsoil',))
nc_depth_W=do.createVariable('depth_W',np.float64,('levWsoil',))
nc_depth_T.units='depth of T in meter'
nc_depth_W.units='depth of W in meter'

#--- initial profiles

##atmospheric variables
nc_uIN=do.createVariable('uIN',np.float64,('nt','lev'))
nc_uIN.description='zonal wind'
nc_uIN.units='m/s'
nc_vIN=do.createVariable('vIN',np.float64,('nt','lev'))
nc_vIN.description='meridional wind'
nc_vIN.units='m/s'
nc_wIN=do.createVariable('wIN',np.float64,('nt','lev'))
nc_wIN.description='veritcal velocity'
nc_wIN.units='m/s'
nc_thIN=do.createVariable('thIN',np.float64,('nt','lev'))
nc_thIN.description='potential temperature'
nc_thIN.units='K'
nc_qvIN=do.createVariable('qvIN',np.float64,('nt','lev'))
nc_qvIN.description='water vapour specific moisture'
nc_qvIN.units='kg/kg'
nc_qcIN=do.createVariable('qcIN',np.float64,('nt','lev'))
nc_qcIN.description='liquid water specific moisture'
nc_qcIN.units='kg/kg'
nc_pIN=do.createVariable('pIN',np.float64,('nt','lev'))
nc_pIN.description='atmospheric pressure'
nc_pIN.units='Pa'
nc_tkeIN=do.createVariable('tkeIN',np.float64,('nt','lev'))
nc_tkeIN.description='turbulence kinetic energy'
nc_tkeIN.units='m^2/s^2'
###values
nc_height[0,:]=hlevels
nc_time[:]=0.0
nc_uIN[:]=u
nc_uIN[:,0]=nc_uIN[:,1]#no gradient of u in stratosphere
nc_vIN[:]=v
nc_vIN[:,0]=nc_vIN[:,1]#no gradient of u in stratosphere
nc_wIN[:]=0.0
nc_tkeIN[:]=tke
nc_thIN[:]=thl
nc_thIN[:,0]=nc_thIN[:,1]+dtdz_str*(nc_height[:,0]-nc_height[:,1])#statospheric gradient of pot. temperature, levels above are extrapolated according to this gradinet in SCM
nc_qvIN[:]=qt
nc_qvIN[:,0]=nc_qvIN[:,1] #no gradient of moisture in stratosphere
nc_qcIN[:]=0.0
nc_pIN[:]=-999.0#should not be used in the simulation
nc_pIN[:,0]=p_surf

##soil variables
nc_T_SO=do.createVariable('T_SO',np.float64,('nt','levTsoil'))
nc_T_SO.description='soil temperature'
nc_T_SO.units='K'
nc_W_SO=do.createVariable('W_SO',np.float64,('nt','levWsoil'))
nc_W_SO.description='water content'
nc_W_SO.units='kg m-2'
#values
nc_T_SO[:]=0.0
nc_W_SO[:]=0.0

##forcings
nc_uGEO=do.createVariable('uGEO',np.float64,('nt','lev'))
nc_uGEO.description='zonal geostrophic wind'
nc_uGEO.units='m/s'
nc_vGEO=do.createVariable('vGEO',np.float64,('nt','lev'))
nc_vGEO.description='meridional geostrophic wind'
nc_vGEO.units='m/s'
nc_wLS=do.createVariable('wLS',np.float64,('nt','lev'))
nc_wLS.description='veritcal velocity'
nc_wLS.units='m/s'
nc_dTadv=do.createVariable('dTadv',np.float64,('nt','lev'))
nc_dTadv.description='potential temperature tendency - advection'
nc_dTadv.units='K/s'
nc_dTrad=do.createVariable('dTrad',np.float64,('nt','lev'))
nc_dTrad.description='potential temperature tendency - radiation'
nc_dTrad.units='K/s'
nc_dQVadv=do.createVariable('dQVadv',np.float64,('nt','lev'))
nc_dQVadv.description='potential temperature tendency - advection'
nc_dQVadv.units='kg/kg/s'
nc_dUadv=do.createVariable('dUadv',np.float64,('nt','lev'))
nc_dUadv.description='zonal wind endency - advection'
nc_dUadv.units='m/s/s'
nc_dVadv=do.createVariable('dVadv',np.float64,('nt','lev'))
nc_dVadv.description='meridional wind tendency - advection'
nc_dVadv.units='m/s/s'
nc_height[:]=hlevels
###values
nc_uGEO[:]=ug
nc_vGEO[:]=vg
nc_wLS[:]=wls
#nc_wLS[:,0]=nc_wLS[:,1]#no gradient in stratosphere
nc_dTadv[:]=thlls
#nc_dTadv[:,0]=nc_dTadv[:,1]#no gradient in stratosphere
nc_dTrad[:]=nc_dTadv[:]*0.0 # setting radiation tendency to 0
nc_dQVadv[:]=qtls
#nc_dQVadv[:,0]=nc_dQVadv[:,1]#no gradient in stratosphere
nc_dUadv[:]=0.0
nc_dVadv[:]=0.0

#--- surface

##pressure
nc_psurf=do.createVariable('psurf',np.float64,('nt'))
nc_psurf.description='Surface atmospheric pressure'
nc_psurf.units='Pa'
##fluxes
nc_sfc_lat_flx=do.createVariable('sfc_lat_flx',np.float64,('nt'))
nc_sfc_lat_flx.description='Latent heat flux'
nc_sfc_lat_flx.units='W/m2'
nc_sfc_sens_flx=do.createVariable('sfc_sens_flx',np.float64,('nt'))
nc_sfc_sens_flx.description='Sensible heat flux'
nc_sfc_sens_flx.units='W/m2'
##temperature and moisture
nc_ts=do.createVariable('ts',np.float64,('nt'))
nc_ts.description='Surface temperature'
nc_ts.units='K'
nc_tg=do.createVariable('tg',np.float64,('nt'))
nc_tg.description='Skin temperature'
nc_tg.units='K'
nc_qvs=do.createVariable('qvs',np.float64,('nt'))
nc_qvs.description='Surface specific vapor humidity'
nc_qvs.units='kg/kg'
## drag coefficients
nc_Cm=do.createVariable('Cm',np.float64,('nt'))
nc_Cm.description='Drag coefficient for heat'
nc_Cm.units=''
nc_Ch=do.createVariable('Ch',np.float64,('nt'))
nc_Ch.description='Drag coefficient for heat'
nc_Ch.units=''
nc_Cq=do.createVariable('Cq',np.float64,('nt'))
nc_Cq.description='Drag coefficient for moisture'
nc_Cq.units=''
## Friction velocity
nc_ustar=do.createVariable('ustar',np.float64,('nt'))
nc_ustar.description='Friction velocity'
nc_ustar.units='m/s'
###values
nc_sfc_lat_flx[:]=0.0
nc_sfc_sens_flx[:]=0.0
nc_psurf[:]=p_surf
nc_ts[:]=ths
nc_tg[:]=SST
nc_qvs[:]=qs
nc_Ch[:]=CH
nc_Cq[:]=CQ
nc_Cm[:]=CM
nc_ustar[:]=0.0

#--- constants

#nudging parameters
nc_dt_relax=do.createVariable('dt_relax',np.float64)
nc_dt_relax.description='relaxation time step for nudging'
nc_dt_relax.units='s'

nc_dt_relax[:]=dt_relax

#external parameters
nc_longitude=do.createVariable('longitude',np.float64)
nc_longitude.description='longitude'
nc_longitude.units='degrees'
nc_latitude=do.createVariable('latitude',np.float64)
nc_latitude.description='latitude'
nc_latitude.units='degrees'
nc_FR_LAND=do.createVariable('FR_LAND',np.float64)
nc_FR_LAND.description='land area fraction'
nc_FR_LAND.units=''
nc_PLCOV_MX=do.createVariable('PLCOV_MX',np.float64)
nc_PLCOV_MX.description='Plant cover maximum due to land use data'
nc_PLCOV_MX.units=''
nc_LAI_MX=do.createVariable('LAI_MX',np.float64)
nc_LAI_MX.description='Leaf Area Index Maximum'
nc_LAI_MX.units=''
nc_ROOTDP=do.createVariable('ROOTDP',np.float64)
nc_ROOTDP.description='Root depth'
nc_ROOTDP.units=''
nc_RSMIN=do.createVariable('RSMIN',np.float64)
nc_RSMIN.description='Minimal stomata resistence'
nc_RSMIN.units=''
nc_SOILTYP=do.createVariable('SOILTYP',np.float64)
nc_SOILTYP.description='soil type derived from FAO Digital Soil Map of the World'
nc_SOILTYP.units=''
nc_Z0=do.createVariable('Z0',np.float64)
nc_Z0.description='Roughness length'
nc_Z0.units=''
nc_EMIS_RAD=do.createVariable('EMIS_RAD',np.float64)
nc_EMIS_RAD.description='longwave surface emissivity'
nc_EMIS_RAD.units=''
nc_TOPO=do.createVariable('TOPO',np.float64)
nc_TOPO.description='geometric height of the earths surface above sea level'
nc_TOPO.units='m'
nc_LU_CLASS_FRACTION=do.createVariable('LU_CLASS_FRACTION',np.float64,('nclass_lu'))
nc_LU_CLASS_FRACTION.description='Fraction of land use classes in target grid element'
nc_LU_CLASS_FRACTION.units=''
nc_LU_CLASS_FRACTION.lctype='GLOBCOVER2009'

nc_longitude[:]=longitude
nc_latitude[:]=latitude
nc_FR_LAND[:]=1.0
nc_PLCOV_MX[:]=0.0
nc_LAI_MX[:]=0.0
nc_ROOTDP[:]=0.0
nc_RSMIN[:]=0.0
nc_SOILTYP[:]=8
nc_Z0[:]=z0
nc_EMIS_RAD[:]=0.98
nc_TOPO[:]=0.0
nc_LU_CLASS_FRACTION[:]=np.zeros([23])
nc_LU_CLASS_FRACTION[5]=1.0#forrest : Closed (>40%) broadleaved deciduous forest (>5m

today = datetime.today()
do.history = "Created " + today.strftime("%d/%m/%y")

do.close()
