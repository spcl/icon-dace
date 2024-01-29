#!/usr/bin/env python3

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# prepare SCM/LES input and forcing data from ICON native netCDF file
#
# runs as: python3 get_SCM_data_ICON.py lat lon   (lon between -180 and 180!!)
#
# work flow SCM from ICON input:
#  - ICON ini:   read_icon_ana_oper_mem1_40km.s
#  - ICON run:   run_ICON_4_SCMini
#  - SCM ini:    get_SCM_data_ICON.py
#  - SCM extpar: create_SCM_extpar_ICON.py
#  - SCM run:    run_SCM_ICONini
#  - plot SCM:   plot-scm-*.py
#
# Martin Koehler, 2020-04 
# (based on COSMO version from Ivan Bastak Duran and Annika Schomburg)
#-------------------------------------------------------------------

import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import sys, getopt, os       # file handling and operators
from datetime import datetime

#-------------------------------------------------------------------
# arguments

print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:      ', str(sys.argv) )
lat_scm = float(sys.argv[1])
lon_scm = float(sys.argv[2])

#-------------------------------------------------------------------
# setup

print ("Beginning of main -------------------------------------------")

# input and output files: ML, height, extpar, output

inidate  = '2021061700'           # initial date

# select location

#lat_scm = 52.209722  # Lindenberg
#lon_scm = 14.118889  #   - " -
#lat_scm = 33.81      # Mediterranean
#lon_scm = 29.43      #   - " -
#lat_scm = -60.0      # Southern Ocean -40.0  -40.0  -60.0  -60.0
#lon_scm =   0.0      #   - " -          0.0   10.0    0.0   10.0

icon_dir = '../../../'
data_dir = 'scripts/scm/init_create/'
ext_dir  = '/hpc/rhome/routfor/routfox/icon/grids/public/edzw/'

file_ml  = '/hpc/uwork/mkoehler/run-icon/experiments/exp_003_'+inidate+'/exp_003_'+inidate+'_DOM01_ML_0001.nc'
file_z   = '/hpc/uwork/mkoehler/run-icon/experiments/exp_003_'+inidate+'/exp_003_'+inidate+'_z_DOM01_ML_0001.nc'
file_ext = ext_dir + 'icon_extpar_0024_R02B06_G_20200917_tiles.nc'
file_out = icon_dir + data_dir + 'init_SCM_data_ICON_'+inidate+'_lat'+str(round(lat_scm, 1))+'_lon'+str(round(lon_scm, 1))+'.nc'

nc_fid = Dataset(file_ml,  'r')   # Dataset is the class behavior to open the file
nc_zid = Dataset(file_z,   'r')   # and create an instance of the ncCDF4 class
nc_eid = Dataset(file_ext, 'r')

nc_fid.set_auto_mask(False)
nc_zid.set_auto_mask(False)
nc_eid.set_auto_mask(False)

# read input NetCDF file

print ("... read from netCDF ML file:", '\n ', file_ml)

# Extract data from NetCDF file
clat = nc_fid.variables['clat'][:]
clon = nc_fid.variables['clon'][:]

print('  clat: ', type(clat), np.shape(clat), np.amin(clat), 'to', np.amax(clat))
print('  clon: ', type(clon), np.shape(clon), np.amin(clon), 'to', np.amax(clon))

dist = (clon - lon_scm/180.*np.pi)**2 + (clat - lat_scm/180.*np.pi)**2
index_near = np.where(dist == dist.min())
index_near = index_near[0][0]          # convert tuple to number
print('  selected lat/lon (rad) ', clat[index_near], clon[index_near], 'and deg',
  clat[index_near]/np.pi*180.0, clon[index_near]/np.pi*180.0  )
print('  target location: (lat/lon in deg) ', lat_scm, lon_scm)

if lon_scm/180.*np.pi < np.amin(clon) or lon_scm/180.*np.pi > np.amax(clon):
  sys.exit('ERROR: longitude range inconsistency')

# get external parameters

print ('... reading external parameter from extpar file:', '\n ', file_ext)
FR_LAND  = nc_eid.variables['FR_LAND'] [index_near]
PLCOV_MX = nc_eid.variables['PLCOV_MX'][index_near]
LAI_MX   = nc_eid.variables['LAI_MX']  [index_near]
ROOTDP   = nc_eid.variables['ROOTDP']  [index_near]
RSMIN    = nc_eid.variables['RSMIN']   [index_near]
SOILTYP  = nc_eid.variables['SOILTYP'] [index_near]
Z0       = nc_eid.variables['Z0']      [index_near]
EMIS_RAD = nc_eid.variables['EMIS_RAD'][index_near]
TOPO     = nc_eid.variables['topography_c'] [index_near]
LU_CLASS_FRACTION = nc_eid.variables['LU_CLASS_FRACTION'] [:,index_near]

nc_eid.close()
print ('... finished reading external parameters ')

print ('... reading atmospheric fields')

timef         = nc_fid.variables['time']      [:]
timef_units   = nc_fid.variables['time'].units
heights       = nc_zid.variables['z_mc']      [:,index_near]
heights_w     = nc_zid.variables['z_ifc']     [:,index_near]

pres_sfc      = nc_fid.variables['pres_sfc']  [:,index_near]
t_g           = nc_fid.variables['t_g']       [:,index_near]
t_s           = nc_fid.variables['t_s']       [:,index_near]
qv_s          = nc_fid.variables['qv_s']      [:,index_near]
accshfl       = nc_fid.variables['accshfl_s'] [:,index_near]
acclhfl       = nc_fid.variables['acclhfl_s'] [:,index_near]

temp          = nc_fid.variables['temp']      [:,:,index_near]
theta_v       = nc_fid.variables['theta_v']   [:,:,index_near]
exner         = nc_fid.variables['exner']     [:,:,index_near]
rho           = nc_fid.variables['rho']       [:,:,index_near]
u             = nc_fid.variables['u']         [:,:,index_near]
v             = nc_fid.variables['v']         [:,:,index_near]
w             = nc_fid.variables['w']         [:,:,index_near]
o3            = nc_fid.variables['o3']        [:,:,index_near]
tke           = nc_fid.variables['tke']       [:,:,index_near]
qv            = nc_fid.variables['qv']        [:,:,index_near]
qc            = nc_fid.variables['qc']        [:,:,index_near]
qr            = nc_fid.variables['qr']        [:,:,index_near]
qi            = nc_fid.variables['qi']        [:,:,index_near]
#qs           = nc_fid.variables['qs']        [:,:,index_near]
#qg           = nc_fid.variables['qg']        [:,:,index_near]

nlev  = np.shape(heights)  [0]
nlev1 = np.shape(heights_w)[0]
nt    = timef.size
print ('  nlev', nlev, ', nlev1', nlev1, ', nt', nt)

# add time dimension to heights
hhh = np.copy(u)
for it in range(nt):
  hhh[it,:] = heights
heights = hhh
print('  heights [m]: ', type(heights), np.shape(heights), np.amin(heights), 'to', np.amax(heights))
print('  time: ', timef, nc_fid.variables['time'].units)

hhh2 = np.copy(w)
for it in range(nt):
  hhh2[it,:] = heights_w
heights_w = hhh2

# pot. temperature
theta = temp / exner

# deaccumulate the fluxes
shfl = accshfl * 0.0
lhfl = acclhfl * 0.0
dt   = ( timef[1] - timef[0] ) * 60.0     # convert from [min] to [s]
for it in range(1,nt):
   shfl[it] = ( accshfl[it] - accshfl[it-1] ) / dt
   lhfl[it] = ( acclhfl[it] - acclhfl[it-1] ) / dt
#print('  shfl ', shfl, '   lhfl ', lhfl)
#print('  qv [kg/kg]: ', type(qv), np.shape(qv), np.amin(qv), 'to', np.amax(qv))

print ("... finishing reading atmospheric fields")

# get soil fields

print ("... reading soil fields")
depth_T       = nc_fid.variables['depth']   [:]
depth_W       = nc_fid.variables['depth_2'] [:]
soil_temp     = nc_fid.variables['t_so']    [:,:,index_near]
soil_moisture = nc_fid.variables['w_so']    [:,:,index_near]
soil_ice      = nc_fid.variables['w_so_ice'][:,:,index_near]

nc_zid.close()
nc_fid.close()
print ("... finishing reading soil fields")

dt_relax = 3 * 3600.0


#-------------------------------------------------------------------
# write data into netcdf file

print ("... write to netCDF file", file_out)
do=Dataset(file_out,mode='w',format='NETCDF4_CLASSIC')
         
do.description='Initial profiles and forcings generated from ICON data'
do.createDimension('lev'      ,size=nlev)
do.createDimension('lev1'     ,size=nlev1)
do.createDimension('nt'       ,size=timef.size)
do.createDimension('levTsoil' ,size=depth_T.size)
do.createDimension('levWsoil' ,size=depth_W.size)
do.createDimension('nclass_lu',size=LU_CLASS_FRACTION.size)


#--- coordinates

#time
nc_time=do.createVariable('time',np.float64,'nt')
nc_time.units=timef_units

#level - full level
lev=do.createVariable('lev',np.float64,'lev')
lev.units='m'
lev.standard_name='height'
lev[:]=np.arange(nlev)

#level - half level
lev1=do.createVariable('lev1',np.float64,'lev1')
lev1.units='m'
lev1.standard_name='height'
lev1[:]=np.arange(nlev1)

#soil levels
levTsoil=do.createVariable('levTsoil',np.float64,'levTsoil')
levTsoil.units='level'
levTsoil[:]=np.arange(depth_T.size)

levWsoil=do.createVariable('levWsoil',np.float64,'levWsoil')
levWsoil.units='level'
levWsoil[:]=np.arange(depth_W.size)

#atmospheric coordinates
nc_height=do.createVariable('height',np.float64,('nt','lev'))
nc_height.units='meter above surface'

nc_height_ifc=do.createVariable('height_ifc',np.float64,('nt','lev1'))
nc_height_ifc.units='meter above surface'

#soil coordinates
nc_depth_T=do.createVariable('depth_T',np.float64,('levTsoil',))
nc_depth_W=do.createVariable('depth_W',np.float64,('levWsoil',))
nc_depth_T.units='depth of T in meter'
nc_depth_W.units='depth of W in meter'
##values
nc_depth_T[:]=depth_T
nc_depth_W[:]=depth_W

#--- initial profiles

##atmospheric variables
nc_uIN              = do.createVariable('uIN',np.float64,('nt','lev'))
nc_uIN.description  = 'zonal wind'
nc_uIN.units        = 'm/s'
nc_vIN              = do.createVariable('vIN',np.float64,('nt','lev'))
nc_vIN.description  = 'meridional wind'
nc_vIN.units        = 'm/s'
nc_wIN              = do.createVariable('wIN',np.float64,('nt','lev1'))
nc_wIN.description  = 'veritcal velocity'
nc_wIN.units        = 'm/s'
nc_thIN             = do.createVariable('thIN',np.float64,('nt','lev'))
nc_thIN.description = 'potential temperature'
nc_thIN.units       = 'K'
nc_thvIN            = do.createVariable('thvIN',np.float64,('nt','lev'))
nc_thvIN.description= 'virtual potential temperature'
nc_thvIN.units      = 'K'
nc_qvIN             = do.createVariable('qvIN',np.float64,('nt','lev'))
nc_qvIN.description = 'water vapour specific moisture'
nc_qvIN.units       = 'kg/kg'
nc_qcIN             = do.createVariable('qcIN',np.float64,('nt','lev'))
nc_qcIN.description = 'liquid water specific moisture'
nc_qcIN.units       = 'kg/kg'
nc_qiIN             = do.createVariable('qiIN',np.float64,('nt','lev'))
nc_qiIN.description = 'cloud ice water specific moisture'
nc_qiIN.units       = 'kg/kg'
nc_o3IN             = do.createVariable('o3IN',np.float64,('nt','lev'))
nc_o3IN.description = 'ozone mass mixing ratio'
nc_o3IN.units       = 'kg/kg'
nc_tkeIN            = do.createVariable('tkeIN',np.float64,('nt','lev1'))
nc_tkeIN.description= 'turbulence kinetic energy'
nc_tkeIN.units      = 'm^2/s^2'
nc_exnerIN            = do.createVariable('exnerIN',np.float64,('nt','lev'))
nc_exnerIN.description= 'exner pressure'
nc_exnerIN.units      = 'Pa'
nc_rhoIN            = do.createVariable('rhoIN',np.float64,('nt','lev'))
nc_rhoIN.description= 'turbulence kinetic energy'
nc_rhoIN.units      = 'kg/m3'
###values
nc_height[:]     = heights
nc_height_ifc[:] = heights_w
nc_time[:]       = timef    # in seconds
nc_uIN[:]        = u
nc_vIN[:]        = v
nc_wIN[:]        = w
nc_tkeIN[:]      = tke
nc_thIN[:]       = theta
nc_thvIN[:]      = theta_v
nc_exnerIN[:]    = exner
nc_rhoIN[:]      = rho
nc_qvIN[:]       = qv
nc_qcIN[:]       = qc
nc_qiIN[:]       = qi
nc_o3IN[:]       = o3

##soil variables
nc_T_SO=do.createVariable('T_SO',np.float64,('nt','levTsoil'))
nc_T_SO.description='soil temperature'
nc_T_SO.units='K'
nc_W_SO=do.createVariable('W_SO',np.float64,('nt','levWsoil'))
nc_W_SO.description='water content'
nc_W_SO.units='kg m-2'
#values
nc_T_SO[:] = soil_temp
nc_W_SO[:] = soil_moisture

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
#values
nc_uGEO[:]     = 0.0
nc_vGEO[:]     = 0.0
nc_wLS[:]      = 0.0
nc_wLS[:,0]    = 0.0 # no gradient in stratosphere
nc_dTadv[:]    = 0.0
nc_dTadv[:,0]  = 0.0 # no gradient in stratosphere
nc_dTrad[:]    = 0.0 # setting radiation tendency to 0
nc_dQVadv[:]   = 0.0
nc_dQVadv[:,0] = 0.0 # no gradient in stratosphere
nc_dUadv[:]    = 0.0
nc_dVadv[:]    = 0.0 


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
nc_sfc_lat_flx[:]  = lhfl
nc_sfc_sens_flx[:] = shfl
nc_psurf[:]        = pres_sfc
nc_ts[:]           = t_s
nc_tg[:]           = t_g
nc_qvs[:]          = qv_s
nc_Ch[:]           = 0.0      #no input
nc_Cq[:]           = 0.0      #no input
nc_Cm[:]           = 0.0      #no input
nc_ustar[:]        = 0.0      #no input


#--- constants

#nudging parameters
nc_dt_relax=do.createVariable('dt_relax',np.float64)
nc_dt_relax.description='relaxation time step for nudging'
nc_dt_relax.units='s'
##values
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
nc_Z0.units='m'
nc_EMIS_RAD=do.createVariable('EMIS_RAD',np.float64)
nc_EMIS_RAD.description='longwave surface emissivity'
nc_EMIS_RAD.units=''
nc_TOPO=do.createVariable('TOPO',np.float64)
nc_TOPO.description='geometric height of the earths surface above sea level'
nc_TOPO.units='m'
nc_LU_CLASS_FRACTION=do.createVariable('LU_CLASS_FRACTION',np.float64,'nclass_lu')
nc_LU_CLASS_FRACTION.description='Fraction of land use classes in target grid element'
nc_LU_CLASS_FRACTION.units=''
nc_LU_CLASS_FRACTION.lctype='GLOBCOVER2009'
##values
nc_longitude[:]  = clon[index_near]/np.pi*180.0
nc_latitude[:]   = clat[index_near]/np.pi*180.0
nc_FR_LAND[:]    = FR_LAND
nc_PLCOV_MX[:]   = PLCOV_MX
nc_LAI_MX[:]     = LAI_MX
nc_ROOTDP[:]     = ROOTDP
nc_RSMIN[:]      = RSMIN
nc_SOILTYP[:]    = SOILTYP
nc_Z0[:]         = Z0
nc_EMIS_RAD[:]   = EMIS_RAD
nc_TOPO[:]       = TOPO
nc_LU_CLASS_FRACTION[:] = LU_CLASS_FRACTION

today = datetime.today()
do.history = "Created " + today.strftime("%d/%m/%y")

do.close()

print ("END -------------------------------------------------------")
