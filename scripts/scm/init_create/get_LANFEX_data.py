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

# prepare SCM/LES input and forcing data for the LANFEX
# (fog case in England) case in netCDF file
#
# Tobias Goecke, 2019
#-------------------------------------------------------------------

import numpy as np
import sys, getopt, os      # file handling and operators
from netCDF4 import Dataset # netCDF4 package for reading files
import datetime

#-------------------------------------------------------------------


print ("Beginning of main")

icon_dir = '../../../'
data_dir = 'data/scm/init_data/'
file_out = icon_dir + data_dir + 'init_SCM_LANFEX.nc'

# read in forcing from txt files
z_in, theta_in, qv_in, u_in, v_in = [], [], [], [], [] 
file_sound = open('LANFEX_create/my_init_profiles.txt')
for line in file_sound:
    dummy = line.split()
    z_in.append(float(dummy[0]))
    theta_in.append(float(dummy[1]))
    qv_in.append(float(dummy[2]))
    u_in.append(float(dummy[3]))
    v_in.append(float(dummy[4]))

start_date = datetime.datetime(2014,11,24,17)
dates, t_delta, seconds, T_surf = [], [], [], []
file_surf = open('LANFEX_create/surf_temp.txt')
for line in file_surf:
    dummy = line.split()
    datobj = datetime.datetime.strptime(dummy[0],"%Y-%m-%dT%H:%M:%SZ")
    dates.append( datobj )
    t_delta.append( dates[-1] - start_date  )
    #seconds.append( dates[-1].hour*3600.0+dates[-1].minute*60.0+dates[-1].second )
    seconds.append( t_delta[-1].seconds )
    print (dates[-1], t_delta[-1],seconds[-1])
    T_surf.append( float(dummy[1]) )
print ( "len(seconds) = ",len(seconds)," len(T_surf) =  ",len(T_surf))

# interpolate to other timesteps

seconds_target = np.arange(0.0,68600.0,600.0)
T_target = np.zeros_like(seconds_target)
for i in range(len(seconds_target)):
    idx = np.argmin( np.abs(seconds[:]-seconds_target[i])  ) #seconds.index(seconds_target[i])
    if seconds[idx]> seconds_target[i]: print ("PROBLEM")
    if idx+2 >= len(seconds) : idx = len(seconds)-2
    m = (T_surf[idx+1]-T_surf[idx])/(seconds[idx+1]-seconds[idx])
    b = T_surf[idx]-m*seconds[idx]
    T_target[i] = m*seconds_target[i] + b
    #print idx, seconds_target[i], seconds[idx], seconds[idx+1], T_surf[idx], T_surf[idx+1],T_target[i]

seconds = seconds_target
T_surf = T_target
#print " target values : "
#for t in range(len(seconds_target)):
#    print t+1,seconds_target[t], T_target[t]


# read and interpolate ozone

z_O3, O3 = [], []
file_O3 = open('LANFEX_create/ozone.txt')
for line in file_O3:
    dummy = line.split()
    z_O3.append( float(dummy[0]) )
    O3.append( float(dummy[1]) )
    #print "z_O3, O3", z_O3[-1], O3[-1]
 
O3_in = np.zeros_like(z_in)
for i, z in enumerate(z_in):
    idx = np.argmin( np.abs(np.asarray(z_O3)-z) )
    if (z_O3[idx] > z) and (idx != 0): idx = idx - 1 
    if idx+2 >= len(z_O3) : idx = len(z_O3) -2
    m = (O3[idx+1]-O3[idx])/(z_O3[idx+1]-z_O3[idx])
    b = O3[idx]-m*z_O3[idx]
    O3_in[i] = m*z_in[i] + b
    if O3_in[i] <= 0.0: O3_in[i] = 0.0
    print (i, z, idx, z_O3[idx], z_O3[idx+1], O3[idx], O3[idx+1],O3_in[i])
    #if i> 100: exit()


# specify model levels
#hlevels = np.array([0.0,300.,520.,700.,1480.,1500.0,2000.,2100.,3000.,5000.,20000.])
#hlevels = hlevels[hlevels.size::-1] # inverse order like model levels ?

seconds = np.asarray(seconds)

# constants
R_d   = 287.06            # dry air gas constant
g     = 9.81              # grav. acceleration



nt = seconds.shape[0]

# reverse orderting of vertical axes
z_in = z_in[::-1]
theta_in = theta_in[::-1]
qv_in = qv_in[::-1]
u_in = u_in[::-1]
v_in = v_in[::-1]
O3_in = O3_in[::-1]

hlevels = np.asarray(z_in)
# variable initialisation
tg    = np.asarray(T_surf)
z     = np.tile(z_in,(nt,1)) #np.asarray(z_in)
thl   = np.tile(theta_in,(nt,1)) #np.asarray(theta_in)
qt    = np.tile(qv_in,(nt,1)) #np.asarray(qv_in)
u     = np.tile(u_in,(nt,1)) #np.asarray(u_in)
v     = np.tile(v_in,(nt,1)) #np.asarray(v_in)
O3    = np.tile(O3_in,(nt,1)) #np.asarray(v_in)
tke   = np.zeros(np.size(z))
ug    = np.zeros(np.size(z))
wls   = np.zeros(np.size(z))
thlls = np.zeros(np.size(z))
qtls  = np.zeros(np.size(z))
print ("z.shape = ", z.shape)



# normalize profiles to SI
qtls  /= 1000.  # from g/kg to kg/kg
wls   /= 100.   # from cm/s to m/s
thlls /= 86400. # from K/d to K/s
qtls  *= 1.e-8

# surface parameters
ustar        = 0.0 #0.28     # m/s    
sfc_lat_flx  = 0.0 #150.0    # W/m2
sfc_sens_flx = 0.0 #9.0      # W/m2
p_surf       = 102350.0 #101500.0 # surface pressure [Pa]
z0           = 0.1 #0.1      # m
dt_relax     = 1e12 #7200.0   # s

#location Cardington UK
longitude = -0.42 #6.417
latitude  = 52.1 #30.00

#----------------------------------------------------------------

# write data into netcdf file

print ("write to netCDF file ", file_out)
do=Dataset(file_out,mode='w',format='NETCDF4_CLASSIC')
         
do.description='Initial profiles and forcings for LANFEX case'
do.createDimension('lev',size=hlevels.size)
do.createDimension('nt',size=nt)
do.createDimension('ndT',size=1)
do.createDimension('ndW',size=1)
do.createDimension('nclass_lu',size=23)

#coordinates
##time
nc_time=do.createVariable('time',np.float64,('nt',))
nc_time.units='seconds since '+datetime.datetime.strftime(start_date,"%Y-%m-%dT%H:%M:%SZ")
##atmospheric coordinates
nc_height=do.createVariable('height',np.float64,('nt','lev'))
nc_height.units='meter above surface'
##soil coordinates
nc_depth_T=do.createVariable('depth_T',np.float64,('ndT',))
nc_depth_W=do.createVariable('depth_W',np.float64,('ndW',))
nc_depth_T.units='depth of T in meter'
nc_depth_W.units='depth of W in meter'

#level
nc_lev=do.createVariable('lev',np.float64,('lev'))
nc_lev.units='m'
nc_lev.standard_name='level'
nc_lev=np.arange(np.size(hlevels))
#nc_lev[:] = lev[::-1]
#print lev[::-1]

#initial profiles
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
nc_pIN.description='atmospheric pressue'
nc_pIN.units='Pa'
nc_tkeIN=do.createVariable('tkeIN',np.float64,('nt','lev'))
nc_tkeIN.description='turbulence kinetic energy'
nc_tkeIN.units='m^2/s^2'
nc_o3IN             = do.createVariable('o3IN',np.float64,('nt','lev'))
nc_o3IN.description = 'ozone mass mixing ratio'
nc_o3IN.units       = 'kg/kg'
###values
nc_height[:] = np.tile(hlevels,(nt,1)) #hlevels
nc_time[:]=seconds[:]
nc_uIN[:]=u
nc_vIN[:]=v
nc_wIN[:]=0.0
nc_tkeIN[:]=tke
nc_thIN[:]=thl
nc_qvIN[:]=qt
nc_qcIN[:]=0.0
nc_pIN[:]=-999.0#should not be used in the simulation
nc_o3IN[:]=O3
#nc_pIN[0]=p_surf

##soil variables
nc_T_SO=do.createVariable('T_SO',np.float64,('nt','ndT'))
nc_T_SO.description='soil temperature'
nc_T_SO.units='K'
nc_W_SO=do.createVariable('W_SO',np.float64,('nt','ndW'))
nc_W_SO.description='water content'
nc_W_SO.units='kg m-2'
###values
nc_T_SO[:]=0.0
nc_W_SO[:]=0.0

#forcings
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
nc_height[0,:]=hlevels
###values
#nc_time[:]=0.0
nc_uGEO[:]=0.0 #ug
nc_vGEO[:]=0.0
nc_wLS[:]=0.0 #wls
nc_dTadv[:]=0.0 #thlls
nc_dTrad[:]=nc_dTadv[:]*0.0 # setting radiation tendency to 0
nc_dQVadv[:]=0.0 #qtls
nc_dUadv[:]=0.0
nc_dVadv[:]=0.0


#surface
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
nc_Ch=do.createVariable('Ch',np.float64,('nt'))
nc_Ch.description='Drag coefficient for heat'
nc_Ch.units=''
nc_Cq=do.createVariable('Cq',np.float64,('nt'))
nc_Cq.description='Drag coefficient for moisture'
nc_Cq.units=''
nc_Cm=do.createVariable('Cm',np.float64,('nt'))
nc_Cm.description='Drag coefficient for momentum'
nc_Cm.units=''
## Friction velocity
nc_ustar=do.createVariable('ustar',np.float64,('nt'))
nc_ustar.description='Friction velocity'
nc_ustar.units='m/s'
###values
nc_sfc_lat_flx[:]=sfc_lat_flx
nc_sfc_sens_flx[:]=sfc_sens_flx
nc_ts[:]=tg[:]
nc_tg[:] = tg[:]
nc_qvs[:]=0.0
nc_Ch[:]=0.0
nc_Cq[:]=0.0
nc_Cm[:]=0.0
nc_ustar[:]=ustar
nc_psurf[:]=p_surf

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

nc_longitude[:]=longitude
nc_latitude[:]=latitude
nc_FR_LAND[:]=1.0
nc_PLCOV_MX[:]=0.0
nc_LAI_MX[:]=2.0
nc_ROOTDP[:]=0.0
nc_RSMIN[:]=0.0
nc_SOILTYP[:]=6.0
nc_Z0[:]=z0
nc_EMIS_RAD[:]=0.98
nc_TOPO[:]=0.0
#nc_LU_CLASS_FRACTION[:]=0.0 
nc_LU_CLASS_FRACTION[:]=np.zeros([23])
nc_LU_CLASS_FRACTION[14]=1.0  #14:sparse herbaceous or grass
nc_LU_CLASS_FRACTION.lctype='GLC2000'  #'GLOBCOVER2009'
do.close()
