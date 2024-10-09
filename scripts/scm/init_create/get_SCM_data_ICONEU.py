#!/usr/bin/env python3

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

# prepare SCM/LES input and forcing data for realistic cases based on ICON EU forecast and analyses 
#
# runs as: python3 get_SCM_data_ICONEU.py
# requires: cdo
#
# Ivan Bastak Duran, 2020-06
#-------------------------------------------------------------------

from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import xarray as xr
import numpy as np
import sys, getopt, os       # file handling and operators
from datetime import datetime

#-------------------------------------------------------------------
# setup
datetime_str='2020060200'
nt = 26      # files (each 1hr) - forcings for this integration time

# Lindenberg
lat   = 52.209722
lon   = 14.118889
delta = 2.0  # size of averaging boxes in degrees

icon_dir      = '../../../'
data_dir      = 'data/scm/init_data/'

#inputs
##constant input
nc_c = 'ICONEU_cnst/ieaf0001010100.nc'  # Path to filename with constants

##forecast input
file_out = icon_dir + data_dir + "init_SCM"+datetime_str+"delta"+str(delta)+".nc"
#inpath='2020060200'#folder with forecast inputs
inpath=datetime_str#folder with forecast inputs
tmppath='tmp'+datetime_str#temporary folder
#initial conditions
nc_init_str='ieaf'
#forecast
nc_fcst_str='iefff'

damph=5000.#uGEO and vGEO will converge tu u and v at the location --experimental
damph_slope=2e-4  #--experimental
hlevels=np.array([0,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3500,4000,4500,5000,5500,6000,7500,8000,8500,9000,9500,10000,20000])
hlevels=hlevels[hlevels.size::-1] #inverse order like model levels ?
dt_relax     = 7200.0   # s- relaxation time scale
#end of setup-------------------------------------------------------------------
if not os.path.exists(tmppath):
   os.makedirs(tmppath)


# constants
R_d   = 287.06            # dry air gas constant
R_v   = 461.56            # water vapour gas constant
L_v   = 2.501e6           # latent heat of vaporization
g     = 9.81              # grav. acceleration
p0    = 100000.           # reference atm. pressure[Pa]
T0    = 273.15            # zero Celsius in K
cpd   = 1004.7            # the heat capacity of dry air at constant pressure
dtdz_str = 0.04           # stratospheric gradient of potential temperature [K/m]
Omega  = 2 * np.pi / 86400.#angular velocity of Earth
REarth = 6371000.785       #radius of Earth
Hmod1 = 20.0#heigth of lowest model level - cheat because we dont have surface pressure
nlev      = hlevels.size
f2=2.0 * Omega * np.cos(2 * np.pi * lat / 360.) 
f3 = 2.0 * Omega * np.sin(2 * np.pi * lat / 360.)
d_y= delta/360. * 2 * np.pi * REarth # calculate domain by using earth's circumference
d_x = d_y * np.cos(lat*np.pi/180.) # attribute for mapping factor (projection of the part-circle in x direction)
hsize=hlevels.size
trange=np.arange(0,nt+1,dtype='int32')
nt=trange.size


#functions
def diag_ustar(z,bflx,w,z0):
    '''
    diag_ustar outputs friction velocity
    Parameters
    ----------
    z : heigth[m]
    w : wind speed at z[m/s]
    bflx: surface buoyancy flux [m^2/s^3]
    z0:roughness height
    Returns
    -------
    ustar : friction velocity [m/s]
    '''
    am   =  6.0 #MO constant Hoegstoerm(1988)
    bm   = 19.3 #MO constant Hoegstoerm(1988)
    kappa= 0.4  #von Karman constant

    ni=5 #number of iterations
 
    lnz   = np.log( z / z0 )
    c1    = np.pi / 2.0 - 3.0*np.log( 2.0)

    ustar =  w*kappa*lnz
    for itime in range(0,ustar.size):
     if (np.abs(bflx[itime]) > 1.e-6):
      for i in range(0,ni):
       lmo   = -ustar[itime]**3 / ( kappa * bflx[itime] )
       zeta  = z/lmo
       if (zeta > 0.):
         if ( zeta > 1.e10):# large zeta
           ustar[itime] = 1e-10
           exit
         else:
           ustar[itime] =  kappa*w[itime]/(lnz + am*zeta)
       else:
         x     = np.sqrt( np.sqrt( 1.0 - bm*zeta ) ) 
         psi1  = 2.*np.log( 1.0+x )+ np.log( 1.0+x*x ) - 2.*np.arctan( x ) + c1
         ustar[itime] = w[itime]*kappa/(lnz - psi1)

    return ustar

def get_indices(lats,lons,radius,LATX = 52.209722, LONX = 14.118889):
    """ get_indices takes latitude latx and longitude lonx and
    returns the indices indexi and indexj of the nearest grid cell in a 2D array"""
    DISTANCE_TO_X = np.sqrt((lats - LATX)**2 + (lons - LONX)**2)
    indices = np.where(DISTANCE_TO_X < radius )
    return indices[0]
    pass

def get_index(lats,lons,LATX = 52.209722, LONX = 14.118889):
    """ get_indices takes latitude latx and longitude lonx and
    returns the indices indexi and indexj of the nearest grid cell in a 2D array"""
    DISTANCE_TO_X = np.sqrt((lats - LATX)**2 + (lons - LONX)**2)
    index = np.where(DISTANCE_TO_X == DISTANCE_TO_X.min())
    return index[0]
    pass
#end of functions


#read constant input
print ("reading constant fields")
nc_c_x=xr.open_dataset(nc_c)
lats_c=nc_c_x.CLAT.values[:]
lons_c=nc_c_x.CLON.values[:]
# get external parameters
key_index_c = get_index(lats_c,lons_c,lat, lon)
print ('key_index ', key_index_c,lats_c[key_index_c],lons_c[key_index_c],lat, lon)
print ("reading external parameters")
FR_LAND = nc_c_x.lsm.values[key_index_c]
PLCOV_MX = nc_c_x.PLCOV_MX.values[key_index_c]
LAI_MX = nc_c_x.LAI_MX.values[key_index_c]
ROOTDP = nc_c_x.ROOTDP.values[key_index_c]
RSMIN = nc_c_x.rsmin.values[key_index_c]
SOILTYP = nc_c_x.SOILTYP.values[key_index_c]
Z0 = nc_c_x.sr.values[key_index_c]
EMIS_RAD = nc_c_x.EMIS_RAD.values[key_index_c]
LU_CLASS_FRACTION = nc_c_x.FR_LUC.values[:,key_index_c]
nc_c_x.close()
print ('finished reading external parameters ')



# get atmospheric fields
print ("reading atmospheric fields")
print ("initial conditions")

grib_init=nc_init_str+datetime_str
nc_init=grib_init+'.nc'
key_index_s=str(key_index_c[0])
os.system('cdo -f nc selgridcell,'+key_index_s+' '+inpath+'/'+grib_init+' '+tmppath+'/'+nc_init)

nc_init_x=xr.open_dataset(tmppath+'/'+nc_init)
u=nc_init_x.u.values[0,:,0]
nk=u.shape[0]
v=nc_init_x.v.values[0,:,0]
w=nc_init_x.wz.values[0,:,0]
temp=nc_init_x.t.values[0,:,0]
tke=nc_init_x.tke.values[0,:,0]
pressure=nc_init_x.pres.values[0,:,0]
mlheigths_w=nc_init_x.HHL.values[0,:,0]
HSURF=nc_init_x.HHL.values[0,nk,0]
spec_hum=nc_init_x.q.values[0,:,0]
cloud_mix=nc_init_x.clwmr.values[0,:,0]
#rain_mix=nc_init_x.rwmr.values[0,:,0]
cloud_ice_mix=nc_init_x.QI.values[0,:,0]
snow_mix=nc_init_x.snmr.values[0,:,0]

#surface
qv_s=nc_init_x.QV_S.values[0,0]
t_g=nc_init_x.T_G.values[0,0]

#soil
depth_T =  nc_init_x.depth_2.values[:]
depth_W =  nc_init_x.depth.values[:]
soil_temp = nc_init_x.T_SO.values[0,:,0]
soil_moisture = nc_init_x.W_SO.values[0,:,0]
if 'W_SO_ICE' in nc_init_x:
   soil_ice = nc_init_x.W_SO_ICE.values[0,:,0]

nc_init_x.close()
os.system('rm '+tmppath+'/'+nc_init)


#interpolate
mlheigths=np.zeros([nk])
for iz in range(0,nk):
  mlheigths[iz]=0.5*(mlheigths_w[iz]+mlheigths_w[iz+1])

u_in_h=np.interp(hlevels+HSURF, mlheigths[nk:0:-1], u[nk:0:-1])
v_in_h=np.interp(hlevels+HSURF, mlheigths[nk:0:-1], v[nk:0:-1])
w_in_h=np.interp(hlevels+HSURF, mlheigths_w[nk+1:0:-1], w[nk+1:0:-1])
temp_in_h=np.interp(hlevels+HSURF, mlheigths[nk:0:-1], temp[nk:0:-1])
tke_in_h=np.interp(hlevels+HSURF, mlheigths_w[nk+1:0:-1], tke[nk+1:0:-1])
pressure_in_h=np.interp(hlevels+HSURF, mlheigths[nk:0:-1], pressure[nk:0:-1])
spec_hum_in_h=np.interp(hlevels+HSURF, mlheigths[nk:0:-1], spec_hum[nk:0:-1])
cloud_mix_in_h=np.interp(hlevels+HSURF, mlheigths[nk:0:-1], cloud_mix[nk:0:-1])
cloud_ice_mix_in_h=np.interp(hlevels+HSURF, mlheigths[nk:0:-1], cloud_ice_mix[nk:0:-1])
snow_mix_in_h=np.interp(hlevels+HSURF, mlheigths[nk:0:-1], snow_mix[nk:0:-1])

qt_in_h = (spec_hum_in_h + cloud_mix_in_h + cloud_ice_mix_in_h)/\
 (1.0+cloud_mix_in_h + cloud_ice_mix_in_h)
qc_in_h=(cloud_mix_in_h)*(1.0-qt_in_h)
qv_in_h=qt_in_h-qc_in_h
r_mix_in_h=spec_hum_in_h/(1.0-qt_in_h)
tv_in_h=temp_in_h*(1.0+0.608*r_mix_in_h-cloud_mix_in_h-cloud_ice_mix_in_h)
rho_in_h=pressure_in_h/(R_d*tv_in_h)
exner_in_h=(pressure_in_h/p0)**(R_d/cpd)
th_in_h=temp_in_h/exner_in_h
thv_in_h=tv_in_h/exner_in_h

#forecast nudging
print ("nudging")
u_n_f=np.zeros([nt,nk])
v_n_f=np.zeros([nt,nk])
w_n_f=np.zeros([nt,nk+1])
temp_n_f=np.zeros([nt,nk])
pressure_n_f=np.zeros([nt,nk])
spec_hum_n_f=np.zeros([nt,nk])
cloud_mix_n_f=np.zeros([nt,nk])
cloud_ice_mix_n_f=np.zeros([nt,nk])
for it in range(0,nt):
 it0=trange[it]//24
 it1=trange[it]%24
 grib_fcst=nc_fcst_str+"{:02d}".format(it0)+"{:02d}".format(it1)+'0000'
 nc_fcst=nc_fcst_str+"{:02d}".format(it0)+"{:02d}".format(it1)+'0000.nc'
 os.system('cdo -f nc selgridcell,'+key_index_s+' '+inpath+'/'+grib_fcst+' '+tmppath+'/'+nc_fcst)
 nc_fcst_x=xr.open_dataset(tmppath+'/'+nc_fcst)

 u_n_f[it,:]=nc_fcst_x.u[0,:,0].values
 v_n_f[it,:]=nc_fcst_x.v[0,:,0].values
 w_n_f[it,:]=nc_fcst_x.wz[0,:,0].values
 temp_n_f[it,:]=nc_fcst_x.t[0,:,0].values
 pressure_n_f[it,:]=nc_fcst_x.pres[0,:,0].values
 spec_hum_n_f[it,:]=nc_fcst_x.q[0,:,0].values
 cloud_mix_n_f[it,:]=nc_fcst_x.clwmr[0,:,0].values
 cloud_ice_mix_n_f[it,:]=nc_fcst_x.QI[0,:,0].values

 nc_fcst_x.close()
 os.system('rm '+tmppath+'/'+nc_fcst)

##compute model levels
print ("compute heigths of model levels")
qt_n_f = (spec_hum_n_f + cloud_mix_n_f + cloud_ice_mix_n_f)/\
 (1.0+cloud_mix_n_f + cloud_ice_mix_n_f)
qc_n_f=(cloud_mix_n_f)*(1.0-qt_n_f)
qv_n_f=qt_n_f-qc_n_f
r_mix_n_f=spec_hum_n_f/(1.0-qt_n_f)
tv_n_f=temp_n_f*(1.0+0.608*r_mix_n_f-cloud_mix_n_f-cloud_ice_mix_n_f)
rho_n_f=pressure_n_f/(R_d*tv_n_f)
exner_n_f=(pressure_n_f/p0)**(R_d/cpd)
th_n_f=temp_n_f/exner_n_f
thv_n_f=tv_n_f/exner_n_f

pressure_n_ha=np.zeros([nt,nk])
for iz in np.arange(nk-1,0,-1):
 pressure_n_ha[:,iz]=(pressure_n_f[:,iz-1]+pressure_n_f[:,iz])*0.5
pressure_n_ha[:,0]=pressure_n_f[:,0]

mlheigths_n_h=np.zeros([nt,nk+1])
mlheigths_n_h[:,nk]=HSURF
mlheigths_n_h[:,nk-1]=HSURF+Hmod1
for iz in np.arange(nk-2,-1,-1):
   mlheigths_n_h[:,iz]=mlheigths_n_h[:,iz+1]-np.log(pressure_n_ha[:,iz]/pressure_n_ha[:,iz+1])*R_d/0.5/(1.0/tv_n_f[:,iz]+1.0/tv_n_f[:,iz+1])/g

mlheigths_n_f=np.zeros([nt,nk])
for iz in np.arange(nk-1,0,-1):
 mlheigths_n_f[:,iz]=(mlheigths_n_h[:,iz+1]+mlheigths_n_h[:,iz])*0.5


##interpolate to output levels
u_n_h=np.zeros([nt,hsize])
v_n_h=np.zeros([nt,hsize])
qv_n_h=np.zeros([nt,hsize])
qc_n_h=np.zeros([nt,hsize])
th_n_h=np.zeros([nt,hsize])
exner_n_h=np.zeros([nt,hsize])
thv_n_h=np.zeros([nt,hsize])
rho_n_h=np.zeros([nt,hsize])
pressure_n_h=np.zeros([nt,hsize])
for it in range(0,nt):
  u_n_h[it,:]=np.interp(hlevels+HSURF, mlheigths_n_f[it,nk:0:-1], u_n_f[it,nk:0:-1])
  v_n_h[it,:]=np.interp(hlevels+HSURF, mlheigths_n_f[it,nk:0:-1], v_n_f[it,nk:0:-1])
  th_n_h[it,:]=np.interp(hlevels+HSURF, mlheigths_n_f[it,nk:0:-1], th_n_f[it,nk:0:-1])
  thv_n_h[it,:]=np.interp(hlevels+HSURF, mlheigths_n_f[it,nk:0:-1], thv_n_f[it,nk:0:-1])
  exner_n_h[it,:]=np.interp(hlevels+HSURF, mlheigths_n_f[it,nk:0:-1], exner_n_f[it,nk:0:-1])
  rho_n_h[it,:]=np.interp(hlevels+HSURF, mlheigths_n_f[it,nk:0:-1], rho_n_f[it,nk:0:-1])
  qv_n_h[it,:]=np.interp(hlevels+HSURF, mlheigths_n_f[it,nk:0:-1], qv_n_f[it,nk:0:-1])
  qc_n_h[it,:]=np.interp(hlevels+HSURF, mlheigths_n_f[it,nk:0:-1], qc_n_f[it,nk:0:-1])
  pressure_n_h[it,:]=np.interp(hlevels+HSURF, mlheigths_n_f[it,nk:0:-1], pressure_n_f[it,nk:0:-1])

#forecast forcings
print ("forcing")

#locations:6-east,2-south,4-west,8-north,5-center
locations=np.arange(0,9)
key_indices=([[] for ipl in locations])

#get indices for averaging
key_indices[5]=get_indices(lats_c,lons_c,0.5*delta,lat, lon)
key_indices[8]=get_indices(lats_c+0.5*delta,lons_c,0.5*delta,lat, lon)
key_indices[2]=get_indices(lats_c-0.5*delta,lons_c,0.5*delta,lat, lon)
key_indices[4]=get_indices(lats_c,lons_c-0.5*delta,0.5*delta,lat, lon)
key_indices[6]=get_indices(lats_c,lons_c+0.5*delta,0.5*delta,lat, lon)

time=[]
u_f=np.zeros([9,nt,nk])#first dimension is location:6-east,2-south,4-west,8-north,5-center
v_f=np.zeros([9,nt,nk])
w_f=np.zeros([9,nt,nk+1])
temp_f=np.zeros([9,nt,nk])
pressure_f=np.zeros([9,nt,nk])
spec_hum_f=np.zeros([9,nt,nk])
cloud_mix_f=np.zeros([9,nt,nk])
cloud_ice_mix_f=np.zeros([9,nt,nk])

for it in range(0,nt):
 it0=trange[it]//24
 it1=trange[it]%24
 grib_fcst=nc_fcst_str+"{:02d}".format(it0)+"{:02d}".format(it1)+'0000'
 nc_fcst=nc_fcst_str+"{:02d}".format(it0)+"{:02d}".format(it1)+'0000.nc'
 #grib_fcst=nc_fcst_str+"{:02d}".format(trange[it])+'0000'
 #nc_fcst=nc_fcst_str+"{:02d}".format(trange[it])+'0000.nc'
 #os.system('cdo -f nc  copy '+inpath+'/'+grib_fcst+' '+tmppath+'/'+nc_fcst)
 #north
 for iloc in [2,4,5,6,8]:
   key_indices_s=",".join(map(str, key_indices[iloc]))
   os.system('cdo -f nc selgridcell,'+key_indices_s+' '+inpath+'/'+grib_fcst+' '+tmppath+'/'+nc_fcst)
   nc_fcst_x=xr.open_dataset(tmppath+'/'+nc_fcst)
   if(iloc==2):
     time.append(nc_fcst_x.time.values)

   u_f[iloc,it,:]=nc_fcst_x.u[0,:,:].mean(axis=1).values
   v_f[iloc,it,:]=nc_fcst_x.v[0,:,:].mean(axis=1).values
   w_f[iloc,it,:]=nc_fcst_x.wz[0,:,:].mean(axis=1).values
   temp_f[iloc,it,:]=nc_fcst_x.t[0,:,:].mean(axis=1).values
   pressure_f[iloc,it,:]=nc_fcst_x.pres[0,:,:].mean(axis=1).values
   spec_hum_f[iloc,it,:]=nc_fcst_x.q[0,:,:].mean(axis=1).values
   cloud_mix_f[iloc,it,:]=nc_fcst_x.clwmr[0,:,:].mean(axis=1).values
   cloud_ice_mix_f[iloc,it,:]=nc_fcst_x.QI[0,:,:].mean(axis=1).values
   #rain_mix_f[iloc]=nc_fcst_x.rwmr[0,:,:].mean(axis=1).values
   #snow_mix_f[iloc]=nc_fcst_x.snmr[0,:,:].mean(axis=1).values
   nc_fcst_x.close()
   os.system('rm '+tmppath+'/'+nc_fcst)

#get heigths from initial files
HSURFs=np.zeros([9])
for iloc in [2,4,5,6,8]:
 grib_init=nc_init_str+datetime_str
 nc_init=grib_init+'.nc'
 key_indices_s=",".join(map(str, key_indices[iloc]))
 os.system('cdo -f nc selgridcell,'+key_indices_s+' '+inpath+'/'+grib_init+' '+tmppath+'/'+nc_init)
 nc_init_x=xr.open_dataset(tmppath+'/'+nc_init)
 HSURFs[iloc]=nc_init_x.HHL[0,nk,:].mean().values
 nc_init_x.close()
 os.system('rm '+tmppath+'/'+nc_init)

print ("compute heigths of model levels")
qt_f = (spec_hum_f + cloud_mix_f + cloud_ice_mix_f)/\
 (1.0+cloud_mix_f + cloud_ice_mix_f)
qc_f=(cloud_mix_f)*(1.0-qt_f)
r_mix_f=spec_hum_f/(1.0-qt_f)
tv_f=temp_f*(1.0+0.608*r_mix_f-cloud_mix_f-cloud_ice_mix_f)
rho_f=pressure_f/(R_d*tv_f)
th_f=temp_f*(pressure_f/p0)**(-R_d/cpd)

pressure_h=np.zeros([9,nt,nk])
for iz in np.arange(nk-1,0,-1):
 pressure_h[:,:,iz]=(pressure_f[:,:,iz-1]+pressure_f[:,:,iz])*0.5
pressure_h[:,:,0]=pressure_f[:,:,0]

mlheigths_h=np.zeros([9,nt,nk+1])
for iloc in [2,4,5,6,8]:
 mlheigths_h[iloc,:,nk]=HSURFs[iloc]
 mlheigths_h[iloc,:,nk-1]=HSURFs[iloc]+Hmod1
 for iz in np.arange(nk-2,-1,-1):
   mlheigths_h[iloc,:,iz]=mlheigths_h[iloc,:,iz+1]-np.log(pressure_h[iloc,:,iz]/pressure_h[iloc,:,iz+1])*R_d/0.5/(1.0/tv_f[iloc,:,iz]+1.0/tv_f[iloc,:,iz+1])/g

mlheigths_f=np.zeros([9,nt,nk])
for iz in np.arange(nk-1,0,-1):
 mlheigths_f[:,:,iz]=(mlheigths_h[:,:,iz+1]+mlheigths_h[:,:,iz])*0.5

#compute geostrophic wind
ug_f=np.zeros([nt,nk])
vg_f=np.zeros([nt,nk])
for it in range(0,nt):
 ug_f[it,:]= ((1.0/(rho_f[5,it,:]))*(pressure_f[8,it,:]-pressure_f[2,it,:])+g*(mlheigths_f[8,it,:]-mlheigths_f[2,it,:]))/d_y/f3
 vg_f[it,:]= -((1.0/(rho_f[5,it,:]))*(pressure_f[6,it,:]-pressure_f[4,it,:])+g*(mlheigths_f[6,it,:]-mlheigths_f[4,it,:]))/d_x/f3

#compute LS tendencies
ddt_qt_f=np.zeros([nt,nk])
ddt_th_f=np.zeros([nt,nk])
ddt_u_f=np.zeros([nt,nk])
ddt_v_f=np.zeros([nt,nk])
for it in range(0,nt):
  ddt_qt_f[it,:]=-(u_f[5,it,:]*(qt_f[6,it,:]-qt_f[4,it,:])/d_x+v_f[5,it,:]*(qt_f[8,it,:]-qt_f[2,it,:])/d_y)
  ddt_th_f[it,:]=-(u_f[5,it,:]*(th_f[6,it,:]-th_f[4,it,:])/d_x+v_f[5,it,:]*(th_f[8,it,:]-th_f[2,it,:])/d_y)
  ddt_u_f[it,:]=-(u_f[5,it,:]*(u_f[6,it,:]-u_f[4,it,:])/d_x+v_f[5,it,:]*(u_f[8,it,:]-u_f[2,it,:])/d_y)
  ddt_v_f[it,:]=-(u_f[5,it,:]*(v_f[6,it,:]-v_f[4,it,:])/d_x+v_f[5,it,:]*(v_f[8,it,:]-v_f[2,it,:])/d_y)

#interpolate to fixed output heigths
u_h=np.zeros([9,nt,hsize])#first dimension is location:6-east,2-south,4-west,8-north,5-center
v_h=np.zeros([9,nt,hsize])
w_h=np.zeros([9,nt,hsize])
temp_h=np.zeros([9,nt,hsize])
pressure_h=np.zeros([9,nt,hsize])
spec_hum_h=np.zeros([9,nt,hsize])
cloud_mix_h=np.zeros([9,nt,hsize])
cloud_ice_mix_h=np.zeros([9,nt,hsize])
v_h=np.zeros([9,nt,hsize])
rho_h=np.zeros([9,nt,hsize])
ug_h=np.zeros([nt,hsize])
vg_h=np.zeros([nt,hsize])
ddt_qt_h=np.zeros([nt,hsize])
ddt_th_h=np.zeros([nt,hsize])
ddt_u_h=np.zeros([nt,hsize])
ddt_v_h=np.zeros([nt,hsize])
for it in range(0,nt):
 for iloc in [2,4,5,6,8]:
  u_h[iloc,it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], u_f[iloc,it,nk:0:-1])
  v_h[iloc,it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], v_f[iloc,it,nk:0:-1])
  w_h[iloc,it,:]=np.interp(hlevels+HSURF, mlheigths_h[iloc,it,nk+1:0:-1], w_f[iloc,it,nk+1:0:-1])
  temp_h[iloc,it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], temp_f[iloc,it,nk:0:-1])
  pressure_h[iloc,it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], pressure_f[iloc,it,nk:0:-1])
  spec_hum_h[iloc,it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], spec_hum_f[iloc,it,nk:0:-1])
  cloud_mix_h[iloc,it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], cloud_mix_f[iloc,it,nk:0:-1])
  cloud_ice_mix_h[iloc,it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], cloud_ice_mix_f[iloc,it,nk:0:-1])
  rho_h[iloc,it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], rho_f[iloc,it,nk:0:-1])
  ug_h[it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], ug_f[it,nk:0:-1])
  vg_h[it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], vg_f[it,nk:0:-1])
  ddt_qt_h[it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], ddt_qt_f[it,nk:0:-1])
  ddt_th_h[it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], ddt_th_f[it,nk:0:-1])
  ddt_u_h[it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], ddt_u_f[it,nk:0:-1])
  ddt_v_h[it,:]=np.interp(hlevels+HSURF, mlheigths_f[iloc,it,nk:0:-1], ddt_v_f[it,nk:0:-1])

#-------------------------------------------------------------------
# write data into netcdf file

print ("write to netCDF file ", file_out)
levs=np.arange(0,hlevels.size,dtype='int32')
do=Dataset(file_out,mode='w',format='NETCDF4_CLASSIC')
         
do.description='Initial profiles and forcings generated from COSMO data'
do.createDimension('lev',size=levs.size)
do.createDimension('lev1'     ,size=levs.size)#could be changed
do.createDimension('nt', size=trange.size)
do.createDimension('levTsoil',size=depth_T.size)
do.createDimension('levWsoil',size=depth_W.size)
do.createDimension('nclass_lu',size=LU_CLASS_FRACTION.size)

#--- coordinates

#time
nc_time=do.createVariable('time',np.float64,'nt')
nc_time.units='s since'

#level
lev=do.createVariable('lev',np.float64,'lev')
lev.units='m'
lev.standard_name='height'
lev[:]=levs
#level - half level
lev1=do.createVariable('lev1',np.float64,'lev1')
lev1.units='m'
lev1.standard_name='height'
lev1[:]=levs#can be changes

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

##soil variables
nc_height[:]=np.tile(hlevels,trange.size)
nc_height_ifc[:] = nc_height[:]#can be changed
nc_time[:]=((time-time[0]).astype('int')/1e9)[:,0]
nc_uIN[0,:]=u_in_h
nc_uIN[1:,:]=u_n_h[1:,:]
nc_uIN[:,0]=nc_uIN[:,1]#no gradient of u in stratosphere
nc_vIN[0,:]=v_in_h
nc_vIN[1:,:]=v_n_h[1:,:]
nc_vIN[:,0]=nc_vIN[:,1]#no gradient of u in stratosphere
nc_wIN[0,:]=w_in_h
nc_tkeIN[0,:]=tke_in_h
nc_thIN[0,:]=th_in_h
nc_thIN[1:,:]=th_n_h[1:,:]
nc_thIN[:,0]=nc_thIN[:,1]+dtdz_str*(nc_height[:,0]-nc_height[:,1])#statospheric gradient of pot. temperature, levels above are extrapolated according to this gradinet in SCM
nc_thvIN[0,:]=thv_in_h
nc_thvIN[1:,:]=thv_n_h[1:,:]
nc_thvIN[:,0]=nc_thvIN[:,1]+dtdz_str*(nc_height[:,0]-nc_height[:,1])#statospheric gradient of pot. temperature, levels above are extrapolated according to this gradinet in SCM
nc_qvIN[0,:]=qv_in_h
nc_qvIN[1:,:]=qv_n_h[1:,:]
nc_qvIN[:,0]=nc_qvIN[:,1] #no gradient of moisture in stratosphere
nc_qcIN[0,:]=qc_in_h
nc_qcIN[1:,:]=qc_n_h[1:,:]
nc_exnerIN[0,:]=exner_in_h
nc_exnerIN[1:,:]=exner_n_h[1:,:]
nc_rhoIN[0,:]=rho_in_h
nc_rhoIN[1:,:]=rho_n_h[1:,:]
nc_o3IN[0,:]=exner_in_h*0 #missing here, but should be available
nc_o3IN[1:,:]=rho_n_h[1:,:]*0#missing here, but should be available
#data for nudging

##soil variables
nc_T_SO=do.createVariable('T_SO',np.float64,('nt','levTsoil'))
nc_T_SO.description='soil temperature'
nc_T_SO.units='K'
nc_W_SO=do.createVariable('W_SO',np.float64,('nt','levWsoil'))
nc_W_SO.description='water content'
nc_W_SO.units='kg m-2'
#values
nc_T_SO[:]=np.tile(soil_temp,trange.size)
nc_W_SO[:]=np.tile(soil_moisture,trange.size)

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
nc_uGEO[:]=ug_h
nc_vGEO[:]=vg_h
nc_wLS[:]=w_h[5,:,:]
nc_wLS[:,0]=nc_wLS[:,1]#no gradient in stratosphere
nc_dTadv[:]=ddt_th_h
nc_dTadv[:,0]=nc_dTadv[:,1]#no gradient in stratosphere
nc_dTrad[:]=nc_dTadv[:,levs]*0.0 # setting radiation tendency to 0
nc_dQVadv[:]=ddt_qt_h
nc_dQVadv[:,0]=nc_dQVadv[:,1]#no gradient in stratosphere
nc_dUadv[:]=ddt_u_h
nc_dVadv[:]=ddt_v_h

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
nc_sfc_lat_flx[:]=0.0#setting to zero for now
nc_sfc_sens_flx[:]=0.0#setting to zero for now
nc_psurf[:]=pressure_n_ha[:,nk-1]
nc_ts[:]=t_g #should there be a different temperature?
nc_tg[:]=t_g
#nc_qvs[:]=q_vh[:,levs[-1],key_index_fi,key_index_fj]##no better input for now - input from surface file ?
nc_qvs[:]=qv_s
nc_Ch[:]=0.0#no input
nc_Cq[:]=0.0#no input
nc_Cm[:]=0.0#no input
nc_ustar[:]=0.0#no input

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
nc_LU_CLASS_FRACTION.lctype='GLOBCOVER2009'#should be addapted when needed

nc_longitude[:]=lon
nc_latitude[:]=lat
nc_FR_LAND[:]=FR_LAND
nc_PLCOV_MX[:]=PLCOV_MX
nc_LAI_MX[:]=LAI_MX
nc_ROOTDP[:]=ROOTDP
nc_RSMIN[:]=RSMIN
nc_SOILTYP[:]=SOILTYP
nc_Z0[:]=Z0
nc_EMIS_RAD[:]=EMIS_RAD
nc_TOPO[:]=HSURF
nc_LU_CLASS_FRACTION[:]=LU_CLASS_FRACTION

today = datetime.today()
do.history = "Created " + today.strftime("%d/%m/%y")

do.close()


os.system('rm -r ' + tmppath)        
