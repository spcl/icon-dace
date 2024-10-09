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

# prepare SCM/LES input and forcing data for the CKD-MIP
# radiation core comparison from a netCDF file
#
# runs as: python3 get_CKD-MIP.py
#
# select location: (on LCE)
#   module load netcdf nco
#   ncks -d nt,0,0 init_SCM_CKD-MIP.nc init_SCM_CKD-MIP-temp.nc
#
# Martin Koehler, 2020-04
#-------------------------------------------------------------------

import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import sys, getopt, os       # file handling and operators
from datetime import datetime

#-------------------------------------------------------------------
# ncdump from https://docs.python.org

def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        try:
            print("\t\ttype:", repr(nc_fid.variables[key].dtype) )
            for ncattr in nc_fid.variables[key].ncattrs():
                print('\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
            print( "\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print( "NetCDF Global Attributes:" )
        for nc_attr in nc_attrs:
            print( '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)) )
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print( "NetCDF dimension information:" )
        for dim in nc_dims:
            print( "\tName:", dim ) 
            print( "\t\tsize:", len(nc_fid.dimensions[dim]) )
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print( "NetCDF variable information:" )
        for var in nc_vars:
            if var not in nc_dims:
                print( '\tName:', var )
                print( "\t\tdimensions:", nc_fid.variables[var].dimensions )
                print( "\t\tsize:", nc_fid.variables[var].size )
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars

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

#-------------------------------------------------------------------
# setup
   
print ("Beginning of main -------------------------------------------")

# constants
R_d   = 287.06            # dry air gas constant
R_v   = 461.56            # water vapour gas constant
L_v   = 2.501e6           # latent heat of vaporization
g     = 9.81              # grav. acceleration
p0    = 100000.           # reference atm. pressure [Pa]
T0    = 273.15            # zero Celsius [K]
cpd   = 1004.7            # the heat capacity of dry air at constant pressure
dtdz_str  = 0.04          # stratospheric gradient of potential temperature [K/m]
M_dry_air = 0.0289654     # molar mass of dry air [kg/mol]
M_h2o     = 0.018016      # molar mass of water (vapour) [kg/mol] 
M_o3      = 0.047997      # molar mass of ozone [kg/mol]

icon_dir = '../../../'
data_dir = 'data/scm/init_data/'
file_in  = icon_dir + data_dir + 'ckdmip_evaluation1_concentrations.nc'
file_out = icon_dir + data_dir + 'init_SCM_CKD-MIP.nc'

# read input NetCDF file

print ("... read from netCDF file \n ", file_in)

nc_fid = Dataset(file_in, 'r')  # Dataset is the class behavior to open the file
                                # and create an instance of the ncCDF4 class
nc_attrs, nc_dims, nc_vars = ncdump(nc_fid,verb=False  )

# Extract data from NetCDF file
latitude  = nc_fid.variables['latitude'][:]
longitude = nc_fid.variables['longitude'][:]
time      = nc_fid.variables['time'][:]
t_units   = nc_fid.variables['time'].units
ntime     = len(time)

# --- profile variables: shape is time, lev

# model levels
plevels   = nc_fid.variables['pressure_fl'][:,:]
nlev      = np.shape(plevels)[1]
#p_surf   = nc_fid.variables['Ps'][:,0,0] # surface pressure [Pa]
p_surf    = nc_fid.variables['pressure_hl'][:,nlev] # surface pressure [Pa]
exns      = (p_surf/p0)**(R_d/cpd)#Exner function at surface

#find only relevant vertical levels
#for ik in range(0,nnlev):
#   if(np.where(p_surf>plevels[ik])[0].size==ntime):
#      ilev=ik

#plevels  = plevels   [0:ilev+1] 
#nlev     = plevels.size
exn       =(plevels/p0)**(R_d/cpd) #Exner function

# temperature
temp      = nc_fid.variables['temperature_fl'][:,:]    # T
thl       = temp/exn    # THETA

# specific humidity
h2o_mole_frac = nc_fid.variables['h2o_mole_fraction_fl'][:,:]
qt            = (M_dry_air/(h2o_mole_frac*M_h2o)-1)**(-1)  # q

# ozone
o3_mole_frac = nc_fid.variables['o3_mole_fraction_fl'][:,:]
o3_mmr       = o3_mole_frac*M_o3/M_dry_air

tempv     = temp*(1+0.61*qt)                       # virtual temperature
thv       = temp/exn*(1+0.61*qt)                   # virtual pot. temperature
dens      = plevels/(R_d*tempv)

# u-wind and v-wind components
#u        = nc_fid.variables['u'][:,0:nlev,0,0]  
#v        = nc_fid.variables['v'][:,0:nlev,0,0]
u         = 0.0 * qt
v         = 0.0 * qt

# large scale vertical velocity
#omega     = nc_fid.variables['omega'][:,0:nlev,0,0]   # [Pa/s]
omega     = 0.0 * qt
wls       = -omega/dens/g   # vertical velocity [m/s]

#horizontal advection of u/v - not solved
#HdivV    = nc_fid.variables['div'][:,0:nlev,0,0]     # horizontal wind divergence
HdivV     = 0.0 * qt                                  # horizontal wind divergence

# TKE
tke       = 0.0 * qt                       # 0 assumed

# geostrophic u-wind
ug        = 0.0 * qt                       # 0 assumed

#T advection [K/s]
# T advection (add horizontal and vertical) # ATTENTION ????, T but not THETA ????
thlls     = 0.0 * qt

#Q advection [kg/kg/s]
# Q advection (add horizontal and vertical) # ATTENTION ????
qtls      = 0.0 * qt

print('... profile variables:  ntime: ',ntime, '  nlev: ', nlev, '  u-shape: ', np.shape(u))


# --- surface variables: shape is time

# surface fluxes
#sfc_lat_flx  = nc_fid.variables['LH'][:]   # W/m2
#sfc_lat_flx  = nc_fid.variables['lhflx'][:,0,0]   # W/m2
#sfc_sens_flx = nc_fid.variables['SH'][:]   # W/m2
#sfc_sens_flx = nc_fid.variables['shflx'][:,0,0]   # W/m2
sfc_lat_flx  = 0.0 * time
sfc_sens_flx = 0.0 * time

print('surface variables:  ntime: ',ntime, '  p_surf-shape: ', np.shape(p_surf))

# normalize profiles to SI
#qtls   = qtls  / 3600. / 1000.  # from g/kg/h to kg/kg/s
#thlls  = thlls / 3600.          # from K/h to K/s
#p_surf = p_surf * 100.0         # from hPa to Pa
#qt     = qt / 1000.0            # from g/kg to kg/kg
#wls   /= 100.          # omega!!        ????

# surface parameters
#ustar       = 0.28     # m/s            ????
z0           = 0.1      # m              ?forest?
dt_relax     = 7200.0   # s              ????

#temps       = nc_fid.variables['Tsair'][:,0,0]  # surface air temperature
temps        = nc_fid.variables['temperature_hl'][:,nlev]  # surface air temperature
tg           = temps                                       # ground temperature
#rvs          = nc_fid.variables["qs"][:,0,0]
rvs          = 0.0 * time
qvs          = rvs/(1.+rvs)
tempvs       = temps*(1.0+0.61*qvs)
thvs         = temps/exns
denss        = p_surf/R_d/tempvs
#ws          = nc_fid.variables["wspd_srf"][:,0,0]
ws           = 0.0 * time

bflx=(sfc_sens_flx/cpd+thvs*0.61*sfc_lat_flx/L_v)/denss*g/thvs
ustar=diag_ustar(10.0,bflx,ws,z0)

hlevels      = qt*0   # ATTENTION: height levels - need to be computed
#hlevels[:,nlev-1]=nc_fid.variables['phis'][0,0]/g-0.5*(temps[:]+tempv[:,nlev-1])*R_d/g*np.log(plevels[nlev-1]/p_surf)
#first level above surface
hlevels[:,nlev-1]=-0.5*(temps[:]+tempv[:,nlev-1])*R_d/g*np.log(plevels[:,nlev-1]/p_surf)
for ik in range(nlev-2,-1,-1):
   hlevels[:,ik]=hlevels[:,ik+1]-0.5*(tempv[:,ik+1]+tempv[:,ik])*R_d/g*np.log(plevels[:,ik]/plevels[:,ik+1])

topo         = 0.0


#-------------------------------------------------------------------
# write data into netcdf file

print ("... write to netCDF file", file_out)
do=Dataset(file_out,mode='w',format='NETCDF4_CLASSIC')
         
do.description='Initial profiles and forcings for CKD-MIP case'
do.createDimension('lev'      ,size=nlev)
do.createDimension('nt'       ,None)
do.createDimension('levTsoil' ,size=1)
do.createDimension('levWsoil' ,size=1)
do.createDimension('nclass_lu',size=23)

#--- coordinates

#time
nc_time=do.createVariable('time',np.float64,'nt')
nc_time.units=t_units

#level
lev=do.createVariable('lev',np.float64,'lev')
lev.units='m'
lev.standard_name='height'
lev[:]=np.arange(nlev)

#soil levels
levTsoil=do.createVariable('levTsoil',np.float64,'levTsoil')
levTsoil.units='level'
levTsoil[:]=1

levWsoil=do.createVariable('levWsoil',np.float64,'levWsoil')
levWsoil.units='level'
levWsoil[:]=1

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
nc_uIN              = do.createVariable('uIN',np.float64,('nt','lev'))
nc_uIN.description  = 'zonal wind'
nc_uIN.units        = 'm/s'
nc_vIN              = do.createVariable('vIN',np.float64,('nt','lev'))
nc_vIN.description  = 'meridional wind'
nc_vIN.units        = 'm/s'
nc_wIN              = do.createVariable('wIN',np.float64,('nt','lev'))
nc_wIN.description  = 'veritcal velocity'
nc_wIN.units        = 'm/s'
nc_thIN             = do.createVariable('thIN',np.float64,('nt','lev'))
nc_thIN.description = 'potential temperature'
nc_thIN.units       = 'K'
nc_qvIN             = do.createVariable('qvIN',np.float64,('nt','lev'))
nc_qvIN.description = 'water vapour specific moisture'
nc_qvIN.units       = 'kg/kg'
nc_qcIN             = do.createVariable('qcIN',np.float64,('nt','lev'))
nc_qcIN.description = 'liquid water specific moisture'
nc_qcIN.units       = 'kg/kg'
nc_o3IN             = do.createVariable('o3IN',np.float64,('nt','lev'))
nc_o3IN.description = 'ozone mass mixing ratio'
nc_o3IN.units       = 'kg/kg'
nc_pIN              = do.createVariable('pIN',np.float64,('nt','lev'))
nc_pIN.description  = 'atmospheric pressure'
nc_pIN.units        = 'Pa'
nc_tkeIN            = do.createVariable('tkeIN',np.float64,('nt','lev'))
nc_tkeIN.description= 'turbulence kinetic energy'
nc_tkeIN.units      = 'm^2/s^2'
###values
nc_height[:] = hlevels
nc_time[:]   = time
nc_uIN[:]    = u
nc_uIN[:,0]  = nc_uIN[:,1]   # no gradient of u in stratosphere
nc_vIN[:]    = v
nc_vIN[:,0]  = nc_vIN[:,1]   # no gradient of u in stratosphere
nc_wIN[:]    = 0.0
nc_tkeIN[:]  = tke
nc_thIN[:]   = thl
nc_thIN[:,0] = nc_thIN[:,1]+dtdz_str*(nc_height[:,0]-nc_height[:,1])#statospheric gradient of pot. temperature, levels above are extrapolated according to this gradinet in SCM
nc_qvIN[:]   = qt
nc_qvIN[:,0] = nc_qvIN[:,1]  # no gradient of moisture in stratosphere
nc_qcIN[:]   = 0.0
nc_pIN[:]    = plevels
#ozone
nc_o3IN[:]   = o3_mmr

##soil variables
nc_T_SO=do.createVariable('T_SO',np.float64,('nt','levTsoil'))
nc_T_SO.description='soil temperature'
nc_T_SO.units='K'
nc_W_SO=do.createVariable('W_SO',np.float64,('nt','levWsoil'))
nc_W_SO.description='water content'
nc_W_SO.units='kg m-2'
#values
nc_T_SO[:] = 0.0
nc_W_SO[:] = 0.0

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
nc_uGEO[:]     = ug
nc_vGEO[:]     = 0.0
nc_wLS[:]      = wls
nc_wLS[:,0]    = nc_wLS[:,1]#no gradient in stratosphere
nc_dTadv[:]    = thlls
nc_dTadv[:,0]  = nc_dTadv[:,1]#no gradient in stratosphere
nc_dTrad[:]    = nc_dTadv[:]*0.0 # setting radiation tendency to 0
nc_dQVadv[:]   = qtls
nc_dQVadv[:,0] = nc_dQVadv[:,1]#no gradient in stratosphere
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
nc_sfc_lat_flx[:]  = sfc_lat_flx
nc_sfc_sens_flx[:] = sfc_sens_flx
nc_psurf[:]        = p_surf
nc_ts[:]           = temps
nc_tg[:]           = tg
nc_qvs[:]          = qvs
nc_Ch[:]           = 0.0
nc_Cq[:]           = 0.0
nc_Cm[:]           = 0.0
nc_ustar[:]        = ustar


#--- constants

#nudging parameters
nc_dt_relax=do.createVariable('dt_relax',np.float64)
nc_dt_relax.description='relaxation time step for nudging'
nc_dt_relax.units='s'
##values
nc_dt_relax[:]=dt_relax

#external parameters
nc_longitude=do.createVariable('longitude',np.float64,('nt'))
nc_longitude.description='longitude'
nc_longitude.units='degrees'
nc_latitude=do.createVariable('latitude',np.float64,('nt'))
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
nc_longitude[:]  = longitude
nc_latitude[:]   = latitude
nc_FR_LAND[:]    = 1.0
nc_PLCOV_MX[:]   = 0.0
nc_LAI_MX[:]     = 0.0
nc_ROOTDP[:]     = 0.0
nc_RSMIN[:]      = 0.0
nc_SOILTYP[:]    = 8
nc_Z0[:]         = z0
nc_EMIS_RAD[:]   = 0.98
nc_TOPO[:]       = topo
nc_LU_CLASS_FRACTION[:] = np.zeros([23])
nc_LU_CLASS_FRACTION[5] = 1.0      # forrest : Closed (>40%) broadleaved deciduous forest (>5m)

today = datetime.today()
do.history = "Created " + today.strftime("%d/%m/%y")

do.close()
