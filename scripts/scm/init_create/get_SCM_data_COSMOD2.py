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

# prepare SCM/LES input and forcing data for the GASS DCP 
# (diurnal cycle of cloud and precipitation) case in netCDF file
#
# runs as: python3 get_GASS_DCP_GoAmazon.py
#
# Martin Koehler, 2019-11
#
# first version: June 2017 by Markus Ernst
# small adaptions: August 2017 by Annika Schomburg
#   09 Sept 2017: use liquid water potential temperature instead of 
#                 dry potential temperature
#   25 Oct 2017: A.S: read also topography_c (height above sea level) 
#                     from external parameter file and write out together with
#                     other constant data
# adapted for COSMO analyzes input: Aug 2018 by Ivan Bastak Duran
# basic usage from the command-line: python3 generate_SCM_data.py
#-------------------------------------------------------------------

import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import sys, getopt, os       # file handling and operators
from datetime import datetime

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

def get_indices(lats,lons,LATX = 52.209722, LONX = 14.118889):
    """ get_indices takes latitude latx and longitude lonx and
    returns the indices indexi and indexj of the nearest grid cell in a 2D array"""
    DISTANCE_TO_X = np.sqrt((lats - LATX)**2 + (lons - LONX)**2)
    minimum_indices = np.where(DISTANCE_TO_X == DISTANCE_TO_X.min())
    indexi, indexj = minimum_indices[0][0], minimum_indices[1][0]
    return indexi, indexj
    pass

def get_avg_fields(field, weigths,lats,lons,ix, jx, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright, latx = 52.209722, lonx = 14.118889, delta = 2):
    """get_avg_fields takes an imported COSMO field, latitude and longitude of a
    central point x as well as a delta in degrees and returns height profiles averaged
    in longitude and latitude for location x and 4 surrounding points"""
    avg_left = np.ma.average(field[:,:,ix_down:ix_up,jx_leftleft:jx], axis=(2,3),weights=weigths[:,:,ix_down:ix_up,jx_leftleft:jx])
    avg_up = np.ma.average(field[:,:,ix:ix_upup,jx_left:jx_right], axis=(2,3),weights=weigths[:,:,ix:ix_upup,jx_left:jx_right])
    avg_right = np.ma.average(field[:,:,ix_down:ix_up,jx:jx_rightright], axis=(2,3),weights=weigths[:,:,ix_down:ix_up,jx:jx_rightright])
    avg_down = np.ma.average(field[:,:,ix_downdown:ix,jx_left:jx_right], axis=(2,3),weights=weigths[:,:,ix_downdown:ix,jx_left:jx_right])
    avg_centre = np.ma.average(field[:,:,ix_down:ix_up,jx_left:jx_right], axis=(2,3),weights=weigths[:,:,ix_down:ix_up,jx_left:jx_right])
    return avg_left, avg_up, avg_right, avg_down, avg_centre
    pass

def get_avg_fields_s(field,lats,lons,ix, jx, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright, latx = 52.209722, lonx = 14.118889, delta = 2):
    """get_avg_fields_s takes an imported COSMO field, latitude and longitude of a
    central point x as well as a delta in degrees and returns height profiles averaged
    in longitude and latitude for location x and 4 surrounding points"""
    avg_left = np.average(field[ix_down:ix_up,jx_leftleft:jx])
    avg_up = np.average(field[ix:ix_upup,jx_left:jx_right])
    avg_right = np.average(field[ix_down:ix_up,jx:jx_rightright])
    avg_down = np.average(field[ix_downdown:ix,jx_left:jx_right])
    avg_centre = np.average(field[ix_down:ix_up,jx_left:jx_right])
    return avg_left, avg_up, avg_right, avg_down, avg_centre
    pass

def get_surrounding_points(lats,lons,latx = 52.209722, lonx = 14.118889, delta = 2):
    """get_surrounding_points takes latitude latx and longitude lonx of a central
    point x and returns the needed indices for x and the 4 surrounding points"""
    ix, jx = get_indices(lats,lons,latx, lonx)
    ix_left, jx_left = get_indices(lats,lons,latx, lonx - (delta/2.))
    ix_up, jx_up =get_indices(lats,lons,latx + (delta/2.), lonx)
    ix_right, jx_right =get_indices(lats,lons,latx, lonx + (delta/2.))
    ix_down, jx_down =get_indices(lats,lons,latx - (delta/2.), lonx)

    ix_leftleft, jx_leftleft = get_indices(lats,lons,latx, lonx - (delta))
    ix_upup, jx_upup =get_indices(lats,lons,latx + (delta), lonx)
    ix_rightright, jx_rightright =get_indices(lats,lons,latx, lonx + (delta))
    ix_downdown, jx_downdown =get_indices(lats,lons,latx - (delta), lonx)

    return ix, jx, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright
    pass

def calculate_liquid_potential_temperature(T,P,qc):
    """calculate_potential_temperature takes temperature fields T and pressure
    fields P to calculate and return potential temperature fields for the 
    whole domain. This is afterwards used to compute the liquid water potential
    temperature"""
    RDCP = 0.286
    Lv = 2.5008*10**6
    Cp = 1004.64
    P0 = np.ones_like(P) * 100000.
    POT_T = T * (P0/P)**(RDCP)
    POT_T_l = POT_T - Lv/Cp*(POT_T/T)*(qc)
    return POT_T_l,POT_T
    pass

def calculate_tendencies(LSGRADIENT_X1, LSGRADIENT_X2, U_LS1, U_LS2, U_LS3, LSGRADIENT_X3 = 0):
    """calculate_tendencies takes two LSGradients east-west, north-south LSGradient_x1, LSGradient_x2 
    and three components of the large scale wind velocity U_LS1, U_LS2, U_LS3 and returns the LS horizonzal
    advection LSA and the KS vertical advection SUB and the applied nudging tendency NUD"""
    LSA = - (U_LS1 * LSGRADIENT_X1 + U_LS2 * LSGRADIENT_X2)
    SUB =  - U_LS3 * LSGRADIENT_X3
    # NUD = - ((AVG(LSQUANTITY)) - LSQUANTITY)/TAU
    return LSA, SUB
    PASS

def interpolate_levels(LEVELS):
    """interpolate_levels takes an k height-level array LEVELS with values on the edges of the grid and interpolates them to the middle
    the function returns a k-1 height level array int_levels """
    int_levels = np.zeros(len(LEVELS) - 1)
    LEVELS = np.append(LEVELS, 0)
    for i in range(len(int_levels)):
        int_levels[i] = (LEVELS[i+1] + LEVELS[i]) / 2.
        #ibd check half/full levels ordering
    return int_levels
    pass

def interpolate_pressure(pressure_profile_at_ij, height_levels_at_ij, surface_pressure_at_ij, final_height_levels):
    """interpolate_pressure takes a pressure profile at a specific point, the height levels at that specific point, the surface pressure at that specific point
    and an arbitrary array of height levels. It then interpolates the pressure to that given height levels and gives them back"""
    # create new array with the logarithm of the pressure profile
    logp = np.log(pressure_profile_at_ij)
    # create new array for the goal pressures
    interpolated_pressure = np.zeros(final_height_levels.shape[0])
    # for loop that goes over all goal_heights
    for i in range(final_height_levels.shape[0]):
        # if height lower than lowest layer use surface pressure
        if final_height_levels[i] < height_levels_at_ij[-1]:
            interpolated_pressure[i] = surface_pressure_at_ij
        else:
            # find closest vertical layer lower than the goal
            height_index = find_lnearest_index(height_levels_at_ij, final_height_levels[i])
            # steigung bestimmen m
            m = (logp[height_index - 1] - logp[height_index])/(height_levels_at_ij[height_index - 1] - height_levels_at_ij[height_index])
            # achsenabschnitt bestimmen, find b
            b = logp[height_index] - m * height_levels_at_ij[height_index]
            # interpolieren
            # take the exponential to convert back
            interpolated_pressure[i] = np.exp(m * final_height_levels[i] + b)
    return interpolated_pressure
    pass

def interpolate_linearly(profile_at_ij, height_levels_at_ij, surface_data_at_ij, final_height_levels):
    """interpolate_pressure takes a pressure profile at a specific point, the height levels at that specific point, the surface pressure at that specific point
    and an arbitrary array of height levels. It then interpolates the pressure to that given height levels and gives them back"""
    # create new array for the goal pressures
    interpolated_profile = np.zeros(final_height_levels.shape[0])
    # for loop that goes over all goal_heights
    for i in range(final_height_levels.shape[0]):
        # if height lower than lowest layer use surface data
        if final_height_levels[i] < height_levels_at_ij[-1]:
            interpolated_profile[i] = surface_data_at_ij
        else:
            # find closest vertical layer lower than the goal
            height_index = find_lnearest_index(height_levels_at_ij, final_height_levels[i])
            # steigung bestimmen m
            m = (profile_at_ij[height_index - 1] - profile_at_ij[height_index])/(height_levels_at_ij[height_index - 1] - height_levels_at_ij[height_index])
            # achsenabschnitt bestimmen, find b
            b = profile_at_ij[height_index] - m * height_levels_at_ij[height_index]
            # interpolieren
            interpolated_profile[i] = m * final_height_levels[i] + b
    return interpolated_profile
    pass


def find_lnearest_index(array,value):
    """find_lnearest_index is a helper function that takes an array and a value as arguments
    and returns the index of the array that produces the lower-nearest value to the given value"""
    idx = (np.abs(array-value)).argmin()
    if array[idx] > value:
        return idx+1
    else:
        return idx
    pass


#-------------------------------------------------------------------
# setup
#cosmo_path="./"
#cosmo_path="/home/bastakdu/Documents/HErZ/Annika/"
icon_dir = '../../../'
data_dir = 'data/scm/init_data/'
file_out = icon_dir + data_dir + "init_SCM2019082600dt2.nc"

##setup
# Lindenberg
lat = 52.209722
lon = 14.118889
damph=5000.#uGEO and vGEO will converge tu u and v at the location --experimental
damph_slope=2e-4  #--experimental
delta = 2.0  #size of averaging boxes in degree
hlevels=np.array([0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3500,4000,4500,5000,5500,6000,7500,8000,8500,9000,9500,10000,20000])
hlevels=hlevels[hlevels.size::-1] #inverse order like model levels ?
dt_relax     = 7200.0   # s- relaxation time scale
   
print ("Beginning of main")
#constants
nlev      = hlevels.size
f2=2.0 * Omega * np.cos(2 * np.pi * lat / 360.) 
f3 = 2.0 * Omega * np.sin(2 * np.pi * lat / 360.)

# get input data
#atmospheric
nc_f = '/home/bastakdu/Documents/HErZ/Annika/data/2019082600_an/3d.nc'  # Path to your filename
#surface
nc_s = '/home/bastakdu/Documents/HErZ/Annika/data/20180401-04/surf201804.nc'  # Path to your filename
#soil
nc_so = '/home/bastakdu/Documents/HErZ/Annika/data/2019082600_an/soil.nc'  # Path to your filename
#constant
nc_c = '/home/bastakdu/Documents/HErZ/Annika/data/constD2/lmboden.nc'  # Path to your filename

#sub-domain for the computation of the forcings
d_y= delta/360. * 2 * np.pi * REarth # calculate domain by using earth's circumference
d_x = d_y * np.cos(lat*np.pi/180.) # attribute for mapping factor (projection of the part-circle in x direction)
print ('dx:',d_x)
print ('dy:',d_y)

nc_cid = Dataset(nc_c, 'r')  

lats_c=nc_cid.variables['lat'][:]
lons_c=nc_cid.variables['lon'][:]


# get external parameters
key_index_ci,key_index_cj = get_indices(lats_c,lons_c,lat, lon)
print ('key_index ', key_index_ci, key_index_cj,lats_c[key_index_ci,key_index_cj],lons_c[key_index_ci,key_index_cj],lat, lon)
print ("reading external parameters")
FR_LAND = nc_cid.variables['FR_LAND'][key_index_ci,key_index_cj]
PLCOV_MX = nc_cid.variables['PLCOV_MX'][key_index_ci,key_index_cj]
LAI_MX = nc_cid.variables['LAI_MX'][key_index_ci,key_index_cj]
ROOTDP = nc_cid.variables['ROOTDP'][key_index_ci,key_index_cj]
RSMIN = nc_cid.variables['RSMIN'][key_index_ci,key_index_cj]
SOILTYP = nc_cid.variables['SOILTYP'][key_index_ci,key_index_cj]
Z0 = nc_cid.variables['Z0'][key_index_ci,key_index_cj]
EMIS_RAD = nc_cid.variables['EMIS_RAD'][key_index_ci,key_index_cj]
LU_CLASS_FRACTION = nc_cid.variables['LU_CLASS_FRACTION'][:,key_index_ci,key_index_cj]
TOPO = nc_cid.variables['HSURF'][key_index_ci,key_index_cj]
nc_cid.close()
print ('finished reading external parameters ')



nc_fid = Dataset(nc_f, 'r')  
nc_cid = Dataset(nc_c, 'r')  
# get atmospheric fields
print ("reading atmospheric fields")
lats_f=nc_fid.variables['lat'][:]
lons_f=nc_fid.variables['lon'][:]
key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright = get_surrounding_points(lats_f,lons_f,lat,lon,delta)
#print ('rightright', jx_leftleft,jx_rightright)
#print ('key_index ', key_index_fi,key_index_fj,lats_f[key_index_fi,key_index_fj],lons_f[key_index_fi,key_index_fj],lat, lon)

#cutting just the region of interest
lats_f=nc_fid.variables['lat'][ix_downdown:ix_upup,jx_leftleft:jx_rightright]
lons_f=nc_fid.variables['lon'][ix_downdown:ix_upup,jx_leftleft:jx_rightright]
timef =  nc_fid.variables['time'][:]
timef_units=nc_fid.variables['time'].units
heights =  nc_fid.variables['height'][:]
heights_w =  nc_fid.variables['height_2'][:]
temp = nc_fid.variables['t'][:,:,ix_downdown:ix_upup,jx_leftleft:jx_rightright]
u=nc_fid.variables['u'][:,:,ix_downdown:ix_upup,jx_leftleft:jx_rightright]
v=nc_fid.variables['v'][:,:,ix_downdown:ix_upup,jx_leftleft:jx_rightright]
w=nc_fid.variables['wz'][:,:,ix_downdown:ix_upup,jx_leftleft:jx_rightright]
pressure=nc_fid.variables['pres'][:,:,ix_downdown:ix_upup,jx_leftleft:jx_rightright]
spec_hum=nc_fid.variables['q'][:,:,ix_downdown:ix_upup,jx_leftleft:jx_rightright]
cloud_mix=nc_fid.variables['clwmr'][:,:,ix_downdown:ix_upup,jx_leftleft:jx_rightright]
rain_mix=nc_fid.variables['rwmr'][:,:,ix_downdown:ix_upup,jx_leftleft:jx_rightright]
cloud_ice_mix=nc_fid.variables['QI'][:,:,ix_downdown:ix_upup,jx_leftleft:jx_rightright]
snow_mix=nc_fid.variables['snmr'][:,:,ix_downdown:ix_upup,jx_leftleft:jx_rightright]
graupel=nc_fid.variables['grle'][:,:,ix_downdown:ix_upup,jx_leftleft:jx_rightright]
HSURF = nc_cid.variables['HSURF'][ix_downdown:ix_upup,jx_leftleft:jx_rightright]
lctype=(nc_cid.rawdata).split(',')[0]

key_index_fi=key_index_fi-ix_downdown
key_index_fj=key_index_fj-jx_leftleft
ix_up=ix_up-ix_downdown
ix_upup=ix_upup-ix_downdown
ix_down=ix_down-ix_downdown
ix_downdown=ix_downdown-ix_downdown
jx_right=jx_right-jx_leftleft
jx_rightright=jx_rightright-jx_leftleft
#print ('rightright', jx_leftleft,jx_rightright)
jx_left=jx_left-jx_leftleft
jx_leftleft=jx_leftleft-jx_leftleft
#print ('rightright', jx_leftleft,jx_rightright)



print ("finishing reading atmospheric fields")
nc_fid.close()
nc_cid.close()

# get surface fields
nc_sid = Dataset(nc_s, 'r')  
print ("reading surface fields")
lats_s=nc_sid.variables['lat'][:]
lons_s=nc_sid.variables['lon'][:]
key_index_si,key_index_sj = get_indices(lats_s,lons_s,lat, lon)
#print ('key_index ', key_index_si,key_index_sj,lats_s[key_index_si,key_index_sj],lons_s[key_index_si,key_index_sj],lat, lon)

#consider cutting just the region of interest
time =  nc_sid.variables['time'][:]
height_1 =  nc_sid.variables['height'][:]
height_2 =  nc_sid.variables['height_2'][:]
#fluxes or values ...to be read
#qvs =  nc_sid.variables['QV_S'][key_index_si,key_index_sj]


print ("finishing reading surface fields")
nc_sid.close()

# get soil fields
nc_soid = Dataset(nc_so, 'r')  
print ("reading soil fields")
lats_so=nc_soid.variables['lat'][:]
lons_so=nc_soid.variables['lon'][:]
key_index_soi,key_index_soj = get_indices(lats_so,lons_so,lat, lon)
#print ('key_index ', key_index_soi,key_index_soj,lats_so[key_index_soi,key_index_soj],lons_so[key_index_soi,key_index_soj],lat, lon)

nc_attrs, nc_dims, nc_vars = ncdump(nc_soid, verb=False)
time =  nc_soid.variables['time'][:]
depth_T =  nc_soid.variables['depth'][:]
depth_W =  nc_soid.variables['depth_2'][:]
soil_temp = nc_soid.variables['T_SO'][:]
soil_moisture = nc_soid.variables['W_SO'][:]
if 'W_SO_ICE' in nc_vars:
   soil_ice = nc_soid.variables['W_SO_ICE'][:]

print ("finishing reading soil fields")
nc_soid.close()

#calculate input variables
print ("calculating input variables")
q_t = (spec_hum + cloud_mix + cloud_ice_mix)/\
 (1.0+cloud_mix + cloud_ice_mix)
print ("qt mean",q_t.mean())
q_c=(cloud_mix)*(1.0-q_t)

r_mix=spec_hum/(1.0-q_t)
print ("r mean",r_mix.mean())

temp_v=temp*(1.0+0.608*r_mix-cloud_mix-cloud_ice_mix)
print ("temp_v mean",temp_v.mean())

print ("compute heigths of model levels")
mlheigths=np.zeros([timef.size,heights.size,ix_upup-ix_downdown,jx_rightright-jx_leftleft])
temp_h=np.zeros([timef.size,hlevels.size,ix_upup-ix_downdown,jx_rightright-jx_leftleft])
temp_vh=np.zeros(temp_h.shape)
r_mixh=np.zeros(temp_h.shape)
cloud_mixh=np.zeros(temp_h.shape)
cloud_ice_mixh=np.zeros(temp_h.shape)
pressure_h=np.zeros(temp_h.shape)
pressure_hp=np.zeros(temp_h.shape)
u_h=np.zeros(temp_h.shape)
v_h=np.zeros(temp_h.shape)
w_h=np.zeros(temp_h.shape)
q_th=np.zeros(temp_h.shape)
q_vh=np.zeros(temp_h.shape)
q_ch=np.zeros(temp_h.shape)
hmaskz=np.ones(temp_h.shape)
hmaskzm=np.ones(temp.shape)

p_LSwestmh=np.zeros([timef.size,hlevels.size])
p_LSeastmh=np.zeros([timef.size,hlevels.size])
p_LSnorthmh=np.zeros([timef.size,hlevels.size])
p_LSsouthmh=np.zeros([timef.size,hlevels.size])
u_g1mh=np.zeros([timef.size,hlevels.size])
u_g2mh=np.zeros([timef.size,hlevels.size])

for it in range(0,timef.size):

    print ('it', it)
    for ix in range(ix_downdown,ix_upup):
      for iy in range(jx_leftleft,jx_rightright):
        mlheigths[it,heights.size-1,ix,iy]=Hmod1+HSURF[ix,iy]
        for iz in np.arange(heights.size-2,-1,-1):
           mlheigths[it,iz,ix,iy]=mlheigths[it,iz+1,ix,iy]-np.log(pressure[it,iz,ix,iy]/pressure[it,iz+1,ix,iy])*R_d*0.5*(temp_v[it,iz,ix,iy]+temp_v[it,iz+1,ix,iy])/g

    p_LSwestm,p_LSnorthm,p_LSeastm,p_LSsouthm,p_LSm = get_avg_fields(pressure,hmaskzm,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
    mlheigths_LSwestm,mlheigths_LSnorthm,mlheigths_LSeastm,mlheigths_LSsouthm,mlheigths_LS = get_avg_fields(mlheigths,hmaskzm,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)



    t_LSwestm,t_LSnorthm,t_LSeastm,t_LSsouthm,t_LSm = get_avg_fields(temp,hmaskzm,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
    r_LSwestm,r_LSnorthm,r_LSeastm,r_LSsouthm,r_LSm = get_avg_fields(r_mix,hmaskzm,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
    rl_LSwestm,rl_LSnorthm,rl_LSeastm,rl_LSsouthm,rl_LSm = get_avg_fields(cloud_mix,hmaskzm,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
    ri_LSwestm,ri_LSnorthm,ri_LSeastm,ri_LSsouthm,ri_LSm = get_avg_fields(cloud_ice_mix,hmaskzm,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)

    tv_LSm=t_LSm*(1.0+0.608*r_LSm-rl_LSm-ri_LSm)
    tv_LSwestm=t_LSwestm*(1.0+0.608*r_LSwestm-rl_LSwestm-ri_LSwestm)
    tv_LSeastm=t_LSeastm*(1.0+0.608*r_LSeastm-rl_LSeastm-ri_LSeastm)
    tv_LSnorthm=t_LSnorthm*(1.0+0.608*r_LSnorthm-rl_LSnorthm-ri_LSnorthm)
    tv_LSsouthm=t_LSsouthm*(1.0+0.608*r_LSsouthm-rl_LSsouthm-ri_LSsouthm)
    rho_LSm=p_LSm/(R_d*tv_LSm)

    HSURFLSwestm,HSURFLSnorthm,HSURFLSeastm,HSURFLSsouthm,HSURFLSm = get_avg_fields_s(HSURF,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)

    mlheigths_LSwestm[it,heights.size-1]=Hmod1+HSURFLSwestm
    for iz in np.arange(heights.size-2,-1,-1):
      mlheigths_LSwestm[it,iz]=mlheigths_LSwestm[it,iz+1]-np.log(p_LSwestm[it,iz]/p_LSwestm[it,iz+1])*R_d*0.5*(tv_LSwestm[it,iz]+tv_LSwestm[it,iz+1])/g
      mlheigths_LSeastm[it,iz]=mlheigths_LSeastm[it,iz+1]-np.log(p_LSeastm[it,iz]/p_LSeastm[it,iz+1])*R_d*0.5*(tv_LSeastm[it,iz]+tv_LSeastm[it,iz+1])/g
      mlheigths_LSnorthm[it,iz]=mlheigths_LSnorthm[it,iz+1]-np.log(p_LSnorthm[it,iz]/p_LSnorthm[it,iz+1])*R_d*0.5*(tv_LSnorthm[it,iz]+tv_LSnorthm[it,iz+1])/g
      mlheigths_LSsouthm[it,iz]=mlheigths_LSsouthm[it,iz+1]-np.log(p_LSsouthm[it,iz]/p_LSsouthm[it,iz+1])*R_d*0.5*(tv_LSsouthm[it,iz]+tv_LSsouthm[it,iz+1])/g

    u_g1m= -((1.0/(rho_LSm))*(p_LSnorthm-p_LSsouthm)+g*(mlheigths_LSnorthm-mlheigths_LSsouthm))/d_y/f3
    u_g2m= ((1.0/(rho_LSm))*(p_LSeastm-p_LSwestm)+g*(mlheigths_LSeastm-mlheigths_LSwestm))/d_x/f3


    u_g1mh[it,:]=np.interp(hlevels, mlheigths_LS[it,heights.size-1:0:-1], u_g1m[it,heights.size-1:0:-1])
    u_g2mh[it,:]=np.interp(hlevels, mlheigths_LS[it,heights.size-1:0:-1], u_g2m[it,heights.size-1:0:-1])

    p_LSwestmh[it,:]=np.interp(hlevels, mlheigths_LSwestm[it,heights.size-1:0:-1], p_LSwestm[it,heights.size-1:0:-1])
    p_LSeastmh[it,:]=np.interp(hlevels, mlheigths_LSeastm[it,heights.size-1:0:-1], p_LSeastm[it,heights.size-1:0:-1])
    p_LSnorthmh[it,:]=np.interp(hlevels, mlheigths_LSnorthm[it,heights.size-1:0:-1], p_LSnorthm[it,heights.size-1:0:-1])
    p_LSsouthmh[it,:]=np.interp(hlevels, mlheigths_LSsouthm[it,heights.size-1:0:-1], p_LSsouthm[it,heights.size-1:0:-1])

    for ix in range(ix_downdown,ix_upup):
      for iy in range(jx_leftleft,jx_rightright):
        temp_vh[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], temp_v[it,heights.size-1:0:-1,ix,iy])
        temp_h[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], temp[it,heights.size-1:0:-1,ix,iy])
        r_mixh[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], r_mix[it,heights.size-1:0:-1,ix,iy])
        cloud_mixh[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], cloud_mix[it,heights.size-1:0:-1,ix,iy])
        cloud_ice_mixh[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], cloud_ice_mix[it,heights.size-1:0:-1,ix,iy])
        pressure_h[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], pressure[it,heights.size-1:0:-1,ix,iy])
        u_h[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], u[it,heights.size-1:0:-1,ix,iy])
        v_h[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], v[it,heights.size-1:0:-1,ix,iy])
        w_h[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights_w.size-1:0:-1,ix,iy], w[it,heights_w.size-1:0:-1,ix,iy])
        #w_h[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], w[it,heights.size-1:0:-1,ix,iy])
        q_th[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], q_t[it,heights.size-1:0:-1,ix,iy])
        q_vh[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], spec_hum[it,heights.size-1:0:-1,ix,iy])
        q_ch[it,:,ix,iy]=np.interp(hlevels, mlheigths[it,heights.size-1:0:-1,ix,iy], q_c[it,heights.size-1:0:-1,ix,iy])
        #pressure_hp=pressure_h
        for iz in np.arange(0,hlevels.size):
          ilps=np.where(hlevels[iz]>mlheigths[it,:,ix,iy])[0]
          if (ilps.size>0):
           ilp=ilps[0]
           #print("ilp:",ilp,"iz:",iz)
           pressure_hp[it,iz,ix,iy]=pressure[it,ilp,ix,iy]*np.exp(9.81*(mlheigths[it,ilp,ix,iy]-hlevels[iz])/R_d/temp_v[it,ilp,ix,iy])
        #pressure_h=pressure_hp

        for iz in range(0,hlevels.size):
           hmask=np.where(hlevels[iz]<mlheigths[it,heights.size-1,:,:])
           hmaskz[it,iz,hmask[0],hmask[1]]=0



print ("compute means")
_,_,_,_,t_LS = get_avg_fields(temp_h,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
_,_,_,_,r_LS = get_avg_fields(r_mixh,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
_,_,_,_,rl_LS = get_avg_fields(cloud_mixh,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
_,_,_,_,ri_LS = get_avg_fields(cloud_ice_mixh,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
print ("tempLS mean",t_LS.mean())
print ("rLS mean",r_LS.mean())
print ("rlLS mean",rl_LS.mean())
print ("riLS mean",ri_LS.mean())



p_LSwest,p_LSnorth,p_LSeast,p_LSsouth,p_LS = get_avg_fields(pressure_h,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
print ("pressure mean",p_LSwest.mean(),p_LSnorth.mean(),p_LSeast.mean(),p_LSsouth.mean(),p_LS.mean())

tv_LS=t_LS*(1.0+0.608*r_LS-rl_LS-ri_LS)
rho_LS=p_LS/(R_d*tv_LS)
print ("tvLS mean",tv_LS.mean())
print ("rhoLS mean",rho_LS.mean())

# pressure gradient
p_LS_eastwestgrad=(p_LSeast-p_LSwest)/d_x
p_LS_northsouthgrad=(p_LSnorth-p_LSsouth)/d_y
#p_LS_eastwestgrad=(p_LSeastmh-p_LSwestmh)/d_x
#p_LS_northsouthgrad=(p_LSnorthmh-p_LSsouthmh)/d_y

# vertical profiles of geostrophic wind profiles 
# u_g,i
#alternative: u_g1= (-1.0/(rho_LS*f3))*p_LS_northsouthgrad
#alternative: u_g2= (1.0/(rho_LS*f3))*p_LS_eastwestgrad
u_g1=u_g1mh
u_g2=u_g2mh

#converge to the u and v profiles
wdamph=np.ones(hlevels.size)
idamph=(np.where(hlevels>damph))[0]
wdamph[idamph]=np.maximum(1.0-(hlevels[idamph]-damph)*damph_slope,0.0)
u_g1=u_g1*wdamph+u_h[:,:,key_index_fi,key_index_fj]*(1-wdamph)
u_g2=u_g2*wdamph+v_h[:,:,key_index_fi,key_index_fj]*(1-wdamph)

w_g=np.sqrt(u_g1**2+u_g2**2)
print ("ug1 mean",u_g1.mean())
print ("ug2 mean",u_g2.mean())
#correction for vertical layers bellow orography
for it in range(0,timef.size):
  for iz in range(0,hlevels.size):
    #if(hmaskz[it,iz,:,:].mean()<0.975):
    if(hmaskz[it,iz,:,:].mean()<0.5):
      u_g1[it,iz]=u_g1[it,iz-1]
      u_g2[it,iz]=u_g2[it,iz-1]

u_1LSwest,u_1LSnorth,u_1LSeast,u_1LSsouth,u_1LS = get_avg_fields(u_h,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
u_2LSwest,u_2LSnorth,u_2LSeast,u_2LSsouth,u_2LS = get_avg_fields(v_h,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
u_3LSwest,u_3LSnorth,u_3LSeast,u_3LSsouth,u_3LS = get_avg_fields(w_h,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
print ("u1LS mean",u_1LS.mean())
print ("u2LS mean",u_2LS.mean())
print ("u3LS mean",u_3LS.mean())
u_1LS_eastwestgrad=(u_1LSeast-u_1LSwest)/d_x
u_1LS_northsouthgrad=(u_1LSnorth-u_1LSsouth)/d_y
u_2LS_eastwestgrad=(u_2LSeast-u_2LSwest)/d_x
u_2LS_northsouthgrad=(u_2LSnorth-u_2LSsouth)/d_y

#theta_l,LS
theta_lh,theta_h = calculate_liquid_potential_temperature(temp_h,pressure_h,cloud_mixh)
theta_lLSwest,theta_lLSnorth,theta_lLSeast,theta_lLSsouth,theta_lLS = get_avg_fields(theta_lh,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
_,_,_,_,theta_LS = get_avg_fields(theta_h,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
# q_t,LS
q_vLSwest,q_vLSnorth,q_vLSeast,q_vLSsouth,qv_LS = get_avg_fields(q_vh,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
_,_,_,_,qc_LS = get_avg_fields(q_ch,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
q_tLSwest,q_tLSnorth,q_tLSeast,q_tLSsouth,q_tLS = get_avg_fields(q_th,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
_,_,_,_,q_cLS = get_avg_fields(cloud_mixh,hmaskz,lats_f,lons_f,key_index_fi,key_index_fj, ix_upup, ix_up, ix_down, ix_downdown, jx_leftleft, jx_left, jx_right, jx_rightright)
# larger-scale gradients (horizontal, vertical) of theta_l,LS, q_t,LS
theta_lLS_eastwestgrad=(theta_lLSeast-theta_lLSwest)/d_x
theta_lLS_northsouthgrad=(theta_lLSnorth-theta_lLSsouth)/d_y
q_tLS_eastwestgrad=(q_tLSeast-q_tLSwest)/d_x
q_tLS_northsouthgrad=(q_tLSnorth-q_tLSsouth)/d_y
q_vLS_eastwestgrad=(q_vLSeast-q_vLSwest)/d_x
q_vLS_northsouthgrad=(q_vLSnorth-q_vLSsouth)/d_y


# larger-scale tendencies (horizontal, vertical) of theta_l,LS, q_t,LS
ddt_theta_LSA, _ = calculate_tendencies(theta_lLS_eastwestgrad, theta_lLS_northsouthgrad, u_1LS, u_2LS, u_3LS)
ddt_q_t_LSA, _ = calculate_tendencies(q_tLS_eastwestgrad, q_tLS_northsouthgrad, u_1LS, u_2LS, u_3LS)
ddt_q_v_LSA, _ = calculate_tendencies(q_vLS_eastwestgrad, q_vLS_northsouthgrad, u_1LS, u_2LS, u_3LS)
ddt_u_1_LSA, _ = calculate_tendencies(u_1LS_eastwestgrad, u_1LS_northsouthgrad, u_1LS, u_2LS, u_3LS)
ddt_u_2_LSA, _ = calculate_tendencies(u_2LS_eastwestgrad, u_2LS_northsouthgrad, u_1LS, u_2LS, u_3LS)
ddt_u_LSP = - (f3 * u_g2)
ddt_v_LSP = (f3 * u_g1)
print ("ddt_theta_LSA",ddt_theta_LSA.mean())
print ("ddt_q_t_LSA",ddt_q_t_LSA.mean())
print ("ddt_q_v_LSA",ddt_q_t_LSA.mean())
print ("ddt_u_LSA",ddt_u_1_LSA.mean())
print ("ddt_v_LSA",ddt_u_2_LSA.mean())
print ("ddt_u_LSP",ddt_u_LSP.mean())
print ("ddt_v_LSP",ddt_v_LSP.mean())


#-------------------------------------------------------------------
# write data into netcdf file
levs=np.where(hmaskz[:].mean(axis=(0,2,3))>0.5)[0] #only levels, which are above the actual station height

print ("write to netCDF file ", file_out)
do=Dataset(file_out,mode='w',format='NETCDF4_CLASSIC')
         
do.description='Initial profiles and forcings generated from COSMO data'
do.createDimension('lev',size=levs.size)
do.createDimension('nt', size=timef.size)
do.createDimension('levTsoil',size=depth_T.size)
do.createDimension('levWsoil',size=depth_W.size)
do.createDimension('nclass_lu',size=LU_CLASS_FRACTION.size)

#--- coordinates

#time
nc_time=do.createVariable('time',np.float64,'nt')
nc_time.units='s'

#level
lev=do.createVariable('lev',np.float64,'lev')
lev.units='m'
lev.standard_name='height'
lev[:]=levs

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
nc_height[:]=np.tile(hlevels[levs],timef.size)
nc_time[:]=timef*3600.0#in seconds
#nc_time[1]=1200.#correction for GOAMAZON time[1]=1119- probably a rounding error in the time variable
nc_uIN[:]=(u_h[:,levs,key_index_fi,key_index_fj])
nc_uIN[:,0]=nc_uIN[:,1]#no gradient of u in stratosphere
nc_vIN[:]=(v_h[:,levs,key_index_fi,key_index_fj])
nc_vIN[:,0]=nc_vIN[:,1]#no gradient of u in stratosphere
nc_wIN[:]=(w_h[:,levs,key_index_fi,key_index_fj])
nc_tkeIN[:]=0.0
nc_thIN[:]=(theta_h[:,levs,key_index_fi,key_index_fj])
nc_thIN[:,0]=nc_thIN[:,1]+dtdz_str*(nc_height[:,0]-nc_height[:,1])#statospheric gradient of pot. temperature, levels above are extrapolated according to this gradinet in SCM
nc_qvIN[:]=(q_vh[:,levs,key_index_fi,key_index_fj])
nc_qvIN[:,0]=nc_qvIN[:,1] #no gradient of moisture in stratosphere
nc_qcIN[:]=(q_ch[:,levs,key_index_fi,key_index_fj])
nc_pIN[:]=pressure_h[:,levs,key_index_fi,key_index_fj]

##soil variables
nc_T_SO=do.createVariable('T_SO',np.float64,('nt','levTsoil'))
nc_T_SO.description='soil temperature'
nc_T_SO.units='K'
nc_W_SO=do.createVariable('W_SO',np.float64,('nt','levWsoil'))
nc_W_SO.description='water content'
nc_W_SO.units='kg m-2'
#values
nc_T_SO[:]=(soil_temp[:,:,key_index_soi,key_index_soj])
nc_W_SO[:]=(soil_moisture[:,:,key_index_soi,key_index_soj])

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
nc_uGEO[:]=(u_g1[:,levs])
nc_vGEO[:]=(u_g2[:,levs])
nc_wLS[:]=(u_3LS[:,levs])
nc_wLS[:,0]=nc_wLS[:,1]#no gradient in stratosphere
nc_dTadv[:]=(ddt_theta_LSA[:,levs])
nc_dTadv[:,0]=nc_dTadv[:,1]#no gradient in stratosphere
nc_dTrad[:]=nc_dTadv[:,levs]*0.0 # setting radiation tendency to 0
nc_dQVadv[:]=(ddt_q_v_LSA[:,levs])
nc_dQVadv[:,0]=nc_dQVadv[:,1]#no gradient in stratosphere
nc_dUadv[:]=(ddt_u_1_LSA[:,levs])
nc_dVadv[:]=(ddt_u_2_LSA[:,levs])

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
nc_psurf[:]=pressure_h[:,levs[-1],key_index_fi,key_index_fj] #direct input from surface file would be better
nc_ts[:]=temp_h[:,levs[-1],key_index_fi,key_index_fj]#no better input for now - should be replaced with input from surface file
nc_tg[:]=nc_ts[:]#no better input for now - input from surface file
nc_qvs[:]=q_vh[:,levs[-1],key_index_fi,key_index_fj]##no better input for now - input from surface file ?
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
nc_TOPO[:]=TOPO
nc_LU_CLASS_FRACTION[:]=LU_CLASS_FRACTION

today = datetime.today()
do.history = "Created " + today.strftime("%d/%m/%y")

do.close()
