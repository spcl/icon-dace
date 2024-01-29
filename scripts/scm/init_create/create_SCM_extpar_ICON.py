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

#
# create extpar file for SCM/LES runs with TERRA
#
# runs as: python3 create_SCM_extpar_ICON.py lat lon   (lon between -180 and 180!!)
#
# alternative: ncks -d cell,74925,74956 icon_extpar_0024_R02B06_G_20180209_tiles.nc test.nc
#
# work flow SCM from ICON input:
#  - ICON ini:   read_icon_ana_oper_mem1_40km.s
#  - ICON run:   run_ICON_4_SCMini
#  - SCM ini:    get_SCM_data_ICON.py
#  - SCM extpar: create_SCM_extpar_ICON.py
#  - SCM run:    run_SCM_ICONini
#  - plot SCM:   plot-scm-*.py
#
# Martin Koehler, 2020-10 
#-------------------------------------------------------------------

import numpy as np
from netCDF4 import Dataset
from subprocess import call
import sys

# arguments

print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:      ', str(sys.argv) )
lat_scm = float(sys.argv[1])
lon_scm = float(sys.argv[2])

# select location
                     # Lindenberg  Ocean
#lat_scm = 52.21      # 52.209722   -60.0
#lon_scm = 14.12      # 14.118889   120.0    (lon between -180 and 180!!)

location='lat'+ "{:.1f}".format(lat_scm) +'_lon'+ "{:.1f}".format(lon_scm)

# number of points

npts = 32            # 32 minimum points for SCM
nt   = 12            # 12 month for climatologies

# Set the input/output directories

icon_dir      = '../../../'
data_dir      = 'scripts/scm/init_create/'
grid_dir      = '/hpc/uwork/mkoehler/scm/data/grid_data/'
ext_dir       = '/hpc/rhome/routfor/routfox/icon/grids/public/edzw/'

file_ext_icon = ext_dir  + 'icon_extpar_0024_R02B06_G_20180209_tiles.nc'
file_grid_scm = icon_dir + grid_dir + 'Torus_Triangles_4x4_2500m.nc'
file_ext_scm  = icon_dir + data_dir + 'extpar_Torus_Triangles_4x4_2500m_'+location+'.nc'

nc_eid = Dataset(file_ext_icon, 'r')

# find location in grid

clat = nc_eid.variables['clat'][:]
clon = nc_eid.variables['clon'][:]

print('  clat: ', type(clat), np.shape(clat), np.amin(clat), 'to', np.amax(clat))
print('  clon: ', type(clon), np.shape(clon), np.amin(clon), 'to', np.amax(clon))

dist = (clon - lon_scm/180.*np.pi)**2 + (clat - lat_scm/180.*np.pi)**2
index_near = np.where(dist == dist.min())
print('  selected index: ', index_near)
index_near = index_near[0][0]          # convert tuple to number
print('  selected lat/lon (rad) ', clat[index_near], clon[index_near], 'and deg',
  clat[index_near]/np.pi*180.0, clon[index_near]/np.pi*180.0  )
print('  target location: (lat/lon in deg) ', lat_scm, lon_scm)

if lon_scm/180.*np.pi < np.amin(clon) or lon_scm/180.*np.pi > np.amax(clon):
  sys.exit('ERROR: longitude range inconsistency')

# get external parameters

print ('... reading external parameter from extpar file: \n ', file_ext_icon)

SOILTYP      = nc_eid.variables['SOILTYP']      [index_near]
FR_LAND      = nc_eid.variables['FR_LAND']      [index_near]
ICE          = nc_eid.variables['ICE']          [index_near]
PLCOV_MX     = nc_eid.variables['PLCOV_MX']     [index_near]
LAI_MX       = nc_eid.variables['LAI_MX']       [index_near]
RSMIN        = nc_eid.variables['RSMIN']        [index_near]
URBAN        = nc_eid.variables['URBAN']        [index_near]
FOR_D        = nc_eid.variables['FOR_D']        [index_near]
FOR_E        = nc_eid.variables['FOR_E']        [index_near]
EMIS_RAD     = nc_eid.variables['EMIS_RAD']     [index_near]
ROOTDP       = nc_eid.variables['ROOTDP']       [index_near]
Z0           = nc_eid.variables['Z0']           [index_near]
lon          = nc_eid.variables['lon']          [index_near]
lat          = nc_eid.variables['lat']          [index_near]
clon         = nc_eid.variables['clon']         [index_near]
clat         = nc_eid.variables['clat']         [index_near]
NDVI_MAX     = nc_eid.variables['NDVI_MAX']     [index_near]
topography_c = nc_eid.variables['topography_c'] [index_near]
SSO_STDH     = nc_eid.variables['SSO_STDH']     [index_near]
SSO_THETA    = nc_eid.variables['SSO_THETA']    [index_near]
SSO_GAMMA    = nc_eid.variables['SSO_GAMMA']    [index_near]
SSO_SIGMA    = nc_eid.variables['SSO_SIGMA']    [index_near]
T_CL         = nc_eid.variables['T_CL']         [index_near]
FR_LAKE      = nc_eid.variables['FR_LAKE']      [index_near]
DEPTH_LK     = nc_eid.variables['DEPTH_LK']     [index_near]
TOPO_CLIM    = nc_eid.variables['TOPO_CLIM']    [index_near]

NDVI         = nc_eid.variables['NDVI']         [:,index_near]
NDVI_MRAT    = nc_eid.variables['NDVI_MRAT']    [:,index_near]
AER_BC       = nc_eid.variables['AER_BC']       [:,index_near]
AER_DUST     = nc_eid.variables['AER_DUST']     [:,index_near]
AER_ORG      = nc_eid.variables['AER_ORG']      [:,index_near]
AER_SO4      = nc_eid.variables['AER_SO4']      [:,index_near]
AER_SS       = nc_eid.variables['AER_SS']       [:,index_near]
ALB          = nc_eid.variables['ALB']          [:,index_near]
ALNID        = nc_eid.variables['ALNID']        [:,index_near]
ALUVD        = nc_eid.variables['ALUVD']        [:,index_near]
T_SEA        = nc_eid.variables['T_SEA']        [:,index_near]
W_SNOW       = nc_eid.variables['W_SNOW']       [:,index_near]
T_2M_CLIM    = nc_eid.variables['T_2M_CLIM']    [:,index_near]
     
LU_CLASS_FRACTION = nc_eid.variables['LU_CLASS_FRACTION'] [:,index_near]


#-------------------------------------------------------------------
# write data into netcdf file

print ("... write to SCM extpar netCDF file \n ", file_ext_scm)
# nc_tid = Dataset(file_temp,mode='w',format='NETCDF4_CLASSIC')
nc_tid = Dataset(file_ext_scm,mode='w',format='NETCDF4_CLASSIC')

nlu  = size=LU_CLASS_FRACTION.size

nc_tid.createDimension('cell'     , size=npts)
nc_tid.createDimension('time'     , size=nt)
nc_tid.createDimension('nclass_lu', size=nlu)
nc_tid.createDimension('level'   , nlu)

var    = nc_tid.createVariable('SOILTYP','i4','cell')
var[:] = np.full(npts,SOILTYP)

var    = nc_tid.createVariable('FR_LAND',np.float64,'cell')
var[:] = np.full(npts,FR_LAND)

var    = nc_tid.createVariable('ICE',np.float64,'cell')
var[:] = np.full(npts,ICE)

var    = nc_tid.createVariable('PLCOV_MX',np.float64,'cell')
var[:] = np.full(npts,PLCOV_MX)

var    = nc_tid.createVariable('LAI_MX',np.float64,'cell')
var[:] = np.full(npts,LAI_MX)

var    = nc_tid.createVariable('RSMIN',np.float64,'cell')
var[:] = np.full(npts,RSMIN)

var    = nc_tid.createVariable('URBAN',np.float64,'cell')
var[:] = np.full(npts,URBAN)

var    = nc_tid.createVariable('FOR_D',np.float64,'cell')
var[:] = np.full(npts,FOR_D)

var    = nc_tid.createVariable('FOR_E',np.float64,'cell')
var[:] = np.full(npts,FOR_E)

var    = nc_tid.createVariable('EMIS_RAD',np.float64,'cell')
var[:] = np.full(npts,EMIS_RAD)

var    = nc_tid.createVariable('ROOTDP',np.float64,'cell')
var[:] = np.full(npts,ROOTDP)

var    = nc_tid.createVariable('Z0',np.float64,'cell')
var[:] = np.full(npts,Z0)

var    = nc_tid.createVariable('lat',np.float64,'cell')
var[:] = np.full(npts,lat)

var    = nc_tid.createVariable('lon',np.float64,'cell')
var[:] = np.full(npts,lon)

var    = nc_tid.createVariable('clat',np.float64,'cell')
var[:] = np.full(npts,clat)

var    = nc_tid.createVariable('clon',np.float64,'cell')
var[:] = np.full(npts,clon)

var    = nc_tid.createVariable('NDVI_MAX',np.float64,'cell')
var[:] = np.full(npts,NDVI_MAX)

var    = nc_tid.createVariable('topography_c',np.float64,'cell')
var[:] = np.full(npts,topography_c)

var    = nc_tid.createVariable('SSO_STDH',np.float64,'cell')
var[:] = np.full(npts,SSO_STDH)

var    = nc_tid.createVariable('SSO_THETA',np.float64,'cell')
var[:] = np.full(npts,SSO_THETA)

var    = nc_tid.createVariable('SSO_GAMMA',np.float64,'cell')
var[:] = np.full(npts,SSO_GAMMA)

var    = nc_tid.createVariable('SSO_SIGMA',np.float64,'cell')
var[:] = np.full(npts,SSO_SIGMA)

var    = nc_tid.createVariable('T_CL',np.float64,'cell')
var[:] = np.full(npts,T_CL)

var    = nc_tid.createVariable('FR_LAKE',np.float64,'cell')
var[:] = np.full(npts,FR_LAKE)

var    = nc_tid.createVariable('DEPTH_LK',np.float64,'cell')
var[:] = np.full(npts,DEPTH_LK)

var    = nc_tid.createVariable('TOPO_CLIM',np.float64,'cell')
var[:] = np.full(npts,TOPO_CLIM)


var    = nc_tid.createVariable('NDVI',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = NDVI

var    = nc_tid.createVariable('NDVI_MRAT',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = NDVI_MRAT

var    = nc_tid.createVariable('AER_BC',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = AER_BC

var    = nc_tid.createVariable('AER_DUST',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = AER_DUST

var    = nc_tid.createVariable('AER_ORG',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = AER_ORG

var    = nc_tid.createVariable('AER_SO4',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = AER_SO4

var    = nc_tid.createVariable('AER_SS',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = AER_SS

var    = nc_tid.createVariable('ALB',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = ALB

var    = nc_tid.createVariable('ALNID',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = ALNID

var    = nc_tid.createVariable('ALUVD',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = ALUVD

var    = nc_tid.createVariable('T_SEA',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = T_SEA

var    = nc_tid.createVariable('W_SNOW',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = W_SNOW

var    = nc_tid.createVariable('T_2M_CLIM',np.float64,('time','cell'))
for n in range(npts):  var[:,n] = T_2M_CLIM


var    = nc_tid.createVariable('LU_CLASS_FRACTION',np.float64,('nclass_lu','cell'))
for n in range(npts):  var[:,n] = LU_CLASS_FRACTION

nc_tid.rawdata = "GLOBCOVER2009, FAO DSMW, GLOBE, Lake Database"

nc_tid.close


# merge grid and extpar temp files into final SCM extpar file

#call('cdo -O merge %s %s %s'%(file_grid_scm, file_temp, file_ext_scm), shell=True)
