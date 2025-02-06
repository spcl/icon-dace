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

'''
Purpose
    Read ICON SCM output and convert to GOAmazon GASS-DCP project 
Details
    netcdf4 in and out
    run as: python3 GOAmazon_reformat_output.py
Author
    Martin Koehler, DWD
Revision history
    20200707 -- Initial version
'''

from netCDF4 import Dataset
import numpy as np
from numpy import dtype
import sys
from GOAmazon_reformat_variables import varmeta


# files

#date='20140215'
date = sys.argv[1]
yyyy = date[:4]

#location = 'GoAmazon'
location = 'SGP'
#yyyy     = '2004'
#expname = 'SCM_GASS_DCP_GoAmazon_2015_test'
expname = 'SCM_GASS_DCP_SGP'

#yyyy = np.genfromtxt('/hpc/uwork/mkoehler/run-icon/scm/' + expname + '/YYYY', dtype=None)
#yyyy = open('/hpc/uwork/mkoehler/run-icon/scm/' + expname + '/YYYY', 'r').read()
print(yyyy)

expdir  = '/hpc/uwork/mkoehler/run-icon/scm/' + expname + '/' + yyyy

#date='20140215'
date = sys.argv[1]

nc_pl_file  =  expdir + '/scm-out/scm_GASS_out_PL_'+date+'T000000Z.nc'
nc_pl_fid   =  Dataset(nc_pl_file, 'r')
nc_ml_file  =  expdir + '/scm-out/scm_GASS_out_ML_'+date+'T000000Z.nc'
nc_ml_fid   =  Dataset(nc_ml_file, 'r')
nc_flx_file =  expdir + '/scm-out/scm_GASS_flx_ML_'+date+'T000000Z.nc'
nc_flx_fid  =  Dataset(nc_flx_file, 'r')

print('open files:')
print(nc_pl_file)
print(nc_ml_file)


# definitions

clat = -3.11
clon = -60.02

nlat = 1
nlon = 1
nlev = 40
ntime= nc_pl_fid.variables['time'][:].shape[0]
print('ntime: ',ntime)


# write file

# open a netCDF file to write
ncout = Dataset(expdir + '/scm-GASS-format/' + location + '_icon_u_nudging.'+date+'.nc', 'w', format='NETCDF4')

# define axis size
ncout.createDimension('time', None)  # unlimited
ncout.createDimension('lat', nlat)
ncout.createDimension('lon', nlon)
ncout.createDimension('lev', nlev)

# metadata
ncout.description = 'simulated and processed by Martin Koehler (martin.koehler@dwd.de)'

# create time axis
time = ncout.createVariable('time', dtype('double').char, ('time',))
time.long_name = 'time'
time.units     = nc_pl_fid.variables['time'].units
time.units     = time.units.replace("minutes", "days")  # convert minutes since to days since

# create vertical axis
lev = ncout.createVariable('lev', dtype('double').char, ('lev'))
lev.long_name = 'pressure'
lev.units     = 'hPa'
lev.positive  = 'down'

# create latitude axis
lat = ncout.createVariable('lat', dtype('double').char, ('lat'))
lat.long_name = 'latitude'
lat.units     = 'degrees_north'

# create longitude axis
lon = ncout.createVariable('lon', dtype('double').char, ('lon'))
lon.long_name = 'longitude'
lon.units     = 'degrees_east'

# fill axis values
time[:] = nc_pl_fid.variables['time'] [:] / 1440.0    # convert minutes to days
lev [:] = nc_pl_fid.variables['plev'] [:]
lon [:] = clon
lat [:] = clat

# multiple variables
icon_pl_names  = ['u', 'v', 'temp', 'tot_qv_dia', 'rh', 'geopot', 'clc', 'tot_qc_dia', 'tot_qi_dia',
  ('ddt_temp_radsw','ddt_temp_radlw','ddt_temp_turb','ddt_temp_drag','ddt_temp_pconv','ddt_temp_gscp'),
  ('ddt_qv_turb','ddt_qv_conv','ddt_qv_gscp'), 'ddt_temp_radsw', 'ddt_temp_radlw',   
  ('ddt_temp_radsw','ddt_temp_radlw','ddt_temp_turb','ddt_temp_drag','ddt_temp_pconv','ddt_temp_gscp','ddt_temp_dyn'),
  'ddt_temp_dyn', 'ddt_temp_gscp',   
  ('ddt_temp_radsw','ddt_temp_radlw','ddt_temp_turb','ddt_temp_drag','ddt_temp_pconv'),
  ('ddt_qv_turb','ddt_qv_conv')]

icon_ml_names  = ['tqv_dia', 't_g',
  'pres_sfc', 'u_10m', 'v_10m', 'tqc_dia', 'tqi_dia', ('tqc_dia','tqi_dia'), ('tqr','tqs'), 'clct' ]

icon_flx_names = ['tot_prec', ('rain_con','snow_gsp'), 'accsod_t', 'accsou_t',
  'accthb_t', ('accsob_s','accsod_s'), 'accsod_s', 'accthu_s', 'accthd_s', 'accshfl_s', 'acclhfl_s']

#icon_pl_names  = []
#icon_ml_names  = []
#icon_flx_names = []
#icon_pl_names  = ['u',('ddt_qv_turb', 'ddt_qv_conv')]
#icon_ml_names  = ['u_10m']
#icon_flx_names = ['tot_prec','accsod_t']

print('..............................')
print('vertical variables: ', icon_pl_names)
print('surface variables:  ', icon_ml_names)
print('flux variables:     ', icon_flx_names)


# PL variables
for icon_name in icon_pl_names:
  print('..............................')
  var       = varmeta(icon_name)                                      # just get metadata
  print('metadata output:                ', var.out_name, var.long_name, var.units)

  # array for many variables
  varmany = np.zeros((ntime, 40, 10))
  print('varmany output shape:           ', varmany.shape)
  ivar = 0

  names = icon_name if isinstance(icon_name, tuple) else [icon_name]  # convert to list
  for n_name in names[:]:
    print('reading data for              ', n_name)
    vardata     = nc_pl_fid.variables[n_name][:]
    print('  variable input shape:         ', vardata.shape)
    vardata_mn  = vardata[:,:,:].mean(axis=2)                         # mean over 32 points
    print('  variable output shape & mean: ', vardata_mn.shape, np.mean(vardata_mn))
    varmany[:,:,ivar] =  vardata_mn
    ivar = ivar + 1

  args      = [icon_name,varmany[:,:,0],varmany[:,:,1],varmany[:,:,2],varmany[:,:,3],varmany[:,:,4],
               varmany[:,:,5],varmany[:,:,6]]
  var       = varmeta(*args)                                          # rescaling of variable data

  # create variable array
  varout    = ncout.createVariable(var.out_name, dtype('float').char, ('time', 'lev', 'lat', 'lon'))
  varout.long_name     = var.long_name
  varout.units         = var.units
  varout.missing_value = -9999

  # write variable data
  print('  final output shape & mean:  ', var.data.shape, np.mean(var.data))
  varout[:] = var.data


# ML variables
for icon_name in icon_ml_names:
  print('..............................')
  var       = varmeta(icon_name)                                      # just get metadata
  print('metadata output:              ', var.out_name, var.long_name, var.units)

  # array for many variables
  varmany = np.zeros((ntime, 10))
  print('varmany output shape:         ', varmany.shape)
  ivar = 0

  names = icon_name if isinstance(icon_name, tuple) else [icon_name]  # convert to list
  for n_name in names[:]:
    print('reading data for              ', n_name)
    vardata     = nc_ml_fid.variables[n_name][:]
    print('  variable input shape:         ', vardata.shape)
    if n_name == 'u_10m' or n_name == 'v_10m':
      vardata_mn  = vardata[:,0,:].mean(axis=1)                       # mean over 32 points
    else:
      vardata_mn  = vardata[:,:].mean(axis=1)                         # mean over 32 points
    print('  variable output shape & mean: ', vardata_mn.shape, np.mean(vardata_mn))
    varmany[:,ivar] =  vardata_mn
    ivar = ivar + 1

  var       = varmeta(icon_name,varmany[:,0],varmany[:,1])            #  rescaling of variable data

  # create variable array
  varout    = ncout.createVariable(var.out_name, dtype('float').char, ('time', 'lat', 'lon'))
  varout.long_name     = var.long_name
  varout.units         = var.units
  varout.missing_value = -9999

  # write variable data
  print('  final output shape & mean:  ', var.data.shape, np.mean(var.data))
  varout[:] = var.data


# FLX variables
for icon_name in icon_flx_names:
  print('..............................')
  var       = varmeta(icon_name)                                      # just get metadata
  print('metadata output:              ', var.out_name, var.long_name, var.units)

  # array for many variables
  varmany = np.zeros((241, 10))
  print('varmany output shape:         ', varmany.shape)
  ivar = 0

  names = icon_name if isinstance(icon_name, tuple) else [icon_name]  # convert to list
  for n_name in names[:]:
    print('reading data for              ', n_name)
    vardata     = nc_flx_fid.variables[n_name][:]
    print('  variable input shape:         ', vardata.shape)
    vardata_mn  = vardata[:,:].mean(axis=1)                           # mean over 32 points
    print('  variable output shape & mean: ', vardata_mn.shape, np.mean(vardata_mn))
    varmany[:,ivar] =  vardata_mn
    ivar = ivar + 1

  var       = varmeta(icon_name,varmany[:,0],varmany[:,1])            #  rescaling of variable data

  # create variable array
  varout    = ncout.createVariable(var.out_name, dtype('float').char, ('time', 'lat', 'lon'))
  varout.long_name     = var.long_name
  varout.units         = var.units
  varout.missing_value = -9999

  # write variable data
  print('  final output shape & mean:  ', var.data.shape, np.mean(var.data))
  varout[:] = var.data[0::2]


# close files
ncout.close()


