# 
# -*- coding: utf-8 -*-
#
# util_inout.py - Utilities for reading and writing netcdf files
# 
# Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#

"""
###############################################################################
## Utility for reading and writing netcdf files - analog zu mo_inout.f90
#
###############################################################################
## Author: Stefan Hagemann (Hereon, Institute of Coastal Systems, GER)
##
"""

from netCDF4 import Dataset
import numpy as np

#
# ***************************************************
class Mapdata:

    def __init__(self):
        self.nlon = 0
        self.nlat = 0
        self.xlon = []
        self.xlat = []
        self.vars = []     # List of variable names
        self.value = []
        self.units = []    # List of variable units
        self.itime = 0     # Time dependence: 0/1 = no/yes

#        
# ***************************************************
    def openw(self, dnam):      # Opens new file for writing, Returns stream = Dataset class

        ofile = Dataset(dnam, 'w', format='NETCDF4')   # or 'NETCDF4_CLASSIC'?
        ofile.description = 'Map data'

        if self.nlon > 0:
            ofile.createDimension('lon', self.nlon) 
            x = ofile.createVariable('lon', np.float32, ('lon')) 
        if self.nlat > 0:
            ofile.createDimension('lat', self.nlat)
            y = ofile.createVariable('lat', np.float32, ('lat')) 
        else:
            y = ofile.createVariable('lat', np.float32, ('lon')) 

        x.standard_name = "longitude"
        x.axis = "X"
        x.long_name = "longitude"
        x.units = "degrees_east" 
        y.standard_name = "latitude"
        y.axis = "Y"
        y.long_name = "latitude"
        y.units = "degrees_north"
        x[:] = self.xlon
        y[:] = self.xlat

        if self.itime == 1:
            self.otime = ofile.createDimension('time', None) 
            if self.nlat > 0:
                vardimtuple = ('time','lat','lon')
            else:
                vardimtuple = ('time', 'lon')
        else:
            if self.nlat > 0:
                vardimtuple = ('lat','lon')
            else:
                vardimtuple = ('lon')

        nvar = len(self.vars)
        if nvar == 1:
            field = ofile.createVariable(str(self.vars[0]), np.float32, vardimtuple) 
            field.standard_name = str(self.vars[0])
            field.units = self.units[0]
        else:
            field = []
            for i in range(nvar):
                fv = ofile.createVariable(self.vars[i], np.float32, vardimtuple) 
                fv.standard_name = self.vars[i]
                fv.units = self.units[i]
                field.append(fv)

        self.varstream = field

        return ofile
#        
# ***************************************************
    def write(self, fdat, ivar=0, yyyymmdd=0):      # Writes array fdat into self.varstream[ivar] of ofile
   
        print('ivar = ', ivar, '  ndim of varstream= ', self.varstream.ndim) 
        if self.varstream.ndim == 2:
            if self.itime == 1:
                self.otime = yyyymmdd
                self.varstream[:,:] = fdat
            else:
                self.varstream[:,:] = fdat
        else:
            if self.itime == 1:
                self.otime = yyyymmdd
                self.varstream[:] = fdat
            else:
                self.varstream[:] = fdat
#        
# ***************************************************
    def close(self, ofile):      # Close ofile
        ofile.close() 
#
# ***************************************************
#        
# ***************************************************
    def openr(self, dnam):      # Opens file dnam for reading, Returns stream = Dataset class

        df = Dataset(dnam, 'r')

        # check coordinate names
        if 'lon' in df.variables.keys():
            clon = 'lon'
            clat = 'lat'
        elif 'clon' in df.variables.keys():
            clon = 'clon'
            clat = 'clat'
        else:
            print('Neither lon nor clon is in variable list of ', dnam)
            print(' --> ERROR --> STP')
            exit()

        ndim = len(df.variables[clat].dimensions)

        if ndim == 1:
            self.nlon = df.variables[clon].size
            self.nlat = df.variables[clat].size
            self.xlon = df.variables[clon][:]
            self.xlat = df.variables[clat][:]
        elif ndim == 2:
            print('Dim of lat: ', len(df.variables[clat].dimensions))
            print('Dim of lon: ', len(df.variables[clon].dimensions))
            self.nlon = df.dimensions[df.variables[clat].dimensions[1]].size
            self.nlat = df.dimensions[df.variables[clat].dimensions[0]].size
            self.xlon = df.variables[clon][:]
            self.xlat = df.variables[clat][:]

        print ('NLON = ', self.nlon, '  NLAT = ', self.nlat)

        return df
#
# ***************************************************
    def varinfo(self, df, ivar=1):      # Provide info on variables in data stream df
                                      # Ivar = 0, all variables, 1 = non coordinate variables 

        print(df.variables.keys())
        nvar = 0
        varlist = []

        coord_list = ['lon', 'time', 'lat', 'level', 'time_bnds', 'clon', 'clon_bnds', 
                      'clat', 'clat_bnds', 'depth_2', 'depth_2_bnds']
     
        for var in df.variables.keys():
            if ivar == 1:
                if var in coord_list:
                    continue 
                if var == 'lon' or var == 'time' or var == 'lat' or var == 'level':
                    continue 
            varlist.append(var)
        nvar = len(varlist)
        print(nvar, ' Variables found') 

        return varlist



                 
        


