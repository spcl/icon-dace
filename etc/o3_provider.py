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

from yac import *
import xarray as xr
import numpy as np
import isodate
from datetime import date
from glob import glob
import f90nml
from os.path import exists

import sys

VERBOSE = 2

iso_data_interval='P1M'
iso_coupling_interval='PT90M' # radiation time step

dataPath = "/pool/data/ECHAM6/input/r0008/T127/ozone/"
fileRoot = "T127_ozone_"

NAMELIST=sys.argv[1]
nml_fname = glob(NAMELIST)[0]
nml = f90nml.read(nml_fname)

irad_o3  = nml['aes_rad_nml']['aes_rad_config'][0]['irad_o3']
print ( "o3_provider: found irad_o3 = ", irad_o3, ' in ',  nml_fname )

try :
    lyr_perp = nml['aes_rad_nml']['aes_rad_config'][0]['lyr_perp']
    print ( "o3_provider: found lyr_perp = ", lyr_perp, " in ",  nml_fname )
except :
    print ( "o3_provider: lyr_perp not set in ",  nml_fname )
    lyr_perp = False

if ( lyr_perp ) :
    yr_perp = nml['aes_rad_nml']['aes_rad_config'][0]['yr_perp']
    print ( "o3_provider: found yr_perp = ", yr_perp, " in ",  nml_fname )

if ( irad_o3 == 6 ) :
    scenario = "picontrol"
if ( irad_o3 == 5 ) :
    scenario = "historical"

if ( irad_o3 != 5 and irad_o3 != 6 ) :
    print ( "o3_provider: irad_o3 =", irad_o3, "is not supported" )
    raise SystemExit(1); exit

try :
    lrad_yac = nml['aes_rad_nml']['aes_rad_config'][0]['lrad_yac']
    print ( "o3_provider: found lrad_yac = ", lrad_yac, " in ",  nml_fname )
except :
    print ( "o3_provider: lrad_yac not set in ",  nml_fname )
    lrad_yac = False

if ( not lrad_yac ) :
    print ( "o3_provider: lrad_yac = .FALSE. cannot be used when running o3_provider " )
    raise SystemExit(1); exit

# scenario = "ssp119"
# scenario = "ssp126"
# scenario = "ssp370"
# scenario = "ssp585"
# scenario = "perpetual" ; scenario_year = "historical_2014"

switch_year = False
scenario_year = 2014

watch_this = 0
watch_other = 0

data_interval = isodate.parse_duration(iso_data_interval)

def filename_dflt ( dataPath , fileRoot, scenario ) :
    filename = dataPath+fileRoot+scenario+".nc"
    return filename

def filename_perp ( dataPath , fileRoot, scenario_year ) :
    filename = dataPath+fileRoot+scenario_year+".nc"
    
    if ( exists(filename) ) :
        return filename
    else :
        print ( filename, " not available!" )
        raise SystemExit(1) ; exit

def filename_year ( dataPath , fileRoot, scenario, year ) :
    if ( scenario[:3] == "ssp" and year > 2014 ) :
        filename = dataPath+fileRoot+scenario+"_"+str(year)+".nc"
    else:
        filename = dataPath+fileRoot+"historical"+"_"+str(year)+".nc"
        
    if ( exists(filename) ) :
        return filename
    else :
        print ( filename, " not available!" )
        raise SystemExit(1); exit

yac = YAC()
def_calendar(Calendar.PROLEPTIC_GREGORIAN)
comp = yac.def_comp(f"o3_provider")

start_date = isodate.parse_datetime(yac.start_datetime)

if ( scenario == "picontrol" ) :
    filename = filename_dflt( dataPath, fileRoot ,scenario )
elif ( scenario == "perpetual" ) :
    filename = filename_perp( dataPath, fileRoot, scenario_year )
elif ( scenario[:3] == "ssp" or scenario == "historical" ) :
    filename = filename_year( dataPath, fileRoot, scenario, start_date.year )
else :
    print ( "scenario ", scenario, " not supported!" )
    raise SystemExit(1); exit

dataset = xr.open_dataset(filename, decode_times=False)
other_dataset = dataset

if ( VERBOSE > 0 ) :
    print ( "o3_provider: Reading ", filename, flush=True )

deg2rad = np.pi/180
lon = deg2rad*dataset["lon"]
lat = deg2rad*dataset["lat"]
dims = ["lat", "lon"]
grid = Reg2dGrid("o3_grid", lon, lat)
points = grid.def_points(Location.CORNER, lon, lat)

o3_field = Field.create("o3", comp, points, dataset.dims["plev"],
                        iso_coupling_interval, TimeUnit.ISO_FORMAT)

yac.enddef()

end_datetime = isodate.parse_datetime(yac.end_datetime)

while isodate.parse_datetime(o3_field.datetime) < end_datetime:

    this_date  = isodate.parse_datetime(o3_field.datetime)
    this_month = this_date.month
    this_year  = this_date.year

    if ( this_month < 12) :
        mid_of_this_month = (date(this_year, this_month+1, 1) -
                             date(this_year, this_month, 1)).days * 43200
    else :
        mid_of_this_month = 31 * 43200

    sec_of_this_month = (this_date.day-1) * 86400 + \
                         this_date.hour   * 3600 +  \
                         this_date.minute * 60 +    \
                         this_date.second  
    
    if ( sec_of_this_month <= mid_of_this_month ) :
        switch = -1
    else :
        switch = 1

    other_date = this_date + switch * data_interval
    other_month = other_date.month
    other_year = other_date.year

    if ( scenario[:3] == "ssp" or scenario == "historical" ) :

        if ( other_year != this_year and not switch_year ) :
            filename = filename_year( dataPath, fileRoot, scenario, other_year )
            other_dataset = xr.open_dataset(filename, decode_times=False)
            switch_year = True
            if ( VERBOSE > 0 ) :
                print ( "o3_provider: Opened dataset for year", other_year, flush=True  )

        if ( other_year == this_year and switch_year ) :
            switch_year = False
            other_dataset = dataset
            if ( VERBOSE > 0 ) :
                print ( "o3_provider: Complete switch to year", other_year, flush=True  )

    if ( other_month < 12 ) :
        mid_of_other_month = (date(this_year, other_month+1, 1) -
                              date(this_year, other_month, 1)).days * 43200
    else :
        mid_of_other_month = 31 * 43200
 
    other_wght = switch * ( sec_of_this_month - mid_of_this_month ) / \
                          ( mid_of_this_month + mid_of_other_month )
    this_wght  = 1 - other_wght

    if ( this_month != watch_this ) :
        if ( watch_this == 0 ) :
            if ( VERBOSE > 1 ) :
                print ( "o3_provider: Reading this month", this_month , flush=True )
            o3_this_month  = dataset["O3"][{"time":this_month-1}].values.reshape((1,-1))
        else :
            if ( VERBOSE > 1 ) :
                print ( "o3_provider: Flipp this month to", this_month, flush=True  )
            o3_this_month = o3_other_month
        watch_this = this_month

    # with an extra copy one could save some more reading from file
    if ( other_month != watch_other ) :
        if ( VERBOSE > 1 ) :
            print ( "o3_provider: Reading other month", other_month, flush=True  )
        o3_other_month = other_dataset["O3"][{"time":other_month-1}].values.reshape((1,-1))
        watch_other = other_month

    if VERBOSE > 1 :
        print ( "o3_provider: ", this_date.isoformat() , ' wght ' , 
                "%8.6f" % this_wght, '*' , "%2i" % this_month , '+' ,
                "%8.6f" % other_wght, '*', "%2i" % other_month, flush=True )

    o3_array = this_wght * o3_this_month + other_wght * o3_other_month

    o3_field.put(o3_array)

