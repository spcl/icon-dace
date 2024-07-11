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

DRYRUN = False

if ( not DRYRUN ) :
    from yac import *

import xarray as xr
import numpy as np
import pandas as pd

import isodate
import datetime as dt
from datetime import datetime

from glob import glob
import f90nml
from os.path import exists

import socket
import sys

from cdo import *
cdo = Cdo(tempdir='/scratch/m/m300083/tmp')

# ------------------------------------------
#
# some early definition changable by user
#
VERBOSE = 2

iso_data_interval='P1M'
PERP_YEAR = 1850

# ------------------------------------------
# get setting from namelist
# ------------------------------------------

NAMELIST=sys.argv[1]
nml_fname = glob(NAMELIST)[0]
nml = f90nml.read(nml_fname)

perpetual_year = PERP_YEAR

try :
    iso_coupling_interval = nml['aes_phy_nml']['aes_phy_config'][0]['dt_rad']
    print ( "o3_provider: found dt_rad =", iso_coupling_interval, "in", nml_fname )
except :
    print ( "o3_provider: dt_rad not set in", nml_fname )
    raise SystemExit(1); exit

try :
    irad_o3  = nml['aes_rad_nml']['aes_rad_config'][0]['irad_o3']
    print ( "o3_provider: found irad_o3 =", irad_o3, "in", nml_fname )
except :
    print ( "o3_provider: irad_o3 not set in",  nml_fname )
    raise SystemExit(1); exit

try :
    lyr_perp = nml['aes_rad_nml']['aes_rad_config'][0]['lyr_perp']
    print ( "o3_provider: found lyr_perp =", lyr_perp, "in", nml_fname )
except :
    print ( "o3_provider: lyr_perp not set in",  nml_fname )
    lyr_perp = False

if ( lyr_perp ) :
    yr_perp = nml['aes_rad_nml']['aes_rad_config'][0]['yr_perp']
    print ( "o3_provider: found yr_perp =", yr_perp, "in", nml_fname )

if ( irad_o3 == 6 ) :
    scenario = "picontrol"

if ( irad_o3 == 5 ) :
    scenario = "historical"
    if ( lyr_perp ) :
        scenario = "perpetual"
        perpetual_year = yr_perp

if ( irad_o3 != 5 and irad_o3 != 6 ) :
    print ( "o3_provider: irad_o3 =", irad_o3, "is not supported" )
    raise SystemExit(1); exit

try :
    coupled_to_o3 = nml['coupling_mode_nml']['coupled_to_o3']
    print ( "o3_provider: found coupled_to_o3 =", coupled_to_o3, "in", nml_fname )
except :
    print ( "o3_provider: coupled_to_o3 not set in", nml_fname )
    coupled_to_o3 = False

if ( not coupled_to_o3 ) :
    print ( "o3_provider: coupled_to_o3 = .FALSE. cannot be used when running o3_provider" )
    raise SystemExit(1); exit

# ------------------------------------------
# definition of fuctions
# ------------------------------------------

def get_hostname () :

    fqdn = socket.getfqdn().split(".")

    if ( "nid" == fqdn[0][:3] ) :
        hostname = "Lumi"
    elif ( "lvt.dkrz.de" == fqdn[1]+"."+fqdn[2]+"."+fqdn[3] ) :
        hostname = "Levante"
    else :
        hostname = "unknown"
        raise ValueError(f'Host cannot be detected')

    return hostname
    
def filename_year ( dataPath, fileRoot, scenario, year ) :
    filename = dataPath+fileRoot+str(year)+".nc"
    return filename

def input_historical ( year ) :

    if ( 1849 < year < 1900 ) : year_range = "185001-189912"
    if ( 1899 < year < 1950 ) : year_range = "190001-194912"
    if ( 1949 < year < 2000 ) : year_range = "195001-199912"
    if ( 1999 < year < 2015 ) : year_range = "200001-201412"

    if ( year < 1849 or 2015 < year):
        raise ValueError(f'input_historical got illegal input {year}')

    if ( get_hostname() == "Levante" ) :
        dataPath = "/work/kd0956/INPUT4MIPS/data/input4MIPs/CMIP6/CMIP/UReading/UReading-CCMI-1-0/atmos/mon/vmro3/gn/v20160711/"
    if ( get_hostname() == "Lumi" ) :
         dataPath = "/appl/local/climatedt/pool/data/ICON/grids/public/mpim/common/CMIP6_ozone"

    fileRoot = "vmro3_input4MIPs_ozone_CMIP_UReading-CCMI-1-0_gn_"
    filename = filename_year( dataPath, fileRoot, "historical", year_range )
 
    return filename

def input_scenario ( scenario, year ) :

    if ( 2014 < year < 2050 ) : year_range = "201501-204912"
    if ( 2049 < year < 2100 ) : year_range = "205001-209912"

    if ( year < 2014 or year > 2100 ):
        raise ValueError(f'input_scenario got illegal input {year}')

    if ( scenario == "ssp119" ) :
        vdate="v20190201"
    else :
        vdate="v20181101"

    if ( get_hostname() == "Levante" ) :
        dataPath = "/work/kd0956/INPUT4MIPS/data/input4MIPs/CMIP6/ScenarioMIP/UReading/UReading-CCMI-"+scenario+"-1-0/atmos/mon/vmro3/gn/"+vdate+"/"
    if ( get_hostname() == "Lumi" ) :
        dataPath = "/appl/local/climatedt/pool/data/ICON/grids/public/mpim/common/CMIP6_ozone"

    fileRoot = "vmro3_input4MIPs_ozone_ScenarioMIP_UReading-CCMI-"+scenario+"-1-0_gn_"
    filename = filename_year( dataPath, fileRoot, scenario, year_range )

    return filename

def input_file ( scenario, year, perpetual_year=PERP_YEAR ) :

    if ( scenario == "picontrol" ) :
        filename = input_historical ( 1850 )

    elif ( scenario == "historical" ) :

        if ( year < 2015 ) :
            filename = input_historical ( year )
            if ( VERBOSE > 0 ) :
                print ( "o3_provider: reading from historical" )
        else :
            filename = input_scenario ( "ssp370", year )
            if ( VERBOSE > 0 ) :
                print ( "o3_provider: reading from ssp370" )

    elif ( scenario[:3] == "ssp") :

        if ( year < 2015 ) :
            filename = input_historical ( year )
            if ( VERBOSE > 0 ) :
                print ( "o3_provider: reading from historical" )
        else :
            filename = input_scenario ( scenario, year )
            if ( VERBOSE > 0 ) :
                print ( "o3_provider: reading from", scenario )

    elif ( scenario == "perpetual" ) :

        if ( 1849 < perpetual_year < 2015 ) :
            filename = input_historical ( perpetual_year )
            if ( VERBOSE > 0 ) :
                print ( "o3_provider: reading from historical" )
        elif ( 2014 < perpetual_year < 2100 ) :
            scenario = "ssp370"
            filename = input_scenario ( scenario, perpetual_year )
            if ( VERBOSE > 0 ) :
                print ( "o3_provider: reading from", scenario )

    if ( exists(filename) ) :
        return filename

    else :
        raise OSError("o3_provider:" +filename + "not available!" )

# ------------------------------------------
# definition of filelist
# ------------------------------------------

def filelist ( scenario, start_year, end_year ) :
    db = {}
    for year in range(start_year,end_year+1):
        db[year] = input_file(scenario, year)

    filenames = []
    filenames.append(db[start_year])
    for year in range(start_year+1,end_year+1):
        yearm1 = year-1
        if ( db[year] != db[yearm1] ) :
            filenames.append(db[year])
    filelist = ' '.join(filenames)
    return filelist

# ------------------------------------------
# yac initialisation
# ------------------------------------------

if ( not DRYRUN ) :
    yac = YAC()
    def_calendar(Calendar.PROLEPTIC_GREGORIAN)
    o3_comp = yac.def_comp(f"o3_provider")

# get coordinates
# ------------------------------------------

filename = input_file("historical", 1850)
dataset = xr.open_dataset(filename, decode_times=False)

deg2rad = np.pi/180

lon = deg2rad*dataset["lon"]
lat = deg2rad*dataset["lat"]

# pressure levels: read, invert and store as string
# ------------------------------------------

plev_string = str(dataset.sizes["plev"])
plev_string += ''.join(["%10.2f" % x for x in dataset["plev"][::-1].values*100 ])

if ( not DRYRUN ) :
    o3_grid   = Reg2dGrid("o3_grid", lon, lat)
    o3_points = o3_grid.def_points(Location.CORNER, lon, lat)

    o3_field  = Field.create("o3", o3_comp, o3_points, dataset.sizes["plev"],
                             iso_coupling_interval, TimeUnit.ISO_FORMAT)

dataset.close()

if ( not DRYRUN ) :
    yac.enddef()

# note that start and end dates are only available after yac.enddef

if ( not DRYRUN ) :
    start_date = isodate.parse_datetime(yac.start_datetime)
    end_date   = isodate.parse_datetime(yac.end_datetime)
else :
    start_date = isodate.parse_datetime('1970-01-01T00:00:00.000')
    end_date   = isodate.parse_datetime('1970-03-01T00:00:00.000')

coupling_interval = isodate.parse_duration(iso_coupling_interval)
data_interval     = isodate.parse_duration(iso_data_interval)

if ( "picontrol" == scenario ) :
    start_year = PERP_YEAR
    end_year   = PERP_YEAR
else:
    start_year = start_date.year-1
    end_year   = end_date.year+1

# ------------------------------------------
# get list of all input files for time loop 
# ------------------------------------------

filenames = filelist ( scenario , start_year, end_year )

# ------------------------------------------
# reading the data
# ------------------------------------------

model_date = start_date
vmro3_date = start_date

if ( "perpetual" == scenario or "picontrol" ==  scenario ) :
    vmro3_date = model_date.replace(year=PERP_YEAR)

if ( VERBOSE > 0 ) :
    print ( "o3_provider: reading from", filenames )

# file read with cdo and calendar conversion
input = '-select,name=vmro3,year={1}/{2} {0}'.format(filenames,start_year,end_year )
ds = (cdo.settunits('days', input = input, options = '-r',
                    returnXDataset = True)).convert_calendar("standard", use_cftime=True)

# ------------------------------------------
# time loop
# ------------------------------------------

while model_date < end_date:

    vmro3_date = model_date

    if ( scenario == "perpetual" or scenario == "picontrol" ) :

        vmro3_date = model_date.replace(year=perpetual_year)

        if ( vmro3_date < dt.datetime( perpetual_year,  1, 16, 12, 0, 0) ) :

            auxil_date   = dt.datetime( perpetual_year, 12, 16, 12, 0, 0 )
            ds_prev_elem = ds.sel(time=auxil_date, method='nearest')
            o3_prev_date = datetime.strptime(str(ds_prev_elem['time'].values),
                                     '%Y-%m-%d %H:%M:%S' ).isoformat()
            o3_prev_date = isodate.parse_datetime(o3_prev_date).replace(year=perpetual_year-1)
            o3_prev_date = datetime.strptime(str(o3_prev_date),
                                                 '%Y-%m-%d %H:%M:%S').isoformat()
        else :

            ds_prev_elem = ds.sel(time=vmro3_date, method='ffill')
            o3_prev_date = datetime.strptime(str(ds_prev_elem['time'].values),
                                     '%Y-%m-%d %H:%M:%S' ).isoformat()

        if ( vmro3_date > dt.datetime( perpetual_year, 12, 16, 12, 0, 0) ) :

            auxil_date   = dt.datetime( perpetual_year,  1, 16, 12, 0, 0 )
            ds_next_elem = ds.sel(time=auxil_date, method='nearest')
            o3_next_date = datetime.strptime(str(ds_next_elem['time'].values),
                                     '%Y-%m-%d %H:%M:%S' ).isoformat()
            o3_next_date = isodate.parse_datetime(o3_next_date).replace(year=perpetual_year+1)
            o3_next_date = datetime.strptime(str(o3_next_date),
                                     '%Y-%m-%d %H:%M:%S').isoformat()

        else :

            print ( vmro3_date )
            ds_next_elem = ds.sel(time=vmro3_date, method='bfill')
            o3_next_date = datetime.strptime(str(ds_next_elem['time'].values),
                                     '%Y-%m-%d %H:%M:%S' ).isoformat()
    else :

        ds_prev_elem = ds.sel(time=vmro3_date, method='ffill')
        o3_prev_date = datetime.strptime(str(ds_prev_elem['time'].values),
                                     '%Y-%m-%d %H:%M:%S' ).isoformat()

        ds_next_elem = ds.sel(time=vmro3_date, method='bfill')
        o3_next_date = datetime.strptime(str(ds_next_elem['time'].values),
                                     '%Y-%m-%d %H:%M:%S' ).isoformat()

    if ( VERBOSE > 2 ) :
        print ( "o3_provider:", o3_prev_date , ":" , model_date , "/" ,
                                vmro3_date , ":" , o3_next_date)

    if ( o3_next_date == o3_prev_date ) :

        # in case no data are available outside the required time interval

        prev_wght = 0.5

    else :

        # default time interpolation weights

        delta = pd.to_datetime(o3_next_date,utc=True) - pd.to_datetime(o3_prev_date,utc=True)

        delta_sec = delta.total_seconds()

        if ( delta_sec <= 0.0 ) :
            raise ValueError("delta must be larger than 0:",
                             delta_sec, o3_next_date, o3_prev_date)

        prev_wght = 1 - (pd.to_datetime(vmro3_date) - \
                         pd.to_datetime(o3_prev_date)).total_seconds()/delta_sec

        if ( prev_wght < 0.0 or prev_wght > 1.0 ) :
            raise ValueError("weight must be in the range [0:1]",
                             prev_wght, vmro3_date, o3_prev_date)

    next_wght  = 1 - prev_wght

    if ( VERBOSE > 2 ) :
        print ( "o3_provider: delta_sec", delta_sec, ":",
               o3_prev_date , " to " , o3_next_date)

    if ( VERBOSE > 1 ) :
        print ("o3_provider:", vmro3_date , "wght" , 
               "%8.6f" % prev_wght, '*' ,
               "%2i" % dt.datetime.fromisoformat(o3_prev_date).month , '+' ,
               "%8.6f" % next_wght, '*' ,
               "%2i" % dt.datetime.fromisoformat(o3_next_date).month)

    # apply weights, revert data on vertical (pressure level) axis and send out

    o3_prev_elem = ds_prev_elem["vmro3"].values
    o3_next_elem = ds_next_elem["vmro3"].values

    o3_array = prev_wght * o3_prev_elem[::-1,:,:] + next_wght * o3_next_elem[::-1,:,:] 

    if ( not DRYRUN ) :
        o3_field.put(o3_array)

    model_date = model_date + coupling_interval

if VERBOSE > 0 :
    print ( "o3_provider: Done \n", flush=True  )

