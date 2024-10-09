#!/bin/bash

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

#############################################
#  This file was written by Harel Muskatel E-mail: muskatelh@ims.gov.il
#  
#  from the website:
#  https://confluence.ecmwf.int/display/ECRAD
#  download 3D monthly aerosol climatology derived from CAMS reanalysis system as described by Bozzo et al. (2020):FTP site containing netCDF file (274 MB)
#  rename the file: aerosol_cams_3d_climatology_2003-2013_orig.nc
#
#  CDO is needed for this script to work (version 2.1.0 and above)
#############################################

BaseName=aerosol_cams_3d_climatology_2003-2013
GridName=C3_DOM01 # an example

# name of the original CAMS climatology file
origFile=${BaseName}_orig.nc
# temporal text file
txtFile=${BaseName}.txt

ncdump $origFile > $txtFile

search="2003-2013_orig"
replace="2003-2013"
sed -i -e "s|$search|$replace|g" $txtFile

search="month = 12 ;"
replace="time = UNLIMITED ; // (12 currently)"
sed -i -e "s|$search|$replace|g" $txtFile

search="short month(month) ;"
replace="float time(time) ;"
sed -i -e "s|$search|$replace|g" $txtFile

search="Month"
replace="unstructured"
sed -i -e "s|$search|$replace|g" $txtFile

search="month:long_name"
replace="time:CDI_grid_type"
sed -i -e "s|$search|$replace|g" $txtFile

search="(month, lev, lat, lon)"
replace="(time, lev, lat, lon)"
sed -i -e "s|$search|$replace|g" $txtFile

search="month = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ;"
replace="time = 1.111011e+07,1.111021e+07,1.111031e+07,1.111041e+07,1.111051e+07,1.111061e+07,1.111071e+07,1.111081e+07,1.111091e+07,1.111101e+07,1.111111e+07,1.111121e+07;"
sed -i -e "s|$search|$replace|g" $txtFile

ncgen -b $txtFile

rm -rf $txtFile

# name of the CAMS climatology file after the above modifications
sourceFile=${BaseName}.nc

# name of the target grid file. File provided by user
TARGETGRID=${GridName}.nc

# name of the output file with the interpulated CAMS climatology on ICON grid
OFILE=icon_cams_clim_${GridName}.nc

startDate='2001-1-15'
startTime='24:00:00'
timeUnit='1months'

GRIDNUM=`cdo sinfov ${TARGETGRID} | grep nvertex=3 | awk '{print $1}'` 
cdo -s -r -P 4 remapbic,"${TARGETGRID}:${GRIDNUM}" -settaxis,${startDate},${startTime},${timeUnit} ${sourceFile} t1.nc
cdo mul -gec,0.0 t1.nc t1.nc t2.nc
cdo add -mulc,0.0 -ltc,0.0 t1.nc t2.nc ${OFILE}
rm -rf t1.nc t2.nc
echo 'Done'

# If you get error like this:
#    Error (cdf_enddef): NetCDF: One or more variable sizes violate format constraints
# It means that you file is too large.

# solution from:  https://stackoverflow.com/questions/17332353/what-is-the-easiest-way-to-convert-netcdf-to-hdf5-on-windows:
# First, see which netCDF format your file is using.
# if you have netCDF installed on your system,
# the easiest way to do this is to use ncdump:

#   > ncdump -k your_file.nc

# If ncdump returns netCDF-4, or netCDF-4  classic  model, then congratulations,
# you already have an HDF file, as netCDF-4 is the netCDF data model implemented using HDF5 as the storage layer.
# These files can be ready by the HDF library version 1.8 or later, and from what I can tell, pytables.
# If ncdump returns classic or 64-bit offset, then you are using netCDF-3 and will need to convert to netCDF-4
# (and thus HDF). The good news is that if you have netCDF installed, you can do the conversion pretty easily
# (...the conversion can become quite time consuming):

#    > nccopy -k 4 your_file.nc your_file_4c.nc
