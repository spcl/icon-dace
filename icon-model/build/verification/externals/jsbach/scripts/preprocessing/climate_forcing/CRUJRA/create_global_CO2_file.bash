#!/bin/bash

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2024, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

#-----------------------------------------------------------------------------
# script generating a global CO2 file from global_co2_ann_1700_2020.txt as used for TRENDY_v10 runs
#------------------------------------------------------------------------------
#############################################################################################################

#all files expected to be located in the workingDir!
workingDir='/work/mj0060/m300316/TRENDY/TRENDY_v10/preprocessed_data/'
outFile='global_CO2_for_trendy_v10_from_global_co2_ann_1700_2020.nc'
patternInFile='/work/mj0060/m300316/CRESCENDO/preprocessed_input/CO2_greenhouse_historical.nc'
co2Data='global_co2_ann_1700_2020.txt' 
co2TemplateFile=${workingDir}/tmp_template.nc

startYear=1700
endYear=2020

firstYearFromFile=1700
thisLine=1 #file starts with 1700 i.e. use line 1

#############################################################################################################
module unload cdo
module load cdo/1.7.0-magicsxx-gcc48


# create a CO2 template netcdf file
co2TemplateFile=${workingDir}/tmp_CO2_template.nc
cdo -s -div -selyear,${firstYearFromFile} -selvar,CO2 ${patternInFile} -selyear,${firstYearFromFile} -selvar,CO2 ${patternInFile} ${co2TemplateFile}

# create a co2 netcdf file for each year
currentYear=${startYear}
while [ ${currentYear} -le ${endYear} ]; do
    echo ${currentYear}
    # get the CO2 value for this year
    if [ ${currentYear} -lt  ${firstYearFromFile} ]; then
      thisCO2Value=${constantCO2}
    else 
      line=$(sed -n ${thisLine}p ${workingDir}/${co2Data})
      thisCO2Value=$(echo ${line} | awk '{print $2}')
      echo "read: ${thisCO2Value}"
      (( thisLine = thisLine + 1 ))
    fi

    # write it into the template and set year correctly
    cdo -s -setyear,${currentYear} -mulc,${thisCO2Value} ${co2TemplateFile} ${workingDir}/tmp_CO2_${currentYear}.nc

    (( currentYear = currentYear + 1 ))
done #while [[ ${currentYear} -le ${last_year} ]]; do
rm -f ${co2TemplateFile}

# merge CO2 files
cdo -s -mergetime ${workingDir}/tmp_CO2_*.nc ${outFile}
rm -f ${workingDir}/tmp_CO2_*.nc

