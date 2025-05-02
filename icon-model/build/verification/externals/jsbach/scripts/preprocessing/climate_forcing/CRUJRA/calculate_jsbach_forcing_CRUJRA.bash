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
# Script to calculate jsbach forcing data from crujra data
#
# The script is based on calculate_jsbach_forcing_with_common_slm_from_remapped_crujra_TRENDY_v9.bash, which were based on 
# previous TRENDY scripts by JN, preprocessing_site_level_forcing_data_JN_140924.tcsh created by Melinda Galfi and adapted by JN
# and preprocessing_change_time_unit_in_forcing_file.ksh by JN
#
# NOTE: new: unit same for both shortwave ("W m-2") and longwave (W/m^2) 
#
#----------------------------------------------------------------------------
# submit as a serial job to mistral -> adapt wall_clock_limit to your needs!
#----------------------------------------------------------------------------
#
###############################################################################
### Batch Queuing System is SLURM
#SBATCH --job-name=calculate_jsbach_forcing_from_remapped_CRUJRA
#SBATCH --output=calculate_jsbach_forcing_from_remapped_CRUJRA.o%j
#SBATCH --error=calculate_jsbach_forcing_from_remapped_CRUJRA.o%j
#SBATCH --partition=prepost
#SBATCH --mem-per-cpu=5120 
#SBATCH --ntasks=1
#SBATCH --account bm0891
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --mail-user=julia.nabel@mpimet.mpg.de # Set your e-mail address

#============================================================================
#============================================================================
# 'User' definitions

# directory where the crujra data is located
workingPath='/work/mj0060/m300316/TRENDY/TRENDY_v10/preprocessed_data/preprocessed_climate/'
inPath='/scratch/m/m300316/trendy_v10/mapped_climate/'

oldPrefix='crujra_v2.2_T63_'
newPrefix='Climate_crujra_v2.2_T63_'

first_year=1921 #1901
last_year=2020 #2020

#---------------------------------------------- constants
newVarName=( longwave qair precip shortwave tmin tmax wspeed ) #wspeed needs to be last item in the list
newUnit=( 'W/m2' 'g/g' 'mm/day' 'W/m2' 'degC' 'degC' 'm/s' )
# old units: 'W/m^2' 'kg/kg' 'mm/6h' 'W m-2' 'Degrees Kelvin' 'Degrees Kelvin' 'm/s' 'm/s' 
varNameInFileName=( dlwrf spfh pre tswrf tmin tmax ugrd vgrd ) # one more entry then newVarName because of windspeed calculation from ugrd and vgrd
varNameInFile=( dlwrf spfh pre tswrf tmin tmax ugrd vgrd )

# list of all variables that should adhere to lower limits (given in the same order in limitsList)
limitVarList=( precip wspeed shortwave longwave qair )
limitsList=( 0.0 0.0 0.0 0.0 0.0 )

kelvin_to_c=273.15
secondsIn6h=21600

fileWithCruncepSLM='/work/mj0060/m300316/TRENDY/TRENDY_v4/input_spinUp/jsbach_T63dCRUNCEP2015_11tiles_5layers_1860.nc'

#ncap2='/sw/aix53/nco-4.0.3/bin/ncap2'

#============================================================================
#============================================================================
# 'Program'

module load nco
module unload cdo
module load cdo/1.7.0-magicsxx-gcc48

currentDate=$(date +'%m/%d/%Y')
currentYear=${first_year}
while [[ ${currentYear} -le ${last_year} ]]; do

    echo 'Currently processing year: '${currentYear}'...'

    # iterate over the required variables
    fileList=''
    varIndex=0
    while [[ ${varIndex} -lt ${#newVarName[@]} ]]; do

        thisNewVarName=${newVarName[varIndex]}
        thisNewUnit=${newUnit[varIndex]}

        thisOldVarFileName=${varNameInFileName[varIndex]}
        thisOldVarNameInFile=${varNameInFile[varIndex]}

        thisOrgFileName=${oldPrefix}${thisOldVarFileName}_${currentYear}.nc
        thisOrgFile=${inPath}/${thisOrgFileName}

        cdo -s -setunit,${thisNewUnit} -chname,${thisOldVarNameInFile},${thisNewVarName} -selvar,${thisOldVarNameInFile} ${thisOrgFile} ${workingPath}/tmp_${thisNewVarName}.nc

        if [[ "${thisNewVarName}" == 'precip' ]]; then
            cdo -s -daysum ${workingPath}/tmp_${thisNewVarName}.nc ${workingPath}/tmp2_${thisNewVarName}.nc
        elif [[ "${thisNewVarName}" == 'tmin' ]]; then
            cdo -s -daymin -subc,${kelvin_to_c} ${workingPath}/tmp_${thisNewVarName}.nc ${workingPath}/tmp2_${thisNewVarName}.nc
        elif [[ "${thisNewVarName}" == 'tmax' ]]; then
            cdo -s -daymax -subc,${kelvin_to_c} ${workingPath}/tmp_${thisNewVarName}.nc ${workingPath}/tmp2_${thisNewVarName}.nc
        elif [[ "${thisNewVarName}" == 'wspeed' ]]; then
            # for wspeed it is expected that thisOldVarFileName contains the second component at the next index 
            (( nextIndex = varIndex + 1 ))
            secondComponentVarName=${varNameInFileName[nextIndex]}
            secondComponentVarNameInFile=${varNameInFile[nextIndex]}
            secondComponentFile=${inPath}/${oldPrefix}${secondComponentVarName}_${currentYear}.nc

            cdo -s -setunit,${thisNewUnit} -chname,${secondComponentVarNameInFile},${thisNewVarName} ${secondComponentFile} ${workingPath}/tmp_${secondComponentVarNameInFile}.nc
            cdo -s -dayavg -sqrt -add -mul ${workingPath}/tmp_${thisNewVarName}.nc ${workingPath}/tmp_${thisNewVarName}.nc \
                -mul ${workingPath}/tmp_${secondComponentVarNameInFile}.nc ${workingPath}/tmp_${secondComponentVarNameInFile}.nc ${workingPath}/tmp2_${thisNewVarName}.nc
        else
            #longwave, shortwave and qair 
            cdo -s -dayavg ${workingPath}/tmp_${thisNewVarName}.nc ${workingPath}/tmp2_${thisNewVarName}.nc
        fi

        #cut to common cruncep slm
        cdo -s -ifthen -selvar,slm ${fileWithCruncepSLM} -selvar,${thisNewVarName} ${workingPath}/tmp2_${thisNewVarName}.nc ${workingPath}/tmp_${thisNewVarName}_${currentYear}.nc

        fileList="${fileList} ${workingPath}/tmp_${thisNewVarName}_${currentYear}.nc"

        (( varIndex = varIndex + 1 ))
    done #while [[ ${varIndex} -lt ${#varInFileName[@]} ]]; do

    # merge and change relative to absolute time axis and convert to cdo netcdf
    cdo -s -a -f nc -merge ${fileList} ${workingPath}/tmp_${newPrefix}${currentYear}.nc
    rm -f ${fileList}

    #--- make each year a leap year 
    #get the 28 of February
    cdo -s -seldate,${currentYear}-02-28 ${workingPath}/tmp_${newPrefix}${currentYear}.nc ${workingPath}/tmp2_${newPrefix}2802${currentYear}.nc
    #set it to 29
    cdo -s -setdate,${currentYear}-02-29 ${workingPath}/tmp2_${newPrefix}2802${currentYear}.nc ${workingPath}/tmp2_${newPrefix}2902${currentYear}.nc

    #merge it 
    cdo -s -mergetime ${workingPath}/tmp_${newPrefix}${currentYear}.nc ${workingPath}/tmp2_${newPrefix}2902${currentYear}.nc ${workingPath}/tmp2_${newPrefix}${currentYear}.nc
    rm -f ${workingPath}/tmp_${newPrefix}${currentYear}.nc
    rm -f ${workingPath}/tmp2_${newPrefix}2902${currentYear}.nc ${workingPath}/tmp2_${newPrefix}2802${currentYear}.nc

    # set the required limits for the different variables
    varIndex=0
    while [[ ${varIndex} -lt ${#limitVarList[@]} ]]; do
        thisVar=${limitVarList[varIndex]}
        echo 'Currently setting the limit for : '${thisVar}'...'

        thisLimit=${limitsList[varIndex]}
        ncap2 -s "where(${thisVar}<${thisLimit}) ${thisVar}=${thisLimit};" ${workingPath}/tmp2_${newPrefix}${currentYear}.nc ${workingPath}/tmp_${newPrefix}${currentYear}.nc
        mv ${workingPath}/tmp_${newPrefix}${currentYear}.nc ${workingPath}/${newPrefix}${currentYear}.nc
        (( varIndex = varIndex + 1 ))
    done #while [[ ${varIndex} -lt ${#varList[@]} ]]; do

    # remove file_name from global attributes
    ncatted -O -a file_name,global,d,, ${workingPath}/${newPrefix}${currentYear}.nc

    # remove several of the global attributes
    ncatted -O -a conventions,global,d,, ${workingPath}/${newPrefix}${currentYear}.nc ${workingPath}/${newPrefix}${currentYear}.nc
    ncatted -O -a id,global,d,, ${workingPath}/${newPrefix}${currentYear}.nc

    # add/replace some
    ncatted -O -a title,global,o,c,"CRUJRA v2.2 based offline forcing for JSBACH" ${workingPath}/${newPrefix}${currentYear}.nc
    ncatted -O -a institution,global,o,c,"Max Planck Institute for Meteorology" ${workingPath}/${newPrefix}${currentYear}.nc

    # replace some
  ncatted -O -a comment,global,o,c,"${currentDate}: Derived from crujra data downloaded 28.07.21 from trendy-v10@trendy.ex.ac.uk. Scripts: remap_CRUJRA.bash and calculate_jsbach_forcing_CRUJRA.bash - julia.nabel@mpimet.mpg.de" ${workingPath}/${newPrefix}${currentYear}.nc -h
  ncatted -O -a contact,global,o,c,"julia.nabel@mpimet.mpg.de" ${workingPath}/${newPrefix}${currentYear}.nc
  ncatted -O -a history,global,d,, ${workingPath}/${newPrefix}${currentYear}.nc

    (( currentYear = currentYear + 1 ))
done #while [[ ${currentYear} -le ${last_year} ]]; do
rm -f ${workingPath}/tmp*.nc

