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
# Script to remap CRUJRA climate data from 0.5deg (r720x360)
#
# The script is based on remap_to_T63_common_weight_TRENDY_v10.bash down to remap_to_T63_TRENDY_v4.bash, which in turn were based on
# Change_grid_from_T63_to_R720x360_140108.tcsh (Author: J. Nabel, 2014/01/07; modified by J. Pongratz, 2014/01/08)
# which had the purpose to remap all files to be handed in for trendy v2 from T63 to 0.5deg x 0.5deg (r720x360)
#
# The script calculates the weights once with conservative remapping
# and then remaps all files with given postfix and creates new files with prefix "crujra_v2.2_T63_"
# ++ in contrast to cruncep no grid creation of shifts are required 
#    (cruncep climate data did not suit the 'r720x360' cdo grid at all, seems as crujra does)
#
#----------------------------------------------------------------------------
# submit as a serial job -> adapt wall_clock_limit to your needs!
#----------------------------------------------------------------------------
#
###############################################################################
### Batch Queuing System is SLURM
#SBATCH --job-name=remap_CRUJRA
#SBATCH --output=remap_CRUJRA.o%j
#SBATCH --error=remap_CRUJRA.o%j
#SBATCH --partition=prepost
#SBATCH --mem-per-cpu=5120
#SBATCH --ntasks=1
#SBATCH --account mj0060
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --mail-user=julia.nabel@mpimet.mpg.de # Set your e-mail address

#============================================================================
#============================================================================
# 'User' definitions

module load cdo

# directory where the crujra data is located
inpath='/scratch/m/m300316/raw_data/climate_forcing/CRUJRA2022/'

# prefix of the original files and for the new file names
orgPrefix='crujra.v2.3.5d.'
orgPostfix='.365d.noc.nc'
newPrefix='crujra_v2.3_R2B4_'

cdoGrid=/pool/data/ICON/grids/public/mpim/0043/icon_grid_0043_R02B04_G.nc

# directory where the new remapped data should be located
outpath='/scratch/m/m300316/mapped_climate/'

# variable names in the file names
varNameInFileName=(  dlwrf pre tmax tmin spfh vgrd ugrd  )
#varNameInFileName=( tswrf )

# and variable names in the files
varNameInFile=( dlwrf pre tmax tmin spfh vgrd ugrd )
#varNameInFile=( tswrf )

# time span for which the remapping should be done
first_year=1901
last_year=2021

#ncap2='/sw/aix53/nco-4.0.3/bin/ncap2'

#============================================================================
#============================================================================
# 'Program'

module load nco


#file for the weights
#remap_weight_file_name='tmp_remap_weights_tswrf.nc'
remap_weight_file_name='tmp_remap_weights.nc'

# iterate over the variables
varIndex=0
while [[ ${varIndex} -lt ${#varNameInFileName[@]} ]]; do

    thisVar=${varNameInFileName[varIndex]}
    thisVarName=${varNameInFile[varIndex]}

    # iterate over the years to be remapped
    currentYear=${first_year}
    while [[ ${currentYear} -le ${last_year} ]]; do

        thisFileName=${orgPrefix}${thisVar}.${currentYear}${orgPostfix}

        thisFile=${inpath}/${thisFileName}
           
        # generate weights if required
        [[ ! -f ${outpath}${remap_weight_file_name} ]] && cdo -s -P 2 gencon,${cdoGrid} ${thisFile} ${outpath}${remap_weight_file_name}

        # remap
        cdo -s -P 2 -remap,${cdoGrid},${outpath}${remap_weight_file_name} ${thisFile} ${outpath}${newPrefix}${thisVar}_${currentYear}.nc

        (( currentYear = currentYear + 1 ))
    done #while [[ ${currentYear} -le ${last_year} ]]; do

    (( varIndex = varIndex + 1 ))
done #while [[ ${varIndex} -lt ${#varInFileName[@]} ]]; do
   
rm -f ${outpath}'tmp_T63_grid.nc'
rm -f ${outpath}'tmp_remap_weights.nc'



