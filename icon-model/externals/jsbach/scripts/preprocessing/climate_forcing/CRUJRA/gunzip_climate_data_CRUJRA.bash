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

#----------------------------------------------------------------------------
# submit as a serial job to mistral -> adapt wall_clock_limit to your needs!
#----------------------------------------------------------------------------
#
###############################################################################
### Batch Queuing System is SLURM
#SBATCH --job-name=gunzip_climate_data_CRUJRA
#SBATCH --output=gunzip_climate_data_CRUJRA.o%j
#SBATCH --error=gunzip_climate_data_CRUJRA.o%j
#SBATCH --partition=prepost
#SBATCH --mem-per-cpu=5120 
#SBATCH --ntasks=1
#SBATCH --account bm0891
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --mail-user=julia.nabel@mpimet.mpg.de # Set your e-mail address

thisUnzipPath=/scratch/m/m300316/trendy_v10/raw_data/

var=$1

# unzip zipped files
#gzfiles=$(ls ${thisUnzipPath}/*.gz)
gzfiles=$(ls ${thisUnzipPath}/*${var}.*gz )

for file in ${gzfiles[@]}; do
    gunzip ${file}
done

