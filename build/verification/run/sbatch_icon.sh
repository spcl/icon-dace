#!/bin/bash

# Take the first argument ($1) from the command line. 
# If $1 is empty, the :? check will crash with the error message.
export XTSTEP=${1:? "Error: Please provide a time step as an argument. Usage: ./submit.sh <XTSTEP> <XICONV>"}
export XICONV=${2:? "Error: Please provide a icon version as an argument. Usage: ./submit.sh <XTSTEP> <XICONV>"}

# Construct the experiment name with the suffix
export EXPNAME="exclaim_ape_R02B04_dt${XTSTEP}_v${XICONV}_0008_R02B05_G"

# Submit to Slurm, overriding the internal #SBATCH headers
sbatch --job-name="${EXPNAME}" \
       --output="/capstor/scratch/cscs/pmazumde/gitspace/icon-dace/build/verification/run/LOG.SAVEME.${EXPNAME}.%j.o" \
       --error="/capstor/scratch/cscs/pmazumde/gitspace/icon-dace/build/verification/run/LOG.SAVEME.${EXPNAME}.%j.o" \
       --export=ALL,XTSTEP=${XTSTEP},XICONV=${XICONV} \
       exp.exclaim_ape_R02B04.run

