#!/bin/sh
#
# hd_put.sh - Moves HD output from scratch to work on levante incl. some gzip compression. 
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
# Puts all data of a HD model run into work and compresses with gzip
# It also transfers some files to another machine
#
set -ex
#
while [[ -n $1 ]] ; do
  case $1 in
    --id | -i         ) exp=$2      ;;
    -*                  ) err_strg="\nERROR: invalid option $1 !\n"
                          [[ "$1" != $( echo "$1" | tr -d = ) ]] && \
                          err_strg=$err_strg"  Please use blank instead of '=' to separate option's argument.\n"
                          usage ;;
    *                   ) export CPL_MOD=$1 ;;
  esac
  shift
done
#
# *** Provide experiment ID via command line 
if [ "$exp" = "" ]; then
   echo "No Exp. ID is set. Please use Option -i <exp ID>"
   echo "  --> e.g. 7055045. --> TERMINATION!"
   exit 
fi
#
#exp=7062007

DRUN=/scratch/g/g260122/hd
DSAV=/work/gg0302/g260122/HD/output/$exp
DHZG=kopie
#
cd ${DRUN}/${exp}/out
set +ex
mkdir $DSAV
set -ex
#
mv *${exp}* $DSAV
#
# Archive Executable
cd ..
EXE=`ls ${exp}*.exe`
if [ -e ${DSAV}/${EXE} ] ; then
  echo 'Executable already archived'
else
  cp -p ${EXE} ${DSAV}
fi
cd $DSAV
gzip *${exp}*.log ${exp}_hd*restart*.nc
#

exit
