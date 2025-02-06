#!/bin/ksh                                        #@ shell = /usr/bin/ksh

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#==========================================================================
#  SLURM Directives:
#    qos        quite: normal or large
#    time       time limit h:m:s
#    job-name   job_name
#    output     output
#    error      error
#    mail-type  notification 
#==========================================================================

#SBATCH --qos=normal
#SBATCH --time=06:00:00
#SBATCH --job-nam=ICON_%TASK%
#SBATCH --output=%SMSJOBOUT%
#SBATCH --error=%SMSJOBOUT%
###SBATCH --mail-type=FAIL


#==========================================================================
#  shell settings
#==========================================================================

set -e
set -u
set -x

SUITE=%SUITE%
TASK=%TASK%
SMSPASS=%SMSPASS%
SMSNODE=%SMSNODE%
SMSNAME=%SMSNAME%
SMSHOME=%SMSHOME%
SMS_PROG=%SMS_PROG%
WSHOST=%WSHOST%
SCHOST=%SCHOST%
USER=%USER%
SMSJOBOUT=%SMSJOBOUT%
LOGDIR=%SMSOUT%

#==========================================================================
#  SMS Definitions
#==========================================================================

LOADL_STEP_ID=${LOADL_STEP_ID:=NOT_SET}
QSUB_REQID=${QSUB_REQID:=NOT_SET}

if [[ $LOADL_STEP_ID != NOT_SET ]] ; then
   SMSRID=$(echo $LOADL_STEP_ID | cut -f2 -d.)
   JOB_ID=$LOADL_STEP_ID
elif [[ $QSUB_REQID != NOT_SET ]] ; then
   SMSRID=$(echo $QSUB_REQID | cut -f1 -d.)
   JOB_ID=$QSUB_REQID
else
  SMSRID=$$
  JOB_ID=$SMSRID
fi

export SMSPASS SMSNODE SMSNAME SMS_PROG SMSHOME SMSJOBOUT

export WSHOST SCHOST USER SUITE

#==========================================================================
#  Tell sms that it started
#==========================================================================

smsinit $SMSRID

#==========================================================================

#==========================================================================
#  Define error handler
#==========================================================================

ERROR() {
  set -x
  set +e
  wait
  smsabort 
  trap 0
  date
  times
  echo "environment was:"
  printenv | sort
  exit
}

#==========================================================================
#  Define trap handler, including "normal" exit before smscomplete
#==========================================================================

export SMS_SIGNAL_LIST="1 2 3 4 5 6 7 8 9 10 11 12 13 15 24"

trap ERROR 0
trap '{ echo "Killed by a signal"; ERROR ; }' \
     $SMS_SIGNAL_LIST

[[ -d $TMPDIR ]] && cd $TMPDIR

date


