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
#  Loadleveler Directives
#==========================================================================
#@ job_name         = ICON_%TASK%
#@ output           = %SMSJOBOUT%
#@ error            = %SMSJOBOUT%
#@ notification     = never
#@ shell            = /usr/bin/ksh
#@ job_type         = %SCJOBTYPE%
#@ class            = %SCJOBCLASS%
#@ total_tasks      = %SCTOTALTASKS%
#@ wall_clock_limit = %SCWALLCLOCKLIMIT%
#@ resources        = ConsumableMemory(%SCMEM%) ConsumableCpus(%SCCPUS%)
#@ queue
#@ core_limit = unlimited


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

