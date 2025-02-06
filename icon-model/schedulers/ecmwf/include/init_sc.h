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
#PBS -S /usr/bin/ksh
#PBS -l EC_job_name=ICON_%TASK%
#PBS -o %SMSJOBOUT%
#PBS -j oe
#PBS -m n
#PBS -q ns

##PBS -l EC_job_name=ICON_%TASK%
##PBS -q nf
##PBS -l EC_total_tasks=24
##PBS -l EC_hyperthreads=2
##PBS -l EC_threads_per_task=1
##PBS -o %SMSJOBOUT%
##PBS -j oe
##PBS -m n

##PBS -l EC_job_name=ICON_%TASK%
##PBS -q np
##PBS -S /usr/bin/ksh
##PBS -o %SMSJOBOUT%
##PBS -j oe
##PBS -m n
##PBS -l walltime=%SCWALLCLOCKLIMIT%
##PBS -l EC_total_tasks=96
##PBS -l EC_threads_per_task=4
##PBS -l EC_hyperthreads=2

##PBS -e %SMSJOBOUT%_err
##PBS -q np
##PBS -l EC_total_tasks=48
##PBS -S /usr/bin/ksh
##PBS -l walltime=%SCWALLCLOCKLIMIT%
##PBS -l EC_threads_per_task=2
##PBS -l EC_total_tasks=%SCTOTALTASKS%
##PBS -q %SCJOBCLASS%

#  Loadleveler Directives
#==========================================================================
#@ core_limit = unlimited
##PBS -N ICON_%TASK% (max 15 characters)
##PBS -q np

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
# To give access to job output file in PBS spool on Cray systems
#==========================================================================
 
if [[ $HOST = @(cc*) ]]; then
  _real_pbs_outputfile=/var/spool/PBS/spool/${PBS_JOBID}.OU
  _pbs_outputfile=/nfs/moms/$HOST${_real_pbs_outputfile}
  _running_output=${SMSJOBOUT}.running
  ln -sf $_pbs_outputfile $_running_output
fi

#==========================================================================
#  Define trap handler, including "normal" exit before smscomplete
#==========================================================================

export SMS_SIGNAL_LIST="1 2 3 4 5 6 7 8 9 10 11 12 13 15 24"

trap ERROR 0
trap '{ echo "Killed by a signal"; ERROR ; }' \
     $SMS_SIGNAL_LIST

[[ -d $TMPDIR ]] && cd $TMPDIR

date

