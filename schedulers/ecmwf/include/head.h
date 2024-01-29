#!/bin/ksh

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

set -e # stop the shell on first error
set -u # fail when using an undefined variable
set -x # echo script lines as they are executed

# Defines the three variables that are needed for any
# communication with SMS

export SMS_PROG=%SMS_PROG%  # SMS Remote Procedure Call number
export SMSNODE=%SMSNODE%    # The name sms that issued this task
export SMSNAME=%SMSNAME%    # The name of this current task
export SMSPASS=%SMSPASS%    # A unique password
export SMSTRYNO=%SMSTRYNO%  # Current try number of the task

# Tell SMS we have stated
# The SMS variable SMSRID will be set to parameter of smsinit
# Here we give the current PID.

smsinit $$ 

# Defined a error handler

ERROR() {
	set +e        # Clear -e flag, so we don't fail
	smsabort      # Notify SMS that something went wrong
	trap 0        # Remove the trap
	exit 0        # End the script
}

# Trap any calls to exit and errors caught by the -e flag

trap ERROR 0

# Trap any signal that may cause the script to fail

trap '{ echo "Killed by a signal"; ERROR ; }' 1 2 3 4 5 6 7 8 10 12 13 15

