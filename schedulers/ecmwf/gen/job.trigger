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

# @ shell    = /bin/ksh
# @ job_name = ICON_TRIGGER
# @ class    = express
##@ class    = ts
# @ job_type = serial
# @ notification = error
# @ notify_user  = Florian.Prill@dwd.de
# @ queue

# submitted with
# ecaccess-job-submit -noDirectives -eventIds 167 -queueName ecgate /home/ms/de/dfi0/ICON_r2B06_10d/def/smsfiles/job.trigger 


set -x

date

# trigger ICON SMS suite ------------------------------
cdp << ENDCDP
set SMS_PROG 903409
myalias
force complete /ifs_trigger/fct_ifs
#force queued  /icon/forecast/check_progress
status -a
exit
ENDCDP

exit
