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

# terminate sms:
cdp << EOF
   myalias
   cancel -fy icon
   halt -y
   terminate -y
EOF

ps -ef | grep sms | grep ${USER}

rm sms.check*
rm icon/*.job*
rm -rf $HOME/sms/*
rm -rf $HOME/sms_server