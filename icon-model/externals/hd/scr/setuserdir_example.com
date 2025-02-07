#!/bin/ksh
#
# setuserdir_example.com - Example for setuserdir.com (sub-script of run_hdmodel.ksh) with setting the main HD directories
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
#
export HDMAIN=/pf/g/g260122/hdmodel               # HD main Directory that includes script, grid and log directories.
export HDDIR=/scratch/g/g260122/hd                # Run Directory
export HDFILE=/work/gg0302/g260122/HD/input       # HD input data Directory, e.g. HD parameter and start files.
export HDFORCING=/work/gg0302/g260122/HD/forcing  # HD directory with forcing data in subdirectories
#
ierrset=0
if [ -e $HDMAIN ];then echo "$HDMAIN exists" ; else
  echo "$HDMAIN does NOT exists!"
  echo "   --> You need to create it"
  ierrset=1
fi
if [ -e $HDFILE ];then echo "$HDFILE exists" ; else
  echo "$HDFILE does NOT exists!"
  echo "   --> You need to create it"
  ierrset=1
fi
if [ -e $HDDIR ];then echo "$HDDIR exists" ; else
  echo "$HDDIR does NOT exists!"
  echo "   --> You need to create it"
  ierrset=1
fi
if [ -e $HDFORCING ];then echo "$HDFORCING exists" ; else
  echo "$HDFORCING does NOT exists!"
  echo "   --> You need to create it"
  ierrset=1
fi
if [ ${ierrset} == 1 ]; then
  echo "Some HD directories does not exist --> TERMINATION"
  export ierrset=1
else
  export ierrset=0
fi

