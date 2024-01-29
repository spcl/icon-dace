#! /bin/ksh

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#---------------------------------------------------------------------
# Read ICON initialized analysis data (0h forecast) for initialization of 
# ICON forecasts at ECMWF
#
# Info:
#   - requires eccert dei2
#   - runs on rcl on command line
#
# Martin Koehler, Nov. 2014
#---------------------------------------------------------------------

set -ex

inidate_start=2021061700
inidate___end=2021061700

res=R02B06

datadir=/hpc/uwork/mkoehler/run-icon/icon_init

#---------------------------------------------------------------------

inidates=$(eval_numlist -d -i 24 ${inidate_start}-${inidate___end})
integer inidate2

# ecaccess-file-mkdir ec:/dei2/icon/icon_input/icon_ana || true

cd ${datadir}

for inidate in ${inidates[*]} ; do

  iconfile=icon2icon_${res}_G_${inidate}.grb
  dm1=`date_calc.py -d ${inidate} -s 1`
  inidate2=dm1+21                              # date - 3h

  cat > sky_icon_request << EOF
    reqColl proc=parallel timeout=0 
    read db=roma cat=icogle_ass_fc_rout letype=101 enum=1 d=${inidate2} s[s]=10800  bin info=countPlus f=${iconfile} lvt1=!100
    read db=roma cat=icogle_ass_an_rout letype=101 enum=1 d=${inidate}  gptype=0    bin info=countPlus f=${iconfile} p=w_so lin=0
    read db=roma cat=icogle_ass_an_rout letype=101 enum=1 d=${inidate}  gptype=0    bin info=countPlus f=${iconfile} p=w_snow,t_snow,w_i,rho_snow
    read db=roma cat=icogle_ass_an_rout letype=101 enum=1 d=${inidate}  gptype=0    bin info=countPlus f=${iconfile} p=freshsnw,h_snow
EOF

  sky -v sky_icon_request
# ecaccess-file-put ${iconfile} ec:/dei2/icon/icon_input/icon_ana/${iconfile}

done


exit

#    read db=numex cat=icogl400l90_ass_an_exp1 exp=${expid} letype=101 enum=1 d=${inidate}  gptype=!201 bin info=countPlus f=${iconfile} p=u,v,t,p,qv
#    read db=numex cat=icogl400l90_ass_an_exp1 exp=${expid} letype=101 enum=1 d=${inidate}  gptype=-1   bin info=countPlus f=${iconfile} p=t_so,fr_ice,h_ice,t_ice
