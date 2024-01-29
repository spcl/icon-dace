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

#------------------------------------------------------------
# convert ICON SCM netcdf output to DEPHY v0 format for plotting
# with DephyDiags tool
#
# 202006 Martin Koehler, DWD
#------------------------------------------------------------
#

module load cdo
module load netcdf
module load nco

dir=/hpc/uwork/mkoehler/icon-run/scm/SCM_ARM/
#dir=/hpc/uwork/mkoehler/icon-run/scm/SCM_ARM_uf/

iconfile=scm_out_ML_19970621T000000Z.nc

cd $dir

# mean of 32 points

cdo fldmean $iconfile scm_out_ensmean.nc

# change names  icon                dephy
#----------------------------------------
#               qv or tot_qv_dia    qv
#               clc                 rneb
#               theta=T/exner       theta
#               z_mc                zf

cdo chname,z_mc,zf        scm_out_ensmean.nc temp1.nc
cdo chname,z_ifc,zh       temp1.nc           temp2.nc

cdo expr,rneb=clc/100.0   scm_out_ensmean.nc temp3.nc
cdo expr,theta=temp/exner scm_out_ensmean.nc temp4.nc
cdo merge        temp2.nc temp3.nc  temp4.nc temp5.nc

ncrename -d height_2,levf temp5.nc           temp6.nc
ncrename -d height_3,levh temp6.nc           icon-scm.nc

# clean

rm -rf temp0.nc temp1.nc temp2.nc temp3.nc temp4.nc temp5.nc

