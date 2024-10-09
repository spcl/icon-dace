# weight_vtk.sh - Script to generate weight files for HD model coupled to YAC
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Author: Ha Ho-Hagemann (Hereon, Germany)
# Contact: ha.hagemann@hereon.de
#_________________________________________
#
YAC_dir="/work/gg0302/g260062/GCOAST_oas5/yac/contrib"

var="RDC2NEMO"
var="RUNOFF_S"
var="RUNOFF_G"

workdir="/scratch/g/g260062/hd/hd_gl05deg"
wgtin="${workdir}/${var}.nc"
wgtout="${workdir}/${var}.vtk"
icongrid="/work/gg0302/g260062/HD/run_hd/input/icon_grid_0020_R02B05_O.nc"  # this was used to generate the weight file
hd_grid="-179.99,89.99,179.99,-89.99,720,360"
#hd_grid="-179.75,89.75,179.75,-89.75,720,360"

echo "icongrid=" ${icongrid}
echo "hd_grid="  ${hd_grid}

#----------------------------------------------------------------------------------
cd ${YAC_dir}

if [ $var == "RDC2NEMO" ];then
 echo "./weights2vtk.x -S g -s ${hd_grid} -T i -t ${icongrid} -w ${wgtin} -o ${wgtout}"
 ./weights2vtk.x -S g -s ${hd_grid} -T i -t ${icongrid} -w ${wgtin} -o ${wgtout}
else
 echo "./weights2vtk.x -T g -t ${hd_grid} -S g -s ${hd_grid} -w ${wgtin} -o ${wgtout}"
 ./weights2vtk.x -T g -t ${hd_grid} -S g -s ${hd_grid} -w ${wgtin} -o ${wgtout}
fi

#------------------------------------------------------------------------------------------------------------
####./weights2vtk.x -S g -T i -s -11,72,69,27,960,540 -t icon_grid_0020_R02B05_O.nc -w RDC2NEMO.nc -o RDC2NEMO.vtk
