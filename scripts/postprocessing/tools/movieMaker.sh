# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

set -x
remoteDir="/work/mh0287/users/stephan/Icon/icon-dev.tst/experiments/xmpiom.r11009.b5.2"
remoteFilePattern="$remoteDir/Output/xm*nc $remoteDir/xm*nc"
localFilePattern=xm*.nc
weights="weights_r2b05_to_r360x180.nc"
plotToolDir="/work/mh0287/users/ram/icon/scripts/postprocessing/tools"
plotTool="icon_plot.ncl"
ls $remoteFilePattern
#== preselect
for file in $remoteFilePattern; do echo "cdo -sellevel,100 -selname,u,v,wet_c $file _$(basename $file); mv _$(basename $file) $(basename $file)";done > sel.sh
prun.rb -d sel.sh 
#== velocity add
for file in $localFilePattern; do echo "cdo -setunit,'m/s'  -expr,'velocity=sqrt(u*u+v*v)' $file veloc_$file; cdo merge $file veloc_$file _$file; mv _$file $file; rm veloc_$file;"; done > vel.sh
prun.rb -d vel.sh 
cdo cat $localFilePattern r2b05.nc
cdo -intntime,5 r2b05.nc _r2b05.nc; mv _r2b05.nc r2b05.nc
cdo -P 20 gennn,r360x180 -seltimestep,1 r2b05.nc $weights
# pre-interpolation for fast plotting
cdo -O remap,r360x180,$weights  r2b05.nc remapnn_r360x180_r2b05.nc
# splitting for fast plotting
mkdir split
cd split/
export CDO_FILE_SUFFIX=.nc
cdo -v splitsel,1 ../remapnn_r360x180_r2b05.nc remapnn_r360x180_r2b05_
cdo -v splitsel,1 ../r2b05.nc r2b05_

# perform plot
module load ncar/ncl-6.1.0
for file in r2*_0*.nc; do echo "nclsh $plotToolDir/$plotTool -altLibDir=$plotToolDir -iFile=$file  -varName=velocity -vecVars=u,v -oFile=$file -oType=png -maskName=wet_c  -mapType=ortho +mapLine -isIcon -resolution=r360x180 -centerLon=-40 -centerLat=0 -selMode=halflog -minVar=0.01 -maxVar=1 -streamLine -rStrg=' ' -bStrg=' ' +withLines; "; done > plot.sh
prun.rb -d plot.sh 

# make movie
/work/mh0287/users/ram/local/bin/mencoder "mf://r2*png" -mf fps=15 -o r2b05-tInt-correctForcing.avi  -profile mpeg4-hq
/work/mh0287/users/ram/local/bin/mencoder "mf://r2*png" -mf fps=25 -o r2b05-tInt-correctForcing-25fps.avi  -profile mpeg4-hq
