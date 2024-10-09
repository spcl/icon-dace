#!/bin/ksh
#   
#
# ****** Prepare nemo mouth mask (e.g. EHYPE) from used runoff climatology  **************
# 
#        Note: NEMO style is using an inverted mask, 0=mouth, 1 otherwise
#
DDIR=/work/gg0302/g260122/HD/nemo
DNIN=GCOAST_runoff_clim.nc
DOUT=/work/gg0302/g260122/HD/input/nemo
DNOUT=gcoast_mouth.nc
DNGRID=${DOUT}/NEMO_grid_GCOAST.nc
#
cd $DDIR
cdo lec,0. -yearmax $DNIN m1.nc
ncwa -O -a time m1.nc m2.nc
ncks -O -x -v time m2.nc m3.nc
ncrename -h -O -d lon,x m3.nc
ncrename -h -O -d lat,y m3.nc
ncap -O -s 'nmhd_msk=int(runoff)' m3.nc m4.nc
ncks -O -v nmhd_msk m4.nc ${DOUT}/${DNOUT}
ncks -h -O -v lat,lon $DNGRID grid.nc                 #extract 2D fields of lan and lat
ncks -h -a -v lon,lat -A grid.nc ${DOUT}/${DNOUT}

##ncks -O -x -v runoff m4.nc ${DOUT}/${DNOUT}
#

