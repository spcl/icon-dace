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

#  Calculate barotropic stream function from icon ocean model output
#
#  Author: Stephan Lorenz, MPIfMet, 06/2012
#
#  Input : averaged icon-ocean standard output file including variables u_vint and wet_c 
#  Output: interpolated input/output/plot-files named u_vint/psi/nclpsi
#
#  Method:
#   - vertical interpolation is done in the model, variable u_vint
#   - horizontal interpolation done by cdo remapcon, conservative remapping
#   - convert to SERVICE format, calculate meridional integral incl. using
#     land-sea-mask, fortran program psiread
#   - convert badk to NETCDF format using grid description file
#   - simple plot using NCL with shell-script nclsh
#          

set -e

# Resolution:
# cdo remapnn and remapdis generate shaky isolines
# cdo remapcon generates missing values at Poles if resolution of regular target grid is too fine
#   R2B04-r180x90 and R2B05-r360x180 match well

### input parameter
avgfile=${1}             #  file including wet_c and average of u_vint
resol=${2:-r180x90}      #  default resolution: 2x2 degrees, 180x90 gridpoints
filestr=${3:-$(basename $1)}         #  pattern index for filenames
weightfile=${4}          #  optional weights for interpolation

# Default name of input file is avg.$filestr
if [ "$filestr" == "$avgfile" ]; then
  xtag=${avgfile##avg.}
  filestr=${xtag%%.nc}
fi

echo "PSI: Input file is '$avgfile'"
echo "PSI: Tag is '$filestr'"
echo "PSI: Resolution is '$resol'"
echo "PSI: Weightsfile is '$weightfile'"
### input parameter end

#resol=r${nlon}x${nlat}
#resol='r180x90'

# file names
inpfile=uvint.${resol}.$filestr
outfile=psi.${resol}.$filestr
plotfile=nclpsi.$filestr

rlon=${resol%x*}
nlon=${rlon#r}
nlat=${resol#*x}

echo "PSI: Working on files: input: $avgfile ; output: $outfile ; plot: $plotfile"

# select wet_c (lsm at surface) and u_vint (vertical integral of u)
#  - do interpolation if neither inpfile nor outfile is available
#  - sellevidx: select first level index at surface
if [ ! -s "$inpfile.nc" ] && [ ! -s "$outfile.nc" ]; then
  echo "PSI -> no result file available, do interpolation"
  if [ -z "$weightfile" ]; then
    echo "PSI -> weightfile is not given! Use remapcon for u_vint and remapnn for wet_c ..."
    #cdo -P 8 remapdis,$resol -selvar,u_vint $avgfile xsrc.${resol}_uint.nc
    cdo -P 8 setmisstoc,0.0 -remapcon,$resol -selvar,u_vint $avgfile xsrc.${resol}_uint.nc
    cdo -P 8 remapnn,$resol  -selvar,wet_c -sellevidx,1 $avgfile xsrc.${resol}_wetc.nc
    #cdo merge xsrc.${resol}_wetc.nc xsrc.${resol}_uint.nc $inpfile.nc
    cdo merge xsrc.${resol}_uint.nc xsrc.${resol}_wetc.nc $inpfile.nc  # -O
    # Memory fault in r2b5?
    #cdo -P 8 merge -remapnn,$resol  -selvar,wet_c -sellevidx,1 $avgfile -selvar,u_vint $avgfile $inpfile.nc
    rm xsrc.${resol}_wetc.nc xsrc.${resol}_uint.nc
  else
    echo "PSI -> weightfile is given:'$weightfile'. Use remap ..."
    cdo -P 8 remap,$resol,$weightfile -selvar,u_vint,wet_c, $avgfile $inpfile.nc
  fi
else
  echo "PSI -> Use input for the psi-computation:'$inpfile.nc'"
fi

# if outfile is available do plotting only:
if [ -s "$outfile.nc" ]; then
  echo "PSI -> run plot program only "
else

# store grid description, convert to service format
echo "PSI -> cdo -v griddes $inpfile.nc > griddes.$resol"
cdo griddes $inpfile.nc > griddes.$resol
echo "PSI -> cdo -v -f srv copy $inpfile.nc $inpfile.srv"
cdo -f srv copy $inpfile.nc $inpfile.srv
#rm $inpfile.nc

cat > scr-psiread.f90 <<EOF
!-------------------------------------------------------------------------  
!
!
!!  Calculation of horizontal stream function
!
!>
!!
!! @par Revision History
!! Developed  by  Stephan Lorenz, MPI-M (2012).
!!  based on code from MPIOM
!   ignore vertical dimension
!
! TODO: diffuse from ocean to land, cut land points
! TODO: implement variable output dimension (1 deg resolution) and smoothing extent
!! 
PROGRAM psiread

IMPLICIT NONE

INTEGER, PARAMETER ::  rho_ref = 1025.022            ! reference density

INTEGER, PARAMETER ::  nlat = $nlat                  ! meridional dimension of regular grid
INTEGER, PARAMETER ::  nlon = $nlon                  ! zonal dimension of regular grid

! smoothing area is 2*jsmth-1 lat/lon areas of 1 deg
INTEGER, PARAMETER ::  jsmth = 3                  
INTEGER            :: jb, jc, i_startidx, i_endidx
INTEGER            :: jlat, jlon, jx, jy
INTEGER            :: isrv(8), isrvu(8)


REAL               :: z_lat_dist, erad, pi
REAL               :: z_uint_reg(nlon,nlat)     ! vertical integral on regular grid
REAL               :: psi_reg(nlon,nlat)        ! horizontal stream function
REAL               :: wet_c(nlon,nlat)          ! slm

!CHARACTER(len=max_char_length), PARAMETER :: routine = ('mo_oce_diagnostics:calc_psi')

!-----------------------------------------------------------------------

psi_reg(:,:)    = 0.0
z_uint_reg(:,:) = 0.0

! (2) read barotropic system: after cdo merge above, now read first uint then wet_c

  open (11,file="$inpfile.srv", form='unformatted')
  open (80,file="$outfile.srv", form='unformatted')

  read (11) isrvu
  write (*,*) isrvu
  read (11) z_uint_reg(:,:)

  read (11) isrv
  read (11) wet_c(:,:)

  write(80) (isrv(jb),jb=1,8)
  write(80) ((wet_c(jlon,jlat),jlon=1,nlon),jlat=1,nlat)


  ! (3) calculate meridional integral on regular grid starting from north pole:
  DO jlat = nlat-1, 1, -1
    z_uint_reg(:,jlat) = z_uint_reg(:,jlat) + z_uint_reg(:,jlat+1)
  END DO

  ! (4) calculate stream function: scale with length of meridional resolution:
  erad = 6.371229e6                 !  earth's radius [m]
  pi   = 3.141592653
  z_lat_dist = pi/real(nlat)*erad   !  z_lat_dist = dlat* pi*R/180 ; dlat=180/nlat

  ! multiply with -1: Gulf/Kuroshio and ACC positive
  psi_reg(:,:) = -z_uint_reg(:,:) * z_lat_dist * wet_c(:,:) * 1.0e-6


  write(80) (isrvu(jb),jb=1,8)
  write(80) ((psi_reg(jlon,jlat),jlon=1,nlon),jlat=1,nlat)


END PROGRAM psiread
!-------------------------------------------------------------------------  
EOF

echo "PSI -> compile and run program scr-psiread.x"
gfortran -o scr-psiread.x scr-psiread.f90
./scr-psiread.x
rm scr-psiread.* $inpfile.srv
if [ "$resol" == "r180x90" ] ; then rm $inpfile.nc; fi

# convert back to netcdf
echo "PSI -> cdo -f nc -g griddes.$resol chvar,var4,psi -chvar,var1,wet_c $outfile.srv $outfile.nc"
cdo -f nc -g griddes.$resol chvar,var4,psi -chvar,var1,wet_c $outfile.srv $outfile.nc
rm $outfile.srv griddes.$resol

fi  #  run fortran program

# plot with nclsh:

echo "PSI -> plot using icon_plot.ncl:"
nclsh /pool/data/ICON/tools/icon_plot.ncl -altLibDir=/pool/data/ICON/tools \
  -iFile=$outfile.nc -oFile=$plotfile -varName=psi -timeStep=0 -oType=png \
  -maskName=wet_c -selMode=manual -minVar=-150 -maxVar=150 -numLevs=15 -tStrg='r11413 dtime=1h' \
  -plotLevs=-150,-100,-75,-50,-30,-20,-10,-5,0,5,10,20,30,50,75,100,150 -withLineLabels > /dev/null

# nclsh $ICONPLOT \
#   -iFile=$outfile.nc -oFile=$plotfile -varName=psi -timeStep=0 -oType=ps \
#   -selMode=manual -minVar=-60 -maxVar=-10 -numLevs=25 -bStrg=' ' -maskName=wet_c \
#   -withLineLabels \

# -mapLLC=-77,20 -mapURC=-70,35 #north-atlantic gyre
#  -mapLLC=120,10 -mapURC=160,40 # japanise WBC
# -plotLevs=-150,-100,-75,-50,-30,-20,-15,-10,-5,0,5,10,15,20,30,50,75,100,150 -withLineLabels
# -plotLevs=-150,-100,-75,-50,-30,-20,-10,-5,0,5,10,20,30,50,75,100,150 -withLineLabels
# -maxView \

# -selMode=manual -minVar=-240 -maxVar=210 -numLevs=15
# -maskName=wet_c -selMode=manual -minVar=-240 -maxVar=210 -numLevs=15 \
# -maskName=wet_c -selMode=manual -minVar=-100 -maxVar=100 -numLevs=20 \

exit

# -selMode=halflog -scaleLimit=1

