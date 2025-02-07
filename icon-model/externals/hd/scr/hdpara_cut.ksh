#!/bin/ksh
#
# hdpara_cut.ksh - Cut out a region from the global 5 Min-HD parameter file 
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
#   Cutout region from global HD parameter file at 5 Min. spatial resolution 
#
#
# Parse command line parameters
#
while [[ -n $1 ]] ; do
  case $1 in
    --file     | -f     ) file=$2                ; shift    ;;
    --help     | -h     ) help=1                 ; break    ;;
    --region   | -r     ) IREG=$2                           ;;
    -*                  ) err_strg="\nERROR: invalid option $1 !\n"
                          [[ "$1" != $( echo "$1" | tr -d = ) ]] && \
                          err_strg=$err_strg"  Please use blank instead of '=' to separate option's argument.\n"
                          usage ;;
    *                   ) export cline=$1 ;;
  esac
  shift
done

echo 'Programm cuts out a predefined region from the global 5 Min HD model parameter file'
if [ "$file" = "" ] || [ "$help" = "1" ]; then
   if [ "$help" != "1" ] ; then echo "Input HD parameter file not set --> This is necessary!" ; fi
   echo " "
   echo "Help-Listing: Run with option -h or --help or without any option"
   echo " "
   echo "You need to use the option -f or --file to specify the global 5 Min HD model parameter file"
   echo " "
   echo "Options for resolution/domain settings with -r/--region. Please use Option -r <No. or name>"
   echo " No.: 2  = euro5min  -->   European 5 Min. domain: -11°W-69°E, 27-72°N"
   echo "      3  = aus5min   --> Australian 5 Min. domain:  112-154°E, 10-45°S"
   echo "      4  = sea5min   --> South East Asian 5 Min. domain:  89-154°E, 37°N-45°S"
   echo "     50  = medwrf    --> Mediterranean WRF domain:  -16°W-50°E, 24°N-64°N"
   exit
fi
#
case ${IREG} in
  2 | euro5min ) cdocom="--no_history selindexbox,2029,2988,217,756"
                 cdoarea="--no_history selindexbox,1,1,217,756" 
                 TAG=euro5min
                 lon0=2028
                 lat0=216  ;;
  3 | aus5min )  TAG=aus5min    # Australia: 112-154 E, 10-45 S
                 cdocom="--no_history selindexbox,3505,4008,1201,1620"
                 cdoarea="--no_history selindexbox,1,1,1201,1620" 
                 lon0=3504
                 lat0=1200  ;;
  4 | sea5min )  TAG=sea5min    # South East Asia plus Australia
                 cdocom="--no_history -sellonlatbox,89,154,-45,37"
                 cdoarea="--no_history selindexbox,1,1,637,1620"    # lon: no dim -> 1:1, Lat: lat0:lat0+Range in degree * 12 
                 lon0=3228      # = (180 + 89) * 12
                 lat0=636   ;;  # = (90 - 37) * 12
 50 | medwrf  )  TAG=medwrf    # Mediterranean WRF domain
                 cdocom="--no_history -sellonlatbox,-16,50,24,64"
                 cdoarea="--no_history selindexbox,1,1,313,792"     # lat: 313:312+480=792
                 lon0=1968      # = (180 -16) * 12
                 lat0=312   ;;  # = (90 - 64) * 12
  * ) echo "Region not defined"  ; exit  ;;
esac
#
# Read paths definitions from file setuserdir.com
#     HDMAIN              # HD main Directory that includes script, grid and log directories.
#     HDDIR               # Run Directory
#     HDFILE              # HD input data Directory, e.g. HD parameter and start files.
#     HDFORCING           # HD directory with forcing data in subdirectories
source setuserdir.com
#
cd $HDDIR
#
DNAM=$file         # Input Netcdf file
dn1=${DNAM##*/}   # File name without path
dnsub=${dn1%.*}   # Remove all characters until the '.' from the end of the string backwards.
DNOUT=${HDFILE}/${dnsub}_${TAG}.nc
#
#     -v is copying, -x excludes specified variables from the copying.
#     -v copys specified variables to out.nc, -A appends these variabes to old file
ncks -h -O --no_abc -x -v AREA,FIBNEW,FILNEW,FDIR,ALF_K,ALF_N,ARF_K,ARF_N,AGF_K $DNAM ex.nc
cdo $cdocom ex.nc $DNOUT
cdo $cdoarea -selvar,AREA $DNAM area.nc
ncks -h --no_abc -v AREA -A area.nc $DNOUT

cdo $cdocom -selvar,FILNEW,FIBNEW,FDIR,ALF_K,ALF_N,ARF_K,ARF_N,AGF_K $DNAM direction.nc

cat > corr_dir.f90 << EOF
  PROGRAM corr_dir
!
!   ******** Programm that corrects target indices for the shift in origin, 
!   ******** and outflows at the border of the cutted region become artificial sinks.
!
!   ******** Vs.1 - Dec. 2017 - HAG
!
    use netcdf

    IMPLICIT NONE
!
!   *** NETCDF variables
    INTEGER :: ierr, ncid, lat_dimid, lon_dimid, dimid, dimids(2), &
                 lon_varid, lat_varid
    INTEGER, DIMENSION(8) :: varid
!
    INTEGER :: nlon, nlat, jl, jb, icor
    INTEGER, PARAMETER :: kshift_n = $lat0
    INTEGER, PARAMETER :: kshift_w = $lon0
    REAL, DIMENSION(:,:), ALLOCATABLE :: fdat
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: fdir           ! RDF
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: filnew         ! Longitude index of Flow Destination according to FDIR
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: fibnew         ! Latitude index of Flow Destination according to FDIR
    REAL, DIMENSION(:,:), ALLOCATABLE :: alf_k
    REAL, DIMENSION(:,:), ALLOCATABLE :: alf_n
    REAL, DIMENSION(:,:), ALLOCATABLE :: arf_k
    REAL, DIMENSION(:,:), ALLOCATABLE :: arf_n
    REAL, DIMENSION(:,:), ALLOCATABLE :: agf_k

    CHARACTER (len=128) :: dn_para = "direction.nc"
!!    CHARACTER (len=128) :: dn_out = "$DNOUT"

!
!   ******** Open and read dimensions of HD parameter file
    ierr = nf90_open(dn_para, NF90_WRITE, ncid)

    ierr = nf90_inq_dimid(ncid,'lon',lon_dimid)
    ierr = nf90_inquire_dimension(ncid,lon_dimid,len=nlon)
    ierr = nf90_inq_dimid(ncid,'lat',lat_dimid)
    ierr = nf90_inquire_dimension(ncid,lat_dimid,len=nlat)
    WRITE(*,*) "Dimensions: NL =" , nlon, "  NB = ", nlat
!  
!     *** Feld Dimensionierung
    ALLOCATE(fdat(nlon,nlat))
    ALLOCATE(fdir(nlon,nlat))
    ALLOCATE(filnew(nlon,nlat))
    ALLOCATE(fibnew(nlon,nlat))
    ALLOCATE(alf_k(nlon,nlat))
    ALLOCATE(alf_n(nlon,nlat))
    ALLOCATE(arf_k(nlon,nlat))
    ALLOCATE(arf_n(nlon,nlat))
    ALLOCATE(agf_k(nlon,nlat))

    ierr = nf90_inq_varid(ncid,'FDIR',varid(1))
    ierr = nf90_get_var(ncid,varid(1),fdat)
    fdir(:,:) = NINT(fdat(:,:))
    ierr = nf90_inq_varid(ncid,'FILNEW',varid(2))
    ierr = nf90_get_var(ncid,varid(2),fdat)
    filnew(:,:) = NINT(fdat(:,:))
    ierr = nf90_inq_varid(ncid,'FIBNEW',varid(3))
    ierr = nf90_get_var(ncid,varid(3),fdat)
    fibnew(:,:) = NINT(fdat(:,:))
!
    ierr = nf90_inq_varid(ncid,'ALF_K',varid(4))
    ierr = nf90_get_var(ncid,varid(4),alf_k)
    ierr = nf90_inq_varid(ncid,'ALF_N',varid(5))
    ierr = nf90_get_var(ncid,varid(5),alf_n)
    ierr = nf90_inq_varid(ncid,'ARF_K',varid(6))
    ierr = nf90_get_var(ncid,varid(6),arf_k)
    ierr = nf90_inq_varid(ncid,'ARF_N',varid(7))
    ierr = nf90_get_var(ncid,varid(7),arf_n)
    ierr = nf90_inq_varid(ncid,'AGF_K',varid(8))
    ierr = nf90_get_var(ncid,varid(8),agf_k)

!
!   *** Correction of Index fields 
    icor = 0
    DO jb=1,nlat
    DO jl=1,nlon
      filnew(jl,jb) = filnew(jl,jb) - kshift_w
      fibnew(jl,jb) = fibnew(jl,jb) - kshift_n
      IF (filnew(jl,jb).EQ.0 .OR. filnew(jl,jb).EQ.nlon+1 .OR. &
          fibnew(jl,jb).EQ.0 .OR. fibnew(jl,jb).EQ.nlat+1) THEN
        icor = icor +1
        filnew(jl,jb) = jl
        fibnew(jl,jb) = jb
        fdir(jl,jb) = 5
        alf_k(jl,jb) = 0.
        alf_n(jl,jb) = 0.
        arf_k(jl,jb) = 0.
        arf_n(jl,jb) = 0.
        agf_k(jl,jb) = 0.
      ENDIF
    ENDDO
    ENDDO
    WRITE(*,*) icor, 'Border points became sinks'   ! Note that over sinks, K and N parameters must be zero

    fdat(:,:) = FLOAT(fdir(:,:))
    ierr = nf90_put_var(ncid,varid(1),fdat)
    fdat(:,:) = FLOAT(filnew(:,:))
    ierr = nf90_put_var(ncid,varid(2),fdat)
    fdat(:,:) = FLOAT(fibnew(:,:))
    ierr = nf90_put_var(ncid,varid(3),fdat)

    ierr = nf90_put_var(ncid,varid(4),alf_k)
    ierr = nf90_put_var(ncid,varid(5),alf_n)
    ierr = nf90_put_var(ncid,varid(6),arf_k)
    ierr = nf90_put_var(ncid,varid(7),arf_n)
    ierr = nf90_put_var(ncid,varid(8),agf_k)

    ierr = nf90_close(ncid)
    WRITE(*,*) TRIM(dn_para)," is closed."
  END PROGRAM corr_dir
EOF
#
#
# *** Compiling - On levante, necessary netcdf include and library statements are obtained via
# *** the calls listed below.
#
# module load netcdf-fortran/4.5.3-intel-oneapi-mpi-2021.5.0-intel-2021.5.0
# module load netcdf-c/4.8.1-intel-oneapi-mpi-2021.5.0-intel-2021.5.0
LPATH="/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/lib" 
NC_INCLUDE="-I/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/include"
NC_LIB1="`/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/bin/nf-config --flibs`"
NC_LIB="$NC_LIB1 -Wl,-rpath,$LPATH"
NC_INC1="`/sw/spack-levante/netcdf-c-4.8.1-7dq6g2/bin/nc-config --cflags`" 
NC_INC2="`/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/bin/nf-config --fflags`" 

echo $NC_LIB
echo $NC_INCLUDE
#
ifort corr_dir.f90 -o corr_dir.x ${NC_INCLUDE} ${NC_LIB} 
#
./corr_dir.x
ncks -h --no_abc -v FDIR,FILNEW,FIBNEW,ALF_K,ALF_N,ARF_K,ARF_N,AGF_K -A direction.nc $DNOUT
#
echo '**** Check whether file name is consistent with usage in HD model run script !!!!!!'
#
exit


