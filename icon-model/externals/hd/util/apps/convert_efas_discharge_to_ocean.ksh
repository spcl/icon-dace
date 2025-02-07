#!/bin/ksh -l
#
# convert_efas_discharge_to_ocean.ksh - Converts discharge on the EFAS grid to a selected NEMO ocean grid
# 
# Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
# *** Script/program converts discharge (EXP=7056110) on the EFAS grid to a selected ocean grid. 
#     The EFAS mouths are directly converted to mouth points on the ocean grid.
#     Program requires a HD coupling file that, e.g. has been generated with generate_hd_couplimg_file.ksh 
#     Note that the discharge variable must be named friv or var0 for iclim=0.
#
#     cdo setname,friv -daymean -shifttime,-3h riv_dis_efas_hist_2018.nc /work/gg0302/g260122/HD/output/7056110/7056110_meanflow_2018.nc
set -e
#
while [[ -n $1 ]] ; do
  case $1 in
    --id        | -i    ) EXP=$2                 ; shift    ;;
    --nemo      | -n    ) INEMO=$2               ; shift    ;;
    --ybeg      | -y    ) YBEG=$2                ; shift    ;;
    --yend      | -z    ) YEND=$2                ;;
    -*                  ) err_strg="\nERROR: invalid option $1 !\n"
                          [[ "$1" != $( echo "$1" | tr -d = ) ]] && \
                          err_strg=$err_strg"  Please use blank insted of '=' to separate option's argument.\n"
                          usage ;;
    *                   ) export CPL_MOD=$1 ;;
  esac
  shift
done
#
# *** Set Experiment ID / List 
if [ "$EXP" = "" ] ; then
   echo "No experiement number is set. Please use the following options"
   echo "     -i <7-digit ID>,   ID = Source exp. no., Default: Last Exp. no. "
   echo "     -n  <NEMO grid>    1 = GCOAST (Def.), 2=North Sea, 3=Nordic (Bal-MFC)"
   echo "     -y  <First Year>   e.g. YYYY, Default: 2018, if 0 use preset values"
   echo "     -z  <Last Year>    e.g. YYYY, Default: <First Year>"
   exit
fi
#
# *** Set NEMO grid: 1 = GCOAST, 2=North Sea, 3=Nordic (Bal-MFC)
if [ "$INEMO" = "" ] ; then echo "NEMO grid not set --> use default: GCOAST" ; INEMO=1 ; fi
#
# *** Set Begin and End year 
if [ "$YBEG" = "" ] ; then echo "First year is not set --> use default: 2018" ; YBEG=2018 ; fi
if [ "$YEND" = "" ] ; then echo "Last year is not set --> use first year" ; YEND=$YBEG ; fi
#
# *********** Data for Program *************************
# *** This needs to be EDITED  *************************

DSRC=/home/g/g260122/hdmodel/util      # Dir. with source code files
DIN=/work/gg0302/g260122/HD/input
DRUN=/scratch/g/g260122/tmp
TAG_IN='efas'        # Tag for discharge data
DATA=/work/gg0302/g260122/HD/output/$EXP   # Discharge data directory
DATA_OUT=$DRUN                             # Data output directory
#
# *** Conversion of discharge or bgc flow in one go? - Standard values are set below 
# *** Do not change if you do not know what they mean! 
ICONV=1  # 0=NO conversion, 1=daily data, 2=monthly climatology, 3=sequence of 5 bgc fields
         # 4 = annual sequence of several bgc variables and time steps, e.g. monthly
ISEPIN=1      # Inflow in separate files per year (No/yes = 0/1)
ISEPOUT=1     # Flow on ocean in separate files per year (No/yes = 0/1)
ICAL=0        # Calendar: 0=normal, 1=360 day years
IMACH=2       # Machine: 1 Mistral, 2 Levante
#
#
# NEMO grid
case $INEMO in 
   1 ) TAG_OUT='nemo-gcoast'      # Tag for ocean model grid
       DN_HDCOUPLE='/work/gg0302/g260122/HD/input/efas/hdcouple_efas_to_nemo-gcoast.nc'  ;;
   2 ) TAG_OUT='nemo-northsea'    # Tag for ocean model grid
       DN_HDCOUPLE='/work/gg0302/g260122/HD/input/efas/hdcouple_efas_to_nemo-ns.nc'  ;;
   3 ) TAG_OUT='nemo-nordic'      # Tag for ocean model grid
       DN_HDCOUPLE='/work/gg0302/g260122/HD/input/efas/hdcouple_efas_to_nemo-nordic.nc'  ;;
   * ) echo "Grid no. $INEMO does not exist --> TERMINATE" ; exit  ;;
esac
#
# *****************************************************************************
#
DNRIV=${DATA}/${EXP}_meanflow_YYYY.nc    # Def. HD name ${EXP}_meanflow_YYYY.nc for ISEPIN=1
DNOUTFLOW=${DATA_OUT}/${TAG_IN}_on_${TAG_OUT}_YYYY.nc
#
cd $DRUN
#
# *** Write namelist
cat > convert_inflow_ctl.nml << EOF
&convert_inflow_ctl
  dn_inflow = "$DNRIV"
  dn_outflow = "$DNOUTFLOW"
  iconv = $ICONV              ! Method of conversion
  isep_in = $ISEPIN           ! Input data in separate files (No/Yes=0/1)
  isep_out = $ISEPOUT         ! Output data in separate files (No/Yes=0/1)
  ybeg = $YBEG                ! Start year of data
  yend = $YEND                ! End year of data
/
EOF


#
# ******* Compiling **********
case $IMACH in
  1 ) echo "Mistral does not exist anymore --> STOP" ;exit  ;;

# Levante: requires module load netcdf-fortran/4.5.3-intel-oneapi-mpi-2021.5.0-intel-2021.5.0   for netcdf   
#          module load netcdf-c/4.8.1-intel-oneapi-mpi-2021.5.0-intel-2021.5.0                  for netcdf-c
  2 ) module load netcdf-fortran/4.5.3-intel-oneapi-mpi-2021.5.0-intel-2021.5.0
      module load netcdf-c/4.8.1-intel-oneapi-mpi-2021.5.0-intel-2021.5.0

      LPATH="/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/lib" 
      LIBS="-lnetcdf"
        NC_INCLUDE="-I/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/include"
      NC_LIB1="`/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/bin/nf-config --flibs`"
      NC_LIB="$NC_LIB1 -Wl,-rpath,$LPATH"
      NC_INC1="`/sw/spack-levante/netcdf-c-4.8.1-7dq6g2/bin/nc-config --cflags`" 
      NC_INC2="`/sw/spack-levante/netcdf-fortran-4.5.3-r5r3ev/bin/nf-config --fflags`"  ;;
esac 
NC_INCLUDE="$NC_INC2 $NC_INC1"
#
F90=ifort
#
# *** Allocate HD mouths
##${F90} ${DSRC}/mo_grid.f90 ${DSRC}/mo_time.f90 ${DSRC}/mo_flow_inout.f90 ${DSRC}/mo_interpol.f90 ${DSRC}/mo_convert.f90 ${DSRC}/convert_discharge.f90 -o convert.exe $NC_INCLUDE $NC_LIB

#
# ******************************************************
#

cat > convert_hdtoocean.f90 << EOF
PROGRAM convert_hdtoocean
!
!   ******** Programm that takes reads the HD river discharge at HD river mouth points 
!            and transfers it to the ocean model grid and the associated coastal mouth points.
!            It reads the HD coupling file 
!
!   *** HD coupling file variables
!   *** HD grid 
!   ***   FMOU_HD_NEMO = Mask with HD mouth boxes that have an associated ocean model mouth point.
!   ***   INDEXX = x-Indices of nearest ocean model mouth points
!   ***   INDEXY = y-Indices of nearest ocean model mouth points
!   *** Ocean model grid 
!   *** FMOU_HD_ON_NEMO = Mask on which with HD mouth boxes are mapped
!
!   ******** Version 1.0 - January 2023
!            Programmierung und Entwicklung: Stefan Hagemann 
!            Based on the conversion part of convert_discharge.f90                  
!
    use netcdf
    use mo_grid,         ONLY: model
    use mo_convert,      ONLY: convert_inflow, open_coupling_file, read_coupling_file

    IMPLICIT NONE
!
    DOUBLE PRECISION, PARAMETER :: xmiss=-9999. 
    DOUBLE PRECISION, PARAMETER :: zeps = 1.E-6 
    CHARACTER (LEN=20)  :: ctitle ="${TAG_IN}_on_${TAG_OUT}"   ! Description of converted data
!
!   *** Define source and target grid
!
!   *** Grid info
    CHARACTER (len=192) :: dn_couple = "${DN_HDCOUPLE}"   ! HD coupling file 
    TYPE(model)        :: mask_src      ! Source grid
    TYPE(model)        :: mask_target   ! Target grid
!
!   *** Arrays on source grid, e.g. HD
!
!   Mask of source mouths with a valid nearest target mouth
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask_src_mapped   ! Stores FMOU_HD_NEMO
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ix_target         ! Stores INDEXX
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: iy_target         ! Stores INDEXY
!   Ocean model (e.g. NEMO) arrays
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask_on_ocean     ! Stores FMOU_HD_ON_NEMO
!
!   *** Parameters to steer the conversion of discharge/bgc inflow: convert_inflow_ctl
    CHARACTER (len=192) :: dn_inflow 
    CHARACTER (len=192) :: dn_outflow       
    INTEGER :: iconv = 1    ! Method of conversion: 0 = no, 1=Daily data (def.),
                            ! 2=monthly clim., 3=5 bgc sequence
                            ! 4=General sequence of bgc variables and timesteps
    INTEGER :: isep_in =1   ! Input data in separate files (No/Yes=0/1)
    INTEGER :: isep_out = 1 ! Output data in separate files (No/Yes=0/1)
    INTEGER :: ybeg         ! Start year of data
    INTEGER :: yend         ! End year of data

    INTEGER :: lu_nml = 20
    LOGICAL :: logque
    NAMELIST /convert_inflow_ctl/ dn_inflow, dn_outflow, iconv, isep_in, isep_out, ybeg, yend
!
!   *** Read Namelist convert_inflow_ctl.nml
    INQUIRE(FILE='convert_inflow_ctl.nml', EXIST=logque)
    IF (logque) THEN
      OPEN (UNIT=lu_nml, FILE='convert_inflow_ctl.nml', STATUS='OLD')
      READ (lu_nml, NML=convert_inflow_ctl)
      CLOSE (lu_nml)
    ENDIF
!
!   *** Read Coupling file and define arrays
    CALL open_coupling_file(dn_couple, mask_src, mask_target)

    ! Allocation of coupling fields
    ALLOCATE(mask_src_mapped(mask_src%nlon,mask_src%nlat))
    ALLOCATE(ix_target(mask_src%nlon,mask_src%nlat))
    ALLOCATE(iy_target(mask_src%nlon,mask_src%nlat))

    WRITE(*,*) "--> allocating for: NX =" , mask_target%nlon, "  NY = ", mask_target%nlat
    ALLOCATE(mask_on_ocean(mask_target%nlon, mask_target%nlat))

    CALL read_coupling_file(dn_couple, mask_src_mapped, ix_target, iy_target, mask_on_ocean)
!
!   *** Conduct conversion
    CALL convert_inflow(dn_inflow, dn_outflow, ctitle, iconv, ybeg, yend, isep_in, isep_out, &
         mask_src, mask_target, xmiss, &
         mask_src_mapped, ix_target, iy_target, mask_on_ocean)      
!
!   ******** Programmende

END PROGRAM convert_hdtoocean
!	  
EOF
#
set +e
rm convert_hdtooc.exe
set -e
ulimit -s unlimited

${F90} ${DSRC}/mo_grid.f90 ${DSRC}/mo_time.f90 ${DSRC}/mo_flow_inout.f90 ${DSRC}/mo_interpol.f90 ${DSRC}/mo_convert.f90 convert_hdtoocean.f90 -o convert_hdtooc.exe $NC_INCLUDE $NC_LIB
#
# Run the Program 
if [ -e discharge.nc ];then rm discharge.nc ; fi

./convert_hdtooc.exe
#

