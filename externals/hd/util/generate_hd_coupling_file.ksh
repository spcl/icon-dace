#!/bin/ksh -l
#   
# generate_hd_coupling.ksh - Generate HD coupling file used in OASIS coupling and discharge conversion 
# 
# Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
#   
# *** Script/Program to generate a transformation matrix that associates each HD
#     mouth point (i.e. source grid) with a corresponding coastal mouth point on the target ocean grid. 
#     Here, the source coordinates are directly linked to the target mouth points, 
#     i.e. no additional remapping/interpolation required if this coupling file is used.
# 
#     Input: HD parameter file - Variable FDIR to idenfiy the HD mouth points
#            Ocean model land sea mask 
#
# Parse command line parameters
#
set -e
#
while [[ -n $1 ]] ; do
  case $1 in
    --id        | -i    ) EXP=$2                 ; shift    ;;    # Name tag or Experiment ID
    --source    | -s    ) dn_hd=$2               ; shift    ;;    # HD model parameter file (source grid)
    --omask     | -o    ) dn_ocean=$2            ; shift    ;;    # Ocean model land sea mask file 
    --bound     | -b    ) SBOUND=$2              ; shift    ;;    # Search radius outside of regional ocean domain 
    --mode      | -m    ) IMODE=$2               ; shift    ;;    # 0: Sea=0, 1: Sea=1 (Def.)
    --rad       | -r    ) SRAD=$2                ; shift    ;;    # Search radius in km for mapping
    --exc       | -x    ) IEXB=$2                ;;               # Exclude boundary 0=no 1=yes (Def.)
    -*                  ) err_strg="\nERROR: invalid option $1 !\n"
                          [[ "$1" != $( echo "$1" | tr -d = ) ]] && \
                          err_strg=$err_strg"  Please use blank instead of '=' to separate option's argument.\n"
                          usage ;;
    *                   ) export CPL_MOD=$1 ;;
  esac
  shift
done
#
# *** Set input files 
if [ "$dn_hd" = "" ] || [ "$dn_ocean" = "" ]; then
   echo "Either HD source file or Ocean model mask are not set. Please use the following options:"
   echo "  "
   echo "                   -i ,               ID = Experiment name tag or ID. Default if not set: hd_to_ocean "
   echo "                   -s <Source>,       = HD model parameter file"
   echo "                   -o <Ocean target>, Ocean model land sea mask file "
   echo "                   -m <Mode No.>,     Ocean mask type:  0: Sea=0, 1: Sea=1 (Def.)"
   echo "                   -r <Radius in km>, Search radius in km for mapping, e.g. 200 for NEMO-Europe (Def.)"
   echo "                   -x <0 or 1>,       Exclude boundary boxes from target mask: 0=no (Def.) 1=yes "
   exit
else
   echo "  Input files: "
   if [ -s $dn_hd ] ; then 
     echo "        HD parameter: " $dn_hd 
   else
     echo "        HD parameter file does not exist!: " $dn_hd 
     exit
   fi
   if [ -s $dn_ocean ] ; then 
     echo " Ocean land sea mask: " $dn_ocean
   else
     echo " Ocean land sea mask file does not exist: " $dn_ocean
     exit
   fi
fi
# *** Set Experiment ID 
if [ "$EXP" = "" ]; then
   echo "Experiment ID/Tag is not set. Please use Option -i <ID/Tag>"
   EXP=hd_to_ocean
   echo "  --> Default setting EXP=" $EXP " is used"
else
   echo "  --> Experiment ID is set to EXP = " $EXP 
fi
#
# *** Set ocean model mask type 
if [ "$IMODE" = "" ]; then
   echo "Ocean model mask type not set. Please use Option -m <IMODE>"
   echo "  --> Default setting IMODE = 1 is used"
   echo "IMODE = 0    Sea points = 0 (land = 1)"
   echo "IMODE = 1    Sea points = 1 (land = 0 or missing value)"
   IMODE=1
else
   case $IMODE in
     0  ) echo "Sea points in ocean mask = 0 " ;;
     1  ) echo "Sea points in ocean mask = 1 " ;;
     *  ) echo "ERROR - Option does not exist: " $IMODE ; exit
   esac
fi
#
# *** Set search radius in km for mapping
if [ "$SRAD" = "" ]; then
  SRAD=200.  # Search radius in km for mapping, e.g. 700 for ICON-global, 200 for NEMO-Europe
fi
echo "Search radius for mapping: $SRAD km"
#
# *** Set search radius in degree for outside ocean model domain
if [ "$SBOUND" = "" ]; then
  SBOUND=0.0833333  # Search radius in degree outside of regional ocean domain (0 for global grids)
fi
echo "Search radius outside ocean model domain: $SBOUND degree"
#
#
# *** Exclude boundary boxes from target mask
if [ "$IEXB" = "" ]; then
  LEXB=False
else
   case $IEXB in
     0  ) LEXB=False ; echo "Boundary boxes are included in target mask " ;;
     1  ) LEXB=True ; echo "Boundary boxes are excluded from target mask" ;;
     *  ) echo "ERROR - Option does not exist for IEXB: " $IEXB ; exit
   esac
fi
#
set -x
#
# *********** Data for Program *************************
# *** This needs to be EDITED  *************************
DSRC=/home/g/g260122/hdmodel/util      # Dir. with src files
##DIN=/work/gg0302/g260122/HD/input
##DOUT=/work/gg0302/g260122/HD/data
DOUT=/scratch/g/g260122/tmp
#
IMACH=2 # Machine: 1 Mistral, 2 Levante
#
#
# ******************************************************
#
# *** Compiling
case $IMACH in
  1 ) source /sw/rhel6-x64/etc/profile.mistral     # Necessary for module command
      [[ -n `whence ifort` ]] && module unload intel
      module load intel/18.0.4
      LPATH="/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/lib"
      LIBS="-lnetcdf"
      NC_INCLUDE="-I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/include"
      NC_LIB="`/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/bin/nf-config --flibs`"
      NC_INC1="`/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-gcc48/bin/nc-config --cflags`" 
      NC_INC2="`/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/bin/nf-config --fflags`"   ;;

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
# ******************************************************
#
#
# *** Write namelist
cat > hdtoocean_ctl.nml << EOF
&hdtoocean_ctl
  cexp = "$EXP"
  irad = 0                ! 1 = Convert coordinates from Radian to degree 
  icon = 0                ! 1 = ICON grid, i.e. only 1 coordinate dimension
  search_radius = $SRAD   ! Search radius in km
  search_boundary = $SBOUND        ! Search radius in deg outside of regional ocean domain (0 for global grids)
  lexbound = .${LEXB}.    ! Exclude boundary boxes from target mask?
  name_src = "HD-$EXP"    ! Name of Source grid
  name_target = "OM-$EXP" ! Name of Target grid
/
EOF
if [ -e $DOUT ] ; then echo 'Output directory exists: ' $DOUT ; else
  mkdir $DOUT
fi
#
set +e
rm -f convert.exe gen_mouth.exe
rm coast_oceangrid.nc
set -e
#
# ******************************************************
# *** Generate mask with coastal ocean points
#
# *** Run program that generates potential mouths as coastal ocean points from ocean model sea mask 
# *** that neighbor land points or that neighbor missing values at the NSEW side.
# ***  --> Output: coast_oceangrid.nc  
#          1 = Coastal Ocean Point, 0 = Other ocean points, 2 = land point

cdo setvar,mask -eqc,$IMODE ${dn_ocean} ocean_mask.nc
${F90} ${DSRC}/mo_grid.f90 ${DSRC}/generate_mouth.f90 -o gen_mouth.exe $NC_INCLUDE $NC_LIB
./gen_mouth.exe
if [ -s coast_oceangrid.nc ] ; then 
   echo "potential river mouths on ocean grid created from ocean mask !"
else
   echo "coast_oceangrid.nc does not exists = Failure in gen_mouth.exe --> TERMINATION of SCRIPT !"
   exit
fi 
#
# ******************************************************
#
# *** Allocate HD mouths

##${F90} ${DSRC}/mo_grid.f90 ${DSRC}/mo_time.f90 ${DSRC}/mo_flow_inout.f90 ${DSRC}/mo_interpol.f90 ${DSRC}/convert_discharge.f90 -o convert.exe $NC_INCLUDE $NC_LIB

##ulimit -s 102400      # Necessary on levante, otherwise segmentation fault may happen.
ulimit -s unlimited      # Necessary on levante, otherwise segmentation fault may happen.

${F90} ${DSRC}/mo_grid.f90  ${DSRC}/mo_interpol.f90 ${DSRC}/hdtoocean.f90 -o hdtoocean.exe $NC_INCLUDE $NC_LIB
cdo setctomiss,0. -eqc,0. -setmissval,-9999 -selvar,FDIR $dn_hd hdmouth_mask.nc
cdo setvar,FMOUTH hdmouth_mask.nc rivmouth_source.nc
# rm hdmouth_mask.nc 

#
# *** Run the Program to create data that will be used for the remapping of mouth points
# ***  --> Creates files hd_to_ocean_mouth.nc, hdmouth_on_oceangrid.nc
./hdtoocean.exe 
DNOUT=`ls hdcouple*${EXP}.nc`
DACT=`pwd`
#
# *** Add metainfo to global attributes
ncatted -O -h -a srad,global,o,c,"$SRAD" $DNOUT
ncatted -O -h -a sbound,global,o,c,"$SBOUND" $DNOUT

echo 'Coupling file: ' $DACT'/'$DNOUT ' generated'
#rm coast_oceangrid.nc
#rm hdtoocean.exe

exit
#

