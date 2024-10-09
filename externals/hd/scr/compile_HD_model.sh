#!/bin/sh
#
# compile_HD_model.sh - Script to compile the HD model on DKRZ HLRE
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann and Ha Ho-Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
# *** Script to compile HD model on Mistral or Levante
#
# Parse command line parameters
#
while [[ -n $1 ]] ; do
  case $1 in
    --compiler | -c     ) compiler=$2                ; shift    ;;
    --help     | -h     ) help=1                     ; break    ;;
    --oasis    | -o     ) oasis_coupling=$2          ; shift    ;;
    --yac      | -y     ) yac_coupling=$2            ; shift    ;;
    --resolution  | -r  ) HDRES=$2                             ;;
    -*                  ) err_strg="\nERROR: invalid option $1 !\n"
                          [[ "$1" != $( echo "$1" | tr -d = ) ]] && \
                          err_strg=$err_strg"  Please use blank instead of '=' to separate option's argument.\n"
                          usage ;;
    *                   ) export cline=$1 ;;
  esac
  shift
done
#
# *** Set Resolution/Domain 
if [ "$HDRES" = "" ] || [ "$help" = "1" ]; then
   if [ "$help" != "1" ] ; then echo "Resolution/Domain is not set --> This is necessary!" ; fi
   echo " "
   echo "Help-Listing: Run with Option -h or --help or without any option"
   echo " "
   echo "Options for resolution/domain settings with -r/--resolution. Please use Option -r <No.>"
   echo " No.: 0  = 05deg     --> Global 0.5 degree domain"
   echo "      1  = 5min      --> Global or regional 5 Min. domain"
   echo " "
   echo "Option for compiler settings with -c/--compiler: Please use Option -c <No.> or <compiler name>"
   echo " No.: 1  = intel           --> Setting for Mistral (Bull/Atos) at DKRZ"
   echo "      2  = bull            --> Old setting for Mistral at DKRZ"
   echo "      3  = intel-strand    --> Setting for Strand at HZG"
   echo "      4  = openmpi         --> Setting for Mistral at DKRZ"
   echo "      5  = bullxmpi        --> New setting for Mistral at DKRZ"
   echo "      6  = intel-levante   --> Setting for Levante (new Bull/Atos) at DKRZ"
   echo "      7  = openmpi-levante --> Setting for Levante at DKRZ (Default)"

   echo " "
   echo "Options for using OASIS coupling with -o/--oasis. Please use Option -o <No.> or ON/OFF"
   echo " No.: 0  = OFF --> No coupling to OASIS  (default)"
   echo "      1  = ON  --> Coupling to OASIS "

   echo " "
   echo "Options for using YAC coupling with -y/--yac. Please use Option -y <No.> or ON/OFF"
   echo " No.: 0  = OFF --> No coupling to YAC  (default)"
   echo "      1  = ON  --> Coupling to YAC "
   exit
fi
#
# *** Compiler set?
if [ "$compiler" = "" ]; then
   echo "Compiler not set. Please use Option -c <No.> or <compiler name>"
   echo " No.: 1  = intel --> Default (if omitted) - Setting for Mistral at DKRZ"
   echo "      2  = bull = BULL = BULLMPI --> Old Setting for Mistral at DKRZ"
   echo "      3  = intel-strand    --> Setting for Strand at HZG"
   echo "      4  = openmpi         --> Setting for Mistral at DKRZ"
   echo "      5  = bullxmpi        --> New Setting for Mistral at DKRZ"
   echo "      6  = intel-levante   --> Setting for Levante"
   echo "      7  = openmpi-levante --> Setting for Levante at DKRZ"
   ICOMP=7           # Default = openmpi-levante
elif [ "$compiler" = "intel" ] || [ "$compiler" = "1" ] ; then
  ICOMP=1
elif [ "$compiler" = "bull" ] || [ "$compiler" = "bullmpi" ] ||  \
     [ "$compiler" = "BULL" ] || [ "$compiler" = "BULLMPI" ] ||  \
     [ "$compiler" = "2" ]  ; then
  ICOMP=2
elif [ "$compiler" = "intel-strand" ] || [ "$compiler" = "3" ] ; then
  ICOMP=3
elif [ "$compiler" = "openmpi" ] || [ "$compiler" = "4" ]  ; then
  ICOMP=4
elif [ "$compiler" = "bullxmpi" ] || [ "$compiler" = "5" ]  ; then
  ICOMP=5
elif [ "$compiler" = "intel-levante" ] || [ "$compiler" = "6" ] ; then
  ICOMP=6
elif [ "$compiler" = "openmpi-levante" ] || [ "$compiler" = "7" ]  ; then
  ICOMP=7
fi 
#
# *** Coupling with OASIS
if [ "$oasis_coupling" = "" ]; then
   echo "Oasis coupling not set. Please use Option -o <No.>"
   echo " No.: 0  = No coupling to OASIS  (default)"
   echo "      1  = Coupling to OASIS "
   export COUP_OAS='no' 	# must be set to yes for COUPLED runs in 3 models
elif [ "$oasis_coupling" = "ON" ] || [ "$oasis_coupling" = "1" ]; then
  export COUP_OAS='yes' 	# Setting for COUPLED runs in 3 models
else
  export COUP_OAS='no' 	
fi

# *** Coupling with YAC
if [ "$yac_coupling" = "" ]; then
   echo "YAC coupling not set. Please use Option -y <No.>"
   echo " No.: 0  = No coupling to YAC  (default)"
   echo "      1  = Coupling to YAC "
   export COUP_YAC='no' 	# must be set to yes for COUPLED runs in 3 models
elif [ "$yac_coupling" = "ON" ] || [ "$yac_coupling" = "1" ]; then
  export COUP_YAC='yes' 	# Setting for COUPLED runs in 3 models
else
  export COUP_YAC='no' 	
fi
 
export NOMPI='no' 	# MPI
#export NOMPI='yes' 	# NOMPI
echo 'ICOMP = ' $ICOMP  ' and HDRES=' $HDRES 'and COUP_OAS=' $COUP_OAS 'and COUP_YAC=' $COUP_YAC

#
case ${HDRES} in
  0 | 05deg    ) export HD_5MIN='no'       # 0.5 degree HD model resolution
                 EXE="hd_05.exe" ;;
  [123] | 5min ) export HD_5MIN='yes'      # needs yes for 5 Min HD model resolution
                 EXE="hd_5min.exe" ;;
  * )  echo 'ERROR: Definition does not exist for HDRES=' $HDRES ; exit ;;
esac

#
CODEDIR=../code

cd $CODEDIR

if [ -s obj ] ; then echo "Object directory exists" 
else
  mkdir obj
fi
#
set +u
. ${MODULESHOME}/init/ksh
set -u

#
# ****** set modules for compiler
case ${ICOMP} in
    1  )  [[ -n `whence ifort` ]] && module unload intel       # Intel auf Mistral
          module load intel/18.0.4
          [[ -n `whence mpiifort` ]] && module unload intelmpi
          module load intelmpi/5.1.3.223
          cp -p ../code/Fopts_Intel ../code/Fopts   ;;
    2  )  [[ -n `whence ifort` ]] && module unload intel
          module load intel/18.0.4
          [[ -n `whence mpif90` ]] && module unload bullxmpi_mlx
          module load mxm/3.4.3082    # mxm/3.3.3002
          module load bullxmpi_mlx/bullxmpi_mlx-1.2.9.2   # bullxmpi_mlx/bullxmpi_mlx-1.2.8.3
          cp -p ../code/Fopts_BULLMPI ../code/Fopts   ;;
    3  )  [[ -n `whence ifort` ]] && module unload intel       # Intel auf Strand
          module load compilers/intel/2018.1.163
          [[ -n `whence mpiifort` ]] && module unload intelmpi
          module load intelmpi/2018.1.163
          cp -p ../code/Fopts_Intel_Strand ../code/Fopts   ;;
    4  )  [[ -n `whence ifort` ]] && module unload intel       # openmpi on Mistral
          module load intel/17.0.6
          [[ -n `whence mpif90` ]] && module unload openmpi
          module load openmpi/2.0.2p2_hpcx-intel14
          module load gcc/6.4.0
          cp -p ../code/Fopts_openmpi_mistral ../code/Fopts   ;;
    5  )  [[ -n `whence ifort` ]] && module unload intel
          module load intel/17.0.6
          [[ -n `whence mpif90` ]] && module unload bullxmpi_mlx
          module load mxm/3.4.3082    # mxm/3.3.3002
          module load gcc/6.4.0
          module load bullxmpi_mlx/bullxmpi_mlx-1.2.9.2   # bullxmpi_mlx/bullxmpi_mlx-1.2.8.3
          cp -p ../code/Fopts_bullxmpi ../code/Fopts   ;;
    6  )  [[ -n `whence mpiifort` ]] && module unload intel-oneapi-mpi
          module load intel-oneapi-mpi/2021.5.0-intel-2021.5.0
          [[ -n `whence gcc` ]] && module unload gcc
          module load intel-oneapi-compilers/2022.0.1-gcc-11.2.0
          cp -p ../code/Fopts_oneapi-mpi ../code/Fopts  ;;
    7  )  [[ -n `whence mpifort` ]] && module unload openmpi
          module load openmpi/4.1.2-intel-2021.5.0
          cp -p ../code/Fopts_openmpi_levante ../code/Fopts  ;;
    *  )  echo 'Compiler ' $ICOMP ' does not exist' ; exit ;;
esac
#
make clean
make
#
# Rename output file
mv hd.exe ../${EXE}
#


