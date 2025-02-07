# compile_HD_model.sh - Script to compile the HD model in a coupled system
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Author: Ha Ho-Hagemann (Hereon, Germany)
# Contact: ha.hagemann@hereon.de
#_________________________________________
#
#
# First, please adapt WORK_DIR & SRC_DIR if neccessary
#
WORK_DIR=$(pwd)
SRC_DIR=/work/gg0302/g260062/GCOAST_oas5

#=================== CHOOSE HERE WHAT TO COMPILE !!! ==============
#computer="mistral"
computer="levante"

#......................
coupling="no"
coupling="yes"

#......................
coupler="OASIS"
coupler="YAC"

#......................
#Initial values, don't change two following lines!!!
compile_OASIS="no"
compile_YAC="no"

#......................
if [ $coupling == "yes" ];then
 if [ $coupler == "OASIS" ];then
### compile_OASIS="yes"
  compile_OASIS="no"		# if already compiled OASIS
 else
###  compile_YAC="yes"
  compile_YAC="no"		# if already compiled YAC
 fi
fi

#......................
compile_HD="yes"
#compile_HD="no"

if [ $coupling == "no" ];then
 version_HD="5.1"
else
 if [ $coupler == "OASIS" ];then
  version_HD="5.1_oas4_new"		# OASIS interface
 else
  version_HD="5.1_yac"		# YAC interface
 fi
fi

md_domain="global0.5deg"
#md_domain="global5min"
#md_domain="euro5min"

#===============================================================
if [ $computer == "mistral" ];then
#mistral
 NETCDFF_DIR=/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14
 NETCDFC_DIR=/sw/rhel6-x64/netcdf/netcdf_c-4.4.0-parallel-openmpi2-intel14
 HDF5_DIR=/sw/rhel6-x64/hdf5/hdf5-1.8.18-parallel-openmpi2-intel14
 GRIBAPI_DIR=/sw/rhel6-x64/grib_api/grib_api-1.15.0-gcc48
 MODULES="'gcc/6.4.0 intel/17.0.6 openmpi/2.0.2p2_hpcx-intel14'"
 GCCLIB="/sw/rhel6-x64/gcc/gcc-6.4.0/lib64"
 ECCODES_DIR="/sw/rhel6-x64/eccodes/eccodes-2.22.0-gcc64"
 PYTHON='/sw/spack-rhel6/miniforge3-4.9.2-3-Linux-x86_64-pwdbqi/bin/python'
else
# levante
 NETCDFF_DIR=/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g
 NETCDFC_DIR=/sw/spack-levante/netcdf-c-4.8.1-2k3cmu
 HDF5_DIR=/sw/spack-levante/hdf5-1.12.1-tvymb5
 SZIP_ROOT=/sw/spack-levante/libaec-1.0.5-gij7yv
 MPIINC=/sw/spack-levante/openmpi-4.1.2-yfwe6t/include
 MPILIB=/sw/spack-levante/openmpi-4.1.2-yfwe6t/lib
 MODULES="'gcc-11.2.0-7jcqrc intel-oneapi-compilers/2022.0.1-gcc-11.2.0 openmpi/4.1.2-intel-2021.5.0'"
 GCCLIB="/sw/spack-levante/gcc-11.2.0-7jcqrc/lib64"
 PYTHON='/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin/python3'
fi
#===============================================================
OASIS_DIR=${SRC_DIR}/oasis3-mct

#..........................................
YAXT_DIR=${SRC_DIR}/yaxt
YAC_DIR=${SRC_DIR}/yac

#..........................................
HD_DIR=${WORK_DIR}

#===============================================================
#echo "WORK_DIR=${WORK_DIR}" > ${WORK_DIR}/paths.inc
#echo "OASIS_DIR=${OASIS_DIR}" >> ${WORK_DIR}/paths.inc
#echo "HD_DIR=${HD_DIR}" >> ${WORK_DIR}/paths.inc

#===============================================================
# Step 1.1: Compile oasis3-mct
if [ ${compile_OASIS} == "yes" ];then
 cd ${OASIS_DIR}/util/make_dir/
 echo "include ${OASIS_DIR}/util/make_dir/make.${computer}.dkrz" > ${OASIS_DIR}/util/make_dir/make.inc

  sed \
      -e s%@{NETCDFF_DIR}%${NETCDFF_DIR}%g \
      -e s%@{NETCDFC_DIR}%${NETCDFC_DIR}%g \
      -e s%@{OASIS_DIR}%${OASIS_DIR}%g \
      <${WORK_DIR}/config/make.oasis.${computer}.dkrz>${OASIS_DIR}/util/make_dir/make.${computer}.dkrz

 ./run_make.sh
 cd ${WORK_DIR}
fi	# compile_OASIS

#===============================================================
# Step 1.2: Compile yaxt & yac
if [ ${compile_YAC} == "yes" ];then
 cd ${YAXT_DIR}
 ./run_compile_yaxt.sh

 cd ${YAC_DIR}
 ./run_compile_yac.sh

 cd ${WORK_DIR}
fi	# compile_YAC

#===============================================================
# Step 2: Compile HD
if [ ${compile_HD} == "yes" ];then
 if [ $computer == "mistral" ];then
  sed \
      -e s%@{NETCDFF_DIR}%${NETCDFF_DIR}%g \
      -e s%@{NETCDFC_DIR}%${NETCDFC_DIR}%g \
      -e s%@{HDF5_DIR}%${HDF5_DIR}%g \
      -e s%@{OASIS_DIR}%${OASIS_DIR}%g \
      <${WORK_DIR}/config/Fopts_HD_openmpi_${computer}>${HD_DIR}/../code/Fopts
 else
  sed \
      -e s%@{NETCDFF_DIR}%${NETCDFF_DIR}%g \
      -e s%@{NETCDFC_DIR}%${NETCDFC_DIR}%g \
      -e s%@{HDF5_DIR}%${HDF5_DIR}%g \
      -e s%@{SZIP_ROOT}%${SZIP_ROOT}%g \
      -e s%@{OASIS_DIR}%${OASIS_DIR}%g \
      -e s%@{MPIINC}%${MPIINC}%g \
      -e s%@{MPILIB}%${MPILIB}%g \
      -e s%@{YAXT_DIR}%${YAXT_DIR}%g \
      -e s%@{YAC_DIR}%${YAC_DIR}%g \
      <${WORK_DIR}/config/Fopts_HD_openmpi_${computer}>${HD_DIR}/../code/Fopts
 fi
 
 BUILD_HD=BUILD/build_HD${version_HD}_${computer}
 if [ ! -d ${HD_DIR}/${BUILD_HD} ];then
  mkdir -p ${HD_DIR}/${BUILD_HD}
 fi

 cd ../code	# HD main code directory

 if [ -s obj ] ; then echo "Object directory exists" 
 else
  mkdir obj
 fi

 export NOMPI='no' 	# MPI
# export NOMPI='yes' 	# NOMPI

 case ${md_domain} in
  global0.5deg ) export HD_5MIN='no'       # 0.5 degree HD model resolution
                 EXE=${HD_DIR}/${BUILD_HD}/"hd_05.exe" ;;
  global5min   ) export HD_5MIN='yes'      # needs yes for 5 Min HD model resolution
                 EXE=${HD_DIR}/${BUILD_HD}/"hd_5min.exe" ;;
  euro5min     ) export HD_5MIN='yes'      # needs yes for 5 Min HD model resolution
                 EXE=${HD_DIR}/${BUILD_HD}/"hd_euro.exe" ;;
  *            )  echo 'ERROR: Definition does not exist for md_domain=' $md_domain ; exit ;;
 esac

 if [ $coupling == "no" ];then
 #----------- STAND-ALONE ---------------------------
  export COUP_OAS='no' 	
  export COUP_YAC='no'
 else
 #----------- COUPLED ---------------------------
  if [ $coupler == "OASIS" ];then
   export COUP_OAS='yes'
   export COUP_YAC='no'
  elif [ $coupler == "YAC" ];then
   export COUP_OAS='no' 	
   export COUP_YAC='yes'
  else
   echo "ONLY support OASIS and YAC coupler"
   exit
  fi
 fi # coupling

 make clean
 make

 echo ""
 echo "Compiling is successful"

 mv hd.exe ${EXE}
 ls -l ${EXE}

 cd ${WORK_DIR}
fi

#================================================================
cd ${WORK_DIR}

