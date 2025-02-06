#! /bin/bash

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2024, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

main_name="preprocess_LUH2.bash"
# Cleaning up
case ${1} in

    "clean")
        echo -n "Removing binaries and machine config..."
        rm -rf ${BIN_DIR}
        rm config
        rm -f user_config
        rm -f *.o*
        rm ${main_name}
        echo " done"
        exit
        ;;
    *)
        echo "To clean executables call ./compile clean."
        ;;

esac

echo "On machine " $HOSTNAME

case $HOSTNAME in

  *"levante"*)
      echo "Settings for levante."
      module purge
            
      echo module load cdo > config
      echo module load nco >> config
      echo "export ncap=ncap2" >> config

      echo spack load intel-oneapi-compilers >> config
      echo spack load netcdf-c %intel@2021.5.0 /f7hh57 >> config
      echo spack load netcdf-fortran %intel@2021.5.0 /pvmcx6 >> config
            
      spack load intel-oneapi-compilers
      spack load netcdf-c %intel@2021.5.0 /f7hh57
      spack load netcdf-fortran %intel@2021.5.0 /pvmcx6
      
      lib_path=/sw/spack-levante
      netcdf_path=$lib_path/netcdf-fortran-4.5.3-pvmcx6
      netcdfc_path=$lib_path/netcdf-c-4.8.1-f7hh57
      hdf5_path=$lib_path/hdf5-1.12.1-4l5deg
      aec_path=$lib_path/libaec-1.0.5-gij7yv
      
      CC=ifort
      FFLAGS="-Wl,--unresolved-symbols=ignore-in-object-files"
      ;;
  *"mlogin"*)
      echo "Settings for mistral"
      
      echo module unload nco > config
      echo module load nco/4.6.7-gcc48 >> config
      echo module unload cdo >> config
      echo module load cdo/1.7.0-magicsxx-gcc48 >> config
      echo "export ncap=ncap" >> config
    
      echo module load nag >> config
      echo module load netcdf_c/4.3.2-gcc48 >> config
      echo module load netcdf-fortran/4.5.3-nag-7.0 >> config
      
      module load nag
      module load netcdf_c/4.3.2-gcc48
      module load netcdf-fortran/4.5.3-nag-7.0
      
      lib_path=/sw/rhel6-x64
      netcdf_path=$lib_path/netcdf/netcdf_fortran-4.4.2-static-nag60
      netcdfc_path=$lib_path/netcdf/netcdf_c-4.3.2-static-gcc48
      hdf5_path=$lib_path/hdf5/hdf5-1.8.14-static-threadsafe-gcc48
      aec_path=$lib_path/sys/libaec-0.3.2-static-gcc48

      CC=nagfor
      FFLAGS="-colour -O3 -w=uep"
      ;;
  *"jessy"*)
    echo "Settings for jessy"

    echo module unload nco > config
    echo module load nco/4.6.7-gcc48 >> config
    echo module unload cdo >> config
    echo module load cdo/1.7.0-magicsxx-gcc48 >> config
    echo "export ncap=ncap" >> config

    module load nag

    lib_path=/sw/jessie-x64/
    netcdf_path=$lib_path/netcdf_fortran-4.4.2-static-nag60
    netcdfc_path=$lib_path/netcdf-4.3.3.1-static-gccsys
    hdf5_path=$lib_path/hdf5-1.8.16-static-gccsys
    szip_path=$lib_path/szip-2.1-static-gccsys

    CC=nagfor
    FFLAGS="-colour -nan -gline -O0 -C=all -w=uep"
    ;;
  *)
    echo "Machine not recognized. Please add your environment. STOP"
    exit
    ;;
esac

# Set library paths
case $HOSTNAME in

    *"jessy"*)
    NETCDF="-I$netcdf_path/include -lnetcdff -L$netcdf_path/lib -lnetcdf -L$netcdfc_path/lib"
    HDF5LIB="-L$hdf5_path/lib -lhdf5_hl -lhdf5 -ldl"
    SZIP="-L$szip_path//lib"
    GCC="-lsz -lz"

    LDFLAGS="$NETCDF $HDF5LIB $SZIP $GCC"
    ;;
    *)
    NETCDF="-I$netcdf_path/include -lnetcdff -L$netcdf_path/lib -lnetcdf -L$netcdfc_path/lib"
    HDF5LIB="-L$hdf5_path/lib"
    AECLIB="-L$aec_path/lib"
    LIBA="$netcdfc_path/lib/libnetcdf.a $hdf5_path/lib/libhdf5_hl.a $hdf5_path/lib/libhdf5.a -lrt $aec_path/lib*/libsz.a $aec_path/lib*/libaec.a"
    GCC="-lz -lm -ldl -lcurl -lpthread"

    LDFLAGS="$NETCDF $HDF5LIB $AECLIB $LIBA $GCC"
    ;;
esac

# Work directories
SRC_DIR=./src
BIN_DIR=./bin
TOOLS_DIR=./tools
# Source files
SOURCE=( $SRC_DIR/* )

if [ ! -d ${BIN_DIR} ]; then
    mkdir ${BIN_DIR}
fi
# Compile fortran programs
for src in ${SOURCE[@]}; do
    TARGET=`basename ${src} .f90`.x

    if [[ ! -e ${BIN_DIR}/${TARGET} ]]; then
        echo "Compiling ${src} -> ${TARGET}"
        ${CC} ${FFLAGS} ${src} -o ${BIN_DIR}/${TARGET} ${LDFLAGS}
    elif [[ ${src} -nt ${BIN_DIR}/${TARGET} ]]; then
        echo "${src} has been modified -> compiling"
        ${CC} ${FFLAGS} ${src} -o ${BIN_DIR}/${TARGET} ${LDFLAGS}
    else
        echo "Nothing to do."
        echo ${TARGET} "up to date." 
    fi
done

# Create main run script
# Get the options
while getopts ":he:p:c:" option; do
    case $option in
        h) # display Help
            Help
            exit;;
        e) # Enter a email address
            user_email=$OPTARG;;
        p) # Enter a number
            project_number=$OPTARG;;
        c) # Enter config file
            user_config=$OPTARG;;
        \?) # Invalid option
            echo "Error: Invalid option"
            exit;;
    esac
done

if [[ ! -e ${main_name} || ${TOOLS_DIR}/MAIN_${main_name} -nt ${main_name} ]]; then
    echo "Creating ${main_name}"
    cp ${TOOLS_DIR}/MAIN_${main_name} ${main_name}
fi

if [[ ! -z ${user_email} ]]; then
    sed -i 's/yourname@YOUR-INSTITUTE-DOMAIN/'"${user_email}"'/g' ${main_name}
fi
if [[ ! -z ${project_number} ]]; then
    sed -i 's/YOUR-PROJECT/'"${project_number}"'/g' ${main_name}
fi
echo ${user_config}
if [[ ! -z ${user_config} ]]; then
    sed -i 's^YOUR-CONFIG-FILE^'"${user_config}"'^g' ${main_name}
fi
