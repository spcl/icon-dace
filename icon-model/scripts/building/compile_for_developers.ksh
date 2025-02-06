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

#  compdev_mpipc.ksh
#  Simple script to quickly compile and link files to update solely the executable
#  "control_model" without compiling sources other than specified

#  Author: Stephan Lorenz - 2010/02
#  Modified - 2011/05

# Compiles a module file for control_model
#   - much faster than "make" in case of simple development of single routines
#   - runs after successful configure/make to update control_model
#   - several fortran files can be chosen as input
#   - without arguments only the linker is invoked
#   - the current directory must be below the ICON main directory (e.g. run or src)

#  - caution: new module binaries are written to src/module.o, i.e. not to modules/module.mod 
#  - caution: compile as much as necessary since no dependicies as in make are known

# Revisions:
#  - now working directly in icon-dev/src - compiling only two sources (2010-06-07)
#  - now reading objects from Makefile (2010-08-02)
#  - now compiling sources from input arguments (2010-11-18)
#  - now running on all system via $arch (2010-11-18)
#  - also on multiple builds (2011-02-09)

# Further work:
#  - compiler options should be read from Makefile (not yet)

set -e
# you must start this script from below the trunk/branch - e.g. from icon-dev/run
# in case of multiple builds start below the specific dir - e.g. icon-dev/nag_nMnO/run
dir=$(pwd -P)
basedir=${dir%/*}
# determine architecture
arch=`ls ${basedir}/build`

# go to basedir
cd ..
# get compiler from Makefile
fccomp=$(grep FC Makefile)
COMP=${fccomp##* }

# go to src directory:
cd src

# objects to link are the same for all compilers but may change with revisions:
#  - get objects from Makefile
OBS=$(sed -n '/OBJS =/,/^$/p' ../build/$arch/src/Makefile | sed 's/OBJS =//' | sed s'/\\//')
objects="-o ../bin/control_model control_model.o $OBS"
#echo objects are: $objects


if [ $COMP == "gfortran" ] ; then

  # compiler and options gfortran:
  compopt="gfortran -I../include -I/sw/lenny-x64/netcdf-4.1.1-static-gcc45/include -I/sw/lenny-x64/hdf5-1.8.5-p1-static/include -I/sw/lenny-x64/szip-2.1-static/include -I/usr/include -J../module -I../module -march=native -O0 -ffast-math -D__LOOP_EXCHANGE -xf95-cpp-input -std=f2003 -fmodule-private -fimplicit-none -fmax-identifier-length=31 -ffree-line-length-99 -Wall -Wcharacter-truncation -Wconversion -Wunderflow -Wunused-parameter -g -fbacktrace -fbounds-check -D__ICON__ -DNOMPI -c ../../../src"

  # loader with objects and options gfortran:
  loadobj="gfortran -march=native -O0 -ffast-math -D__LOOP_EXCHANGE $objects \
          -L../lib -lsupport -L../lib -llapack -lblas -L/sw/lenny-x64/netcdf-4.1.1-static-gcc45/lib -lnetcdff -lnetcdf -L/sw/lenny-x64/hdf5-1.8.5-p1-static/lib -lhdf5_hl -lhdf5 -L/sw/lenny-x64/szip-2.1-static/lib -lsz -L/usr/lib -lz"

elif [ $COMP == "nagfor" ] ; then

  # compiler and options nag:
  compopt="nagfor -I../include -I/sw/lenny-x64/netcdf-4.1.1-static-nag52/include -I/sw/lenny-x64/hdf5-1.8.5-p1-static/include -I/sw/lenny-x64/szip-2.1-static/include -I/usr/include    -mdir ../module -I../module -float-store -colour -nan -maxcontin=99 -fpp -f2003 -gline -g -C=all -wmismatch=mpi_send,mpi_isend,mpi_recv,mpi_irecv,mpi_bcast,mpi_allreduce,mpi_gather,mpi_gatherv,mpi_scatter,mpi_scatterv,nf_put_att_double,nf_def_var,nf_get_att_int,nf_put_att_int -w=uep  -D__ICON__ -DNOMPI  -c ../../../src"

  # loader with objects and options nag:
  # loadobj="nagfor -I../module -float-store -gline -g -maxcontin=99 -fpp -f2003 -C=all -wmismatch=mpi_send,mpi_isend,mpi_recv,mpi_irecv,mpi_bcast,mpi_allreduce,mpi_gather,mpi_gatherv,mpi_scatter,mpi_scatterv,nf_put_att_double,nf_def_var,nf_get_att_int -w=uep $objects -L../lib -lsupport  -L../lib -llapack -lblas -L/sw/etch-ia32/netcdf-3.6.3/lib -lnetcdf /zmaw/sw/sw/etch-ia32/gcc-4.3.3/bin/../lib/gcc/$arch/4.3.3/libgcc.a    -L/sw/etch-ia32/mpich2-1.2.1-nag52/lib -lmpichf90 -lmpich -lpthread -lrt"
  loadobj="nagfor  -I../module -float-store -colour -nan -maxcontin=99 -fpp -f2003 -gline -g -C=all -wmismatch=mpi_send,mpi_isend,mpi_recv,mpi_irecv,mpi_bcast,mpi_allreduce,mpi_gather,mpi_gatherv,mpi_scatter,mpi_scatterv,nf_put_att_double,nf_def_var,nf_get_att_int,nf_put_att_int -w=uep  -D__ICON__ -DNOMPI $objects  -L../lib -lsupport -L../lib -llapack -lblas -L/sw/lenny-x64/netcdf-4.1.1-static-nag52/lib -lnetcdff -lnetcdf -L/sw/lenny-x64/hdf5-1.8.5-p1-static/lib -lhdf5_hl -lhdf5 -L/sw/lenny-x64/szip-2.1-static/lib -lsz -L/usr/lib -lz"

#elif [ $COMP == "ifort" ] ; then

  # compiler and options intel:
  # compopt="ifort -I../include -I/sw/etch-ia32/netcdf-3.6.3/include    -I/sw/etch-ia32/mpich2-1.2.1-intel11/include   -module ../module -I../module -O3 -msse2 -mieee-fp -fpe0 -pc64 -fpp -traceback -DNOMPI -c ../../../src" 
  # loader with objects and options intel:
  # loadobj="ifort -I../module -O3 -msse2 -mieee-fp -fpe0 -pc64 -fpp -traceback $objects -L/sw/etch-ia32/mpich2-1.2.1-intel11/lib -lmpichf90 -lmpich -lpthread -lrt"

#elif [ $COMP == "pgi" ] ; then

  # compiler and options pgi:
  # compopt="pgf95 -I../include -I/sw/etch-ia32/netcdf-3.6.3/include -I/sw/etch-ia32/mpich2-1.2.1-pgi9/include -module ../module -I../module -fast -fastsse -g -Mpreprocess -Mchkstk  -Mrecursive -c ../../../src"

  # loader with objects and options pgi:
  # loadobj="pgf95 -I../module -fast -fastsse -g -Mpreprocess -Mchkstk  -Mrecursive $objects /zmaw/sw/sw/etch-ia32/gcc-4.3.3/bin/../lib/gcc/$arch/4.3.3/libgcc.a    -L/sw/etch-ia32/mpich2-1.2.1-pgi9/lib -lmpichf90 -lmpich -lpthread -lrt"

else
  echo " compiler not setup yet"
  exit
fi

# change to build directory - starting from below the trunk/branch - e.g. from run
cd ../build/$arch/src

# do compiling and loading
#  - files to compile in input
for i do
  echo " *** *** *** compiler $COMP compiles file $i"
  $compopt/$i
  shift
done

echo " *** *** *** linking objects to control_model"
$loadobj

