# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

# Makefile for ww_main
#
PROGRAM = ww_main

VPATH = ../../../src/atm_phy_nwp

FC = ifort

FFDEBUG = # -debug -ftrapuv -warn -fp-stack-check -fstack-protector-all -fpe-all=0 -check arg_temp_created,format,output_conversion,uninit,pointers,stack,bounds
FFLAGS  = -I. ${GRIB_API_INCLUDE} -O3 -fpp1 -traceback -fpe0 -fp-speculation safe -fp-model source -sox $(FFDEBUG) -list -Winline -DONLYWW

LD  = $(FC)

LDFLAGS  = ${GRIB_API_LIB_FTN} # -Wl,"--print-map"

OBJS  = mo_wwonly.o mo_nwp_ww.o ww_main.o

.SUFFIXES: .f90

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

$(PROGRAM): $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $(PROGRAM)

mo_wwonly.o: mo_wwonly.f90

mo_nwp_ww.o: mo_nwp_ww.f90 mo_wwonly.o

ww_main.o: ww_main.f90 mo_nwp_ww.o mo_wwonly.o

clean:
	 rm $(OBJS) mo_nwp_ww.mod mo_wwonly.mod
