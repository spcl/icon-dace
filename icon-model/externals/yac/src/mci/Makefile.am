# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

lib_LIBRARIES = libyac_mci.a

include_HEADERS = \
    yac.h

if FCMODUC
yac_mod = YAC.$(FCMODEXT)
else
yac_mod = yac.$(FCMODEXT)
endif

nodist_include_HEADERS = \
    $(yac_mod)

libyac_mci_a_SOURCES = \
    component.c \
    component.h \
    config_yaml.c \
    config_yaml.h \
    couple_config.c \
    couple_config.h \
    dist_merge.c \
    dist_merge.h \
    event.c \
    event.h \
    fields.c \
    fields.h \
    instance.c \
    instance.h \
    mpi_handshake.c \
    mpi_handshake.h \
    utils_mci.h \
    version.h \
    yac.c \
    yac.h \
    yac_global.F90 \
    yac_module.F90

AM_CPPFLAGS = \
    -I$(top_srcdir)/src/core

AM_CFLAGS = \
    $(PIC_CFLAGS) \
    $(OPENMP_CFLAGS) \
    $(MPI_CFLAGS) \
    $(YAXT_CFLAGS) \
    $(NETCDF_CFLAGS) \
    $(FYAML_CFLAGS) \
    $(MTIME_CFLAGS)

AM_FCFLAGS = $(PIC_FCFLAGS)

yac_global.$(OBJEXT): $(yac_mod)

$(yac_mod): yac_module.$(OBJEXT)
	@if test ! -f $@; then rm -f $<; $(MAKE) $(AM_MAKEFLAGS) $<; fi

all_local =

if ENABLE_DEPRECATED

include_HEADERS += \
    yac_interface.h

if FCMODUC
mo_yac_finterface_mod = MO_YAC_FINTERFACE.$(FCMODEXT)
else
mo_yac_finterface_mod = mo_yac_finterface.$(FCMODEXT)
endif

nodist_include_HEADERS += \
    $(mo_yac_finterface_mod)

libyac_mci_a_SOURCES += yac_mci_deprecated.f90

yac_mci_deprecated.$(OBJEXT): $(yac_mod)

$(mo_yac_finterface_mod): yac_mci_deprecated.$(OBJEXT)
	@if test ! -f $@; then rm -f $<; $(MAKE) $(AM_MAKEFLAGS) $<; fi

all_local += manifest.txt

manifest.txt: $(libyac_mci_a_OBJECTS)
	$(AM_V_GEN):;{ \
	  for obj in $(libyac_mci_a_OBJECTS); do \
	    echo "$(abs_builddir)/$$obj"; \
	  done; \
	} >$@

MOSTLYCLEANFILES = manifest.txt

endif ENABLE_DEPRECATED

all-local: $(all_local)

mostlyclean-local:
	-rm -f *.$(FCMODEXT)

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
FCDEFS = $(FCPPDEF)HAVE_CONFIG_H
DEFAULT_FCINCLUDES = $(FCPPINC).@yac_fcisrc@ $(FCPPINC)$(top_builddir)/src
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_FCINCLUDES) $(FCINCLUDES) \
	$(AM_FCFLAGS) $(FCFLAGS)
