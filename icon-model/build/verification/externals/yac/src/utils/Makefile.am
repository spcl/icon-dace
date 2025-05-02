# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

if FCMODUC
yac_utils_mod = YAC_UTILS.$(FCMODEXT)
else
yac_utils_mod = yac_utils.$(FCMODEXT)
endif

if ENABLE_UTILS
lib_LIBRARIES = libyac_utils.a
nodist_include_HEADERS = \
    $(yac_utils_mod) \
    yac_utils.h
else !ENABLE_UTILS
check_LIBRARIES = libyac_utils.a
endif !ENABLE_UTILS

libyac_utils_a_SOURCES = \
    duplicate_stencils.c \
    duplicate_stencils.h \
    generate_cubed_sphere.c \
    generate_cubed_sphere.h \
    generate_reg2d.c \
    generate_reg2d.h \
    grid2vtk.c \
    grid2vtk.h \
    read_cube_csv_grid.c \
    read_cube_csv_grid.h \
    read_fesom_grid.c \
    read_fesom_grid.h \
    read_grid.c \
    read_grid.h \
    read_icon_grid.c \
    read_icon_grid.h \
    read_mpiom_grid.c \
    read_mpiom_grid.h \
    read_scrip_grid.c \
    read_scrip_grid.h \
    read_woa_data.c \
    read_woa_data.h \
    test_function.c \
    test_function.h \
    vtk_output.c \
    vtk_output.h \
    yac_utils.F90

yac_utils_h_headers = \
    generate_cubed_sphere.h \
    generate_reg2d.h \
    grid2vtk.h \
    read_cube_csv_grid.h \
    test_function.h \
    vtk_output.h \
    read_grid.h \
    read_fesom_grid.h \
    read_icon_grid.h \
    read_mpiom_grid.h \
    read_scrip_grid.h \
    read_woa_data.h \
    duplicate_stencils.h

AM_CPPFLAGS = \
    -I$(top_builddir)/src \
    -I$(top_srcdir)/src/core

AM_CFLAGS = \
    $(PIC_CFLAGS) \
    $(MPI_CFLAGS) \
    $(YAXT_CFLAGS) \
    $(NETCDF_CFLAGS)

AM_FCFLAGS = $(PIC_FCFLAGS)

$(yac_utils_mod): yac_utils.$(OBJEXT)
	@if test ! -f $@; then rm -f $<; $(MAKE) $(AM_MAKEFLAGS) $<; fi

all_local =

if ENABLE_UTILS

if ENABLE_DEPRECATED

if FCMODUC
mo_yac_utils_mod = MO_YAC_UTILS.$(FCMODEXT)
else
mo_yac_utils_mod = mo_yac_utils.$(FCMODEXT)
endif

nodist_include_HEADERS += \
    $(mo_yac_utils_mod)

libyac_utils_a_SOURCES += yac_utils_deprecated.f90

yac_utils_deprecated.$(OBJEXT): $(yac_utils_mod)

$(mo_yac_utils_mod): yac_utils_deprecated.$(OBJEXT)
	@if test ! -f $@; then rm -f $<; $(MAKE) $(AM_MAKEFLAGS) $<; fi

endif ENABLE_DEPRECATED

all_local += yac_utils.h

# REUSE-IgnoreStart
yac_utils.h: $(yac_utils_h_headers)
	$(AM_V_GEN):;{ \
	echo '// Copyright (c) 2024 The YAC Authors'; \
	echo '//'; \
	echo '// SPDX-License-Identifier: BSD-3-Clause'; \
	echo ''; \
	echo '#ifndef YAC_UTILS_H'; \
	echo '#define YAC_UTILS_H'; \
	echo ''; \
	echo '#include "yac_core.h"'; \
	echo ''; \
	for f in $(yac_utils_h_headers); do \
	  echo "// $$f"; \
	  sed '1,\|^// YAC PUBLIC HEADER START$$|d;\|^// YAC PUBLIC HEADER STOP$$|,$$d' `test -f "$$f" || echo '$(srcdir)/'`"$$f"; \
	  echo ''; \
	done; \
	echo '#endif // YAC_UTILS_H'; \
	} >$@
# REUSE-IgnoreEnd

MOSTLYCLEANFILES = yac_utils.h

endif ENABLE_UTILS

all-local: $(all_local)

mostlyclean-local:
	-rm -f *.$(FCMODEXT)

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
FCDEFS = $(FCPPDEF)HAVE_CONFIG_H
DEFAULT_FCINCLUDES = $(FCPPINC).@yac_fcisrc@ $(FCPPINC)$(top_builddir)/src
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_FCINCLUDES) $(FCINCLUDES) \
	$(AM_FCFLAGS) $(FCFLAGS)
