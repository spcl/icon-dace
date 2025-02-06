# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

noinst_PROGRAMS =

if ENABLE_UTILS
if HAVE_NETCDF
noinst_PROGRAMS += \
    weights2vtk.x \
    debug_grid2vtk.x
endif
endif ENABLE_UTILS

AM_CPPFLAGS = \
    -I$(top_builddir)/src/utils \
    -I$(top_builddir)/src/core

AM_CFLAGS = \
    $(MPI_CFLAGS) \
    $(YAXT_CFLAGS) \
    $(NETCDF_CFLAGS)

weights2vtk_x_LDADD = \
    $(top_builddir)/src/utils/libyac_utils.a \
    $(top_builddir)/src/core/libyac_core.a \
    $(NETCDF_CLIBS) \
    $(YAXT_CLIBS) \
    $(MPI_CLIBS) \
    $(RPATH_CFLAGS)

debug_grid2vtk_x_LDADD = \
    $(top_builddir)/src/utils/libyac_utils.a \
    $(top_builddir)/src/core/libyac_core.a \
    $(NETCDF_CLIBS) \
    $(YAXT_CLIBS) \
    $(MPI_CLIBS) \
    $(RPATH_CFLAGS)

clean-local:
	-rm -rf *.x.dSYM
