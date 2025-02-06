# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

noinst_PROGRAMS =

if ENABLE_MCI
if ENABLE_UTILS
if HAVE_NETCDF
noinst_PROGRAMS += \
    OASIS3_MCT_example/OASIS3_MCT_example.x \
    interpolation_benchmark/toy_scrip.x \
    perf_toy/perf_toy_icon.x \
    toy_atm_ocn/toy_icon_atm.x \
    toy_atm_ocn/toy_icon_ocn.x \
    toy_atm_ocn/toy_mpiom_ocn.x \
    toy_atm_ocn/toy_woa_ocn.x \
    toy_callback/toy_icon_callback.x \
    toy_icon_runoff/toy_icon_runoff_atm.x \
    toy_icon_runoff/toy_icon_runoff_ocn.x \
    toy_multi/toy_multi_unstruct.x \
    toy_multi/toy_multi_unstruct_2.x

# The following executables can be built without NetCDF but do not make sense
# without the executables from the list above:
noinst_PROGRAMS += \
    perf_toy/perf_toy_cube.x \
    toy_callback/toy_cube_callback.x \
    toy_multi/toy_multi_cube.x \
    toy_multi/toy_multi_reg2d.x
endif
endif ENABLE_UTILS
endif ENABLE_MCI

if ENABLE_MCI
if ENABLE_UTILS
noinst_PROGRAMS += \
    toy_atm_ocn/toy_reg2d_atm.x \
    toy_atm_ocn/toy_reg2d_ocn.x \
    toy_coupling/toy_coupling.x
endif
endif

if ENABLE_MCI
if HAVE_NETCDF
noinst_PROGRAMS += \
    toy_output/toy_output.x
endif
endif ENABLE_MCI

if ENABLE_MCI
noinst_PROGRAMS += \
    dynamic_config/dynamic_config.x \
    toy_dummy/dummy_atmosphere.x \
    toy_dummy/dummy_atmosphere_c.x \
    toy_dummy/dummy_io.x \
    toy_dummy/dummy_io_c.x \
    toy_dummy/dummy_ocean.x \
    toy_dummy/dummy_ocean_c.x
endif ENABLE_MCI

if ENABLE_UTILS
if HAVE_NETCDF
noinst_PROGRAMS += \
    access_core/access_core.x
endif
endif ENABLE_UTILS

if HAVE_NETCDF
noinst_PROGRAMS += \
    interpolation_benchmark/generate_OASIS_mask_file.x
endif

AM_CPPFLAGS = \
    -I$(top_srcdir)/src/mci \
    -I$(top_builddir)/src/core \
    -I$(top_builddir)/src/utils

AM_CFLAGS = \
    $(OPENMP_CFLAGS) \
    $(MPI_CFLAGS) \
    $(YAXT_CFLAGS) \
    $(NETCDF_CFLAGS) \
    $(FYAML_CFLAGS) \
    $(MTIME_CFLAGS)

AM_FCFLAGS = \
    $(FCMODINC)$(top_builddir)/src/mci \
    $(FCMODINC)$(top_builddir)/src/core \
    $(FCMODINC)$(top_builddir)/src/utils \
    $(MPI_FCFLAGS)

LDADD_COMMON =

if ENABLE_UTILS
LDADD_COMMON += \
    $(top_builddir)/src/utils/libyac_utils.a
endif ENABLE_UTILS

if ENABLE_MCI
LDADD_COMMON += \
    $(top_builddir)/src/mci/libyac_mci.a
endif ENABLE_MCI

LDADD_COMMON += \
    $(top_builddir)/src/core/libyac_core.a \
    $(LAPACK_CLIBS) \
    $(MTIME_CLIBS) \
    $(FYAML_CLIBS) \
    $(NETCDF_CLIBS) \
    $(YAXT_CLIBS)

LDADD = $(LDADD_COMMON) $(MPI_CLIBS) $(RPATH_CFLAGS)
FCLDADD = $(LDADD_COMMON) $(OPENMP_FCLIBS) $(MPI_FCLIBS) $(RPATH_FCFLAGS)

AM_DEFAULT_SOURCE_EXT = .c
EXTRA_DIST = toy_atm_ocn/toy_common.h

access_core_access_core_x_LDADD = $(FCLDADD)
access_core_access_core_x_SOURCES = access_core/access_core.F90

toy_dummy_dummy_atmosphere_x_LDADD = $(FCLDADD)
toy_dummy_dummy_atmosphere_x_SOURCES = toy_dummy/dummy_atmosphere.f90

toy_dummy_dummy_io_x_LDADD = $(FCLDADD)
toy_dummy_dummy_io_x_SOURCES = toy_dummy/dummy_io.f90

toy_dummy_dummy_ocean_x_LDADD = $(FCLDADD)
toy_dummy_dummy_ocean_x_SOURCES = toy_dummy/dummy_ocean.f90

toy_multi_toy_multi_cube_x_SOURCES = \
    toy_multi/toy_multi_common.c \
    toy_multi/toy_multi_common.h \
    toy_multi/toy_multi_cube.c

toy_multi_toy_multi_reg2d_x_SOURCES = \
    toy_multi/toy_multi_common.c \
    toy_multi/toy_multi_common.h \
    toy_multi/toy_multi_reg2d.c

toy_multi_toy_multi_unstruct_2_x_SOURCES = \
    toy_multi/toy_multi_common.c \
    toy_multi/toy_multi_common.h \
    toy_multi/toy_multi_unstruct_2.c

toy_multi_toy_multi_unstruct_x_SOURCES = \
    toy_multi/toy_multi_common.c \
    toy_multi/toy_multi_common.h \
    toy_multi/toy_multi_unstruct.c

toy_output_toy_output_x_SOURCES = toy_output/toy_output.c

clean-local:
	-rm -f */*.vtk

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
FCDEFS = $(FCPPDEF)HAVE_CONFIG_H
DEFAULT_FCINCLUDES = $(FCPPINC).@yac_fcisrc@ $(FCPPINC)$(top_builddir)/src
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_FCINCLUDES) $(FCINCLUDES) \
	$(AM_FCFLAGS) $(FCFLAGS)
