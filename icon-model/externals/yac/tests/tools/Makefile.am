# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

AM_TESTSUITE_SUMMARY_HEADER = ' for $(PACKAGE_STRING) ($(subdir))'
if ENABLE_OPENMP
AM_TESTS_ENVIRONMENT = OMP_NUM_THREADS=$${OMP_NUM_THREADS-1}; export OMP_NUM_THREADS;
endif

TEST_EXTENSIONS = .x .sh

test_PROGRAMS_ =

xfail_test_PROGRAMS_ =

# Build test programs:
check_PROGRAMS = $(test_PROGRAMS_) $(xfail_test_PROGRAMS_)

test_SCRIPTS_ = \
    test_weights2vtk.sh \
    test_debug_grid2vtk.sh

xfail_test_SCRIPTS_ =

# Build programs associated with test scripts:

if ENABLE_UTILS
if HAVE_NETCDF
check_PROGRAMS += \
    test_weights2vtk.x \
    test_debug_grid2vtk.x
endif
endif ENABLE_UTILS

TESTS = \
    $(test_PROGRAMS_) \
    $(test_SCRIPTS_) \
    $(xfail_test_PROGRAMS_) \
    $(xfail_test_SCRIPTS_)

XFAIL_TESTS = \
    $(xfail_test_PROGRAMS_) \
    $(xfail_test_SCRIPTS_)

if !ENABLE_CONCURRENT_MPI_TESTS
# Run MPI tests one by one:
#   The following dependencies are generated with
#   $(top_srcdir)/config/mpi_test_deps.sh.
endif

AM_CPPFLAGS = \
    -I$(top_builddir)/src/core \
    -I$(top_srcdir)/tests/common

AM_CFLAGS = \
    $(MPI_CFLAGS) \
    $(YAXT_CFLAGS) \
    $(NETCDF_CFLAGS)

LDADD = \
    $(top_builddir)/tests/common/libtest.a \
    $(top_builddir)/src/core/libyac_core.a \
    $(NETCDF_CLIBS) \
    $(YAXT_CLIBS) \
    $(MPI_CLIBS) \
    $(RPATH_CFLAGS)
