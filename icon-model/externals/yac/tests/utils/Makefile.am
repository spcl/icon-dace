# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

AM_TESTSUITE_SUMMARY_HEADER = ' for $(PACKAGE_STRING) ($(subdir))'
if ENABLE_OPENMP
AM_TESTS_ENVIRONMENT = OMP_NUM_THREADS=$${OMP_NUM_THREADS-1}; export OMP_NUM_THREADS;
endif

TEST_EXTENSIONS = .x .sh

test_PROGRAMS_ = \
    test_grid2vtk.x \
    test_vtk_output.x

xfail_test_PROGRAMS_ =

# Build test programs:
check_PROGRAMS = $(test_PROGRAMS_) $(xfail_test_PROGRAMS_)

test_SCRIPTS_ = \
    test_duplicate_stencils_parallel.sh \
    test_generate_cubed_sphere.sh \
    test_read_cube_csv.sh \
    test_read_fesom.sh \
    test_read_icon.sh \
    test_read_icon_parallel.sh \
    test_read_mpiom.sh \
    test_read_scrip.sh \
    test_read_scrip_parallel.sh \
    test_read_woa_data.sh

xfail_test_SCRIPTS_ =

# Build programs associated with test scripts:
check_PROGRAMS += \
    test_read_cube_csv.x

if HAVE_NETCDF
check_PROGRAMS += \
    test_generate_cubed_sphere.x \
    test_read_fesom.x \
    test_read_icon.x \
    test_read_mpiom.x \
    test_read_scrip.x \
    test_read_woa_data.x
endif

if TEST_MPI

check_PROGRAMS += \
    test_duplicate_stencils_parallel.x

if HAVE_NETCDF
check_PROGRAMS += \
    test_read_icon_parallel.x \
    test_read_icon_parallel.x \
    test_read_scrip_parallel.x
endif

endif

TESTS = \
    $(test_PROGRAMS_) \
    $(test_SCRIPTS_) \
    $(xfail_test_PROGRAMS_) \
    $(xfail_test_SCRIPTS_)

XFAIL_TESTS = \
    $(xfail_test_PROGRAMS_) \
    $(xfail_test_SCRIPTS_)

if TEST_MPI
if !ENABLE_CONCURRENT_MPI_TESTS
# Run MPI tests one by one:
#   The following dependencies are generated with
#   $(top_srcdir)/config/mpi_test_deps.sh.
test_read_scrip_parallel.log: test_read_icon_parallel.log
test_read_icon_parallel.log: test_duplicate_stencils_parallel.log
endif
endif

AM_CPPFLAGS = \
    -I$(top_srcdir)/src/core \
    -I$(top_srcdir)/src/utils \
    -I$(top_srcdir)/tests/common

AM_CFLAGS = \
    $(OPENMP_CFLAGS) \
    $(MPI_CFLAGS) \
    $(YAXT_CFLAGS) \
    $(NETCDF_CFLAGS)

LDADD = \
    $(top_builddir)/src/utils/libyac_utils.a \
    $(top_builddir)/tests/common/libtest.a \
    $(top_builddir)/src/core/libyac_core.a \
    $(LAPACK_CLIBS) \
    $(FYAML_CLIBS) \
    $(NETCDF_CLIBS) \
    $(YAXT_CLIBS) \
    $(MPI_CLIBS) \
    $(RPATH_CFLAGS)
