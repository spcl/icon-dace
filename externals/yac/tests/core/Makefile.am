# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

AM_TESTSUITE_SUMMARY_HEADER = ' for $(PACKAGE_STRING) ($(subdir))'
if ENABLE_OPENMP
AM_TESTS_ENVIRONMENT = OMP_NUM_THREADS=$${OMP_NUM_THREADS-1}; export OMP_NUM_THREADS;
endif

TEST_EXTENSIONS = .x .sh

test_PROGRAMS_ = \
    test_angle.x \
    test_area.x \
    test_basic_grid.x \
    test_basic_grid_data.x \
    test_bnd_sphere_part.x \
    test_cell_bnd_circle.x \
    test_circle.x \
    test_clipping.x \
    test_compute_overlap_area.x \
    test_dist_grid_utils.x \
    test_gcxgc.x \
    test_gcxlatc.x \
    test_geometry.x \
    test_interp_stack_config.x \
    test_interval_tree.x \
    test_lat_clipping.x \
    test_latcxlatc.x \
    test_loncxlatc.x \
    test_loncxlonc.x \
    test_mergesort.x \
    test_partial_areas.x \
    test_point_in_cell.x \
    test_point_sphere_part.x \
    test_pxgc.x \
    test_quicksort.x

xfail_test_SCRIPTS_ = \
    test_mpi_error.sh

# Build test programs:
check_PROGRAMS = $(test_PROGRAMS_) $(xfail_test_PROGRAMS_)

test_SCRIPTS_ = \
    test_basic_grid_parallel.sh \
    test_dist_grid_pair_parallel.sh \
    test_fortran_api.sh \
    test_group_comm.sh \
    test_interp_grid_parallel.sh \
    test_interp_method_avg_parallel.sh \
    test_interp_method_callback_parallel.sh \
    test_interp_method_conserv_parallel.sh \
    test_interp_method_creep_parallel.sh \
    test_interp_method_file_parallel.sh \
    test_interp_method_hcsbb_parallel.sh \
    test_interp_method_ncc_parallel.sh \
    test_interp_method_nnn_parallel.sh \
    test_interp_method_parallel.sh \
    test_interp_method_rbf_parallel.sh \
    test_interp_method_spmap_parallel.sh \
    test_interp_weights_parallel.sh \
    test_interpolation_exchange.sh \
    test_interpolation_parallel1.sh \
    test_interpolation_parallel2.sh \
    test_interpolation_parallel3.sh \
    test_interpolation_parallel4.sh \
    test_interpolation_parallel5.sh \
    test_interpolation_parallel6.sh \
    test_io_config.sh \
    test_openmp.sh \
    test_proc_sphere_part_parallel.sh

# Build programs associated with test scripts:

check_PROGRAMS += \
    test_openmp.x

if TEST_MPI

check_PROGRAMS += \
    test_group_comm.x \
    test_interp_method_conserv_parallel.x \
    test_interpolation_exchange.x \
    test_interpolation_parallel1.x \
    test_interpolation_parallel1_c.x \
    test_interpolation_parallel2.x \
    test_interpolation_parallel6.x \
    test_io_config.x \
    test_mpi_error.x \
    test_proc_sphere_part_parallel.x

if HAVE_NETCDF
check_PROGRAMS += \
    test_basic_grid_parallel.x \
    test_dist_grid_pair_parallel.x \
    test_fortran_api.x \
    test_interp_grid_parallel.x \
    test_interp_method_avg_parallel.x \
    test_interp_method_callback_parallel.x \
    test_interp_method_creep_parallel.x \
    test_interp_method_file_parallel.x \
    test_interp_method_hcsbb_parallel.x \
    test_interp_method_ncc_parallel.x \
    test_interp_method_nnn_parallel.x \
    test_interp_method_parallel.x \
    test_interp_method_rbf_parallel.x \
    test_interp_method_spmap_parallel.x \
    test_interp_weights_parallel.x \
    test_interpolation_parallel3.x \
    test_interpolation_parallel4.x \
    test_interpolation_parallel5.x
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
test_proc_sphere_part_parallel.log: test_mpi_error.log
test_mpi_error.log: test_io_config.log
test_io_config.log: test_interpolation_parallel6.log
test_interpolation_parallel6.log: test_interpolation_parallel5.log
test_interpolation_parallel5.log: test_interpolation_parallel4.log
test_interpolation_parallel4.log: test_interpolation_parallel3.log
test_interpolation_parallel3.log: test_interpolation_parallel2.log
test_interpolation_parallel2.log: test_interpolation_parallel1.log
test_interpolation_parallel1.log: test_interpolation_exchange.log
test_interpolation_exchange.log: test_interp_weights_parallel.log
test_interp_weights_parallel.log: test_interp_method_spmap_parallel.log
test_interp_method_spmap_parallel.log: test_interp_method_rbf_parallel.log
test_interp_method_rbf_parallel.log: test_interp_method_parallel.log
test_interp_method_parallel.log: test_interp_method_nnn_parallel.log
test_interp_method_nnn_parallel.log: test_interp_method_ncc_parallel.log
test_interp_method_ncc_parallel.log: test_interp_method_hcsbb_parallel.log
test_interp_method_hcsbb_parallel.log: test_interp_method_file_parallel.log
test_interp_method_file_parallel.log: test_interp_method_creep_parallel.log
test_interp_method_creep_parallel.log: test_interp_method_conserv_parallel.log
test_interp_method_conserv_parallel.log: test_interp_method_callback_parallel.log
test_interp_method_callback_parallel.log: test_interp_method_avg_parallel.log
test_interp_method_avg_parallel.log: test_interp_grid_parallel.log
test_interp_grid_parallel.log: test_group_comm.log
test_group_comm.log: test_dist_grid_pair_parallel.log
test_dist_grid_pair_parallel.log: test_basic_grid_parallel.log
test_basic_grid_parallel.log: test_fortran_api.log
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

AM_FCFLAGS = \
    $(FCPPINC)$(top_srcdir)/tests/common \
    $(FCMODINC)$(top_builddir)/tests/common \
    $(FCMODINC)$(top_builddir)/src/core \
    $(MPI_FCFLAGS)

LDADD_COMMON = \
    $(top_builddir)/tests/common/libtest.a \
    $(top_builddir)/src/core/libyac_core.a \
    $(LAPACK_CLIBS) \
    $(NETCDF_CLIBS)

LDADD = $(LDADD_COMMON) $(YAXT_CLIBS) $(MPI_CLIBS) $(RPATH_CFLAGS)
FCLDADD = $(LDADD_COMMON) $(OPENMP_FCLIBS) $(YAXT_CLIBS) $(MPI_FCLIBS) $(RPATH_FCFLAGS)

test_dist_grid_pair_parallel_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(LDADD)

test_fortran_api_x_LDADD = $(FCLDADD)
test_fortran_api_x_SOURCES = test_fortran_api.F90
test_fortran_api.$(OBJEXT): $(utest_FCDEPS)

test_gcxgc_x_SOURCES = test_gcxgc.c test_cxc.c test_cxc.h

test_gcxlatc_x_SOURCES = test_gcxlatc.c test_cxc.c test_cxc.h

test_interp_grid_parallel_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(LDADD)

test_interp_method_parallel_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(LDADD)

test_interp_weights_parallel_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(LDADD)

test_interpolation_parallel1_x_LDADD = $(FCLDADD)
test_interpolation_parallel1_x_SOURCES = test_interpolation_parallel1.F90
test_interpolation_parallel1.$(OBJEXT): $(utest_FCDEPS)

test_latcxlatc_x_SOURCES = test_latcxlatc.c test_cxc.c test_cxc.h

test_loncxlatc_x_SOURCES = test_loncxlatc.c test_cxc.c test_cxc.h

test_loncxlonc_x_SOURCES = test_loncxlonc.c test_cxc.c test_cxc.h

test_pxgc_x_SOURCES = test_pxgc.c test_cxc.c test_cxc.h

test_basic_grid_parallel_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(LDADD)

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
FCDEFS = $(FCPPDEF)HAVE_CONFIG_H
DEFAULT_FCINCLUDES = $(FCPPINC).@yac_fcisrc@ $(FCPPINC)$(top_builddir)/src
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_FCINCLUDES) $(FCINCLUDES) \
	$(AM_FCFLAGS) $(FCFLAGS)
