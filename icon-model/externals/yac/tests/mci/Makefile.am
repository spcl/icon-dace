# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

AM_TESTSUITE_SUMMARY_HEADER = ' for $(PACKAGE_STRING) ($(subdir))'
if ENABLE_OPENMP
AM_TESTS_ENVIRONMENT = OMP_NUM_THREADS=$${OMP_NUM_THREADS-1}; export OMP_NUM_THREADS;
endif

TEST_EXTENSIONS = .x .sh

test_PROGRAMS_ = \
    test_events.x

xfail_test_PROGRAMS_ =

# Build test programs:
check_PROGRAMS = $(test_PROGRAMS_) $(xfail_test_PROGRAMS_)

test_SCRIPTS_ = \
    test_abort.sh \
    test_component_config.sh \
    test_couple_config.sh \
    test_def_comps.sh \
    test_def_datetime.sh \
    test_def_grid.sh \
    test_def_mask.sh \
    test_def_points.sh \
    test_dummy_coupling.sh \
    test_dummy_coupling2_c.sh \
    test_dummy_coupling3.sh \
    test_dummy_coupling3_c.sh \
    test_dummy_coupling4_c.sh \
    test_dummy_coupling5.sh \
    test_dummy_coupling6.sh \
    test_dummy_coupling6_c.sh \
    test_dummy_coupling7.sh \
    test_dummy_coupling7_c.sh \
    test_dummy_coupling8_c.sh \
    test_dummy_coupling9.sh \
    test_dynamic_config.sh \
    test_init_comm_final.sh \
    test_init_final.sh \
    test_instance_parallel1.sh \
    test_instance_parallel2.sh \
    test_instance_parallel3.sh \
    test_instance_parallel4.sh \
    test_mpi_handshake.sh \
    test_mpi_handshake_c.sh \
    test_multithreading.sh \
    test_query_routines.sh \
    test_query_routines_c.sh \
    test_restart.sh \
    test_restart2.sh \
    test_version.sh


# Build programs associated with test scripts:
check_PROGRAMS += \
    test_abort.x \
    test_abort_c.x \
    test_def_comps.x \
    test_def_datetime.x \
    test_def_grid.x \
    test_def_mask.x \
    test_def_points.x \
    test_init_comm_final.x \
    test_init_final.x \
    test_version.x

if TEST_MPI

check_PROGRAMS += \
    test_component_config.x \
    test_couple_config.x \
    test_dummy_coupling4_c.x \
    test_dummy_coupling5_c.x \
    test_dummy_coupling5_dble.x \
    test_dummy_coupling5_real.x \
    test_dummy_coupling6_c.x \
    test_dummy_coupling7_c.x \
    test_dummy_coupling8_c.x \
    test_dummy_coupling9.x \
    test_dummy_coupling9_c.x \
    test_dummy_coupling_c.x \
    test_dummy_coupling_dble.x \
    test_dummy_coupling_real.x \
    test_instance_parallel2.x \
    test_instance_parallel3.x \
    test_instance_parallel4.x \
    test_mpi_handshake_c.x \
    test_query_routines_c.x \
    test_restart2.x

if HAVE_NETCDF
check_PROGRAMS += \
    test_dummy_coupling2_c.x \
    test_dummy_coupling3_c.x \
    test_dynamic_config.x \
    test_dynamic_config_c.x \
    test_instance_parallel1.x \
    test_restart_c.x \
    test_restart_dble.x
endif

if TEST_PTHREAD
check_PROGRAMS += \
    test_multithreading.x
endif

if TEST_YAXT_FC

check_PROGRAMS += \
    test_dummy_coupling6_dble.x \
    test_dummy_coupling6_real.x \
    test_dummy_coupling7_dble.x \
    test_dummy_coupling7_real.x \
    test_mpi_handshake.x \
    test_query_routines.x

if HAVE_NETCDF
check_PROGRAMS += \
    test_dummy_coupling3_dble.x \
    test_dummy_coupling3_real.x
endif

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
test_restart2.log: test_restart.log
test_restart.log: test_query_routines_c.log
test_query_routines_c.log: test_query_routines.log
test_query_routines.log: test_multithreading.log
test_multithreading.log: test_mpi_handshake_c.log
test_mpi_handshake_c.log: test_mpi_handshake.log
test_mpi_handshake.log: test_instance_parallel4.log
test_instance_parallel4.log: test_instance_parallel3.log
test_instance_parallel3.log: test_instance_parallel2.log
test_instance_parallel2.log: test_instance_parallel1.log
test_instance_parallel1.log: test_init_final.log
test_init_final.log: test_init_comm_final.log
test_init_comm_final.log: test_dynamic_config.log
test_dynamic_config.log: test_dummy_coupling9.log
test_dummy_coupling9.log: test_dummy_coupling8_c.log
test_dummy_coupling8_c.log: test_dummy_coupling7_c.log
test_dummy_coupling7_c.log: test_dummy_coupling7.log
test_dummy_coupling7.log: test_dummy_coupling6_c.log
test_dummy_coupling6_c.log: test_dummy_coupling6.log
test_dummy_coupling6.log: test_dummy_coupling5.log
test_dummy_coupling5.log: test_dummy_coupling4_c.log
test_dummy_coupling4_c.log: test_dummy_coupling3_c.log
test_dummy_coupling3_c.log: test_dummy_coupling3.log
test_dummy_coupling3.log: test_dummy_coupling2_c.log
test_dummy_coupling2_c.log: test_dummy_coupling.log
test_dummy_coupling.log: test_def_points.log
test_def_points.log: test_def_mask.log
test_def_mask.log: test_def_grid.log
test_def_grid.log: test_def_datetime.log
test_def_datetime.log: test_def_comps.log
test_def_comps.log: test_couple_config.log
test_couple_config.log: test_component_config.log
test_component_config.log: test_abort.log
endif
endif

AM_CPPFLAGS = \
    -I$(top_srcdir)/src/mci \
    -I$(top_srcdir)/src/core \
    -I$(top_srcdir)/src/utils \
    -I$(top_srcdir)/tests/common

AM_CFLAGS = \
    $(OPENMP_CFLAGS) \
    $(MPI_CFLAGS) \
    $(YAXT_CFLAGS) \
    $(NETCDF_CFLAGS) \
    $(FYAML_CFLAGS) \
    $(MTIME_CFLAGS)

AM_FCFLAGS = \
    $(FCPPINC)$(top_srcdir)/tests/common \
    $(FCMODINC)$(top_builddir)/tests/common \
    $(FCMODINC)$(top_builddir)/src/mci \
    $(MPI_FCFLAGS) \
    $(YAXT_FCFLAGS)

LDADD_COMMON = \
    $(top_builddir)/tests/common/libtest.a \
    $(top_builddir)/src/mci/libyac_mci.a \
    $(top_builddir)/src/core/libyac_core.a \
    $(LAPACK_CLIBS) \
    $(MTIME_CLIBS) \
    $(FYAML_CLIBS) \
    $(NETCDF_CLIBS)

LDADD = $(LDADD_COMMON) $(YAXT_CLIBS) $(MPI_CLIBS) $(RPATH_CFLAGS)
FCLDADD = $(LDADD_COMMON) $(OPENMP_FCLIBS)
if TEST_YAXT_FC
FCLDADD += $(YAXT_FCLIBS)
else
FCLDADD += $(YAXT_CLIBS)
endif
FCLDADD += $(MPI_FCLIBS) $(RPATH_FCFLAGS)

test_def_comps_x_LDADD = $(FCLDADD)
test_def_comps_x_SOURCES = test_def_comps.F90
test_def_comps.$(OBJEXT): $(utest_FCDEPS)

test_def_datetime_x_LDADD = $(FCLDADD)
test_def_datetime_x_SOURCES = test_def_datetime.F90
test_def_datetime.$(OBJEXT): $(utest_FCDEPS)

test_def_grid_x_LDADD = $(FCLDADD)
test_def_grid_x_SOURCES = test_def_grid.F90
test_def_grid.$(OBJEXT): $(utest_FCDEPS)

test_def_mask_x_LDADD = $(FCLDADD)
test_def_mask_x_SOURCES = test_def_mask.F90
test_def_mask.$(OBJEXT): $(utest_FCDEPS)

test_def_points_x_LDADD = $(FCLDADD)
test_def_points_x_SOURCES = test_def_points.F90
test_def_points.$(OBJEXT): $(utest_FCDEPS)

test_dummy_coupling2_c_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(LDADD)

test_dummy_coupling3_c_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(LDADD)

test_dummy_coupling9_x_LDADD = $(FCLDADD)
test_dummy_coupling9_x_SOURCES = test_dummy_coupling9.F90
test_dummy_coupling9.$(OBJEXT): $(utest_FCDEPS)

test_dummy_coupling_dble_x_LDADD = $(FCLDADD)
test_dummy_coupling_dble_x_SOURCES = test_dummy_coupling_dble.F90 test_dummy_coupling.inc
test_dummy_coupling_dble.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling.inc

test_dummy_coupling_real_x_LDADD = $(FCLDADD)
test_dummy_coupling_real_x_SOURCES = test_dummy_coupling_real.F90 test_dummy_coupling.inc
test_dummy_coupling_real.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling.inc

test_dummy_coupling3_dble_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(FCLDADD)
test_dummy_coupling3_dble_x_SOURCES = test_dummy_coupling3_dble.F90 test_dummy_coupling3.inc
test_dummy_coupling3_dble.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling3.inc

test_dummy_coupling3_real_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(FCLDADD)
test_dummy_coupling3_real_x_SOURCES = test_dummy_coupling3_real.F90 test_dummy_coupling3.inc
test_dummy_coupling3_real.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling3.inc

test_dummy_coupling5_dble_x_LDADD = $(FCLDADD)
test_dummy_coupling5_dble_x_SOURCES = test_dummy_coupling5_dble.F90 test_dummy_coupling5.inc
test_dummy_coupling5_dble.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling5.inc

test_dummy_coupling5_real_x_LDADD = $(FCLDADD)
test_dummy_coupling5_real_x_SOURCES = test_dummy_coupling5_real.F90 test_dummy_coupling5.inc
test_dummy_coupling5_real.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling5.inc

test_dummy_coupling6_dble_x_LDADD = $(FCLDADD)
test_dummy_coupling6_dble_x_SOURCES = test_dummy_coupling6_dble.F90 test_dummy_coupling6.inc
test_dummy_coupling6_dble.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling6.inc

test_dummy_coupling6_real_x_LDADD = $(FCLDADD)
test_dummy_coupling6_real_x_SOURCES = test_dummy_coupling6_real.F90 test_dummy_coupling6.inc
test_dummy_coupling6_real.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling6.inc

test_dummy_coupling7_dble_x_LDADD = $(FCLDADD)
test_dummy_coupling7_dble_x_SOURCES = test_dummy_coupling7_dble.F90 test_dummy_coupling7.inc
test_dummy_coupling7_dble.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling7.inc

test_dummy_coupling7_real_x_LDADD = $(FCLDADD)
test_dummy_coupling7_real_x_SOURCES = test_dummy_coupling7_real.F90 test_dummy_coupling7.inc
test_dummy_coupling7_real.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling7.inc

test_restart_c_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(LDADD)

test_restart_dble_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(FCLDADD)
test_restart_dble_x_SOURCES = test_restart_dble.F90
test_restart_dble.$(OBJEXT): $(utest_FCDEPS)

test_restart2_x_LDADD = $(top_builddir)/src/utils/libyac_utils.a $(LDADD)

test_mpi_handshake_x_LDADD = $(FCLDADD)
test_mpi_handshake_x_SOURCES = test_mpi_handshake.F90
test_mpi_handshake.$(OBJEXT): $(utest_FCDEPS)

test_init_final_x_LDADD = $(FCLDADD)
test_init_final_x_SOURCES = test_init_final.F90
test_init_final.$(OBJEXT): $(utest_FCDEPS)

test_init_comm_final_x_LDADD = $(FCLDADD)
test_init_comm_final_x_SOURCES = test_init_comm_final.F90
test_init_comm_final.$(OBJEXT): $(utest_FCDEPS)

test_multithreading_x_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
test_multithreading_x_LDADD = $(PTHREAD_LIBS) $(LDADD)

test_version_x_LDADD = $(FCLDADD)
test_version_x_SOURCES = test_version.F90
test_version.$(OBJEXT): $(utest_FCDEPS)

test_abort_x_LDADD = $(FCLDADD)
test_abort_x_SOURCES = test_abort.F90

test_dynamic_config_x_LDADD = $(FCLDADD)
test_dynamic_config_x_SOURCES = test_dynamic_config.F90
test_dynamic_config.$(OBJEXT): $(utest_FCDEPS)

test_query_routines_x_LDADD = $(FCLDADD)
test_query_routines_x_SOURCES = test_query_routines.F90
test_query_routines.$(OBJEXT): $(utest_FCDEPS)

mostlyclean-local:
	-rm -f *.$(FCMODEXT)

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
FCDEFS = $(FCPPDEF)HAVE_CONFIG_H
DEFAULT_FCINCLUDES = $(FCPPINC).@yac_fcisrc@ $(FCPPINC)$(top_builddir)/src
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_FCINCLUDES) $(FCINCLUDES) \
	$(AM_FCFLAGS) $(FCFLAGS)
