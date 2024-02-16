AM_TESTS_ENVIRONMENT = $(SHLIB_ENVIRONMENT_MAKE)

TEST_EXTENSIONS = .x .sh

test_PROGRAMS_ =                                    \
        test_angle.x                                \
        test_area.x                                 \
        test_cell_bnd_circle.x                      \
        test_circle.x                               \
        test_clipping.x                             \
        test_lat_clipping.x                         \
        test_compute_overlap_area.x                 \
        test_dist_grid_utils.x                      \
        test_events.x                               \
        test_pxgc.x                                 \
        test_gcxgc.x                                \
        test_gcxlatc.x                              \
        test_generate_cubed_sphere.x                \
        test_geometry.x                             \
        test_grid.x                                 \
        test_grid2vtk.x                             \
        test_interval_tree.x                        \
        test_latcxlatc.x                            \
        test_loncxlatc.x                            \
        test_loncxlonc.x                            \
        test_mergesort.x                            \
        test_partial_areas.x                        \
        test_point_in_cell.x                        \
        test_bnd_sphere_part.x                      \
        test_point_sphere_part.x                    \
        test_quicksort.x                            \
        test_read_cube_csv.x                        \
        test_vtk_output.x                           \
        test_interp_stack_config.x

xfail_test_PROGRAMS_ =

# Build test programs:
check_PROGRAMS = $(test_PROGRAMS_) $(xfail_test_PROGRAMS_)

test_SCRIPTS_ =                                        \
        test_abort.sh                                  \
        test_component_config.sh                       \
        test_couple_config.sh                          \
        test_instance_parallel1.sh                     \
        test_instance_parallel2.sh                     \
        test_instance_parallel3.sh                     \
        test_instance_parallel4.sh                     \
        test_def_comps.sh                              \
        test_def_datetime.sh                           \
        test_def_grid.sh                               \
        test_def_mask.sh                               \
        test_def_points.sh                             \
        test_dist_grid_pair_parallel.sh                \
        test_group_comm.sh                             \
        test_interp_grid_parallel.sh                   \
        test_interp_method_parallel.sh                 \
        test_interp_method_avg_parallel.sh             \
        test_interp_method_callback_parallel.sh        \
        test_interp_method_conserv_parallel.sh         \
        test_interp_method_creep_parallel.sh           \
        test_interp_method_nnn_parallel.sh             \
        test_interp_method_nnn_parallel2.sh            \
        test_interp_method_rbf_parallel.sh             \
        test_interp_method_file_parallel.sh            \
        test_interp_method_hcsbb_parallel.sh           \
        test_interp_method_spmap_parallel.sh           \
        test_interp_weights_parallel.sh                \
        test_mpi_handshake.sh                          \
        test_mpi_handshake_c.sh                        \
        test_dummy_coupling.sh                         \
        test_dummy_coupling2_c.sh                      \
        test_dummy_coupling3.sh                        \
        test_dummy_coupling3_c.sh                      \
        test_dummy_coupling4_c.sh                      \
        test_dummy_coupling5.sh                        \
        test_dummy_coupling6.sh                        \
        test_dummy_coupling6_c.sh                      \
        test_dummy_coupling7.sh                        \
        test_dummy_coupling7_c.sh                      \
        test_dummy_coupling8_c.sh                      \
        test_dummy_coupling9.sh                        \
        test_interpolation_exchange.sh                 \
        test_interpolation_parallel1.sh                \
        test_interpolation_parallel2.sh                \
        test_interpolation_parallel3.sh                \
        test_interpolation_parallel4.sh                \
        test_interpolation_parallel5.sh                \
        test_interpolation_parallel6.sh                \
        test_init_final.sh                             \
        test_init_comm_final.sh                        \
        test_proc_sphere_part_parallel.sh              \
        test_read_fesom.sh                             \
        test_read_icon.sh                              \
        test_read_icon_parallel.sh                     \
        test_read_mpiom.sh                             \
        test_read_scrip.sh                             \
        test_read_woa_data.sh                          \
        test_redirstdout.sh                            \
        test_restart.sh                                \
        test_restart2.sh                               \
        test_version.sh                                \
        test_weights2vtk.sh                            \
        test_dynamic_config.sh                         \
        test_query_routines.sh                         \
        test_multithreading.sh                         \
        test_io_config.sh

xfail_test_SCRIPTS_ =                                  \
        test_mpi_error.sh

# Build programs associated with test scripts:
check_PROGRAMS +=                      \
        test_abort_c.x                 \
        test_abort.x                   \
        test_def_comps.x               \
        test_def_datetime.x            \
        test_def_grid.x                \
        test_def_mask.x                \
        test_def_points.x              \
        test_init_final.x              \
        test_init_comm_final.x         \
        test_version.x

if TEST_NETCDF
check_PROGRAMS +=                        \
        test_read_fesom.x                \
        test_read_icon.x                 \
        test_read_scrip.x                \
        test_read_mpiom.x                \
        test_read_woa_data.x             \
        test_weights2vtk.x
endif

if TEST_MPI
check_PROGRAMS +=                        \
        test_component_config.x          \
        test_couple_config.x             \
        test_group_comm.x                \
        test_mpi_handshake_c.x           \
        test_interpolation_exchange.x    \
        test_interpolation_parallel1.x   \
        test_interpolation_parallel2.x   \
        test_interpolation_parallel6.x   \
        test_instance_parallel2.x        \
        test_instance_parallel3.x        \
        test_instance_parallel4.x        \
        test_dummy_coupling_dble.x       \
        test_dummy_coupling_real.x       \
        test_dummy_coupling_c.x          \
        test_dummy_coupling4_c.x         \
        test_dummy_coupling5_dble.x      \
        test_dummy_coupling5_real.x      \
        test_dummy_coupling5_c.x         \
        test_dummy_coupling6_c.x         \
        test_dummy_coupling7_c.x         \
        test_dummy_coupling8_c.x         \
        test_dummy_coupling9.x           \
        test_dummy_coupling9_c.x         \
        test_mpi_error.x                 \
        test_proc_sphere_part_parallel.x \
        test_redirstdout.x               \
        test_redirstdout_c.x             \
        test_restart2.x                  \
        test_io_config.x

if TEST_PTHREAD
check_PROGRAMS += \
        test_multithreading.x
endif

if TEST_YAXT_FC
check_PROGRAMS +=                   \
        test_dummy_coupling6_dble.x \
        test_dummy_coupling6_real.x \
        test_dummy_coupling7_dble.x \
        test_dummy_coupling7_real.x \
        test_mpi_handshake.x        \
        test_query_routines_c.x     \
        test_query_routines.x
endif

if TEST_NETCDF
check_PROGRAMS +=                             \
        test_instance_parallel1.x             \
        test_dummy_coupling2_c.x              \
        test_dummy_coupling3_c.x              \
        test_interpolation_parallel3.x        \
        test_interpolation_parallel4.x        \
        test_interpolation_parallel5.x        \
        test_read_icon_parallel.x             \
        test_restart_c.x                      \
        test_restart_dble.x                   \
        test_dist_grid_pair_parallel.x        \
        test_interp_grid_parallel.x           \
        test_interp_method_parallel.x         \
        test_interp_method_avg_parallel.x     \
        test_interp_method_nnn_parallel.x     \
        test_interp_method_nnn_parallel2.x    \
        test_interp_method_rbf_parallel.x     \
        test_interp_method_file_parallel.x    \
        test_interp_method_hcsbb_parallel.x   \
        test_interp_method_spmap_parallel.x   \
        test_interp_weights_parallel.x        \
        test_interp_method_conserv_parallel.x \
        test_interp_method_callback_parallel.x \
        test_interp_method_creep_parallel.x    \
        test_dynamic_config_c.x

if TEST_YAXT_FC
check_PROGRAMS +=                   \
        test_dummy_coupling3_dble.x \
        test_dummy_coupling3_real.x \
        test_dynamic_config.x
endif
endif
endif

TESTS =                         \
        $(test_PROGRAMS_)       \
        $(test_SCRIPTS_)        \
        $(xfail_test_PROGRAMS_) \
        $(xfail_test_SCRIPTS_)

XFAIL_TESTS =                   \
        $(xfail_test_PROGRAMS_) \
        $(xfail_test_SCRIPTS_)

if !ENABLE_CONCURRENT_MPI_TESTS
# Run MPI tests ony by one:
test_def_comps.log: test_component_config.log
test_def_datetime.log: test_def_comps.log
test_def_grid.log: test_def_datetime.log
test_def_mask.log: test_def_grid.log
test_def_points.log: test_def_mask.log
test_dist_grid_pair_parallel.log: test_def_points.log
test_dummy_coupling.log: test_dist_grid_pair_parallel.log
test_dummy_coupling2_c.log: test_dummy_coupling.log
test_dummy_coupling3.log: test_dummy_coupling2_c.log
test_dummy_coupling3_c.log: test_dummy_coupling3.log
test_dummy_coupling4_c.log: test_dummy_coupling3_c.log
test_dummy_coupling5.log: test_dummy_coupling4_c.log
test_dummy_coupling6.log: test_dummy_coupling5.log
test_dummy_coupling6_c.log: test_dummy_coupling6.log
test_dummy_coupling7.log: test_dummy_coupling6_c.log
test_dummy_coupling7_c.log: test_dummy_coupling7.log
test_dummy_coupling8_c.log: test_dummy_coupling7_c.log
test_dummy_coupling9.log: test_dummy_coupling8_c.log
test_init_comm_final.log: test_dummy_coupling9.log
test_init_final.log: test_init_comm_final.log
test_mpi_handshake.log: test_init_final.log
test_mpi_handshake_c.log: test_mpi_handshake.log
test_instance_parallel1.log: test_mpi_handshake_c.log
test_instance_parallel2.log: test_instance_parallel1.log
test_instance_parallel3.log: test_instance_parallel2.log
test_instance_parallel4.log: test_instance_parallel3.log
test_interp_grid_parallel.log: test_instance_parallel4.log
test_interp_method_avg_parallel.log: test_interp_grid_parallel.log
test_interp_method_conserv_parallel.log: test_interp_method_avg_parallel.log
test_interp_method_creep_parallel.log: test_interp_method_conserv_parallel.log
test_interp_method_file_parallel.log: test_interp_method_creep_parallel.log
test_interp_method_hcsbb_parallel.log: test_interp_method_file_parallel.log
test_interp_method_nnn_parallel.log: test_interp_method_hcsbb_parallel.log
test_interp_method_nnn_parallel2.log: test_interp_method_nnn_parallel.log
test_interp_method_parallel.log: test_interp_method_nnn_parallel2.log
test_interp_method_rbf_parallel.log: test_interp_method_parallel.log
test_interp_method_spmap_parallel.log: test_interp_method_rbf_parallel.log
test_interp_weights_parallel.log: test_interp_method_spmap_parallel.log
test_interpolation_exchange.log: test_interp_weights_parallel.log
test_interpolation_parallel1.log: test_interpolation_exchange.log
test_interpolation_parallel2.log: test_interpolation_parallel1.log
test_interpolation_parallel3.log: test_interpolation_parallel2.log
test_interpolation_parallel4.log: test_interpolation_parallel3.log
test_interpolation_parallel5.log: test_interpolation_parallel4.log
test_interpolation_parallel6.log: test_interpolation_parallel5.log
test_proc_sphere_part_parallel.log: test_interpolation_parallel6.log
test_read_icon_parallel.log: test_proc_sphere_part_parallel.log
test_redirstdout.log: test_read_icon_parallel.log
test_restart.log: test_redirstdout.log
test_restart2.log: test_restart.log
test_dynamic_config.log:test_restart2.log
test_query_routines.log:test_dynamic_config.log
test_multithreading.log:test_query_routines.log
test_query_routines.log:test_io_config.log
endif

AM_CPPFLAGS = -DVERBOSE
AM_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/contrib \
        $(MPI_CFLAGS) \
        $(YAXT_CFLAGS) \
        $(NETCDF_CFLAGS) \
        $(FYAML_CFLAGS) \
        $(MTIME_CFLAGS)

AM_FCFLAGS = $(FCPPDEF)VERBOSE \
        $(FCMODINC)$(top_builddir)/src \
        $(FCPPINC)$(top_builddir)/src \
        $(MPI_FCFLAGS) \
        $(YAXT_FCFLAGS)

LDADD_COMMON = $(top_builddir)/src/libyac.a \
        $(LAPACK_CLIBS) \
        $(MTIME_CLIBS) \
        $(FYAML_CLIBS) \
        $(NETCDF_CLIBS)

LDADD = $(LDADD_COMMON) $(YAXT_CLIBS) $(MPI_CLIBS)
FCLDADD = $(LDADD_COMMON)
if TEST_YAXT_FC
FCLDADD += $(YAXT_FCLIBS)
else
FCLDADD += $(YAXT_CLIBS)
endif
FCLDADD += $(MPI_FCLIBS)

if FCMODUC
utest_mod = UTEST.$(FCMODEXT)
else
utest_mod = utest.$(FCMODEXT)
endif

EXTRA_DIST = utest.F90 test_macros.inc
utest_FCDEPS = $(utest_mod) test_macros.inc
FCLDADD += utest.$(OBJEXT)

$(utest_mod): utest.$(OBJEXT)
	@if test ! -f $@; then rm -f $<; $(MAKE) $<; fi

test_angle_x_SOURCES = test_angle.c tests.c

test_area_x_SOURCES = test_area.c tests.c test_common.c test_common.h

test_cell_bnd_circle_x_SOURCES = test_cell_bnd_circle.c tests.c test_common.c test_common.h

test_circle_x_SOURCES = test_circle.c tests.c test_common.c test_common.h

test_clipping_x_SOURCES = test_clipping.c tests.c test_common.c test_common.h

test_lat_clipping_x_SOURCES = test_lat_clipping.c tests.c test_common.c test_common.h

test_compute_overlap_area_x_SOURCES = test_compute_overlap_area.c tests.c test_common.c test_common.h

test_couple_config_x_SOURCES = test_couple_config.c tests.c

test_component_config_x_SOURCES = test_component_config.c tests.c test_common.c test_common.h

test_instance_parallel1_x_SOURCES = test_instance_parallel1.c weight_file_common.c weight_file_common.h  \
                                    tests.c test_common.c test_common.h

test_instance_parallel2_x_SOURCES = test_instance_parallel2.c tests.c \
                                    test_common.c test_common.h

test_instance_parallel3_x_SOURCES = test_instance_parallel3.c tests.c \
                                    test_common.c test_common.h

test_instance_parallel4_x_SOURCES = test_instance_parallel4.c tests.c \
                                    test_common.c test_common.h

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

test_dummy_coupling_c_x_SOURCES = test_dummy_coupling_c.c tests.c test_common.c test_common.h

test_dummy_coupling2_c_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_dummy_coupling2_c_x_SOURCES = test_dummy_coupling2_c.c tests.c test_common.c test_common.h

test_dummy_coupling3_c_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_dummy_coupling3_c_x_SOURCES = test_dummy_coupling3_c.c tests.c test_common.c test_common.h

test_dummy_coupling4_c_x_SOURCES = test_dummy_coupling4_c.c tests.c test_common.c test_common.h

test_dummy_coupling5_c_x_SOURCES = test_dummy_coupling5_c.c tests.c test_common.c test_common.h

test_dummy_coupling6_c_x_SOURCES = test_dummy_coupling6_c.c tests.c test_common.c test_common.h

test_dummy_coupling7_c_x_SOURCES = test_dummy_coupling7_c.c tests.c test_common.c test_common.h

test_dummy_coupling8_c_x_SOURCES = test_dummy_coupling8_c.c tests.c test_common.c test_common.h

test_dummy_coupling9_x_LDADD = $(FCLDADD)
test_dummy_coupling9_x_SOURCES = test_dummy_coupling9.F90
test_dummy_coupling9.$(OBJEXT): $(utest_FCDEPS)

test_dummy_coupling9_c_x_SOURCES = test_dummy_coupling9_c.c tests.c test_common.c test_common.h

test_dummy_coupling_dble_x_LDADD = $(FCLDADD)
test_dummy_coupling_dble_x_SOURCES = test_dummy_coupling_dble.F90 test_dummy_coupling.inc
test_dummy_coupling_dble.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling.inc

test_dummy_coupling_real_x_LDADD = $(FCLDADD)
test_dummy_coupling_real_x_SOURCES = test_dummy_coupling_real.F90 test_dummy_coupling.inc
test_dummy_coupling_real.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling.inc

test_dummy_coupling3_dble_x_LDADD = $(top_builddir)/contrib/libgridio.a $(FCLDADD)
test_dummy_coupling3_dble_x_SOURCES = test_dummy_coupling3_dble.F90 test_dummy_coupling3.inc
test_dummy_coupling3_dble.$(OBJEXT): $(utest_FCDEPS) test_dummy_coupling3.inc

test_dummy_coupling3_real_x_LDADD = $(top_builddir)/contrib/libgridio.a $(FCLDADD)
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

test_redirstdout_c_x_SOURCES = test_redirstdout_c.c tests.c test_common.c test_common.h

test_redirstdout_x_LDADD = $(FCLDADD)
test_redirstdout_x_SOURCES = test_redirstdout.F90
test_redirstdout.$(OBJEXT): $(utest_FCDEPS)

test_restart_c_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_restart_c_x_SOURCES = test_restart_c.c tests.c test_common.c test_common.h

test_restart_dble_x_LDADD = $(top_builddir)/contrib/libgridio.a $(FCLDADD)
test_restart_dble_x_SOURCES = test_restart_dble.F90
test_restart_dble.$(OBJEXT): $(utest_FCDEPS)

test_restart2_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_restart2_x_SOURCES = test_restart2.c tests.c test_common.c test_common.h

test_events_x_SOURCES = test_events.c tests.c tests.h

test_mpi_handshake_c_x_SOURCES = test_mpi_handshake_c.c tests.c

test_mpi_handshake_x_LDADD = $(FCLDADD)
test_mpi_handshake_x_SOURCES = test_mpi_handshake.F90
test_mpi_handshake.$(OBJEXT): $(utest_FCDEPS)

test_interpolation_exchange_x_SOURCES = test_interpolation_exchange.c tests.c test_common.c test_common.h

test_interpolation_parallel1_x_SOURCES = test_interpolation_parallel1.c tests.c test_common.c test_common.h

test_interpolation_parallel2_x_SOURCES = test_interpolation_parallel2.c tests.c test_common.c test_common.h

test_interpolation_parallel3_x_SOURCES = test_interpolation_parallel3.c weight_file_common.c weight_file_common.h tests.c test_common.c test_common.h

test_interpolation_parallel4_x_SOURCES = test_interpolation_parallel4.c weight_file_common.c weight_file_common.h tests.c

test_interpolation_parallel5_x_SOURCES = test_interpolation_parallel5.c weight_file_common.c weight_file_common.h tests.c test_common.c test_common.h

test_interpolation_parallel6_x_SOURCES = test_interpolation_parallel6.c tests.c test_common.c test_common.h

test_pxgc_x_SOURCES = test_pxgc.c test_cxc.c tests.c test_cxc.h

test_gcxgc_x_SOURCES = test_gcxgc.c test_cxc.c tests.c test_cxc.h

test_gcxlatc_x_SOURCES = test_gcxlatc.c test_cxc.c tests.c test_cxc.h

test_generate_cubed_sphere_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_generate_cubed_sphere_x_SOURCES = test_generate_cubed_sphere.c tests.c

test_geometry_x_SOURCES = test_geometry.c tests.c

test_grid_x_SOURCES = test_grid.c tests.c test_common.c test_common.h

test_grid2vtk_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_grid2vtk_x_SOURCES = test_grid2vtk.c tests.c test_common.c test_common.h

test_init_final_x_LDADD = $(FCLDADD)
test_init_final_x_SOURCES = test_init_final.F90
test_init_final.$(OBJEXT): $(utest_FCDEPS)

test_init_comm_final_x_LDADD = $(FCLDADD)
test_init_comm_final_x_SOURCES = test_init_comm_final.F90
test_init_comm_final.$(OBJEXT): $(utest_FCDEPS)

test_interval_tree_x_SOURCES = test_interval_tree.c tests.c

test_latcxlatc_x_SOURCES = test_latcxlatc.c test_cxc.c tests.c test_cxc.h

test_loncxlatc_x_SOURCES = test_loncxlatc.c test_cxc.c tests.c test_cxc.h

test_loncxlonc_x_SOURCES = test_loncxlonc.c test_cxc.c tests.c test_cxc.h

test_mergesort_x_SOURCES = test_mergesort.c tests.c

test_partial_areas_x_SOURCES = test_partial_areas.c tests.c test_common.c test_common.h

test_point_in_cell_x_SOURCES = test_point_in_cell.c tests.c test_common.c test_common.h

test_bnd_sphere_part_x_SOURCES = test_bnd_sphere_part.c tests.c test_common.c test_common.h

test_point_sphere_part_x_SOURCES = test_point_sphere_part.c tests.c test_common.c test_common.h

test_proc_sphere_part_parallel_x_SOURCES = test_proc_sphere_part_parallel.c tests.c test_common.c test_common.h

test_quicksort_x_SOURCES = test_quicksort.c tests.c

test_read_cube_csv_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_read_cube_csv_x_SOURCES = test_read_cube_csv.c tests.c

test_vtk_output_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_vtk_output_x_SOURCES = test_vtk_output.c tests.c

test_read_fesom_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_read_fesom_x_SOURCES = test_read_fesom.c tests.c

test_read_icon_parallel_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_read_icon_parallel_x_SOURCES = test_read_icon_parallel.c tests.c test_common.c test_common.h

test_read_icon_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_read_icon_x_SOURCES = test_read_icon.c tests.c test_common.c test_common.h

test_read_mpiom_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_read_mpiom_x_SOURCES = test_read_mpiom.c tests.c test_common.c test_common.h

test_read_scrip_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_read_scrip_x_SOURCES = test_read_scrip.c tests.c test_common.c test_common.h

test_read_woa_data_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_read_woa_data_x_SOURCES = test_read_woa_data.c tests.c test_common.c test_common.h

test_weights2vtk_x_SOURCES = test_weights2vtk.c tests.c test_common.c test_common.h weight_file_common.c weight_file_common.h

test_dist_grid_utils_x_SOURCES = test_dist_grid_utils.c tests.c

test_dist_grid_pair_parallel_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_dist_grid_pair_parallel_x_SOURCES = test_dist_grid_pair_parallel.c tests.c test_common.c test_common.h

test_interp_grid_parallel_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_interp_grid_parallel_x_SOURCES = test_interp_grid_parallel.c tests.c test_common.c test_common.h

test_interp_method_parallel_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_interp_method_parallel_x_SOURCES = test_interp_method_parallel.c tests.c test_common.c test_common.h

test_interp_method_avg_parallel_x_SOURCES = test_interp_method_avg_parallel.c tests.c test_common.c test_common.h

test_interp_method_conserv_parallel_x_SOURCES = test_interp_method_conserv_parallel.c tests.c test_common.c test_common.h

test_interp_method_creep_parallel_x_SOURCES = test_interp_method_creep_parallel.c tests.c test_common.c test_common.h

test_interp_method_nnn_parallel_x_SOURCES = test_interp_method_nnn_parallel.c tests.c test_common.c test_common.h

test_interp_method_nnn_parallel2_x_SOURCES = test_interp_method_nnn_parallel2.c tests.c test_common.c test_common.h

test_interp_method_rbf_parallel_x_SOURCES = test_interp_method_rbf_parallel.c tests.c test_common.c test_common.h

test_interp_method_callback_parallel_x_SOURCES = test_interp_method_callback_parallel.c weight_file_common.h weight_file_common.c tests.c test_common.c test_common.h

test_interp_method_file_parallel_x_SOURCES = test_interp_method_file_parallel.c weight_file_common.h weight_file_common.c tests.c test_common.c test_common.h

test_interp_method_hcsbb_parallel_x_SOURCES = test_interp_method_hcsbb_parallel.c tests.c test_common.c test_common.h

test_interp_method_spmap_parallel_x_SOURCES = test_interp_method_spmap_parallel.c tests.c test_common.c test_common.h

test_interp_weights_parallel_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
test_interp_weights_parallel_x_SOURCES = test_interp_weights_parallel.c tests.c test_common.c test_common.h

test_group_comm_x_SOURCES = test_group_comm.c tests.c test_common.c test_common.h

test_multithreading_x_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
test_multithreading_x_LDADD = $(PTHREAD_LIBS) $(LDADD)
test_multithreading_x_SOURCES = test_multithreading.c

test_version_x_LDADD = $(FCLDADD)
test_version_x_SOURCES = test_version.F90
test_version.$(OBJEXT): $(utest_FCDEPS)

test_abort_x_LDADD = $(FCLDADD)
test_abort_x_SOURCES = test_abort.F90

test_abort_c_x_SOURCES = test_abort_c.c  tests.c test_common.c test_common.h

test_dynamic_config_c_x_SOURCES = test_dynamic_config_c.c tests.c test_common.c test_common.h

test_dynamic_config_x_LDADD = $(FCLDADD)
test_dynamic_config_x_SOURCES = test_dynamic_config.F90
test_dynamic_config.$(OBJEXT): $(utest_FCDEPS)

test_query_routines_c_x_SOURCES = test_query_routines_c.c tests.c

test_query_routines_x_LDADD = $(FCLDADD)
test_query_routines_x_SOURCES = test_query_routines.F90
test_query_routines.$(OBJEXT): $(utest_FCDEPS)

test_interp_stack_config_x_SOURCES = test_interp_stack_config.c tests.c

test_io_config_x_SOURCES = test_io_config.c tests.c test_common.c test_common.h

AUTOMAKE_OPTIONS = color-tests

CLEANFILES = weights_additiona_b.nc

clean-local:
	-rm -rf *.0 *.1 *.2 *.err *.log *.vtk

mostlyclean-local:
	-rm -f *.$(FCMODEXT) i.*.L *.i
	-rm -rf *.x.dSYM

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
FCDEFS = $(FCPPDEF)HAVE_CONFIG_H
DEFAULT_FCINCLUDES = $(FCPPINC).@yac_fcisrc@
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_FCINCLUDES) $(FCINCLUDES) \
	$(AM_FCFLAGS) $(FCFLAGS)
