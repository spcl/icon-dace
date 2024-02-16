noinst_PROGRAMS =              \
        OASIS3_MCT_example/OASIS3_MCT_example.x   \
        toy_dummy/dummy_atmosphere.x     \
        toy_dummy/dummy_atmosphere_c.x   \
        toy_dummy/dummy_io.x             \
        toy_dummy/dummy_io_c.x           \
        toy_dummy/dummy_ocean.x          \
        toy_dummy/dummy_ocean_c.x        \
        dynamic_config/dynamic_config.x \
        perf_toy/perf_toy_cube.x \
        perf_toy/perf_toy_icon.x \
        toy_coupling/toy_coupling.x \
        toy_callback/toy_icon_callback.x    \
        toy_callback/toy_cube_callback.x    \
        toy_icon_runoff/toy_icon_runoff_atm.x  \
        toy_icon_runoff/toy_icon_runoff_ocn.x  \
        toy_multi/toy_multi_cube.x       \
        toy_multi/toy_multi_reg2d.x      \
        toy_multi/toy_multi_unstruct.x   \
        toy_multi/toy_multi_unstruct_2.x \
        interpolation_benchmark/toy_scrip.x \
        toy_atm_ocn/toy_icon_atm.x         \
        toy_atm_ocn/toy_icon_ocn.x         \
        toy_atm_ocn/toy_reg2d_atm.x        \
        toy_atm_ocn/toy_reg2d_ocn.x        \
        toy_atm_ocn/toy_mpiom_ocn.x        \
        toy_atm_ocn/toy_woa_ocn.x

if HAVE_NETCDF
noinst_PROGRAMS +=                 \
        interpolation_benchmark/generate_OASIS_mask_file.x          \
        toy_output/toy_output.x
endif

AM_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/contrib \
        $(MPI_CFLAGS) \
        $(YAXT_CFLAGS) \
        $(NETCDF_CFLAGS) \
        $(FYAML_CFLAGS) \
        $(MTIME_CFLAGS)

AM_FCFLAGS = $(FCMODINC)$(top_builddir)/src $(FCPPINC)$(top_builddir)/src $(MPI_FCFLAGS)

LDADD_COMMON = $(top_builddir)/src/libyac.a \
        $(LAPACK_CLIBS) \
        $(MTIME_CLIBS) \
        $(FYAML_CLIBS) \
        $(NETCDF_CLIBS) \
        $(YAXT_CLIBS)

LDADD = $(LDADD_COMMON) $(MPI_CLIBS)
FCLDADD = $(LDADD_COMMON) $(MPI_FCLIBS)

AM_DEFAULT_SOURCE_EXT = .c
EXTRA_DIST = toy_atm_ocn/toy_common.h

OASIS3_MCT_example_OASIS3_MCT_example_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)

toy_dummy_dummy_atmosphere_x_LDADD = $(FCLDADD)
toy_dummy_dummy_atmosphere_x_SOURCES = toy_dummy/dummy_atmosphere.f90
toy_dummy_dummy_io_x_LDADD = $(FCLDADD)
toy_dummy_dummy_io_x_SOURCES = toy_dummy/dummy_io.f90
toy_dummy_dummy_ocean_x_LDADD = $(FCLDADD)
toy_dummy_dummy_ocean_x_SOURCES = toy_dummy/dummy_ocean.f90

toy_callback_toy_cube_callback_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
toy_callback_toy_icon_callback_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)

perf_toy_perf_toy_cube_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
perf_toy_perf_toy_icon_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)

toy_coupling_toy_coupling_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)

toy_atm_ocn_toy_icon_atm_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
toy_atm_ocn_toy_icon_ocn_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
toy_atm_ocn_toy_reg2d_atm_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
toy_atm_ocn_toy_reg2d_ocn_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
toy_atm_ocn_toy_mpiom_ocn_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
toy_atm_ocn_toy_woa_ocn_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)

toy_icon_runoff_toy_icon_runoff_atm_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
toy_icon_runoff_toy_icon_runoff_ocn_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)

toy_multi_toy_multi_cube_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
toy_multi_toy_multi_cube_x_SOURCES = toy_multi/toy_multi_cube.c toy_multi/toy_multi_common.c toy_multi/toy_multi_common.h
toy_multi_toy_multi_reg2d_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
toy_multi_toy_multi_reg2d_x_SOURCES = toy_multi/toy_multi_reg2d.c toy_multi/toy_multi_common.c toy_multi/toy_multi_common.h
toy_multi_toy_multi_unstruct_2_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
toy_multi_toy_multi_unstruct_2_x_SOURCES = toy_multi/toy_multi_unstruct_2.c toy_multi/toy_multi_common.c toy_multi/toy_multi_common.h
toy_multi_toy_multi_unstruct_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)
toy_multi_toy_multi_unstruct_x_SOURCES = toy_multi/toy_multi_unstruct.c toy_multi/toy_multi_common.c toy_multi/toy_multi_common.h

interpolation_benchmark_toy_scrip_x_LDADD = $(top_builddir)/contrib/libgridio.a $(LDADD)

toy_output_toy_output_x_SOURCES = toy_output/toy_output.c

clean-local:
	-rm -rf *.x.dSYM *.vtk *.0 *.i
