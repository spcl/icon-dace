check_LIBRARIES = libgridio.a

if !ENABLE_LIB_ONLY

noinst_LIBRARIES = libgridio.a
if HAVE_NETCDF
noinst_PROGRAMS = weights2vtk.x
endif

endif

libgridio_a_SOURCES =           \
        generate_cubed_sphere.c \
        generate_cubed_sphere.h \
        generate_reg2d.c        \
        generate_reg2d.h        \
        grid2vtk.c              \
        grid2vtk.h              \
        read_cube_csv_grid.c    \
        read_cube_csv_grid.h    \
        test_function.c         \
        test_function.h         \
        vtk_output.c            \
        vtk_output.h

if HAVE_NETCDF
libgridio_a_SOURCES +=          \
        read_grid.c             \
        read_grid.h             \
        read_fesom_grid.c       \
        read_fesom_grid.h       \
        read_icon_grid.c        \
        read_icon_grid.h        \
        read_mpiom_grid.c       \
        read_mpiom_grid.h       \
        read_scrip_grid.c       \
        read_scrip_grid.h       \
        read_woa_data.c         \
        read_woa_data.h
endif

AM_CFLAGS = \
	-I$(top_srcdir)/src \
	$(MPI_CFLAGS) \
	$(YAXT_CFLAGS) \
	$(NETCDF_CFLAGS)
weights2vtk_x_LDADD = \
	libgridio.a \
	$(top_builddir)/src/libyac.a \
	$(NETCDF_CLIBS) \
	$(MPI_CLIBS) \
	$(YAXT_CLIBS)

clean-local:
	-rm -rf *.x.dSYM
