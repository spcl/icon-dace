TEST_EXTENSIONS = @EXEEXT@ .run

# Serial tests:
TESTS = \
  calendar_test1 \
  cksum_verify \
  pio_write.run \
  pio_write_deco2d.run \
  test_byteswap.run \
  test_cdf_const.run \
  test_cdf_transformation.run \
  test_cksum_ext.run \
  test_cksum_grb.run \
  test_cksum_grb2.run \
  test_cksum_ieg.run \
  test_cksum_nc.run \
  test_cksum_nc2.run \
  test_cksum_nc4.run \
  test_cksum_nc_chunk.run \
  test_cksum_srv.run \
  test_f2003.run \
  test_grib.run \
  test_month_adjust \
  test_resource_copy.run \
  test_table.run

# Update test_cksum_base before running tests that source it:
test_cksum_TEST_LOGS_ = \
  test_cksum_ext.log \
  test_cksum_grb.log \
  test_cksum_grb2.log \
  test_cksum_ieg.log \
  test_cksum_nc.log \
  test_cksum_nc2.log \
  test_cksum_nc4.log \
  test_cksum_srv.log
$(test_cksum_TEST_LOGS_): test_cksum_base

# Parallel tests:
TESTS += \
  pio_cksum_asynch.run \
  pio_cksum_fpguard.run \
  pio_cksum_grb2.run \
  pio_cksum_mpi_fw_at_all.run \
  pio_cksum_mpi_fw_at_reblock.run \
  pio_cksum_mpi_fw_ordered.run \
  pio_cksum_mpinonb.run \
  pio_cksum_nc.run \
  pio_cksum_nc2.run \
  pio_cksum_nc4.run \
  pio_cksum_writer.run \
  pio_write.parallel.run \
  pio_write_deco2d.parallel.run \
  pio_write_dist_array.run \
  test_resource_copy.parallel.run

# Update pio_write_dist_array.run before running tests that source it:
pio_write_dist_array_TEST_LOGS_ = \
  pio_cksum_asynch.log \
  pio_cksum_fpguard.log \
  pio_cksum_grb2.log \
  pio_cksum_mpi_fw_at_all.log \
  pio_cksum_mpi_fw_at_reblock.log \
  pio_cksum_mpi_fw_ordered.log \
  pio_cksum_mpinonb.log \
  pio_cksum_nc.log \
  pio_cksum_nc2.log \
  pio_cksum_nc4.log \
  pio_cksum_writer.log
$(pio_write_dist_array_TEST_LOGS_): pio_write_dist_array.run

# Update pio_write.run before running tests that source it:
pio_write_TEST_LOGS_ = \
  $(pio_write_dist_array_TEST_LOGS_) \
  pio_write.parallel.log \
  pio_write_deco2d.log \
  pio_write_deco2d.parallel.log \
  pio_write_dist_array.log
$(pio_write_TEST_LOGS_): pio_write.run

if !with_concurrent_mpi_checks
pio_cksum_fpguard.log: pio_cksum_asynch.log
pio_cksum_grb2.log: pio_cksum_fpguard.log
pio_cksum_mpi_fw_at_all.log: pio_cksum_grb2.log
pio_cksum_mpi_fw_at_reblock.log: pio_cksum_mpi_fw_at_all.log
pio_cksum_mpi_fw_ordered.log: pio_cksum_mpi_fw_at_reblock.log
pio_cksum_mpinonb.log: pio_cksum_mpi_fw_ordered.log
pio_cksum_nc.log: pio_cksum_mpinonb.log
pio_cksum_nc2.log: pio_cksum_nc.log
pio_cksum_nc4.log: pio_cksum_nc2.log
pio_cksum_writer.log: pio_cksum_nc4.log
pio_write.parallel.log: pio_cksum_writer.log
pio_write_deco2d.parallel.log: pio_write.parallel.log
pio_write_dist_array.log: pio_write_deco2d.parallel.log
test_resource_copy.parallel.log: pio_write_dist_array.log
endif !with_concurrent_mpi_checks

test_PROGRAMS_ = \
  calendar_test1 \
  cksum_read \
  cksum_verify \
  cksum_write \
  cksum_write_chunk \
  pio_write \
  pio_write_deco2d \
  test_byteswap \
  test_grib \
  test_month_adjust \
  test_resource_copy \
  test_table

if ENABLE_NETCDF
test_PROGRAMS_ += \
  test_cdf_read \
  test_cdf_write
endif

if ENABLE_MPI
test_PROGRAMS_ += \
  pio_write.parallel \
  pio_write_deco2d.parallel \
  test_resource_copy.parallel
endif ENABLE_MPI

if with_on_demand_check_programs
check_PROGRAMS = $(test_PROGRAMS_)
else !with_on_demand_check_programs
noinst_PROGRAMS = $(test_PROGRAMS_)
endif !with_on_demand_check_programs

AM_CPPFLAGS = -I$(top_srcdir)/src
AM_CFLAGS = $(PPM_CORE_C_INCLUDE) $(YAXT_C_INCLUDE) $(MPI_C_INCLUDE)

AM_LDFLAGS =
if ENABLE_ALL_STATIC
AM_LDFLAGS += -all-static
endif

LDADD = $(top_builddir)/src/libcdi.la

calendar_test1_SOURCES = calendar_test1.c

cksum_read_SOURCES = \
  cksum_read.c \
  ensure_array_size.c \
  ensure_array_size.h \
  stream_cksum.c \
  stream_cksum.h \
  var_cksum.c \
  var_cksum.h

cksum_verify_SOURCES = cksum_verify.c

cksum_write_SOURCES = \
  cksum_write.c \
  simple_model_helper.c \
  simple_model_helper.h

cksum_write_chunk_SOURCES = \
  cksum_write_chunk.c \
  simple_model_helper.c \
  simple_model_helper.h

pio_write_SOURCES = \
  pio_write.c \
  pio_write.h \
  pio_write_setup_grid.c \
  pio_write_setup_grid.h \
  simple_model.c \
  simple_model_helper.c \
  simple_model_helper.h
nodist_pio_write_parallel_SOURCES = $(pio_write_SOURCES:.c=.parallel.c)
pio_write_parallel_LDADD = $(top_builddir)/src/libcdipio.la $(PPM_CORE_C_LIB) $(YAXT_C_LIB) $(LDADD)

pio_write_deco2d_SOURCES = \
  deco2d_model.c \
  pio_write.c \
  pio_write.h \
  pio_write_setup_grid.c \
  pio_write_setup_grid.h \
  simple_model_helper.c \
  simple_model_helper.h
nodist_pio_write_deco2d_parallel_SOURCES = $(pio_write_deco2d_SOURCES:.c=.parallel.c)
pio_write_deco2d_parallel_LDADD = $(top_builddir)/src/libcdipio.la $(PPM_CORE_C_LIB) $(YAXT_C_LIB) $(LDADD)

test_byteswap_SOURCES = test_byteswap.c

test_grib_SOURCES = test_grib.c

test_month_adjust_SOURCES = test_month_adjust.c

test_resource_copy_SOURCES = test_resource_copy.c
test_resource_copy_LDADD = $(top_builddir)/src/libcdiresunpack.la
nodist_test_resource_copy_parallel_SOURCES = $(test_resource_copy_SOURCES:.c=.parallel.c)
test_resource_copy_parallel_LDADD = $(top_builddir)/src/libcdipio.la $(LDADD)

test_table_SOURCES = test_table.c

test_cdf_read_SOURCES = test_cdf_read.c

test_cdf_write_SOURCES = \
  simple_model_helper.c \
  simple_model_helper.h \
  test_cdf_write.c

clean-local:
	-rm -f *.parallel.c
	-rm -f *.grb *.grb2 *.nc *.nc2 *.nc4 *.srv *.ext *.cksum *.ieg *.stderr *.stdout
	-rm -rf *.dSYM

AUTOMAKE_OPTIONS = color-tests

# Compile programs using Libtool:
COMPILE = $(LTCOMPILE)
FCCOMPILE = $(LTFCCOMPILE)
PPFCCOMPILE = $(LTPPFCCOMPILE)
CXXCOMPILE = $(LTCXXCOMPILE)

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
LTPPFCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=compile $(FC) $(AM_FCFLAGS) $(FCFLAGS)

# Avoid compiling twice by running Libtool with '-static'
#   This is safe because all Libtool libraries in this directory are
#   static-only and contain PIC objects:
AM_CFLAGS += -static
AM_FCFLAGS = -static
AM_CXXFLAGS = -static

# Generate source files for parallel versions of the programs. We do this to
# avoid target-specific flags, which result into rules that do not use the
# compilation commands above.
SUFFIXES = .parallel.c

.c.parallel.c:
	$(AM_V_GEN):;{ \
	  echo '#define USE_MPI 1'; \
	  echo '#include "$<"'; \
	  } >$@
