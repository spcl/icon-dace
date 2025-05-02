SUBDIRS = pio

noinst_PROGRAMS =
check_PROGRAMS =
CLEANFILES =

EXTRA_PROGRAMS = \
  cdi_append \
  cdi_copy \
  cdi_read \
  cdi_write \
  cdi_write_ens \
  cdi_write_hybrid \
  cdi_write_local \
  cdi_write_relativ

# Examples that are part of the test suite:
test_PROGRAMS_ =

if ENABLE_NETCDF
test_PROGRAMS_ += cdi_write_const
endif

if ENABLE_ISOC_INTERFACE
test_PROGRAMS_ += cdi_read_f2003 cdi_write_f2003
endif ENABLE_ISOC_INTERFACE

if with_on_demand_check_programs
check_PROGRAMS += $(test_PROGRAMS_)
else !with_on_demand_check_programs
noinst_PROGRAMS += $(test_PROGRAMS_)
endif !with_on_demand_check_programs

if with_example_programs
noinst_PROGRAMS += $(EXTRA_PROGRAMS)
else !with_example_programs
CLEANFILES += $(EXTRA_PROGRAMS)
endif !with_example_programs
examples-local: $(test_PROGRAMS_) $(EXTRA_PROGRAMS)

AM_CPPFLAGS = -I$(top_srcdir)/src

AM_FCFLAGS =
if ENABLE_ISOC_INTERFACE
AM_FCFLAGS += $(FC_MOD_FLAG)$(top_builddir)/src
endif ENABLE_ISOC_INTERFACE

AM_LDFLAGS =
if ENABLE_ALL_STATIC
AM_LDFLAGS += -all-static
endif

LDADD = $(top_builddir)/src/libcdi.la

cdi_append_SOURCES = cdi_append.c

cdi_copy_SOURCES = cdi_copy.c

cdi_read_SOURCES = cdi_read.c

cdi_write_SOURCES = cdi_write.c

cdi_write_ens_SOURCES = cdi_write_ens.c

cdi_write_hybrid_SOURCES = cdi_write_hybrid.c

cdi_write_local_SOURCES = cdi_write_local.c

cdi_write_relativ_SOURCES = cdi_write_relativ.c

cdi_write_const_SOURCES = cdi_write_const.c

cdi_read_f2003_SOURCES = cdi_read_f2003.f90
cdi_read_f2003_LDADD = $(top_builddir)/src/libcdi_f2003.la $(LDADD)

cdi_write_f2003_SOURCES = cdi_write_f2003.f90
cdi_write_f2003_LDADD = $(top_builddir)/src/libcdi_f2003.la $(LDADD)

if FC_MOD_UPPERCASE
mo_cdi_mod = MO_CDI.$(FCMODEXT)
else !FC_MOD_UPPERCASE
mo_cdi_mod = mo_cdi.$(FCMODEXT)
endif !FC_MOD_UPPERCASE

cdi_read_f2003.$(OBJEXT) cdi_write_f2003.$(OBJEXT): $(top_builddir)/src/$(mo_cdi_mod)

clean-local:
	-rm -f *.grb *.grb2 *.nc
	-rm -rf *.dSYM

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
AM_CFLAGS = -static
AM_FCFLAGS += -static
AM_CXXFLAGS = -static
