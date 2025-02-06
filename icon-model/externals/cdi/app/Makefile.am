bin_PROGRAMS    =
noinst_PROGRAMS =
check_PROGRAMS  =

if ENABLE_CDI_APP
if ENABLE_CDI_LIB
bin_PROGRAMS += cdi
else !ENABLE_CDI_LIB
noinst_PROGRAMS += cdi
endif !ENABLE_CDI_LIB
noinst_PROGRAMS += createtable
else !ENABLE_CDI_APP
if with_on_demand_check_programs
check_PROGRAMS += cdi
else !with_on_demand_check_programs
noinst_PROGRAMS += cdi
endif !with_on_demand_check_programs
endif !ENABLE_CDI_APP

AM_CPPFLAGS = -I$(top_srcdir)/src

AM_LDFLAGS =
if ENABLE_ALL_STATIC
AM_LDFLAGS += -all-static
endif

LDADD = $(top_builddir)/src/libcdi.la

cdi_SOURCES = \
  cdi.c \
  printinfo.c \
  printinfo.h

createtable_SOURCES = createtable.c

clean-local:
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
AM_FCFLAGS = -static
AM_CXXFLAGS = -static
