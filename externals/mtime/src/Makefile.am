# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
noinst_LTLIBRARIES = libmtime_c.la

lib_LTLIBRARIES = libmtime.la
libmtime_la_LDFLAGS = -version-number $(MTIME_VERSION_MAJOR):$(MTIME_VERSION_MINOR):$(MTIME_VERSION_PATCH)
libmtime_la_LIBADD = libmtime_c.la

AM_CPPFLAGS = -I$(top_srcdir)/include

libmtime_c_la_SOURCES = \
    kepler.c \
    mtime_calendar.c \
    mtime_calendar360day.c \
    mtime_calendar365day.c \
    mtime_calendarGregorian.c \
    mtime_date.c \
    mtime_datetime.c \
    mtime_eventHandling.c \
    mtime_eventList.c \
    mtime_iso8601.c \
    mtime_julianDay.c \
    mtime_time.c \
    mtime_timedelta.c \
    mtime_utilities.c \
    orbit.c \
    vsop87.c

libmtime_la_SOURCES = \
    libmtime.f90 \
    mtime_c_bindings.f90 \
    mtime_constants.f90 \
    mtime_error_handling.f90

if ENABLE_FORTRAN_HL
libmtime_la_SOURCES += \
    libmtime_hl.f90 \
    mtime_t_datetime.inc \
    mtime_t_event.inc \
    mtime_t_juliandelta.inc \
    mtime_t_timedelta.inc
endif ENABLE_FORTRAN_HL

if FCMODUC
mtime_c_bindings_mod = MTIME_C_BINDINGS.$(FCMODEXT)
mtime_constants_mod = MTIME_CONSTANTS.$(FCMODEXT)
mtime_error_handling_mod = MTIME_ERROR_HANDLING.$(FCMODEXT)
mtime_mod = MTIME.$(FCMODEXT)
mtime_other_mods =                \
    MTIME_CALENDAR.$(FCMODEXT)    \
    MTIME_DATE.$(FCMODEXT)        \
    MTIME_DATETIME.$(FCMODEXT)    \
    MTIME_EVENTGROUPS.$(FCMODEXT) \
    MTIME_EVENTS.$(FCMODEXT)      \
    MTIME_JULIANDAY.$(FCMODEXT)   \
    MTIME_JULIANDELTA.$(FCMODEXT) \
    MTIME_TIME.$(FCMODEXT)        \
    MTIME_TIMEDELTA.$(FCMODEXT)   \
    MTIME_UTILITIES.$(FCMODEXT)
mtime_print_by_callback_mod = MTIME_PRINT_BY_CALLBACK.$(FCMODEXT)
if ENABLE_FORTRAN_HL
mtime_hl_mod = MTIME_HL.$(FCMODEXT)
endif ENABLE_FORTRAN_HL
else !FCMODUC
mtime_c_bindings_mod = mtime_c_bindings.$(FCMODEXT)
mtime_constants_mod = mtime_constants.$(FCMODEXT)
mtime_error_handling_mod = mtime_error_handling.$(FCMODEXT)
mtime_mod = mtime.$(FCMODEXT)
mtime_other_mods =                \
    mtime_calendar.$(FCMODEXT)    \
    mtime_date.$(FCMODEXT)        \
    mtime_datetime.$(FCMODEXT)    \
    mtime_eventgroups.$(FCMODEXT) \
    mtime_events.$(FCMODEXT)      \
    mtime_julianday.$(FCMODEXT)   \
    mtime_juliandelta.$(FCMODEXT) \
    mtime_time.$(FCMODEXT)        \
    mtime_timedelta.$(FCMODEXT)   \
    mtime_utilities.$(FCMODEXT)
mtime_print_by_callback_mod = mtime_print_by_callback.$(FCMODEXT)
if ENABLE_FORTRAN_HL
mtime_hl_mod = mtime_hl.$(FCMODEXT)
endif ENABLE_FORTRAN_HL
endif !FCMODUC

all_mods =                         \
    $(mtime_c_bindings_mod)        \
    $(mtime_constants_mod)         \
    $(mtime_error_handling_mod)    \
    $(mtime_mod)                   \
    $(mtime_other_mods)            \
    $(mtime_print_by_callback_mod)

if ENABLE_FORTRAN_HL
all_mods += \
    $(mtime_hl_mod)
endif ENABLE_FORTRAN_HL

nodist_include_HEADERS = $(all_mods)

$(mtime_c_bindings_mod): mtime_c_bindings.lo
$(mtime_constants_mod): mtime_constants.lo
$(mtime_error_handling_mod) $(mtime_print_by_callback_mod): mtime_error_handling.lo
$(mtime_mod) $(mtime_other_mods): libmtime.lo

if ENABLE_FORTRAN_HL
$(mtime_hl_mod): libmtime_hl.lo
endif ENABLE_FORTRAN_HL

$(all_mods):
	@if test -z $<; then \
	  echo "Cannot find Fortran source file providing module '$(basename $(@F))'." >&2; \
	elif test ! -f $@; then \
	  rm -f $<; $(MAKE) $(AM_MAKEFLAGS) $<; \
	fi

libmtime.lo:                       \
    $(mtime_c_bindings_mod)        \
    $(mtime_constants_mod)         \
    $(mtime_error_handling_mod)

mtime_c_bindings.lo:            \
    $(mtime_constants_mod)      \
    $(mtime_error_handling_mod)

mtime_iso8601.c: @MAINTAINER_MODE_TRUE@ mtime_iso8601.rl
	$(AM_V_GEN)$(RAGEL) -o $@ mtime_iso8601.rl

if ENABLE_FORTRAN_HL
libmtime_hl.lo: \
    $(mtime_c_bindings_mod) \
    $(mtime_constants_mod) \
    $(mtime_error_handling_mod) \
    $(mtime_mod) \
    mtime_t_datetime.inc \
    mtime_t_event.inc \
    mtime_t_juliandelta.inc \
    mtime_t_timedelta.inc
endif ENABLE_FORTRAN_HL

MOSTLYCLEANFILES = $(all_mods)

EXTRA_DIST = mtime_iso8601.rl

# overrides to make regular compilations also use libtool
COMPILE = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(CC) -static $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)
FCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(FC) -static $(AM_FCFLAGS) $(FCFLAGS)
