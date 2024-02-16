# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
lib_LTLIBRARIES = libmtime.la

AM_CPPFLAGS = -I$(top_srcdir)/include
AM_FCFLAGS =

if ENABLE_OPENMP
# We cannot simply append $(OPENMP_FCFLAG) to AM_FCFLAGS because this breaks
# linking with NAG compiler:
#   a) whenever libtool is called for linking, it duplicates the openmp flag,
#      which is considered as an error by NAG;
#   b) libmtime.la gets an entry inherited_linker_flags=' -openmp' and it
#      becomes impossible to link a C application to libmtime with libtool if
#      the C compiler is gcc, which misinterpretes the flag '-openmp'.
# Therefore, we prepend $(OPENMP_FCFLAG) with '-Xcompiler', which solves both of
# the aforementioned problems: the openmp flag is used when linking the shared
# version of libmtime but is not duplicated, additionaly, it is not saved to
# libmtime.la as an inherited linker flag.
# However, due to a bug in the argument parsing procedure, libtool fails if the
# first argument for the compiler is '-Xcompiler <flag>' or '-Wc,<flag>'.
# Therefore, we additionaly prepend AM_FCFLAGS with a dummy flag '$(FCMODINC).'.
AM_FCFLAGS += $(FCMODINC). -Xcompiler $(OPENMP_FCFLAG)
endif ENABLE_OPENMP

libmtime_la_SOURCES =         \
    kepler.c                  \
    libmtime.f90              \
    libmtime_hl.f90           \
    mtime_c_bindings.f90      \
    mtime_calendar.c          \
    mtime_calendar360day.c    \
    mtime_calendar365day.c    \
    mtime_calendarGregorian.c \
    mtime_constants.f90       \
    mtime_date.c              \
    mtime_datetime.c          \
    mtime_error_handling.f90  \
    mtime_eventHandling.c     \
    mtime_eventList.c         \
    mtime_iso8601.c           \
    mtime_julianDay.c         \
    mtime_t_datetime.inc      \
    mtime_t_event.inc         \
    mtime_t_juliandelta.inc   \
    mtime_t_timedelta.inc     \
    mtime_time.c              \
    mtime_timedelta.c         \
    mtime_utilities.c         \
    orbit.c                   \
    vsop87.c

if FCMODUC
mtime_c_bindings_mod = MTIME_C_BINDINGS.$(FCMODEXT)
mtime_constants_mod = MTIME_CONSTANTS.$(FCMODEXT)
mtime_error_handling_mod = MTIME_ERROR_HANDLING.$(FCMODEXT)
mtime_hl_mod = MTIME_HL.$(FCMODEXT)
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
else !FCMODUC
mtime_c_bindings_mod = mtime_c_bindings.$(FCMODEXT)
mtime_constants_mod = mtime_constants.$(FCMODEXT)
mtime_error_handling_mod = mtime_error_handling.$(FCMODEXT)
mtime_hl_mod = mtime_hl.$(FCMODEXT)
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
endif !FCMODUC

all_mods =                         \
    $(mtime_c_bindings_mod)        \
    $(mtime_constants_mod)         \
    $(mtime_error_handling_mod)    \
    $(mtime_hl_mod)                \
    $(mtime_mod)                   \
    $(mtime_other_mods)            \
    $(mtime_print_by_callback_mod)

nodist_include_HEADERS = $(all_mods)

$(mtime_c_bindings_mod): mtime_c_bindings.lo
$(mtime_constants_mod): mtime_constants.lo
$(mtime_error_handling_mod) $(mtime_print_by_callback_mod): mtime_error_handling.lo
$(mtime_hl_mod): libmtime_hl.lo
$(mtime_mod) $(mtime_other_mods): libmtime.lo

$(all_mods):
	@if test -z $<; then \
	  echo "Cannot find Fortran source file providing module '$(basename $(@F))'." >&2; \
	elif test ! -f $@; then \
	  rm -f $<; $(MAKE) $(AM_MAKEFLAGS) $<; \
	fi

libmtime.lo:                       \
    $(mtime_c_bindings_mod)        \
    $(mtime_constants_mod)         \
    $(mtime_error_handling_mod)    \
    $(mtime_print_by_callback_mod)

libmtime_hl.lo:                 \
    $(mtime_c_bindings_mod)     \
    $(mtime_constants_mod)      \
    $(mtime_error_handling_mod) \
    $(mtime_mod)                \
    mtime_t_datetime.inc        \
    mtime_t_event.inc           \
    mtime_t_juliandelta.inc     \
    mtime_t_timedelta.inc

mtime_c_bindings.lo:            \
    $(mtime_constants_mod)      \
    $(mtime_error_handling_mod)

mtime_error_handling.lo: $(mtime_constants_mod)

mtime_iso8601.c: @MAINTAINER_MODE_TRUE@ mtime_iso8601.rl
	$(AM_V_GEN)$(RAGEL) -o $@ mtime_iso8601.rl

MOSTLYCLEANFILES = $(all_mods)

mostlyclean-local:
	rm -f *.i

install-data-hook:
	sed -i "s,SOFILE = '.*',SOFILE = '$(libdir)/libmtime.so'," $(DESTDIR)$(pythondir)/libmtime.py

python_PYTHON = libmtime.py mtime.py

EXTRA_DIST = mtime_iso8601.rl

