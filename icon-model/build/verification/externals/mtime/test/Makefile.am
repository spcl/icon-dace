# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
check_PROGRAMS = test_runner
TESTS = test_runner test_python.sh
TEST_EXTENSIONS = .sh

if ENABLE_CHECK
check_LTLIBRARIES = libmtime_check.la

AM_CPPFLAGS = -I$(top_srcdir)/include
AM_CFLAGS = $(CHECK_CFLAGS)
test_runner_LDADD = libmtime_check.la $(top_builddir)/src/libmtime.la $(CHECK_CLIBS)

libmtime_check_la_SOURCES = \
    mtime_calendar_test.c   \
    mtime_calendar_test.h   \
    mtime_date_test.c       \
    mtime_date_test.h       \
    mtime_datetime_test.c   \
    mtime_datetime_test.h   \
    mtime_julianDay_test.c  \
    mtime_julianDay_test.h  \
    mtime_time_test.c       \
    mtime_time_test.h       \
    mtime_timedelta_test.c  \
    mtime_timedelta_test.h
endif ENABLE_CHECK

# overrides to make regular compilations also use libtool
COMPILE = $(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(CC) -static $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)
FCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(FC) -static $(AM_FCFLAGS) $(FCFLAGS)
