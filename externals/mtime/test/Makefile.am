# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
EXTRA_LTLIBRARIES = libmtime_check.la
check_PROGRAMS = test_runner
TESTS = test_runner

AM_CPPFLAGS = -I$(top_srcdir)/include

if ENABLE_CHECK
AM_CFLAGS = $(CHECK_CFLAGS)
test_runner_LDADD = libmtime_check.la $(top_builddir)/src/libmtime.la $(CHECK_CLIBS)
endif ENABLE_CHECK

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

clean-local:
	rm -rf *.dSYM

MOSTLYCLEANFILES = $(EXTRA_LTLIBRARIES)

