# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

lib_LIBRARIES = libyac_mtime.a

libyac_mtime_a_SOURCES = \
    include/mtime_calendar.h \
    include/mtime_calendar360day.h \
    include/mtime_calendar365day.h \
    include/mtime_calendarGregorian.h \
    include/mtime_date.h \
    include/mtime_datetime.h \
    include/mtime_eventHandling.h \
    include/mtime_eventList.h \
    include/mtime_iso8601.h \
    include/mtime_julianDay.h \
    include/mtime_time.h \
    include/mtime_timedelta.h \
    include/mtime_utilities.h \
    src/mtime_calendar.c \
    src/mtime_calendar360day.c \
    src/mtime_calendar365day.c \
    src/mtime_calendarGregorian.c \
    src/mtime_date.c \
    src/mtime_datetime.c \
    src/mtime_eventHandling.c \
    src/mtime_eventList.c \
    src/mtime_iso8601.c \
    src/mtime_julianDay.c \
    src/mtime_time.c \
    src/mtime_timedelta.c

AM_CPPFLAGS = -I$(top_srcdir)/mtime/include
AM_CFLAGS = $(PIC_CFLAGS)
