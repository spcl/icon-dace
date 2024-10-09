# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and
# DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
"""\
Object wrapper for libmtime Python bindings

This should eventually replace the original, very low-level bindings
"""

import ctypes

from . import _mtime as libmtime

CALENDAR_TYPE = libmtime.CALENDAR_TYPE


def setCalendar(calendar_type):
    """(Re-)initialize calendar used for date and time calculations"""
    if libmtime.getCalendarType().cType != CALENDAR_TYPE.calendar_not_set:
        libmtime.resetCalendar()
    libmtime.setCalendar(calendar_type)


def calendarToString():
    return libmtime.calendarToString().decode()


class MTime(object):
    """\
    Base class for all other data classes

    If the calendar has not been set before object creation,
    it defaults to Proleptic Gregorian
    """

    def __init__(self):
        if "_lib" not in self.__dict__:
            self._lib = libmtime
            if (
                libmtime.getCalendarType().cType
                == CALENDAR_TYPE.calendar_not_set
            ):
                setCalendar(CALENDAR_TYPE.proleptic_gregorian)


class Date(MTime):
    def __init__(self, spec):
        super(Date, self).__init__()
        self._lib = libmtime
        self._my = libmtime.newDate(str(spec).encode())
        if not self._my:
            raise ValueError("invalid date spec '{0}'".format(spec))

    def __del__(self):
        self._lib.deallocateDate(self._my)

    def __str__(self):
        return libmtime.dateToString(self._my).decode()


class Time(MTime):
    def __init__(self, spec):
        super(Time, self).__init__()
        self._lib = libmtime
        self._my = libmtime.newTime(str(spec).encode())
        if not self._my:
            raise ValueError("invalid time spec '{0}'".format(spec))

    def __del__(self):
        self._lib.deallocateTime(self._my)

    def __str__(self):
        return libmtime.timeToString(self._my).decode()


class DateTime(MTime):
    def __init__(self, spec):
        super(DateTime, self).__init__()
        self._my = libmtime.newDateTime(str(spec).encode())
        if not self._my:
            raise ValueError("invalid date time spec '{0}'".format(spec))

    def __del__(self):
        self._lib.deallocateDateTime(self._my)

    def __str__(self):
        return libmtime.datetimeToString(self._my).decode()

    def __eq__(self, other):
        return libmtime.compareDatetime(self._my, other._my) == 0

    def __ne__(self, other):
        return libmtime.compareDatetime(self._my, other._my) != 0

    def __lt__(self, other):
        return libmtime.compareDatetime(self._my, other._my) < 0

    def __le__(self, other):
        return libmtime.compareDatetime(self._my, other._my) <= 0

    def __gt__(self, other):
        return libmtime.compareDatetime(self._my, other._my) > 0

    def __ge__(self, other):
        return libmtime.compareDatetime(self._my, other._my) >= 0

    def __add__(self, delta):
        result = DateTime(self)
        result._my = libmtime.addTimeDeltaToDateTime(
            result._my, delta._my, result._my
        )
        return result

    def data(self):
        return self._my

    @property
    def date(self):
        return Date(
            libmtime.dateToString(ctypes.byref(self._my.contents.date)).decode()
        )

    @property
    def time(self):
        return Time(
            libmtime.timeToString(ctypes.byref(self._my.contents.time)).decode()
        )


class TimeDelta(MTime):
    def __init__(self, spec):
        super(TimeDelta, self).__init__()
        self._lib = libmtime
        self._my = libmtime.newTimeDelta(str(spec).encode())
        if not self._my:
            raise ValueError("invalid time delta spec '{0}'".format(spec))

    def __del__(self):
        self._lib.deallocateTimeDelta(self._my)

    def __str__(self):
        return libmtime.timedeltaToString(self._my).decode()

    def __rfloordiv__(self, datetimes):
        result = 0
        try:
            (datetime_a, datetime_b) = datetimes
            result_buffer = libmtime._divisionquotienttimespan()
            result = libmtime.divideDatetimeDifferenceInSeconds(
                datetime_a.data(),
                datetime_b.data(),
                self._my,
                ctypes.byref(result_buffer),
            )
        except Exception:
            return NotImplemented
        if not result:
            raise ValueError(
                "time delta '{0}' not valid for division ".format(self)
            )
        return result.contents.quotient

    def data(self):
        return self._my

    def items(self):
        for name, ctype in self._my.contents._fields_:
            yield (name, self._my.contents.__getattribute__(name))
