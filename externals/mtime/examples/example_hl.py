# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
from mtime import *


print('\n\nsetup calendar:')

setCalendar(CALENDAR_TYPE.proleptic_gregorian)
# setCalendar(CALENDAR_TYPE.year_of_365_days)
# setCalendar(CALENDAR_TYPE.year_of_360_days)


# Check the current calendar.

calendar_in_use = calendarToString()
print(calendar_in_use)


print('\n\ntest date interface:')

test_date = Date("2012-09-01")
print(test_date)


print('\n\ntest time interface:')

test_time = Time("12:13:49.654")
print(test_time)


print('\n\ntest datetime and timeddelta interface. Simulating a dummy timer.')

start_date = DateTime("2012-09-01T02:10:00.000")
print(start_date)

stop_date = DateTime("2012-09-10T14:00:00.000")
print(stop_date)

time_step = TimeDelta("PT12H")
print(time_step)

current_date = DateTime("2012-09-01T02:10:00.000")
print(current_date)

print('\n')
while current_date < stop_date:
    print('Model time loop  : ' + str(current_date))
    current_date += time_step
print('\n\n')
