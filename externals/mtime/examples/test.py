# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
from ctypes import *
from libmtime import *

####################################################################


# TEST for calendar.


name = create_string_buffer(" ",max_calendar_str_len);
setCalendar(CALENDAR_TYPE.proleptic_gregorian);
CalType = getCalendarType();
print CalType.cType;


calendarToString(name);
print name.value;
'''
        CalType = getCalendarType();
print CalType.cType;

print 'Hello World'

freeCalendar();
calendarToString(name);
print name.value;
'''

################################################



# Test for julian day.

jd = newJulianDay(100,100000);


print jd.contents.day;
print jd.contents.ms;

print '--------------------'
name = create_string_buffer(" ",max_julianday_str_len);
val = juliandayToString(jd,name);
print name.value;
print val;

deallocateJulianDay(jd);

print jd.contents.day;
print jd.contents.ms;


#####################################################################            

# Test for date.


d = newDate("2099-07-01");
name = create_string_buffer(" ",max_date_str_len);
print 'HELLLLLLLOO'
dateToString(d,name);
print name.value;
'''
dd = newRawDate(1990,11,11);
name_ = create_string_buffer(" ",100);
dateToString(dd,name_);
print name_.value;



ddd = constructAndCopyDate(d);
name__ = create_string_buffer(" ",100);
dateToString(ddd,name__);
print name__.value;

c = compareDate(d,dd);
print c;

dd = replaceDate(d,dd);
name____ = create_string_buffer(" ",100);
dateToString(dd,name____);
print name____.value;

deallocateDate(ddd);
name___ = create_string_buffer(" ",100);
dateToString(ddd,name___);
print name___.value;
'''


#####################################################################

# Test for Time.



'''
t = newTime("01:01:01.009");
name_ = create_string_buffer(" ",100);
timeToString(t,name_);
print name_.value;

tt = newRawTime(12,12,12,999);
name = create_string_buffer(" ",100);
timeToString(tt,name);
print name.value;

ttt = constructAndCopyTime(tt);
name__ = create_string_buffer(" ",100);
timeToString(ttt,name__);
print name__.value;

deallocateTime(ttt);
name___ = create_string_buffer(" ",100);
timeToString(ttt,name___);
print name___.value;


t = replaceTime(tt,t);
name____ = create_string_buffer(" ",100);
timeToString(t,name____);
print name____.value;
'''
#####################


#Test for DateTime




dt = newdatetimefromstring("0-01-01T12:00:00.000");
name = create_string_buffer(" ",100);
datetimeToString(dt,name);
print name.value;

dtt = newdatetimefromraw(2020,12,12,10,01,01,07);
name_ = create_string_buffer(" ",100);
datetimeToString(dtt,name_);
print name_.value;
print 'qqqqqqqqqqqqq'
print datetime_gt(dtt,dtt);


dttt = newdatetimefromconstructandcopy(dt);
name__ = create_string_buffer(" ",100);
datetimeToString(dttt,name__);
print name__.value;

deallocateDateTime(dttt);
name___ = create_string_buffer(" ",100);
datetimeToString(dttt,name___);
print name___.value;

c = compareDatetime(dtt,dtt);
print c;

cc = getNoOfDaysInMonthDateTime(dtt);
print cc;


cc = getNoOfDaysInYearDateTime(dtt);
print cc;

cc = getDayOfYearFromDateTime(dtt);


cc = getDayOfYearFromDateTime(dtt);
print cc;

jd = newJulianDay(0,0);
jd = getJulianDayFromDateTime(dt,jd);
name_____ = create_string_buffer(" ",100);
juliandayToString(jd,name_____);
print name_____.value;

print dt == dt;


#replaceDatetime(dtt,dt);
dt = dtt;
name____ = create_string_buffer(" ",100);
datetimeToString(dt,name____);
print name____.value;




####################################################################

#Test timedelta.



'''
td = newTimeDelta("PT10S");
name = create_string_buffer(" ",100);
timedeltaToString(td,name);
print name.value;


dt = newDateTime("0-01-01T12:00:00.000");
dt_ = newDateTime("0-01-01T00:00:00.000");
dt_ = addTimeDeltaToDateTime(dt,td,dt_);
name_ = create_string_buffer(" ",100);
datetimeToString(dt_,name_);
print name_.value;

deallocateTimeDelta(td);
name__ = create_string_buffer(" ",100);
timedeltaToString(td,name__);
print name__.value;
'''


####################################################################


# Test Event.

e = newEvent("eventCleanup1", "2000-01-01T00:00:00", "2010-01-01T00:00:00", "2013-01-01T00:00:00", "P01M");
name = create_string_buffer(" ",100);
eventToString(e,name);
print name.value;
'''
deallocateEvent(e);
name_ = create_string_buffer(" ",100);
eventToString(e,name_);
print name_.value;
'''


#####################################################################

#Test  eventgroup.


'''
eg = newEventGroup("CleanupEventGroup");
c = addNewEventToEventGroup(e,eg);
print c;
cc = removeEventFromEventGroup("eventCleanup1", eg);
print cc;

deallocateEventGroup(eg);
'''


resetCalendar();


