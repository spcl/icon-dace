# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
## @file libmtime.py
#
# @brief Providing the Python language bindings for libmtime
#
# @author  Luis Kornblueh, Max Planck Institute for Meteorology
# @author  Rahul Sinha, Max Planck Institute for Meteorology
#
# @defgroup PythonBindings libmtime Python languange bindings
# @{
#
#___________________________________________________________________________________________________________


from ctypes import *

#CAREFUL where the python script is invoked from. This must be modified based on where the .so is located.
SOFILE = '.libs/libmtime.so'
#load lib
libmtime = cdll.LoadLibrary(SOFILE)




## @package mtime_calendar
#  @brief Provides the calendar to the users, abstracting the different calendars available.
#
# @details
#
# Three calendar types are provided:
#
#   - a proleptic Gregorian calendar
#   - a calendar with 365 days per year without leap years
#   - a calendar with 360 days per year and each month having 30 days
#
# To use a specific calendar a call to setCalendar with the respective selector must be
# done. The implementation is based on a singleton concept meaning that only one calendar
# can be active at a time. To release a calendar a call to resetCalendar has to be done.
#_____________________________________________________________________________________________


## @brief Enum CALENDAR_TYPE.
#
# USAGE example: setCalendar(CALENDAR_TYPE.proleptic_gregorian);
#
class CALENDAR_TYPE:
  #The calendar types supported.
  calendar_not_set 	= c_int(0).value; 
  proleptic_gregorian 	= c_int(1).value;
  year_of_365_days 	= c_int(2).value;
  year_of_360_days 	= c_int(3).value;

  #when an object of type CALENDAR_TYPE is created, obj.cType stores a value between 0 and 3 
  #denoting calendar_not_set, proleptic_gregorian, year_of_365_days, year_of_360_days.
  def __init__(self,val):
    if val in range(0,4):
      self.cType = val;

  #Allow CALENDAR_TYPE to be passed as an argument to setCalendar();
  @classmethod
  def from_param(cls, obj):
     if obj in range(0,4):
        return obj;
     else:
        return 0;
  

## Provides a string length for toString.
max_calendar_str_len = c_int(32).value;


## 
# @brief Initialize a new calendar.
#
# setCalendar is done at the very begining to select one of the
# provided calendar libraries. It intializes the calendar to one of:
#
# - proleptic_gregorian
# - year_of_365_days
# - year_of_360_days
# 
# The calendar type and hence it's behaviour (Calendar to Julian
# conversion and vice versa) is fixed for the lifetime of the
# selected calendar.  Attempts to change the calendar type on the
# fly is discouraged. The lib has built-in checks to reject
# change attempts at run time.  However, a calendar can be
# "re-initialized" after calling resetCalendar(), but this is not
# advised. 
#
# MANTRA: Know what you are doing before you do it and do it
# right the first time.
#
# USAGE example: setCalendar(CALENDAR_TYPE.proleptic_gregorian);
#
# @param[in]  ct      the calendar type
#		      ct shoud be of class CALENDAR_TYPE.
#
setCalendar                     = libmtime.initCalendar;
setCalendar.argtypes            = [CALENDAR_TYPE];
setCalendar.restype             = None;


## @brief called to discard the selected calendar type
#
# Usage example: resetCalendar();
resetCalendar                    = libmtime.freeCalendar;
resetCalendar.argtype            = [None];
resetCalendar.restype            = None;


## @brief to get an idea, which calendar type is selected
#
# USAGE example: ct = getCalendarType();
#		 ct is an object of class CALENDAR_TYPE.
#		 ct.cType can be
#		 		 - 0 	(calendar_not_set)
#				 - 1	(proleptic_gregorian)
#				 - 2    (year_of_365_days)
#				 - 3    (year_of_360_days)
#
#  @return an integer denoting the calendar in use
getCalendarType                 = libmtime.getCalendarType;
getCalendarType.argtype         = [None];
getCalendarType.restype         = CALENDAR_TYPE;


##
# @brief convert the calendar identifier into a human readable string
#
# USAGE example: out = create_string_buffer(b" ",max_calendar_str_len);
#		 calendarToString(out);
#		 (out will now contain the Calendar type in human readable form)
#
# @param[out]  string      the calendar type verbose
my_calendarToString                = libmtime.calendarToString;
my_calendarToString.argtype        = [c_char_p];
my_calendarToString.restype        = None;

def calendarToString():
    name = create_string_buffer(b" ",max_calendar_str_len);
    my_calendarToString(name);
    return name.value;


##################################################################################################################################


## @package mtime_julianday
#  @brief Julian Day Calendar and some operations supported on julian dates.
#
#  @details


## provides a string length for toString 
max_julianday_str_len = c_int(32).value;

class _julianday(Structure):
  _fields_ = [
		("day",c_int), 	##< the actual Julian day
        	("ms", c_long)	##< the milisecond on that particular day
	     ];


## @brief construct a new Julian date 
#
#  USAGE example: jd = newJulianDay(day,ms);
#		  jd is a pointer to an object of type _julianday.
#
# @param[in] day            the Julian day
# @param[in] ms             an integer denoting the actual milli seconds of a day 
# @return    ret_julianday  a pointer of type(julianday)
newJulianDay                    = libmtime.newJulianDay;
newJulianDay.argtypes           = [c_long,c_long];
newJulianDay.restype            = POINTER(_julianday);



## @brief destructor for a Julian date 
#
#  USAGE example: deallocateJulianDay(jd);
#
#  @param     my_julianday   a pointer of type(julianday)
deallocateJulianDay             = libmtime.deallocateJulianDay;
deallocateJulianDay.argtype     = [POINTER(_julianday)];
deallocateJulianDay.restype     = None;



## @brief get Julian day as a string.
#
#  USAGE example: 	str = create_string_buffer(b" ",max_julianday_str_len);
#			str = juliandayToString(jd,str);
#
#  @param[in]  jd   	a pointer to type(julianday). The Julian day to be converted to a string
#  @param[in]  str   	string where Julian day verbose is written
#  @param[out] str      the Julian day verbose
my_juliandayToString              = libmtime.juliandayToString;
my_juliandayToString.argtype      = [POINTER(_julianday), c_char_p];
my_juliandayToString.restype      = None;


def juliandayToString(jd):
    name = create_string_buffer(b" ",max_julianday_str_len);
    my_juliandayToString(jd,name);
    return name.value;

#################################################################################################################

## @package mtime_date  
#  @brief Date and some operations supported on Date
#
# @details
#
# @author  Luis Kornblueh, Max Planck Institute for Meteorology
# @author  Rahul Sinha, Max Planck Institute for Meteorology



## provides a string length for toString 
max_date_str_len = c_int(32).value;

class _date(Structure):
  _fields_ = [
		("year", c_long),
                ("month", c_int),
                ("day", c_int)	
	     ];

## @brief construct a new date 
#
# USAGE example: d = newDate("1999-07-01");
#
# @param[in] string         an ISO 8601 conforming date string
# @return    ret_date       a pointer of type(_date)
newDate              = libmtime.newDate;
newDate.argtype      = [c_char_p];
newDate.restype      = POINTER(_date);



## @brief destructor for a date 
#
# USAGE example: deallocateDate(d);
#
# @param[in] my_date        a pointer of type(_date)
deallocateDate              = libmtime.deallocateDate;
deallocateDate.argtype      = [POINTER(_date)];
deallocateDate.restype      = None;



## @brief get Date as a string.
#
#  USAGE example: 	str = create_string_buffer(b" ",max_date_str_len);
#			str = dateToString(d,str);
#
#  @param[in]  d        a pointer to type(_date). The Date to be converted to a string
#  @param[in]  str      string where Date verbose is written
#  @param[out] str      the Date verbose
my_dateToString              = libmtime.dateToString;
my_dateToString.argtype      = [POINTER(_date), c_char_p];
my_dateToString.restype      = None;


def dateToString(d):
    name = create_string_buffer(b" ",max_date_str_len);
    my_dateToString(d,name);
    return name.value;

#########################################################################


## @package mtime_time
# @brief Time and some operations supported on Time
#
# @details
#
# @author  Luis Kornblueh, Max Planck Institute for Meteorology
# @author  Rahul Sinha, Max Planck Institute for Meteorology


## provides a string length for toString 
max_time_str_len = c_int(32).value;

class _time(Structure):
  _fields_ = [
                ("hour", c_int),
                ("minute", c_int),
                ("second", c_int),
		("ms", c_int)
             ];


newTime              = libmtime.newTime;
newTime.argtype      = [c_char_p];
newTime.restype      = POINTER(_time);


deallocateTime              = libmtime.deallocateTime;
deallocateTime.argtype      = [POINTER(_time)];
deallocateTime.restype      = None;


my_timeToString              = libmtime.timeToString;
my_timeToString.argtype      = [POINTER(_time), c_char_p];
my_timeToString.restype      = None;


def timeToString(t):
    name = create_string_buffer(b" ",max_time_str_len);
    my_timeToString(t,name);
    return name.value;


###########################################################################

## @package mtime_datetime
# @brief DateTime and some operations supported on DateTime.
#
# @details
#
# @author  Luis Kornblueh, Max Planck Institute for Meteorology
# @author  Rahul Sinha, Max Planck Institute for Meteorology
#

## provides a string length for toString
max_datetime_str_len = c_int(32).value;




class _datetime(Structure):
    _fields_ = [
                ("date", _date),
                ("time",_time),
               ];
    '''
    #Operator overloading.
    #TODO: Couldn't hack the assigment operator. Looks like there is no good way to do this in python :-/

    # a > b.
    def __gt__(self, x):
                ret_val = compareDatetime(self,x);
                if (ret_val == 1):
                  return True;
                else:
                  return False;
    # a < b.
    def __lt__(self, x):
                ret_val = compareDatetime(self,x);
                if (ret_val == -1):
                  return True;
                else:
                  return False;
    # a <= b.
    def __le__(self, x):
                ret_val = compareDatetime(self,x);
                if (ret_val <= 0):
                  return True;
                else:
                  return False;
    # a >= b.
    def __ge__(self, x):
                ret_val = compareDatetime(self,x);
                if (ret_val >= 0):
                  return True;
                else:
                  return False;
    # a == b.
    def __eq__(self, x):
                ret_val = compareDatetime(self,x);
                if (ret_val == 0):
                  return True;
                else:
                  return False;
    # a != b.
    def __ne__(self, x):
                ret_val = compareDatetime(self,x);
                if (ret_val != 0):
                  return True;
                else:
                  return False;
    '''

newDateTime              = libmtime.newDateTime;
newDateTime.argtype      = [c_char_p];
newDateTime.restype      = POINTER(_datetime);


newdatetimefromraw             = libmtime.newRawDateTime;
newdatetimefromraw.argtype     = [c_long, c_int, c_int, c_int, c_int, c_int, c_int];
newdatetimefromraw.restype     = POINTER(_datetime)


newdatetimefromconstructandcopy              = libmtime.constructAndCopyDateTime;
newdatetimefromconstructandcopy.argtype      = [POINTER(_datetime)];
newdatetimefromconstructandcopy.restype      = POINTER(_datetime);


deallocateDateTime              = libmtime.deallocateDateTime;
deallocateDateTime.argtype      = [POINTER(_datetime)];
deallocateDateTime.restype      = None;


my_datetimeToString              = libmtime.datetimeToString;
my_datetimeToString.argtype      = [POINTER(_datetime), c_char_p];
my_datetimeToString.restype      = None;


def datetimeToString(dt):
    name = create_string_buffer(b" ",max_datetime_str_len);
    my_datetimeToString(dt,name);
    return name.value;


replaceDatetime             = libmtime.replaceDatetime;
replaceDatetime.argtype     = [POINTER(_datetime),POINTER(_datetime)];
replaceDatetime.restype     = POINTER(_datetime);


compareDatetime              = libmtime.compareDatetime;
compareDatetime.argtype      = [POINTER(_datetime), POINTER(_datetime)];
compareDatetime.restype      = c_int;


getNoOfDaysInMonthDateTime              = libmtime.getNoOfDaysInMonthDateTime;
getNoOfDaysInMonthDateTime.argtype      = [POINTER(_datetime)];
getNoOfDaysInMonthDateTime.restype      = c_int


getNoOfDaysInYearDateTime              = libmtime.getNoOfDaysInYearDateTime;
getNoOfDaysInYearDateTime.argtype      = [POINTER(_datetime)];
getNoOfDaysInYearDateTime.restype      = c_int


getDayOfYearFromDateTime                = libmtime.getDayOfYearFromDateTime;
getDayOfYearFromDateTime.argtype        = [POINTER(_datetime)];
getDayOfYearFromDateTime.restype        = c_int


getJulianDayFromDateTime                = libmtime.getJulianDayFromDateTime;
getJulianDayFromDateTime.argtype        = [POINTER(_datetime), POINTER(_julianday)];
getJulianDayFromDateTime.restype        = POINTER(_julianday);

############################################################################

## @package mtime_timedelta
# 
# @brief TimeDelta and some operations supported on TimeDelta.
#
# @details
#
# @author  Luis Kornblueh, Max Planck Institute for Meteorology
# @author  Rahul Sinha, Max Planck Institute for Meteorology

## provides a string length for toString
max_timedelta_str_len = c_int(32).value;

class _timedelta(Structure):
    _fields_ = [
                ("sign", c_char),
                ("year", c_long),
                ("month", c_int),
                ("day", c_int),
                ("hour", c_int),
                ("minute", c_int),
                ("second", c_int),
                ("ms", c_int)
               ];

newTimeDelta              = libmtime.newTimeDelta;
newTimeDelta.argtype      = [c_char_p];
newTimeDelta.restype      = POINTER(_timedelta);


deallocateTimeDelta              = libmtime.deallocateTimeDelta;
deallocateTimeDelta.argtype      = [POINTER(_timedelta)];
deallocateTimeDelta.restype      = None;


addTimeDeltaToDateTime              = libmtime.addTimeDeltaToDateTime;
addTimeDeltaToDateTime.argtype      = [POINTER(_datetime), POINTER(_timedelta), POINTER(_datetime)];
addTimeDeltaToDateTime.restype      = POINTER(_datetime);


my_timedeltaToString              = libmtime.timedeltaToString;
my_timedeltaToString.argtype      = [POINTER(_timedelta), c_char_p];
my_timedeltaToString.restype      = c_char_p


def timedeltaToString(td):
    name = create_string_buffer(b" ",max_timedelta_str_len);
    my_timedeltaToString(td,name);
    return name.value;

class _divisionquotienttimespan(Structure):
    _fields_ = [
                ("quotient", c_long),
                ("remainder_in_ms", c_long)
               ]

divideDatetimeDifferenceInSeconds = libmtime.divideDatetimeDifferenceInSeconds
divideDatetimeDifferenceInSeconds.argtype = [
    POINTER(_datetime),
    POINTER(_datetime),
    POINTER(_timedelta),
    POINTER(_divisionquotienttimespan)
]
divideDatetimeDifferenceInSeconds.restype = POINTER(_divisionquotienttimespan)


############################################################################

## @package mtime_event
#
# @brief Definition of the basic event type and its methods.
#
# @details
# 
# @author  Luis Kornblueh, Max Planck Institute for Meteorology
# @author  Rahul Sinha, Max Planck Institute for Meteorology



## provides a string length for toString.
max_event_str_len = c_int(132).value;

class _event(Structure):
   _fields_ = [
		("nextEventInGroup", c_void_p),
		("eventid", c_long),
		("eventname", c_char_p),		

		("eventreferencedatetime", POINTER(_datetime)),
		("eventfirstdatetime", POINTER(_datetime)),
		("eventlastdatetime", POINTER(_datetime)),

		("eventinterval", POINTER(_datetime)),

		("triggercurrentevent", c_bool),

		("nextEventIsFirst", c_bool),

		("eventisFirstInDay", c_bool),
		("eventisFirstInMonth", c_bool),
		("eventisFirstInYear", c_bool),
		("eventisLastInDay", c_bool),
		("eventisLastInMonth", c_bool),
		("eventisLastInYear", c_bool),

		("triggernexteventdatetime", POINTER(_datetime)),
		("triggeredpreviouseventdatetime", POINTER(_datetime)),
	       ];



newEvent              = libmtime.newEvent;
newEvent.argtype      = [c_char_p, c_char_p, c_char_p, c_char_p, c_char_p];
newEvent.restype      = POINTER(_event);


deallocateEvent         = libmtime.deallocateEvent;
deallocateEvent.argtype = [POINTER(_event)];
deallocateEvent.restype = None;


isCurrentEventActive         = libmtime.isCurrentEventActive;
isCurrentEventActive.argtype = [POINTER(_event), POINTER(_datetime)];
isCurrentEventActive.restype = c_bool;

my_eventToString              = libmtime.eventToString;
my_eventToString.argtype      = [POINTER(_event), c_char_p];
my_eventToString.restype      = c_char_p



def eventToString(e):
    name = create_string_buffer(b" ",max_event_str_len);
    my_eventToString(e,name);
    return name.value;


############################################################################################

## @package mtime_eventgroup
#
# @brief Event-groups which contains a list of events.
#
# @details
# 
# @author  Luis Kornblueh, Max Planck Institute for Meteorology
# @author  Rahul Sinha, Max Planck Institute for Meteorology


## provides a string length for toString
max_groupname_str_len = 132

class _eventgroup(Structure):
    _fields_ = [
                ("eventgroupid", c_int),
                ("eventgroupname", c_char_p),
                ("rootEvent", POINTER(_event) ),
               ];


newEventGroup              = libmtime.newEventGroup;
newEventGroup.argtype      = [c_char_p];
newEventGroup.restype      = POINTER(_eventgroup);

deallocateEventGroup         = libmtime.deallocateEventGroup;
deallocateEventGroup.argtype = [POINTER(_eventgroup)];
deallocateEventGroup.restype = None;

addNewEventToEventGroup              = libmtime.addNewEventToEventGroup;
addNewEventToEventGroup.argType      = [POINTER(_event), POINTER(_eventgroup) ];
addNewEventToEventGroup.resType      = c_bool;

removeEventFromEventGroup              = libmtime.removeEventFromEventGroup;
removeEventFromEventGroup.argType      = [c_char_p, POINTER(_eventgroup) ];
removeEventFromEventGroup.resType      = c_bool;

##
# @}
