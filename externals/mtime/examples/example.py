# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#
from libmtime import *
from libmtime import _event


print '\n\nsetup calendar:'
# setup calendar.
setCalendar(CALENDAR_TYPE.proleptic_gregorian);
#setCalendar(CALENDAR_TYPE.year_of_365_days);
# Check the current calendar.
CalType = getCalendarType();
#cType is 0-3 denoting calendar type. They can be compared with CALENDAR_TYPE.calendar_not_set, CALENDAR_TYPE.proleptic_gregorian etc..
if (CalType.cType == CALENDAR_TYPE.calendar_not_set):
	print "CURRENT CAL TYPE = calendar_not_set"
elif (CalType.cType == CALENDAR_TYPE.proleptic_gregorian):
	print "CURRENT CAL TYPE = proleptic_gregorian";
elif (CalType.cType == CALENDAR_TYPE.year_of_365_days):
        print "CURRENT CAL TYPE = year_of_365_days";
elif (CalType.cType == CALENDAR_TYPE.year_of_360_days):
        print "CURRENT CAL TYPE = year_of_360_days";


calendar_in_use = calendarToString();
print calendar_in_use;




print '\n\ntest date interface:'
# test date interface
test_date = newDate("2012-09-01");
if test_date is not None:
	print 'allocated test_date';
	test_date_string = dateToString(test_date);
	print test_date_string
	#Deallocate.
	test_date = deallocateDate(test_date);
	if test_date is None:
		print 'deallocated test_time'




print '\n\n test time interface:'
# test time interface
test_time = newTime("12:13:49.654");
if test_time is not None:
	print 'allocated test_time'
	test_time_string = timeToString(test_time);
	print test_time_string;
	test_time = deallocateTime(test_time);
	if test_time is None:
		print 'deallocated test_time'


print '\n\n#test datetime and timeddelta interface. Simulating a dummy timer.'
#test datetime and timeddelta interface. Simulating a dummy timer.
start_date = newDateTime("2012-09-01T02:10:00.000");
if start_date is not None:
	print 'allocated start_date'
	start_date_string = datetimeToString(start_date);
	print start_date_string;

stop_date = newDateTime("2012-09-10T14:00:00.000");
if stop_date is not None:
        print 'allocated stop_date'
        stop_date_string = datetimeToString(stop_date);
        print stop_date_string;


time_step = newTimeDelta("PT12H");
if time_step is not None:
	print 'allocated time_step'
	time_step_string = timedeltaToString(time_step);
	print time_step_string;

current_date = newDateTime("2012-09-01T02:10:00.000");
if current_date is not None:
        print 'allocated current_date'
        current_date_string = datetimeToString(current_date);
        print current_date_string;


print '\n'
while  (compareDatetime(current_date,stop_date) < 0):
	print 'Model time loop  : ' + datetimeToString(current_date);
	current_date = addTimeDeltaToDateTime(current_date,time_step,current_date);
print '\n\n';


deallocateDateTime(start_date)
deallocateDateTime(stop_date)
deallocateDateTime(current_date)
deallocateTimeDelta(time_step)





print '\n\n#Test Event interface.'
#Test Event interface.

#Crate an event group.
outputEventGroup =  newEventGroup('output driver');

#Create events and add them to event group.
outputEvent = newEvent('output', '2000-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'PT06H');
lret = addNewEventToEventGroup(outputEvent, outputEventGroup);

checkpointEvent = newEvent('checkpoint', '2000-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'P05D');
lret = addNewEventToEventGroup(checkpointEvent, outputEventGroup);

restartEvent = newEvent('restart', '2000-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'P01M');
lret = addNewEventToEventGroup(restartEvent, outputEventGroup);



#List all events in a group.
print  'Events in outputEventGroup:';
#Start at root of Group.
event_ptr = outputEventGroup.contents.rootEvent;
while (1):
	if (event_ptr != None):
		#Cast the void pointer to POINTER(_event). This is a HACK beause nextEventInGroup had to be defied as c_void_p (the author could 
		#not figure out how to create a member element whose type is a pointer to the still being defined class. If you know how to do this,
		#please contact the authors).
		event_ptr = cast(event_ptr, POINTER(_event));
	else:
		#No more events in group.
		break;

	#Do something with the event.
        print eventToString(event_ptr);
        #Go to next node in the group.
        event_ptr = event_ptr.contents.nextEventInGroup;

#Notice that removing checkpointEvent from the group will automatically free the event. You MUST NOT try to free the event using deallocateEvent(.).
removeEventFromEventGroup("checkpoint", outputEventGroup);

#Deallocate event group.
deallocateEventGroup(outputEventGroup);



#Check Events active?
print '\n\nCheck if any event is currently active.'

#Create an event group and add events to it.
SomeEventGroup =  newEventGroup('SomeEventGroup');

Newton = newEvent('Fisher', '2000-01-01T00:00:00', '2010-01-01T00:00:00', '2015-01-01T00:00:00', 'P06M');
addNewEventToEventGroup(Newton, SomeEventGroup);

Einstein = newEvent('Pearson', '2000-01-01T00:00:00', '2010-01-01T00:00:00', '2014-01-01T00:00:00', 'P03M');
addNewEventToEventGroup(Einstein, SomeEventGroup);

Feynman = newEvent('Cox', '2000-01-01T00:00:00', '2010-01-01T00:00:00', '2013-01-01T00:00:00', 'P01Y');
addNewEventToEventGroup(Feynman, SomeEventGroup);


#Simulate a clock.
dummy_current_time_step = newTimeDelta("P01M");
dummy_current_date = newDateTime("2009-09-01T00:00:00.000");
dummy_stop_date = newDateTime("2016-10-01T00:00:00.000");

#Timer.
while  (compareDatetime(dummy_current_date,dummy_stop_date) < 0):

	#Start at root of Group.
	event_ptr = SomeEventGroup.contents.rootEvent;
	#Loop over all events in an event group.
	while (1):
        	if (event_ptr != None):
                	event_ptr = cast(event_ptr, POINTER(_event));
        	else:
                	break;

        	#If this event is active, do something about it.
		if (isCurrentEventActive(event_ptr,dummy_current_date)):
	        	print eventToString(event_ptr) + ' is active';
        	#Go to next node in the group.
	        event_ptr = event_ptr.contents.nextEventInGroup;



	print 'Model time loop  : ' + datetimeToString(dummy_current_date);
        dummy_current_date = addTimeDeltaToDateTime(dummy_current_date,dummy_current_time_step,dummy_current_date);


'''
#TODO: Could not do this. Follow up.
#Test the overloaded operators.
print dt > dtt;
print dt < dtt;
print dt <= dt;
print dt >= dtt; 
print dt == dt;
print dtt != dt;
print dtt != dtt;
'''


#Do before exit.
#TODO: Clenup data structures.
resetCalendar();

