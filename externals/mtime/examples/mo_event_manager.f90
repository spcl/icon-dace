!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
module mo_event_manager

  use mtime

  implicit none

  private

  public :: initEventManager 
  public :: getModelReferenceDate
  public :: addEventGroup
  public :: getEventGroup 
  public :: printEventGroup 
  public :: getEventComponents

  type event_group_list
    type(eventgroup), pointer :: group
  end type event_group_list

  type(event_group_list), allocatable :: model_event_groups(:)

  integer :: model_event_groups_list_member
  integer :: model_event_groups_list_size = 16

  type(datetime), pointer :: model_reference_date => null()

  logical :: linitialized = .false.

contains

  subroutine initEventManager(referenceDate)
    type(datetime), pointer, intent(in) :: referenceDate

    model_reference_date => newDatetime(referenceDate)

    allocate(model_event_groups(model_event_groups_list_size)) 
    model_event_groups_list_member = 0

    linitialized = .true.

  end subroutine initEventManager

  function getModelReferenceDate() result(r)
    type(datetime), pointer :: r

    r => null()
    if (linitialized) then
      r => model_reference_date
    endif

  end function getModelReferenceDate

  function addEventGroup(group) result(handle)
    integer :: handle
    character(len=*), intent(in) :: group
    type(event_group_list), allocatable :: tmp(:)
    character(len=max_groupname_str_len) :: gstring    

    if (.not. linitialized) then
      print *, 'ERROR: event manager not initialized.' 
      stop
    endif

    if (model_event_groups_list_member == model_event_groups_list_size) then
      print *, 'INFO: reallocating event group list.' 
      allocate(tmp(model_event_groups_list_size)) 
      tmp(:) = model_event_groups(:)
      deallocate(model_event_groups)
      allocate(model_event_groups(2*model_event_groups_list_size))
      model_event_groups(:model_event_groups_list_size) = tmp(:)
      deallocate(tmp)
      model_event_groups_list_size = 2*model_event_groups_list_size
      print *, 'INFO: new evcent group list size: ', model_event_groups_list_size 
    endif

    model_event_groups_list_member = model_event_groups_list_member + 1

    model_event_groups(model_event_groups_list_member)%group => newEventGroup(trim(group))
    call getEventGroupName(model_event_groups(model_event_groups_list_member)%group, gstring)
    print *, 'INFO: Added event group: ', trim(gstring)

    handle = model_event_groups_list_member

  end function addEventGroup

  function getEventGroup(handle) result(eventGroupListMember)
    type(eventGroup), pointer :: eventGroupListMember
    integer, intent(in) :: handle
    if (handle > model_event_groups_list_member) then
       eventGroupListMember => null()
     else
       eventGroupListMember =>  model_event_groups(handle)%group
     endif
  end function getEventGroup

  subroutine printEventGroup(handle)
    integer, intent(in) :: handle
    type(eventGroup), pointer :: currentEventGroup
    type(event), pointer :: currentEvent
    character(len=max_eventname_str_len) :: estring
    character(len=max_groupname_str_len) :: egstring

    currentEventGroup => getEventGroup(handle)
    call getEventGroupName(currentEventGroup, egstring)

    currentEvent => getFirstEventFromEventGroup(model_event_groups(handle)%group)

    print *, 'Event list: ', trim(egstring)
    do while (associated(currentEvent))
      call eventToString(currentEvent, estring)
      print *,'   event ', trim(estring)
      currentEvent => getNextEventFromEventGroup(currentEvent)
    enddo
  end subroutine printEventGroup

  subroutine getEventComponents(eventString, referenceDate, timeInterval, startDate, endDate)
    character(len=max_repetition_str_len), intent(in) :: eventString 
    type(datetime),  pointer :: referenceDate
    type(timedelta), pointer :: timeInterval
    type(datetime),  pointer :: startDate
    type(datetime),  pointer :: endDate
    
    character(len=max_repetition_str_len) :: r, s, e, d    
    logical :: lr, ls, le, ld    

    call splitRepetitionString(eventString, r, s, e, d, lr, ls, le, ld)

    if (lr) then
      if (getRepetitions(r) /= -1) then
        print *, 'WARNING: event setup should not have explicit repeat count.'
      endif
    endif
    
    if (ls) then
      startDate => newDatetime(trim(s))
    endif
    
    if (le) then
      endDate => newDatetime(trim(e))
    endif
    
    if (ld) then
      timeInterval => newTimeDelta(trim(d))
    else
      print *, 'ERROR: time interval should be given.'
      stop
    endif

  end subroutine getEventComponents

end module mo_event_manager

