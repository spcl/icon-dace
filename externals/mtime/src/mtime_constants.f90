!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
module mtime_constants

  implicit none
  public 

  !> provides a string length for toString
  integer, parameter :: max_calendar_str_len = 32
  
  !> provides a string length for toString
  integer, parameter :: max_date_str_len = 32

  !> provides a string length for toString
  integer, parameter :: max_datetime_str_len = 32

  !> provides a string length for toString
  integer, parameter :: max_time_str_len = 32

  !> provides a string length for toString
  integer, parameter :: max_julianday_str_len = 32

  !> provides a string length for toString
  integer, parameter :: max_timedelta_str_len = 32

  !> provides a string length for the maximum error string length
  integer, parameter :: max_mtime_error_str_len = 132

  !> provides a string length for toString
  integer, parameter :: max_eventname_str_len = 132
  integer, parameter :: max_event_str_len = 512

  !> provides a string length for toString
  integer, parameter :: max_groupname_str_len = 132

contains
  
end module mtime_constants
