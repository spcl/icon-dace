!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
MODULE mtime_constants

  IMPLICIT NONE
  PUBLIC

  !> provides a string length for toString
  INTEGER, PARAMETER :: max_calendar_str_len = 32

  !> provides a string length for toString
  INTEGER, PARAMETER :: max_date_str_len = 32

  !> provides a string length for toString
  INTEGER, PARAMETER :: max_datetime_str_len = 32

  !> provides a string length for toString
  INTEGER, PARAMETER :: max_time_str_len = 32

  !> provides a string length for toString
  INTEGER, PARAMETER :: max_julianday_str_len = 32

  !> provides a string length for toString
  INTEGER, PARAMETER :: max_timedelta_str_len = 32

  !> provides a string length for the maximum error string length
  INTEGER, PARAMETER :: max_mtime_error_str_len = 132

  !> provides a string length for toString
  INTEGER, PARAMETER :: max_eventname_str_len = 132
  INTEGER, PARAMETER :: max_event_str_len = 512

  !> provides a string length for toString
  INTEGER, PARAMETER :: max_groupname_str_len = 132

CONTAINS

END MODULE mtime_constants
