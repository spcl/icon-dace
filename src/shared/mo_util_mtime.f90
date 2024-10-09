! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! This module contains auxiliary routines for mtime date/time objects.

MODULE mo_util_mtime

  USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t, c_char
  USE mo_kind,                     ONLY: dp, wp, i8
  USE mo_exception,                ONLY: finish
  USE mo_impl_constants,           ONLY: MAX_CHAR_LENGTH
  USE mtime,                       ONLY: datetime, newDatetime, timedelta, newTimeDelta,   &
    &                                    OPERATOR(-), OPERATOR(+), event,      &
    &                                    juliandelta,                          &
    &                                    timeDeltaToJulianDelta,               &
    &                                    deallocateTimeDelta, deallocateDatetime,          &
    &                                    getTimeDeltaFromDateTime,                         &
    &                                    getTotalMillisecondsTimedelta,                    &
    &                                    isCurrentEventActive
  USE mo_time_config,              ONLY: time_config
  USE mo_util_string,              ONLY: t_keyword_list,                   &
                                         associate_keyword, with_keywords, &
                                         int2string
  USE mo_mpi,                      ONLY: p_pe, p_io, p_comm_work, p_bcast

  IMPLICIT NONE

  PUBLIC :: t_mtime_utils, mtime_utils
  PUBLIC ::  FMT_DDHHMMSS_ANNOTATED, FMT_DDDHHMMSS_ANNOTATED, &
    &        FMT_DDHHMMSS, FMT_DDDHHMMSS, FMT_DDDHH, FMT_HHH, FMT_DDHHMMSS_DAYSEP
  PUBLIC :: assumePrevMidnight
  PUBLIC :: assumeNextMidnight
  PUBLIC :: getElapsedSimTimeInSeconds
  PUBLIC :: dummyDateTime
  PUBLIC :: t_datetime_ptr
  PUBLIC :: is_event_active
  PUBLIC :: mtime_convert_netcdf_units
  PUBLIC :: mtime_divide_timedelta
  PUBLIC :: mtime_timedelta_from_fseconds
  PUBLIC :: mtime_timedelta_to_seconds

  TYPE t_datetime_ptr
    TYPE(datetime), POINTER :: ptr => NULL()
  END TYPE t_datetime_ptr

  PRIVATE

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_mtime'

  !> a collection of relative time formats
  CHARACTER(LEN=*), PARAMETER :: FMT_DDHHMMSS_ANNOTATED  = "<dd>D<hh>H<mm>M<ss><ms>S"
  CHARACTER(LEN=*), PARAMETER :: FMT_DDHHMMSS_DAYSEP     = "<dd>D <hh>H<mm>M<ss><ms>S"
  CHARACTER(LEN=*), PARAMETER :: FMT_DDDHHMMSS_ANNOTATED = "<ddd>D<hh>H<mm>M<ss><ms>S"
  CHARACTER(LEN=*), PARAMETER :: FMT_DDHHMMSS            = "<dd><hh><mm><ss><ms>"
  CHARACTER(LEN=*), PARAMETER :: FMT_DDDHHMMSS           = "<ddd><hh><mm><ss><ms>"
  CHARACTER(LEN=*), PARAMETER :: FMT_DDDHH               = "<ddd><hh>"
  CHARACTER(LEN=*), PARAMETER :: FMT_HHH                 = "<hhh>"

  TYPE t_mtime_utils
  CONTAINS
    PROCEDURE, NOPASS :: ddhhmmss => t_mtime_utils_ddhhmmss
  END TYPE t_mtime_utils

  ! global object serving as a collection of static utility routines
  TYPE(t_mtime_utils) :: mtime_utils
  
CONTAINS 

  !> @return DDHHMMSS string of a given time span.
  !
  !  Note: We return the forecast time delta as an ISO 8601-like
  !        string, where the exact format is given by the "format" argument.
  !
  FUNCTION t_mtime_utils_ddhhmmss(start_date, end_date, fmt)  RESULT(result_str)
    CHARACTER(len=MAX_CHAR_LENGTH)  :: result_str
    TYPE(datetime),    INTENT(IN)   :: start_date, end_date
    CHARACTER(LEN=*),  INTENT(IN)   :: fmt
    
    TYPE(timedelta) :: time_delta
    TYPE(juliandelta) :: julian_delta
    TYPE(t_keyword_list), POINTER  :: keywords => NULL()

    ! compute time delta:
    time_delta = end_date - start_date
    ! compute the time delta in a Julian calendar, which allows more
    ! than 30 days:
    CALL timeDeltaToJulianDelta(time_delta,start_date,julian_delta)

    ! now collect the feasible "building-blocks"
    IF (INT(julian_delta%day) < 100) THEN
      CALL associate_keyword("<dd>",  TRIM(int2string(INT(julian_delta%day), '(i2.2)')), keywords)
    ELSE
      CALL associate_keyword("<dd>",  TRIM(int2string(INT(julian_delta%day), '(i0)')),   keywords)
    END IF
    CALL associate_keyword("<ddd>", TRIM(int2string(INT(julian_delta%day), '(i3.3)')),   keywords)
    CALL associate_keyword("<hhh>", TRIM(int2string(INT(julian_delta%day)*24+time_delta%hour, '(i3.3)')),   keywords)
    CALL associate_keyword("<hh>",  TRIM(int2string(time_delta%hour,       '(i2.2)')),   keywords)
    CALL associate_keyword("<mm>",  TRIM(int2string(time_delta%minute,     '(i2.2)')),   keywords)
    CALL associate_keyword("<ss>",  TRIM(int2string(time_delta%second,     '(i2.2)')),   keywords)
    IF (time_delta%ms /= 0) THEN
      CALL associate_keyword("<ms>",  '.'//TRIM(int2string(time_delta%second,'(i3.3)')),   keywords)
    ELSE
      CALL associate_keyword("<ms>",  "", keywords)
    END IF
    ! replace keywords in "ddhhmmss"
    result_str = with_keywords(keywords, fmt)

  END FUNCTION t_mtime_utils_ddhhmmss


  !>
  !! Returns previous 'midnight' datetime for given datetime
  !!
  !! Returns previous 'midnight' datetime for given datetime in datetime-format
  !!
  TYPE(datetime) FUNCTION assumePrevMidnight(current_date)
    TYPE(datetime), INTENT(IN) :: current_date

    assumePrevMidnight = current_date
    !
    assumePrevMidnight%time%hour   = 0 
    assumePrevMidnight%time%minute = 0
    assumePrevMidnight%time%second = 0
    assumePrevMidnight%time%ms     = 0

  END FUNCTION assumePrevMidnight

  !>
  !! Returns nest 'midnight' datetime for given datetime
  !!
  !! Returns next 'midnight' datetime for given datetime in datetime-format
  !!
  TYPE(datetime) FUNCTION assumeNextMidnight(current_date)
    TYPE(datetime), INTENT(IN) :: current_date

    ! local
    TYPE(timedelta), POINTER :: td_1day  => NULL()
    !---------------------------------------------------------

    td_1day =>newTimedelta("P01D")
    assumeNextMidnight = current_date + td_1day
    !
    assumeNextMidnight%time%hour   = 0 
    assumeNextMidnight%time%minute = 0
    assumeNextMidnight%time%second = 0
    assumeNextMidnight%time%ms     = 0

    CALL deallocateTimedelta(td_1day)

  END FUNCTION assumeNextMidnight


  !>
  !! Returns elapsed simulation time in seconds
  !!
  !! Elapsed simulation time is computed as the timedelta 
  !! between the current datetime and the anchor  
  !! datetime. Unless specified otherwise, the anchor date  
  !! is set to the experiment start date time_config%tc_exp_startdate.
  !!
  REAL(wp) FUNCTION getElapsedSimTimeInSeconds(datetime_current, anchor_datetime) RESULT(sim_time)
    !
    TYPE(datetime)          , INTENT(IN) :: datetime_current
    TYPE(datetime), OPTIONAL, INTENT(IN) :: anchor_datetime   ! anchor date
    !
    ! local
    TYPE(timedelta) :: time_diff
    TYPE(datetime)  :: anchor                        ! anchor date for
                                                     ! timeDelta computation
    !---------------------------------------------------------

    If (PRESENT(anchor_datetime)) THEN
      anchor = anchor_datetime
    ELSE
      ! set anchor date to experiment start date
      anchor = time_config%tc_exp_startdate
    ENDIF

    time_diff  =  getTimeDeltaFromDateTime(datetime_current, anchor)
    sim_time   =  getTotalMillisecondsTimedelta(time_diff, anchor)*1.e-3_wp

  END FUNCTION getElapsedSimTimeInSeconds



  !>
  !! Creates a datetime object with a defined dummy date
  !!
  !! Creates a datetime object with a defined dummy date
  !! Can be used for initializing a newly defined datetime 
  !! object and for checking, whether the datetime object 
  !! at hand has already been touched and filed with a meaningful value.
  !!
  TYPE(datetime) FUNCTION dummyDateTime()
    !
    ! local
    TYPE(datetime), POINTER :: dummy_date_ptr
    LOGICAL, SAVE :: dummy_date_is_set = .FALSE.
    TYPE(datetime), SAVE :: dummy_date
    !---------------------------------------------------------

    IF (.NOT. dummy_date_is_set) THEN
      dummy_date_ptr => newDatetime("0-01-01T00:00:00.000")
      dummy_date = dummy_date_ptr
      CALL deallocateDatetime(dummy_date_ptr)
      dummy_date_is_set = .TRUE.
    END IF
    dummydatetime = dummy_date

  END FUNCTION dummyDateTime



  !>
  !! Wrapper for mtime function isCurrentEventActive in order to 
  !! encapsulate the vector-host offloading
  !! needed on the NEC Aurora
  !!
  LOGICAL FUNCTION is_event_active (in_event, mtime_current, offload_mode, plus_slack, opt_lasync)

    TYPE(event),     POINTER, INTENT(INOUT)           :: in_event       !< mtime event to be checked
    TYPE(datetime),  POINTER, INTENT(IN   )           :: mtime_current  !< current_datetime
    LOGICAL,                  INTENT(IN   )           :: offload_mode   !< if TRUE, PE0 is in offloading mode
    TYPE(timedelta), POINTER, INTENT(IN   ), OPTIONAL :: plus_slack
    LOGICAL,                  INTENT(IN   ), OPTIONAL :: opt_lasync     !< if TRUE, broadcast is done by caller

    LOGICAL :: lasync
    LOGICAL :: is_active
  !-----------------------------------------------------------------

    IF (PRESENT(opt_lasync)) THEN
     ! If TRUE, broadcast is done by caller
     lasync = opt_lasync
    ELSE
     ! broadcast is done directly
     lasync = .FALSE.
    ENDIF

    ! If PE0 is detached, execute isCurrentEventActive only on PE0 and broadcast result afterwards
    IF (.NOT. offload_mode .OR. p_pe == p_io) THEN
      IF (PRESENT(plus_slack)) THEN
        is_active = isCurrentEventActive(in_event, mtime_current, plus_slack=plus_slack)
      ELSE
        is_active = isCurrentEventActive(in_event, mtime_current)
      ENDIF
    ENDIF
    IF ( .NOT.lasync .AND. offload_mode ) CALL p_bcast(is_active, p_io, p_comm_work)

    is_event_active = is_active

  END FUNCTION is_event_active



  ! via seconds... should be solved by mtime, could be done without
  ! conversion to a prespecified time unit.
  SUBROUTINE mtime_divide_timedelta (dividend, divisor, quotient)
    TYPE(timedelta),          INTENT(in   ) :: dividend
    TYPE(timedelta),          INTENT(in   ) :: divisor
    REAL(wp),                 INTENT(  out) :: quotient
    INTEGER(i8)                             :: dd, ds

    CALL mtime_timedelta_to_seconds(dividend, dd)
    CALL mtime_timedelta_to_seconds(divisor,  ds)

    quotient = REAL(dd,wp) / REAL(ds,wp)

  END SUBROUTINE mtime_divide_timedelta

  SUBROUTINE mtime_timedelta_to_seconds (mtime_dt, seconds_dt)
    TYPE(timedelta),          INTENT(in   ) :: mtime_dt
    INTEGER(i8),              INTENT(  out) :: seconds_dt
    CHARACTER(*),                 PARAMETER :: routine = &
      & modname//"::mtime_timedelta_to_seconds"

    IF (mtime_dt%year > 0 .OR. mtime_dt%month > 0) THEN
      CALL finish(routine, "year or month nonzero, result is undefined")
    ENDIF

    seconds_dt =   86400 * INT(mtime_dt%day)    &
      &          + 3600  * INT(mtime_dt%hour)   &
      &          + 60    * INT(mtime_dt%minute) &
      &          +         INT(mtime_dt%second)

  END SUBROUTINE mtime_timedelta_to_seconds

  ! convert (fractional) seconds into an mtime timedelta object
  SUBROUTINE mtime_timedelta_from_fseconds(dt, base_dt, td)
    REAL(dp), INTENT(in) :: dt
    TYPE(datetime), INTENT(in) :: base_dt
    TYPE(timedelta), INTENT(out) :: td
    INTEGER :: change_sign
    REAL(dp) :: dtl
    INTERFACE
      SUBROUTINE julianDeltaToTimeDelta(jd, base_dt, td) BIND(c, name='julianDeltaToTimeDelta')
        IMPORT :: juliandelta, datetime, timedelta
        TYPE(juliandelta), INTENT(in) :: jd
        TYPE(datetime), INTENT(in) :: base_dt
        TYPE(timedelta), INTENT(out) :: td
      END SUBROUTINE julianDeltaToTimeDelta
    END INTERFACE
    TYPE(juliandelta) :: jdelta

    jdelta%sign = MERGE(c_char_'+', c_char_'-', dt >= 0.0_dp)
    change_sign = MERGE(1, -1, dt >= 0.0_dp)    
    dtl = ABS(dt)
    jdelta%ms = change_sign*INT(ABS(MOD(dtl, 86400.0_dp)) * 1000.0_dp, c_int64_t)
    jdelta%day = change_sign*INT(dtl/86400.0_dp, c_int64_t)
    CALL juliandeltatotimedelta(jdelta, base_dt, td)
  END SUBROUTINE mtime_timedelta_from_fseconds

  ! the expected input looks similar to this "seconds since 2013-04-24T00:00:00"

  SUBROUTINE mtime_convert_netcdf_units (unit_string, starttime, timeincrement)
    CHARACTER(len=*),         INTENT(in   ) :: unit_string
    TYPE(datetime),  POINTER, INTENT(  out) :: starttime
    TYPE(timeDelta), POINTER, INTENT(  out) :: timeincrement

    INTEGER           :: errno
    INTEGER           :: idx0, idx1, idx2, idx3, idx4, idx5, idx6
    INTEGER           :: month, day
    INTEGER(i8)       :: year
    CHARACTER(len=32) :: iso8601_string
    CHARACTER(len=32) :: timedelta_string

    idx1 = SCAN(unit_string,          " ")
    idx2 = SCAN(unit_string(idx1+1:), " ")
    idx2 = idx1+idx2
    idx3 = SCAN(unit_string(idx2+1:), " ")
    idx3 = idx2+idx3

    idx0 = 1
    idx4 = LEN_TRIM(unit_string)
    idx5 = SCAN(unit_string(idx2+1:idx3-1), "-")
    idx6 = SCAN(unit_string(idx2+idx5+1:idx3-1), "-")
    idx6 = idx5+idx6

    READ(unit_string(idx2+1:idx2+idx5-1),*) year
    READ(unit_string(idx2+idx5+1:idx2+idx6-1),*) month
    READ(unit_string(idx2+idx6+1:idx3-1),*) day

    WRITE(iso8601_string,'(i0,a,i2.2,a,i2.2,a,a)') &
         year, '-', month, '-', day, 'T', unit_string(idx3+1:idx4)

    IF (unit_string(idx1+1:idx2) /= "since") THEN
      CALL finish(modname, "Unit string is in unknown format"//unit_string(idx1+1:idx2))
    ENDIF

    SELECT CASE ( unit_string(:idx1-1) )
    CASE ("seconds")
      timedelta_string = "PT01S"
    CASE ("minutes")
      timedelta_string = "PT01M"
    CASE ("hours")
      timedelta_string = "PT01H"
    CASE ("days")
      timedelta_string = "P01D"
    CASE default
      CALL finish(modname, "Unknown time increment: "//unit_string(:idx1-1))
    END SELECT

    starttime     => newdatetime(TRIM(iso8601_string), errno)
    IF (errno /= 0)  CALL finish(modname, "Error in conversion of date: "//TRIM(iso8601_string))

    timeincrement => newTimeDelta(timedelta_string, errno)
    IF (errno /= 0)  CALL finish(modname, "Error in conversion of unit: "//timedelta_string)

  END SUBROUTINE mtime_convert_netcdf_units

END MODULE mo_util_mtime
