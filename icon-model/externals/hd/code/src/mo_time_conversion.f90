! mo_time_conversion.f90 - Convert between and compare different date/time formats 
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_time_conversion
  !+
  !
  ! mo_time_conversion [module]
  !    convert between and compare different date/time formats
  !
  ! Version:  E5/R1.07+ 04-February-2003
  !
  ! Authors:
  !   I. Kirchner, MPI, May 2000
  !   I. Kirchner, MPI, Juli 2000 revision
  !   I. Kirchner, MPI, November 2000, add year/month length functions
  !   I. Kirchner, FUB, February 2003, revision/code review
  !
  !-

  !+
  !
  ! external modules
  !   mo_exception
  !   mo_time_base
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.
  !
  !-
  USE mo_kind,      ONLY: dp
  USE mo_exception, ONLY: finish, message
  USE mo_time_base, ONLY: IDAYLEN, get_calendar_type, JULIAN, CYL360,     &
                          julian_date, Set_JulianCalendar, Set_JulianDay, &
                          Get_JulianYearLen, Get_JulianMonLen,            & 
                          Get_JulianYearDay,                              &
                          ly360_date, Set_Ly360Calendar, Set_Ly360Day,    &
                          Get_Ly360YearLen, Get_Ly360MonLen,              & 
                          Get_Ly360YearDay,                               &
                          sec2frac, frac2sec

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: STR_LEN_A = 20

  CHARACTER(len=256) :: message_text

  CHARACTER(len=3), PARAMETER :: CMONTHS(12) = &
       (/ 'Jan','Feb','Mar','Apr','May','Jun',&
          'Jul','Aug','Sep','Oct','Nov','Dec' /)

  TYPE, PUBLIC :: time_days        ! relative calendar date and time format
    !+
    !
    ! time_days [structure]
    !   day    [integer]  (day in calendar, -2147483648 ..... 2147483647
    !                      approx. +/-5.8 Mio. years)
    !   second [integer]  (seconds of day, 0,...,86399)
    !
    !-
    PRIVATE
    LOGICAL :: init   = .FALSE.
    INTEGER :: day    = 0
    INTEGER :: second = 0
  END TYPE time_days

  TYPE, PUBLIC :: time_intern      ! internal date and time format
    !+
    !
    ! time_intern [structure]
    !   ymd [integer]  (YYYYYYMMDD year/month/day, -2147480101,....,2147481231)
    !   hms [integer]  (HHMMSS hour/minute/second, 000000,...,235959)
    !
    !-
    PRIVATE
    LOGICAL :: init = .FALSE.
    INTEGER :: ymd  = 0
    INTEGER :: hms  = 0
  END TYPE time_intern

  TYPE, PUBLIC :: time_native      ! normal splitted date and time format
    !+
    !
    ! time_native [structure]
    !   year   [integer]  -2147480101,....,2147481231)
    !   month  [integer] (1,2,...,12)
    !   day    [integer] (1,2,...,MAX=28,29,30,31)
    !   hour   [integer] (00,...,23)
    !   minute [integer] (00,...,59)
    !   second [integer] (00,...,59)
    !
    !-
    PRIVATE
    LOGICAL :: init   = .FALSE.
    INTEGER :: year   = 0
    INTEGER :: month  = 0
    INTEGER :: day    = 0
    INTEGER :: hour   = 0
    INTEGER :: minute = 0
    INTEGER :: second = 0
  END TYPE time_native


  ! ******* INTERFACE SUBROUTINES AND FUNCTIONS

  ! conversion between different formats
  !

  PUBLIC :: TC_convert
  !+
  !
  ! TC_convert [subroutine, interface]
  !   convert date/time formats
  !   (
  !   date1   [time_days|time_intern|time_native] input
  !   date2   [time_days|time_intern|time_native] output
  !   )
  !
  !-
  INTERFACE TC_convert
    MODULE PROCEDURE convert_day2int ! (day, s) --> (ymd, hms)
    MODULE PROCEDURE convert_int2day ! (ymd, hms) --> (day, s)
    MODULE PROCEDURE convert_nat2int ! (yr, mo, dy, hr, m, s) --> (ymd, hms)
    MODULE PROCEDURE convert_int2nat ! (ymd, hms) --> (yr, mo, dy, hr, m, s)
    MODULE PROCEDURE convert_nat2day ! (yr, mo, dy, hr, m, s) --> (day, s)
    MODULE PROCEDURE convert_day2nat ! (day, s) --> (yr, mo, dy, hr, m, s)
  END INTERFACE

  PUBLIC :: TC_set
  !+
  !
  ! TC_set [subroutine, interface]
  !   fill date/time with values, check date/time information
  !   (
  !   i1, i2, ... [integer] input, format dependent (values for date/time)
  !     day, second                            -> time_days
  !     yyyymmdd, hhmmss                       -> time_intern
  !     year, month, day, hour, minute, second -> time_native
  !   date        [time_days|time_intern|time_native] output (new date/time)
  !   )
  !
  !-
  INTERFACE TC_set
    MODULE PROCEDURE set_days   ! fill without check
    MODULE PROCEDURE set_intern ! time check always, date check optional
    MODULE PROCEDURE set_native ! time and date always to be checked
  END INTERFACE

  PUBLIC :: TC_get
  !+
  !
  ! TC_get [subroutine, interface]
  !    get components of a date/time structure
  !    (
  !    date [time_days|time_intern|time_native] input (date/time structure)
  !    i1, i2, ... [integer] output, format dependent 
  !    (structure component values)
  !     day, second                            -> time_days
  !     yyyymmdd, hhmmss                       -> time_intern
  !     year, month, day, hour, minute, second -> time_native
  !    )
  !
  !-
  INTERFACE TC_get
    MODULE PROCEDURE get_days
    MODULE PROCEDURE get_intern
    MODULE PROCEDURE get_native
  END INTERFACE

  !+
  !
  ! <  [operator, logical]
  ! == [operator, logical]
  ! >  [operator, logical]
  !   comparision between different formats, result is a logical
  !   (
  !   date1  [time_days|time_intern|time_native] (left side)
  !   date2  [time_days|time_intern|time_native] (right side)
  !   )
  !
  !-
  PUBLIC :: OPERATOR(<)
  INTERFACE OPERATOR(<)
    MODULE PROCEDURE L_comp_days_lt
    MODULE PROCEDURE L_comp_intern_lt
    MODULE PROCEDURE L_comp_native_lt
  END INTERFACE

  PUBLIC :: OPERATOR(>)
  INTERFACE OPERATOR(>)
    MODULE PROCEDURE L_comp_days_gt
    MODULE PROCEDURE L_comp_intern_gt
    MODULE PROCEDURE L_comp_native_gt
  END INTERFACE

  PUBLIC :: OPERATOR(==)
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE L_comp_days_eq
    MODULE PROCEDURE L_comp_intern_eq
    MODULE PROCEDURE L_comp_native_eq
  END INTERFACE

  PUBLIC :: add_date
  !+
  !
  ! add_date [subroutine, interface]
  !    add a number of days and seconds to a given date and time
  !    (
  !    day     [integer] input (day offset)
  !    second  [integer] input (offset in seconds)
  !    date    [time_days|time_intern|time_native] input, output
  !    )
  !
  !-
  INTERFACE add_date
    MODULE PROCEDURE add2days
    MODULE PROCEDURE add2intern
    MODULE PROCEDURE add2native
  END INTERFACE

  PUBLIC :: print_date
  !+
  !
  ! print_date [subroutine, interface]
  !    print date and time information to standard output
  !    (
  !    date [time_days|time_intern|time_native] (inp)
  !    outtype [character, optional] (inp) define output format (native, days, intern)
  !    mess    [character, optional] (out) store message in string instead of print out
  !    )
  !
  !-
  INTERFACE print_date
    MODULE PROCEDURE print_days
    MODULE PROCEDURE print_intern
    MODULE PROCEDURE print_native
  END INTERFACE
  CHARACTER(len=*), PARAMETER, PUBLIC :: TC_PRN_DAYS   = 'days'
  CHARACTER(len=*), PARAMETER, PUBLIC :: TC_PRN_INTERN = 'intern'
  CHARACTER(len=*), PARAMETER, PUBLIC :: TC_PRN_NATIVE = 'native'

  PUBLIC :: year_len
  !+
  !
  ! year_len [function, interface]
  !    get length of year in days
  !    with input parameter result in full days
  !    without input parameter result real mean length of year
  !    (
  !    date [time_days|time_intern|time_native] input
  !    )
  !
  !-
  INTERFACE year_len
    MODULE PROCEDURE year_len_days
    MODULE PROCEDURE year_len_intern
    MODULE PROCEDURE year_len_native
    MODULE PROCEDURE year_len_mean
  END INTERFACE

  PUBLIC :: month_len
  !+
  !
  ! month_len [function, interface]
  !    get length of month in days
  !    (
  !    date [time_days|time_intern|time_native] input
  !    )
  !
  !-
  INTERFACE month_len
    MODULE PROCEDURE mon_len_days
    MODULE PROCEDURE mon_len_intern
    MODULE PROCEDURE mon_len_native
  END INTERFACE

  PUBLIC :: day_len

  PUBLIC :: day_in_year
  !+
  !
  ! day_in_year [function, interface]
  !    get number of day in the present year
  !    (
  !    date [time_days|time_intern|time_native] input
  !    )
  !    
  !-
  INTERFACE day_in_year
    MODULE PROCEDURE day_in_year_days
    MODULE PROCEDURE day_in_year_intern
    MODULE PROCEDURE day_in_year_native
  END INTERFACE

  !--------------------------------------------------
  ! internal conversions

  PUBLIC  :: IMerge_YMD
  PRIVATE :: Split_YMD

  PUBLIC  :: IMerge_HMS
  PRIVATE :: Split_HMS
  PRIVATE :: Split_Seconds

  PRIVATE :: IConv_Sec2HMS
  PUBLIC  :: IMerge_HMS2Sec

CONTAINS

  !-------------------------------------------------------
  ! conversion between formats

  SUBROUTINE convert_day2int (i_day, o_day)
    TYPE (time_days),   INTENT(in)  :: i_day  ! day, second
    TYPE (time_intern), INTENT(out) :: o_day  ! ymd, hms

    TYPE (julian_date) :: julday
    TYPE (ly360_date) :: ly360day
    INTEGER            :: iy, im, id, isecs, iday
    REAL (dp)          :: zsecs

    IF (.NOT.i_day%init) CALL finish &
         ('mo_time_conversion:convert_day2int','date not initialised')

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      ! remapping of julian day (adjustment at 12UTC) 
      ! and day with adjustment at 00UTC
      !
      zsecs = sec2frac(i_day %second)
      iday  =          i_day %day
      julday%fraction = zsecs - 0.5_dp
      
      IF ((i_day %day < 0).AND. (zsecs > 0.5_dp)) THEN
        iday             = iday + 1
        julday %fraction = zsecs - 1.5_dp
      ELSE IF ((i_day %day > 0) .AND. (zsecs < 0.5_dp)) THEN
        iday             = iday - 1
        julday %fraction = zsecs + 0.5_dp
      END IF
      IF (iday < 0) THEN
        julday %day = AINT(REAL(iday,dp)-0.0001_dp)
      ELSE
        julday %day = AINT(REAL(iday,dp)+0.0001_dp)
      END IF
      CALL Set_JulianCalendar (julday, iy, im, id, isecs)
      
      o_day %ymd  = IMerge_YMD    (iy, im, id)
      o_day %hms  = IConv_Sec2HMS (isecs)
      o_day %init = .TRUE.
    CASE (CYL360)
      ly360day%fraction = sec2frac(i_day %second)
      ly360day%day      = i_day %day
      CALL Set_Ly360Calendar (ly360day, iy, im, id, isecs)

      o_day %ymd  = IMerge_YMD    (iy, im, id)
      o_day %hms  = IConv_Sec2HMS (isecs)
      o_day %init = .TRUE.
    END SELECT

  END SUBROUTINE convert_day2int


  SUBROUTINE convert_int2day (i_day, o_day)
    TYPE (time_intern), INTENT(in)  :: i_day  ! ymd, hms
    TYPE (time_days),   INTENT(out) :: o_day  ! day, second

    TYPE (julian_date) :: julday
    TYPE (ly360_date)  :: ly360day
    INTEGER            :: iy, im, id, ihr, imn, isec, kday
    REAL (dp)          :: zsec

    IF (.NOT.i_day%init) CALL finish &
         ('mo_time_conversion:convert_int2day','date not initialised')

    CALL Split_YMD (i_day %ymd, iy, im, id)
    CALL Split_HMS (i_day %hms, ihr, imn, isec)
    
    o_day %second = IMerge_HMS2Sec (ihr, imn, isec)

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      CALL Set_JulianDay (iy, im, id, o_day%second, julday)
      
      IF (julday %day < 0.0_dp) THEN
        kday = INT(julday %day-0.0001_dp)
      ELSE
        kday = INT(julday %day+0.0001_dp)
      END IF
      
      IF (julday %day < 0.0_dp) THEN
        IF (julday %fraction < -0.5_dp)       kday = kday - 1
        
      ELSE IF (julday %fraction > 0.0_dp) THEN
        IF (.NOT.(julday %fraction < 0.5_dp)) kday = kday + 1
        
      ELSE 
        IF (julday%fraction < -0.5_dp) THEN
          kday = kday - 1
        ELSE IF (.NOT.(julday%fraction < 0.5_dp)) THEN
          kday = kday + 1
        END IF
      END IF
      
      IF      (julday%fraction < -0.5_dp) THEN
        zsec = julday%fraction + 1.5_dp
      ELSE IF (julday%fraction < 0.5_dp) THEN
        zsec = julday%fraction + 0.5_dp
      ELSE
        zsec = julday%fraction - 0.5_dp
      END IF
      isec = frac2sec(zsec)

      o_day%day  = kday

    CASE (CYL360)
      CALL Set_Ly360Day (iy, im, id, o_day%second, ly360day)

      o_day%day = ly360day%day

    END SELECT

    o_day%init = .TRUE.

  END SUBROUTINE convert_int2day


  SUBROUTINE convert_nat2int (i_day, o_day)
    TYPE (time_native), INTENT(in)  :: i_day ! yr, mo, dy, hr, mn, se
    TYPE (time_intern), INTENT(out) :: o_day ! ymd, hms

    INTEGER :: yr, mo, dy, hr, mi, se

    CALL TC_get (i_day, yr, mo, dy, hr, mi, se)
    CALL TC_set (IMerge_YMD(yr, mo, dy), IMerge_HMS(hr, mi, se), o_day)

  END SUBROUTINE convert_nat2int


  SUBROUTINE convert_int2nat (i_day, o_day)
    TYPE (time_intern), INTENT(in)  :: i_day  ! ymd, hms
    TYPE (time_native), INTENT(out) :: o_day  ! yr, mo, dy, hr, mn, se

    INTEGER :: yr, mo, dy, hr, mi, se

    IF (.NOT.i_day%init) CALL finish &
         ('mo_time_conversion:convert_int2nat','date not initialised')

    CALL Split_YMD (i_day %ymd, yr, mo, dy)
    CALL Split_HMS (i_day %hms, hr, mi, se)

    CALL TC_set (yr, mo, dy, hr, mi, se, o_day)

  END SUBROUTINE convert_int2nat


  SUBROUTINE convert_nat2day (i_day, o_day)
    TYPE (time_native), INTENT(in)  :: i_day  ! yr, mo, dy, hr, mn, se
    TYPE (time_days),   INTENT(out) :: o_day  ! day, second

    TYPE (time_intern) :: ymdhms

    IF (.NOT.i_day%init) CALL finish &
         ('mo_time_conversion:convert_nat2day','date not initialised')

    CALL convert_nat2int (i_day, ymdhms)
    CALL convert_int2day (ymdhms, o_day)

  END SUBROUTINE convert_nat2day


  SUBROUTINE convert_day2nat (i_day, o_day)
    TYPE (time_days),   INTENT(in)  :: i_day ! day, second
    TYPE (time_native), INTENT(out) :: o_day ! yr, mo, dy, hr, mn, se

    TYPE (time_intern) :: ymdhms
    INTEGER            :: yr, mo, dy, hr, mi, se

    IF (.NOT.i_day%init) CALL finish &
         ('mo_time_conversion:convert_day2nat','date not initialised')

    CALL convert_day2int (i_day, ymdhms)

    CALL Split_YMD (ymdhms %ymd, yr, mo, dy)
    CALL Split_HMS (ymdhms %hms, hr, mi, se)

    CALL TC_set (yr, mo, dy, hr, mi, se, o_day)

  END SUBROUTINE convert_day2nat



  !-------------------------------------------------------------
  ! comparisons between formats
  ! both input dates must have the structure


  ! using unviversal calendar format (day,second)

  LOGICAL FUNCTION L_comp_days_lt (day1, day2)
    TYPE (time_days), INTENT(in) :: day1, day2

    IF ((.NOT.day1%init) .OR. (.NOT.day2%init)) CALL finish &
         ('mo_time_conversion:L_comp_days_lt','dates not initialised')

    IF ( (day1%day == day2%day .AND. day1%second >= day2%second) &
         .OR. (day1%day > day2%day) ) THEN
      L_comp_days_lt = .FALSE.
    ELSE
      L_comp_days_lt = .TRUE.
    END IF

  END FUNCTION L_comp_days_lt


  LOGICAL FUNCTION L_comp_days_eq (day1, day2)
    TYPE (time_days), INTENT(in) :: day1, day2

    IF ((.NOT.day1%init) .OR. (.NOT.day2%init)) CALL finish &
         ('mo_time_conversion:L_comp_days_eq','dates not initialised')

    IF ( day1%day == day2%day .AND. day1%second == day2%second ) THEN
      L_comp_days_eq = .TRUE.
    ELSE
      L_comp_days_eq = .FALSE.
    END IF

  END FUNCTION L_comp_days_eq


  LOGICAL FUNCTION L_comp_days_gt (day1, day2)
    TYPE (time_days), INTENT(in) :: day1, day2

    IF ((.NOT.day1%init) .OR. (.NOT.day2%init)) CALL finish &
         ('mo_time_conversion:L_comp_days_gt','dates not initialised')

    IF ( (day1%day == day2%day .AND. day1%second <= day2%second) &
         .OR. (day1%day < day2%day) ) THEN
      L_comp_days_gt = .FALSE.
    ELSE
      L_comp_days_gt = .TRUE.
    END IF

  END FUNCTION L_comp_days_gt


  ! using internal date format (ymd, hms)

  LOGICAL FUNCTION L_comp_intern_lt (ymdhms1, ymdhms2)
    TYPE (time_intern), INTENT(in) :: ymdhms1, ymdhms2

    IF ((.NOT.ymdhms1%init) .OR. (.NOT.ymdhms2%init)) CALL finish &
         ('mo_time_conversion:L_comp_intern_lt','dates not initialised')

    IF ( (ymdhms1%ymd == ymdhms2%ymd .AND. ymdhms1%hms >= ymdhms2%hms) &
         .OR. (ymdhms1%ymd > ymdhms2%ymd) ) THEN
      L_comp_intern_lt = .FALSE.
    ELSE
      L_comp_intern_lt = .TRUE.
    END IF

  END FUNCTION L_comp_intern_lt


  LOGICAL FUNCTION L_comp_intern_eq (ymdhms1, ymdhms2)
    TYPE (time_intern), INTENT(in) :: ymdhms1, ymdhms2

    IF ((.NOT.ymdhms1%init) .OR. (.NOT.ymdhms2%init)) CALL finish &
         ('mo_time_conversion:L_comp_intern_eq','dates not initialised')

    IF (ymdhms1%ymd == ymdhms2%ymd .AND. ymdhms1%hms == ymdhms2%hms) THEN
      L_comp_intern_eq = .TRUE.
    ELSE
      L_comp_intern_eq = .FALSE.
    END IF

  END FUNCTION L_comp_intern_eq


  LOGICAL FUNCTION L_comp_intern_gt (ymdhms1, ymdhms2)
    TYPE (time_intern), INTENT(in) :: ymdhms1, ymdhms2

    IF ((.NOT.ymdhms1%init) .OR. (.NOT.ymdhms2%init)) CALL finish &
         ('mo_time_conversion:L_comp_intern_gt','dates not initialised')

    IF ( (ymdhms1%ymd == ymdhms2%ymd .AND. ymdhms1%hms <= ymdhms2%hms) &
         .OR. (ymdhms1%ymd < ymdhms2%ymd) ) THEN
      L_comp_intern_gt = .FALSE.
    ELSE
      L_comp_intern_gt = .TRUE.
    END IF

  END FUNCTION L_comp_intern_gt


  ! using native format (yr, mo, dy, hr, mn, se)

  LOGICAL FUNCTION L_comp_native_lt (date1, date2)
    TYPE (time_native), INTENT(in) :: date1, date2

    TYPE (time_intern) :: ymdhms1, ymdhms2

    CALL convert_nat2int (date1, ymdhms1)
    CALL convert_nat2int (date2, ymdhms2)

    L_comp_native_lt = ymdhms1 < ymdhms2

  END FUNCTION L_comp_native_lt


  LOGICAL FUNCTION L_comp_native_eq (date1, date2)
    TYPE (time_native), INTENT(in) :: date1, date2

    TYPE (time_intern) :: ymdhms1, ymdhms2

    CALL convert_nat2int (date1, ymdhms1)
    CALL convert_nat2int (date2, ymdhms2)

    L_comp_native_eq = ymdhms1 == ymdhms2

  END FUNCTION L_comp_native_eq


  LOGICAL FUNCTION L_comp_native_gt(date1, date2)
    TYPE (time_native), INTENT(in) :: date1, date2

    TYPE (time_intern) :: ymdhms1, ymdhms2

    CALL convert_nat2int (date1, ymdhms1)
    CALL convert_nat2int (date2, ymdhms2)

    L_comp_native_gt = ymdhms1 > ymdhms2

  END FUNCTION L_comp_native_gt


  ! add a number of days and seconds to a given date and time

  SUBROUTINE add2days (days, seconds, my_day)
    INTEGER,          INTENT(in)    :: days, seconds
    TYPE (time_days), INTENT(inout) :: my_day

    REAL(dp)  :: zdayl
    INTEGER   :: idays, isecs

    IF (.NOT.my_day%init) CALL finish &
         ('mo_time_conversion:add2days','date not initialised')

    zdayl = day_len()
    isecs = seconds + my_day%second
    IF (isecs < 0) THEN
       idays = INT((REAL(isecs,dp)-0.001_dp)/zdayl)
    ELSE
       idays = INT((REAL(isecs,dp)+0.001_dp)/zdayl)
    END IF
    isecs = isecs - idays*IDAYLEN
    idays = my_day%day + days + idays

    IF (isecs < 0) THEN
       isecs = IDAYLEN + isecs
       idays = idays - 1
    END IF

    my_day %day    = idays
    my_day %second = isecs

  END SUBROUTINE add2days


  SUBROUTINE add2intern (days, seconds, my_intern)
    INTEGER,            INTENT(in)    :: days, seconds
    TYPE (time_intern), INTENT(inout) :: my_intern

    TYPE (time_days) :: my_day

    CALL convert_int2day (my_intern, my_day)
    CALL add2days    (days, seconds, my_day)
    CALL convert_day2int            (my_day, my_intern)

  END SUBROUTINE add2intern


  SUBROUTINE add2native (days, seconds, my_native)
    INTEGER,            INTENT(in)    :: days, seconds
    TYPE (time_native), INTENT(inout) :: my_native

    TYPE (time_days) :: my_day

    CALL convert_nat2day (my_native, my_day)
    CALL add2days    (days, seconds, my_day)
    CALL convert_day2nat            (my_day, my_native)

  END SUBROUTINE add2native

  !-------------------------------------------------------------
  ! get year and month length in days

  FUNCTION year_len_days (date) RESULT (ix)
    TYPE (time_days), INTENT(in) :: date
    INTEGER                      :: ix

    TYPE (time_native) :: date_nat

    CALL TC_convert  (date, date_nat)

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      ix = Get_JulianYearLen (date_nat%year)
    CASE (CYL360)
      ix = Get_LY360YearLen ()
    END SELECT

  END FUNCTION year_len_days


  FUNCTION year_len_intern (date) RESULT (ix)
    TYPE (time_intern), INTENT(in) :: date
    INTEGER                        :: ix

    TYPE (time_native) :: date_nat

    CALL TC_convert  (date, date_nat)

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      ix = Get_JulianYearLen (date_nat%year)
    CASE (CYL360)
      ix = Get_LY360YearLen ()
    END SELECT

  END FUNCTION year_len_intern


  FUNCTION year_len_native (date) RESULT (ix)
    TYPE (time_native), INTENT(in) :: date
    INTEGER                        :: ix

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      ix = Get_JulianYearLen (date%year)
    CASE (CYL360)
      ix = Get_LY360YearLen ()
    END SELECT

  END FUNCTION year_len_native


  REAL(dp) FUNCTION year_len_mean ()

    ! get mean tropical year length, PCMDI (as used in echam4.5)
    !
    ! year_len = 2450814. / 6710. = 365.25
    ! year_len_mean = 365.25

    IF (get_calendar_type() ==  CYL360) THEN
      year_len_mean = 360.0_dp 
    ELSE
      year_len_mean = 365.2422_dp 
    ENDIF

  END FUNCTION year_len_mean


  FUNCTION mon_len_days (date) RESULT (ix)
    TYPE (time_days), INTENT(in) :: date
    INTEGER                      :: ix

    TYPE (time_native) :: date_nat

    CALL TC_convert  (date,date_nat)

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      ix = Get_JulianMonLen (date_nat%year, date_nat%month)
    CASE (CYL360)
      ix = Get_Ly360MonLen ()
    END SELECT

  END FUNCTION mon_len_days


  FUNCTION mon_len_intern (date) RESULT (ix)
    TYPE (time_intern), INTENT(in) :: date
    INTEGER                        :: ix

    TYPE (time_native) :: date_nat

    CALL TC_convert  (date,date_nat)

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      ix = Get_JulianMonLen (date_nat%year, date_nat%month)
    CASE (CYL360)
      ix = Get_Ly360MonLen ()
    END SELECT

  END FUNCTION mon_len_intern


  FUNCTION mon_len_native (date) RESULT (ix)
    TYPE (time_native), INTENT(in) :: date
    INTEGER                        :: ix

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      ix = Get_JulianMonLen (date%year,date%month)
    CASE (CYL360)
      ix = Get_Ly360MonLen ()
    END SELECT

  END FUNCTION mon_len_native


  FUNCTION day_len () RESULT (ix)
    REAL(dp) :: ix

    ix = REAL(IDAYLEN,dp)

  END FUNCTION day_len

  !------------------------------------------------------------
  ! get day in the year

  FUNCTION day_in_year_days (date) RESULT (ix)
    TYPE (time_days), INTENT(in) :: date
    INTEGER                      :: ix

    TYPE (time_native) :: date_nat

    CALL TC_convert   (date, date_nat)
    ix = day_in_year_native (date_nat)

  END FUNCTION day_in_year_days

  FUNCTION day_in_year_intern (date) RESULT (ix)
    TYPE (time_intern), INTENT(in) :: date
    INTEGER                        :: ix

    TYPE (time_native) :: date_nat

    CALL TC_convert   (date, date_nat)
    ix = day_in_year_native (date_nat)

  END FUNCTION day_in_year_intern


  FUNCTION day_in_year_native (date) RESULT (ix)
    TYPE (time_native), INTENT(in) :: date
    INTEGER                        :: ix

    TYPE (time_days)   :: date_day
    TYPE (julian_date) :: julday
    TYPE (ly360_date)  :: ly360day
    INTEGER            :: yr, mo, dy, hr, mn, se, iday

    CALL TC_get     (date, yr, mo, dy, hr, mn, se)
    CALL TC_convert (date, date_day)
    CALL TC_get           (date_day, iday, se)

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      CALL Set_JulianDay    (yr, mo, dy,     se, julday)
      CALL Get_JulianYearDay                    (julday, yr, dy, se)
    CASE (CYL360)
      CALL Set_Ly360Day    (yr, mo, dy,     se, ly360day)
      CALL Get_Ly360YearDay                    (ly360day, yr, dy, se)
    END SELECT

    ix = dy

  END FUNCTION day_in_year_native

  !-------------------------------------------------------------
  ! print out the date/time information

  SUBROUTINE print_days (day, ptype, mess)
    TYPE (time_days),           INTENT(in)  :: day
    CHARACTER(len=*), OPTIONAL, INTENT(in)  :: ptype
    CHARACTER(len=*), OPTIONAL, INTENT(out) :: mess

    TYPE (time_native)       :: day_nat
    TYPE (time_intern)       :: day_int
    CHARACTER(len=STR_LEN_A) :: ctype
    LOGICAL                  :: lmess

    IF (.NOT.day%init) CALL finish &
         ('mo_time_conversion:print_days','date not initialised')

    ctype = ''      ; IF (PRESENT(ptype)) ctype = TRIM(ptype)
    lmess = .FALSE. ; IF (PRESENT(mess))  lmess = .TRUE.

    SELECT CASE(ctype)
    CASE(TC_PRN_NATIVE)
      CALL TC_convert (day, day_nat)
      IF (lmess) THEN
        CALL print_native  (day_nat, mess=mess)
      ELSE
        CALL print_native  (day_nat)
      END IF

    CASE(TC_PRN_INTERN)
      CALL TC_convert (day, day_int)
      IF (lmess) THEN
        CALL print_intern  (day_int, mess=mess)
      ELSE
        CALL print_intern  (day_int)
      END IF

    CASE default
      SELECT CASE (get_calendar_type())
      CASE (JULIAN)
        WRITE(message_text,'(a,i8,a,i8)') &
             'modified Julian day (00 UT adjusted): ', &
             day %day,' seconds: ',day %second
      CASE (CYL360)
        WRITE(message_text,'(a,i8,a,i8)') &
             '360 day year day (00 UT based): ', &
             day %day,' seconds: ',day %second
      END SELECT
      IF (lmess) THEN
        mess = TRIM(message_text)
      ELSE
        CALL message('',message_text)
      END IF

    END SELECT

  END SUBROUTINE print_days


  SUBROUTINE print_intern (ymdhms, ptype, mess)
    TYPE (time_intern),         INTENT(in)  :: ymdhms
    CHARACTER(len=*), OPTIONAL, INTENT(in)  :: ptype
    CHARACTER(len=*), OPTIONAL, INTENT(out) :: mess

    TYPE (time_days)         :: day_day
    TYPE (time_native)       :: day_nat
    CHARACTER(len=STR_LEN_A) :: ctype
    LOGICAL                  :: lmess

    IF (.NOT.ymdhms%init) CALL finish &
         ('mo_time_conversion:print_intern','date not initialised')

    ctype = ''      ; IF (PRESENT(ptype)) ctype = TRIM(ptype)
    lmess = .FALSE. ; IF (PRESENT(mess))  lmess = .TRUE.

    SELECT CASE(ctype)
    CASE(TC_PRN_DAYS)
      CALL TC_convert (ymdhms, day_day)
      IF (lmess) THEN
        CALL print_days       (day_day, mess=mess)
      ELSE
        CALL print_days       (day_day)
      END IF

    CASE(TC_PRN_NATIVE)
      CALL TC_convert (ymdhms, day_nat)
      IF (lmess) THEN
        CALL print_native     (day_nat,mess=mess)
      ELSE
        CALL print_native     (day_nat)
      END IF

    CASE default
      WRITE(message_text,'(a,i16,a,i6.6)') &
           'date (YYYYMMDD): ',ymdhms%ymd,' time (HHMMSS): ',ymdhms%hms
      IF (lmess) THEN
        mess = TRIM(message_text)
      ELSE
        CALL message('',message_text)
      END IF

    END SELECT

  END SUBROUTINE print_intern


  SUBROUTINE print_native (date, ptype, mess)
    TYPE (time_native),         INTENT(in)  :: date
    CHARACTER(len=*), OPTIONAL, INTENT(in)  :: ptype
    CHARACTER(len=*), OPTIONAL, INTENT(out) :: mess

    TYPE (time_days)         :: day_day
    TYPE (time_intern)       :: day_int
    CHARACTER(len=STR_LEN_A) :: ctype
    LOGICAL                  :: lmess

    IF (.NOT.date%init) CALL finish &
         ('mo_time_conversion:print_native','date not initialised')

    ctype = ''      ; IF (PRESENT(ptype)) ctype = TRIM(ptype)
    lmess = .FALSE. ; IF (PRESENT(mess))  lmess = .TRUE.

    SELECT CASE(ctype)
    CASE(TC_PRN_DAYS)
      CALL TC_convert (date, day_day)
      IF (lmess) THEN
        CALL print_days     (day_day,mess=mess)
      ELSE
        CALL print_days     (day_day)
      END IF

    CASE(TC_PRN_INTERN)
      CALL TC_convert (date, day_int)
      IF (lmess) THEN
        CALL print_intern   (day_int,mess=mess)
      ELSE
        CALL print_intern   (day_int)
      END IF

    CASE default
      WRITE(message_text,'(i2,a,a3,a,i8,3(a,i2.2))') &
               date %day,'. ',CMONTHS(date %month),' ',date %year,&
           ' ',date %hour,':',date %minute,':',date %second
      IF (lmess) THEN
        mess = TRIM(message_text)
      ELSE
        CALL message('',message_text)
      END IF
    END SELECT

  END SUBROUTINE print_native



  !--------------------------------------------------
  ! check and store data into structures

  SUBROUTINE set_days (kday, ktime, my_day)
    INTEGER,          INTENT(in)  :: kday, ktime
    TYPE (time_days), INTENT(out) :: my_day

    IF (ktime < 0 .OR. (IDAYLEN-1) < ktime) CALL finish &
         ('mo_time_conversion:set_days','daytime outside valid range')

    my_day %day    = kday
    my_day %second = ktime
    my_day %init   = .TRUE.

  END SUBROUTINE set_days


  SUBROUTINE set_intern (kymd, khms, my_ymdhms)
    INTEGER,            INTENT(in)  :: kymd, khms
    TYPE (time_intern), INTENT(out) :: my_ymdhms

    TYPE (time_days) :: my_day
    INTEGER          :: ihr, imn, ise

    CALL Split_HMS (khms, ihr, imn, ise)

    IF (ihr < 0 .OR. 23 < ihr) THEN
      WRITE(message_text,*) 'hour [',ihr,'] outside valid range'
      CALL finish('mo_time_conversion:set_intern',message_text)

    ELSE IF (imn < 0 .OR. 59 < imn) THEN
      WRITE(message_text,*) 'minute [',imn,'] outside valid range'
      CALL finish('mo_time_conversion:set_intern',message_text)

    ELSE IF (ise < 0 .OR. 59 < ise) THEN
      WRITE(message_text,*) 'seconds [',ise,'] outside valid range'
      CALL finish('mo_time_conversion:set_intern',message_text)

    END IF

    my_ymdhms %hms  = IMerge_HMS (ihr, imn, ise)
    my_ymdhms %ymd  = kymd
    my_ymdhms %init = .TRUE.

    ! optional check of correct date/time for the selected calendar model
    CALL convert_int2day (my_ymdhms, my_day)

  END SUBROUTINE set_intern


  SUBROUTINE set_native (kyr, kmo, kdy, khr, kmn, kse, my_date)
    INTEGER,            INTENT(in)  :: kyr, kmo, kdy, khr, kmn, kse
    TYPE (time_native), INTENT(out) :: my_date

    TYPE (time_intern) :: my_ymdhms

    IF (kmo < 1 .OR. 12 < kmo) THEN
      WRITE(message_text,*) 'month number [',kmo,'] outside valid range'
      CALL finish('mo_time_conversion:set_native',message_text)

    ELSE IF (kdy < 1 .OR. 31 < kdy) THEN
      WRITE(message_text,*) 'day number [',kdy,'] outside valid range'
      CALL finish('mo_time_conversion:set_native',message_text)

    END IF

    ! check the calendar consistency
    CALL set_intern ( &
         IMerge_YMD (kyr, kmo, kdy), &
         IMerge_HMS (khr, kmn, kse), my_ymdhms)

    my_date %year   = kyr
    my_date %month  = kmo
    my_date %day    = kdy
    my_date %hour   = khr
    my_date %minute = kmn
    my_date %second = kse
    my_date %init   = .TRUE.

  END SUBROUTINE set_native

  !------------------------------------------------------------
  ! get components of structures

  SUBROUTINE get_days (date, day, second)
    TYPE (time_days),  INTENT(in)  :: date
    INTEGER, OPTIONAL, INTENT(out) :: day, second

    IF (.NOT.date%init) CALL finish &
         ('mo_time_conversion:get_days','date not initialised')

    IF (PRESENT(day))    day    = date%day
    IF (PRESENT(second)) second = date%second

  END SUBROUTINE get_days


  SUBROUTINE get_intern (date, ymd, hms)
    TYPE (time_intern), INTENT(in)  :: date
    INTEGER, OPTIONAL,  INTENT(out) :: ymd, hms

    IF (.NOT.date%init) CALL finish &
         ('mo_time_conversion:get_intern','date not initialised')

    IF (PRESENT(ymd)) ymd = date %ymd
    IF (PRESENT(hms)) hms = date %hms

  END SUBROUTINE get_intern


  SUBROUTINE get_native (date, &
       year, month, day, hour, minute, second)
    TYPE (time_native), INTENT(in)  :: date
    INTEGER, OPTIONAL,  INTENT(out) :: &
         year, month, day, hour, minute, second

    IF (.NOT.date%init) CALL finish &
         ('mo_time_conversion:get_native','date not initialised')

    IF (PRESENT(year))   year   = date %year
    IF (PRESENT(month))  month  = date %month
    IF (PRESENT(day))    day    = date %day
    IF (PRESENT(hour))   hour   = date %hour
    IF (PRESENT(minute)) minute = date %minute
    IF (PRESENT(second)) second = date %second

  END SUBROUTINE get_native

  !-------------------------------------------------------------
  ! simple conversions, this part should be only private
  ! independend from the calendar


  FUNCTION IMerge_YMD (ky, km, kd) RESULT (ix)
    !+
    !
    ! IMerge_YMD [function, integer]
    !    merge year, month and day
    !    (
    !    year  [integer] input (year of calendar)
    !    month [integer] input (month of year)
    !    day   [integer] input (day of month)
    !    )
    !
    !-
    INTEGER, INTENT(in) :: ky, km, kd
    INTEGER             :: ix

    INTEGER :: isign, iyr

    isign = SIGN(1,INT(ky))
    iyr   = isign*ky

    ix    = isign*(iyr*10000 + km*100 + kd)

  END FUNCTION IMerge_YMD



  SUBROUTINE Split_YMD (ymd, ky, km, kd)
    !+
    !
    ! Split_YMD [subroutine] (private)
    !   split internal date into year, month and day
    !   (
    !   yyyymmdd  [integer] input (internal date)
    !   year      [integer] output (year of calendar)
    !   month     [integer] output (month of year)
    !   day       [integer] output (day of month)
    !   )
    !
    !-
    INTEGER, INTENT(in)  :: ymd
    INTEGER, INTENT(out) :: ky, km, kd

    INTEGER:: isign, iymd

    isign = SIGN(1,INT(ymd))
    iymd  = isign*ymd

    ky = isign*iymd/10000
    km = MOD(iymd,10000)/100
    kd = MOD(iymd,100)

  END SUBROUTINE Split_YMD



  FUNCTION IMerge_HMS (ih, im, is) RESULT (ix)
    !+
    !
    ! IMerge_HMS [function, integer]
    !   merge hour, minute and second to internal time
    !   (
    !   hour   [integer] input (hour of day, 0...23)
    !   minute [integer] input (minute of hour, 0...59)
    !   second [integer] input (second of minute, 0...59)
    !   )
    !
    !-
    INTEGER, INTENT(in) :: ih, im, is
    INTEGER             :: ix

    ix = ih*10000 + im*100 + is

  END FUNCTION IMerge_HMS



  SUBROUTINE Split_HMS (hms, kh, km, ks)
    !+
    !
    ! Split_HMS [subroutine] (private)
    !   split internal time into hour, minute and second
    !   (
    !   hhmmss [integer] input (internal time)
    !   hour   [integer] output (hour of day, 0...23)
    !   minute [integer] output (minute of hour, 0...59)
    !   second [integer] output (second of minute, 0...59)
    !   )
    !
    !-
    INTEGER, INTENT(in)  :: hms
    INTEGER, INTENT(out) :: kh, km, ks

    kh = hms/10000
    km = MOD(hms,10000)/100
    ks = MOD(hms,100)

  END SUBROUTINE Split_HMS



  SUBROUTINE Split_Seconds (ksecs, kh, km, ks)
    !+
    !
    ! Split_Seconds [subroutine] (private)
    !   split seconds of day into hour, minute and second
    !   (
    !   seconds  [integer] input (seconds of day, 0...86399)
    !   hour     [integer] output (hour of day, 0...23)
    !   minute   [integer] output (minute of hour, 0...59)
    !   second   [integer] output (second of minute, 0,...59)
    !   )
    !
    !-
    INTEGER, INTENT(in)  :: ksecs
    INTEGER, INTENT(out) :: kh, km, ks

    kh = ksecs/3600
    km = MOD(ksecs,3600)/60
    ks = MOD(ksecs,60)

  END SUBROUTINE Split_Seconds



  FUNCTION IConv_Sec2HMS (ksecs) RESULT (ix)
    !+
    !
    ! IConv_Sec2HMS [function, integer] (private)
    !   convert seconds of day into internal time
    !   (
    !   seconds [integer] input (seconds of day)
    !   )
    !
    !-
    INTEGER, INTENT(in)  :: ksecs
    INTEGER              :: ix

    INTEGER :: ih, im, is

    CALL Split_Seconds (ksecs, ih, im, is)
    ix = IMerge_HMS           (ih, im, is)

  END FUNCTION IConv_Sec2HMS



  FUNCTION IMerge_HMS2Sec (kh, km, ks) RESULT (ix)
    !+
    !
    ! IMerge_HMS2Sec [function, integer]
    !    merge hour, minute and second to internal time
    !    (
    !    hour   [integer] input (hour of day, 0...23)
    !    minute [integer] input (minute of hour, 0...59)
    !    second [integer] input (second of minute, 0...59)
    !    )
    !
    !-
    INTEGER, INTENT(in) :: kh, km, ks
    INTEGER             :: ix

    ix = kh*3600 + km*60 + ks

  END FUNCTION IMerge_HMS2Sec

END MODULE mo_time_conversion
