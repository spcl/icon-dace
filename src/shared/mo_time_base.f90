! mo_time_base [module]
!    routines to handle Julian dates and times
!
!
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

MODULE mo_time_base



!!!#define IAU
  !
  ! References:
  !
  ! Montenbruck, Oliver, "Practical Ephmeris Calculations", Ch. 2, pp 33-40.
  ! The 1992 Astronomical Almanac, page B4.
  !
  ! The Julian Date is defined as the the number of days which have elapsed
  ! since the 1st of January of the year 4713 BC 12:00 Universal Time.
  !
  ! Up to 4th October 1582 AD the Julian Calendar was in force with leap
  ! years every four years, but from then on the Gregorian Calendar carried
  ! on from 15th October 1582 with the leap years defined by the rule:
  !
  ! Leap year is every year whose yearly number is divisible by four, but
  ! not by a hundred, or is divisible by four hundred.
  !
  ! At midday on 4th October 1582 2299160.5 Julian days had elapsed.
  ! The Gregorian Calendar then carried on at this point from 15th October
  ! 1582, so its beginning occured on the Julian date 2299160.5.
  !
  ! Note: the astronomical year -4712 corresponds to the year 4713 BC, the
  !       year 0 to the year 1 BC; thereafter the astronomical year match
  !       the year AD, e.g. year 1 = year 1 AD.
  !
  !       This routines work for the years -5877402 BC until 5868098 AD.
  !       However, there are not covered dates by the current GRIB edition 1,
  !       which can handle dates between 1 AD and 25599 AD.
  !
  !       The "Modified Julian Date" is the Julian Date - 2400000.5 and so
  !       has a zero point of 17th November 1858 AD 00:00 Universal Time.
  !
  !       Due to the small time span coverable by the GRIB output there is
  !       no need to use a "Modified Julian Date".
  !
  ! The Julian day number is stored as two doubles to guarantee sufficient
  ! precision on all machines. The complete value is (day+fraction)
  ! although this addition will sometimes lose precision. Note that day
  ! might not be an integer value and fraction might be greater than one.
  !
  !----------------------------------------------------------------------

  USE mo_kind,      ONLY: wp, i8
  USE mo_exception, ONLY: finish, message

  IMPLICIT NONE

  PRIVATE

  !----------------------------------------------------------------------
  ! Subroutines/Functions
  !
  ! Note: negative years are valid
  !
  PUBLIC :: set_calendar_type
  PUBLIC :: get_calendar_type
  !
  PUBLIC :: set_julian_day      ! convert calendar date/time into Julian day
  PUBLIC :: set_julian_calendar ! convert Julian day into calendar date/time
  PUBLIC :: get_julian_yearday  ! select the day number of the year
  PUBLIC :: get_julian_yearlen  ! get the length of the year
  PUBLIC :: get_julian_monlen   ! get the length of the months
  PUBLIC :: print_julian_day    ! print out the date information
  !
  PUBLIC :: set_ly360_day       ! convert calendar date/time into Ly360 day
  PUBLIC :: set_ly360_calendar  ! convert Ly360 day into calendar date/time
  PUBLIC :: get_ly360_yearday   ! select the day number of the year
  PUBLIC :: get_ly360_yearlen   ! get the length of the year
  PUBLIC :: get_ly360_monlen    ! get the length of the months
  PUBLIC :: print_ly360_day     ! print out the date information
  !
  PUBLIC :: sec2frac, frac2sec  ! convert sec to fraction and vice versa
  PUBLIC :: add_time            ! add day, second to Julian date
  !
  !----------------------------------------------------------------------
  ! Structure declarations, variables
  !
  TYPE, PUBLIC :: t_julian_date

    ! julian_date [structure]
    !     day      [real_wp]  (day in calendar)
    !     fraction [real_wp]  (fraction of the day)
    !

    REAL(wp) :: day       ! day since 1. January 4713 BC, 12UTC
    REAL(wp) :: fraction  ! fraction of day
  END TYPE t_julian_date

  TYPE, PUBLIC :: t_ly360_date

    ! t_ly360_date [structure]
    !     day      [integer]  (day in calendar)
    !     fraction [real_wp]  (fraction of the day)
    !

    INTEGER(i8) :: day       ! days since start of experiment
    REAL(wp)    :: fraction  ! fraction of day
  END TYPE t_ly360_date

  ! idaylen [integer, constant]
  !    length of a day in seconds
  !

  INTEGER,  PUBLIC, PARAMETER :: idaylen = 86400
!!$  REAL(wp), PUBLIC, PARAMETER :: rdaylen = REAL(idaylen,wp)  ! cannot be used with NEC compiler
  REAL(wp), PUBLIC, PARAMETER :: rdaylen = 86400._wp

  CHARACTER(len=256) :: message_text

  INTEGER, PUBLIC, PARAMETER :: julian = 0   ! Julian day based full fledged
  INTEGER, PUBLIC, PARAMETER :: cyl360 = 1   ! constant year length of 360 days

  ! used calendar predefined as julian

  INTEGER, SAVE :: used_calendar = julian

CONTAINS

  SUBROUTINE set_calendar_type (calendar)

    INTEGER, INTENT(in) :: calendar

    used_calendar = calendar

  END SUBROUTINE set_calendar_type

  INTEGER FUNCTION get_calendar_type ()

    get_calendar_type = used_calendar

  END FUNCTION get_calendar_type

  !----------------------------------------------------------------------
  !

  SUBROUTINE set_julian_day(ky, km, kd, ksec, my_day)

    ! set_julian__day  [subroutine]
    !    convert year, month, day, seconds into Julian calendar day
    !    (
    !    year   [integer] input (calendar year)
    !    month  [integer] input (month of the year)
    !    day    [integer] input (day of the month)
    !    second [integer] input (seconds of the day)
    !    date   [julian_date] output (Julian day)
    !    )
    !

    INTEGER, INTENT(IN) :: ky
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: kd
    INTEGER, INTENT(IN) :: ksec
    !
    ! for reference: 1. January 1998 00 UTC === Julian Day 2450814.5
    !
    TYPE (t_julian_date), INTENT(out) :: my_day

    INTEGER :: ib, iy, im, idmax
    REAL(wp) :: zd, zsec

    IF ( ksec < 0 .OR. ksec > 86399) THEN
      CALL finish ('mo_time_base:set_julian_day', 'invalid number of seconds')
    ENDIF
    zsec = sec2frac(ksec)

    IF (km <= 2) THEN
      iy = ky-1
      im = km+12
    ELSE
      iy = ky
      im = km
    ENDIF

    IF (iy < 0) THEN
       ib = INT((iy+1)/400)-INT((iy+1)/100)
    ELSE
       ib = INT(iy/400)-INT(iy/100)
    END IF


    ! check the length of the month
    idmax = get_julian_monlen (ky, km)

    IF (kd < 1 .OR. idmax < kd) &
         CALL finish ('mo_time_base:set_julian_day' ,'day in months invalid')

    zd = REAL(FLOOR(365.25_wp*REAL(iy,wp)),wp)+REAL(INT(30.6001_wp*REAL(im+1,wp)),wp) &
        +REAL(ib,wp)+1720996.5_wp+REAL(kd,wp)+zsec

    my_day %day      = AINT(zd,wp)
    my_day %fraction = zd-AINT(zd,wp)

  END SUBROUTINE set_julian_day

  SUBROUTINE set_julian_calendar(my_day, ky, km, kd, ksec)

    ! set_julian_calendar [subroutine]
    !    convert Julian date into year, months, day, seconds of day
    !    (
    !    date   [julian_date]  input (Julian day)
    !    year   [integer] output (calendar year)
    !    month  [integer] output (month of the year)
    !    day    [integer] output (day of the month)
    !    second [integer] output (seconds of the day)
    !    )
    !

    TYPE (t_julian_date), INTENT(IN) :: my_day
    INTEGER, INTENT(OUT)           :: ky
    INTEGER, INTENT(OUT)           :: km
    INTEGER, INTENT(OUT)           :: kd
    INTEGER, INTENT(OUT)           :: ksec

    REAL(wp) :: za, zb, zc, zd, ze, zf, zday, zfrac, zsec

    zday  = my_day%day
    zfrac = my_day%fraction

    za = REAL(FLOOR(zday+zfrac+0.5_wp),wp)

    zb = REAL(FLOOR((za-1867216.25_wp)/36524.25_wp),wp)
    zc = za+zb-REAL(FLOOR(zb/4._wp),wp)+1525._wp


    zd = REAL(FLOOR((zc-122.1_wp)/365.25_wp),wp)
    ze = REAL(FLOOR(365.25_wp*zd),wp)
    zf = REAL(FLOOR((zc-ze)/30.6001_wp),wp)

    kd = INT(zc-ze- REAL(FLOOR(30.6001_wp*zf),wp))
    km = INT(zf-REAL(1+12*FLOOR(zf/REAL(14,wp)),wp))
    ky = INT(zd-4715._wp-REAL((7+km)/10,wp))

    ! convert fraction in seconds of the day
    IF (zfrac < -0.5_wp) THEN
       zsec = 1.5_wp + zfrac
    ELSE IF (zfrac < 0.5_wp) THEN
       zsec = 0.5_wp + zfrac
    ELSE
       zsec = zfrac - 0.5_wp
    END IF
    ksec = frac2sec(zsec)

  END SUBROUTINE set_julian_calendar

  SUBROUTINE get_julian_yearday(my_day, kjy, kd, ksec)

    ! get_julian_yearday [subroutine]
    !    get Julian year, day of year and seconds of day from Julian date
    !    (
    !    date   [julian_date]  input (Julian day)
    !    year   [integer] output (Julian calendar year)
    !    day    [integer] output (day of the year)
    !    second [integer] output (seconds of the day)
    !    )
    !
    ! Day of Astronomical Year (approximated)
    !
    ! it is used for the correct annual cycle
    ! remember the mismatch between calendar and
    ! astronomical year accumulated until 1582
    !
    ! this version gives adjustment until 1583
    !
    TYPE (t_julian_date), INTENT(IN) :: my_day
    INTEGER, INTENT(OUT)           :: kjy
    INTEGER, INTENT(OUT)           :: kd
    INTEGER, INTENT(out)           :: ksec

    TYPE (t_julian_date) :: first_jan, present_day
    INTEGER :: iy, id, im

    ! AD 1 is year 4713 in Julian calendar, 31st of December 1997 the Julian
    ! day 2450814.0 is elapsed at 12 UTC and the Julian year is 6711, so:
    ! year_len = 2450814. / 6710. = 365.25

    ! get the date of the present julian day
    CALL set_julian_calendar(my_day, iy, im, id, ksec)

    ! find the first of January 00UTC
    CALL set_julian_day(iy,  1,  1, 0, first_jan)

    ! find the present Julian Day at 00UTC
    CALL set_julian_day(iy, im, id, 0, present_day)

    ! the day in the year is given
    kd = INT(present_day%day-first_jan%day)+1

    ! get the Julian year
    kjy = iy+4712

  END SUBROUTINE get_julian_yearday

  INTEGER FUNCTION get_julian_yearlen(ky)

    ! get_julian_yearlen  [function, integer]
    !    get the length of a Julian year in days
    !    (
    !    year [integer] input (Calendar year)
    !    )
    !

    INTEGER, INTENT(in) :: ky

    INTEGER :: ylen

    IF ( (MOD(ky,4)==0 .AND. MOD(ky,100)/=0) .OR. MOD(ky,400)==0 ) THEN
       ylen = 366
    ELSE
       ylen = 365
    END IF
    get_julian_yearlen = ylen

  END FUNCTION get_julian_yearlen

  INTEGER FUNCTION get_julian_monlen(ky, km)

    ! get_julian_monlen [function, integer]
    !    get the length of a months in a Julian year
    !    (
    !    year  [integer] input (Calendar year)
    !    month [integer] input (month of the year)
    !    )
    !

    INTEGER, INTENT(in) :: km, ky

    INTEGER :: idmax

    SELECT CASE(km)
    CASE(1,3,5,7,8,10,12)
      idmax = 31
    CASE(4,6,9,11)
      idmax = 30
    CASE(2)
      IF ( (MOD(ky,4)==0 .AND. MOD(ky,100)/=0) .OR. MOD(ky,400)==0 ) THEN
        ! leap year found
        idmax = 29
      ELSE
        idmax = 28
      END IF

    CASE default
      idmax = 0
      CALL finish('mo_time_base:get_julian_monlen', 'month invalid')

    END SELECT
    get_julian_monlen = idmax

  END FUNCTION get_julian_monlen

  SUBROUTINE print_julian_day(my_day)

    ! print_julian_day [subroutine] interface [print]
    !   print Julian/model day on standard output
    !   (
    !   date [julian_date] input (model/Julian calendar day)
    !   )
    !

    TYPE (t_julian_date), INTENT(in) :: my_day

    WRITE(message_text,*) my_day%day+my_day%fraction
    CALL message('print_julian_day', message_text)

  END SUBROUTINE print_julian_day

  SUBROUTINE set_ly360_day(ky, km, kd, ksec, my_day)

    ! set_ly360_day  [subroutine]
    !    convert year, month, day, seconds into Ly360 calendar day
    !    (
    !    year   [integer] input (calendar year)
    !    month  [integer] input (month of the year)
    !    day    [integer] input (day of the month)
    !    second [integer] input (seconds of the day)
    !    date   [t_ly360_date] output (Ly360 day)
    !    )
    !

    INTEGER, INTENT(IN) :: ky
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: kd
    INTEGER, INTENT(IN) :: ksec
    !
    TYPE (t_ly360_date), INTENT(out) :: my_day

    IF ( ksec < 0 .OR. ksec > 86399) THEN
      CALL finish ('mo_time_base:set_ly360_day', 'invalid number of seconds')
    ENDIF

    my_day %day      = INT(360*ky+(km-1)*30+(kd-1),i8)
    my_day %fraction = sec2frac(ksec)


  END SUBROUTINE set_ly360_day

  SUBROUTINE set_ly360_calendar(my_day, ky, km, kd, ksec)

    ! set_ly360_calendar [subroutine]
    !    convert Ly360 date into year, months, day, seconds of day
    !    (
    !    date   [t_ly360_date]  input (Ly360 day)
    !    year   [integer] output (calendar year)
    !    month  [integer] output (month of the year)
    !    day    [integer] output (day of the month)
    !    second [integer] output (seconds of the day)
    !    )
    !

    TYPE (t_ly360_date), INTENT(IN) :: my_day
    INTEGER, INTENT(OUT)            :: ky
    INTEGER, INTENT(OUT)            :: km
    INTEGER, INTENT(OUT)            :: kd
    INTEGER, INTENT(OUT)            :: ksec

    INTEGER :: ir

    IF (my_day%day < 0_i8) THEN
      ky = INT(my_day%day/360_i8)-1
      ir  = INT(MOD(my_day%day, 360_i8))
      km = 12-ir/30
      kd = 30-(MOD(ir,30)+1)
    ELSE
      ky = INT(my_day%day/360_i8)
      ir  = INT(MOD(my_day%day, 360_i8))
      km = ir/30+1
      kd = MOD(ir,30)+1
    ENDIF

    ksec = frac2sec(my_day%fraction)


  END SUBROUTINE set_ly360_calendar

  SUBROUTINE get_ly360_yearday(my_day, ky, kd, ksec)

    ! get_ly360_yearday [subroutine]
    !    get Ly360 year, day of year and seconds of day from Ly360 date
    !    (
    !    date   [t_ly360_date]  input (Ly360 day)
    !    year   [integer] output (Ly360 calendar year)
    !    day    [integer] output (day of the year)
    !    second [integer] output (seconds of the day)
    !    )
    !
    ! Day of Astronomical Year (approximated)
    !
    ! it is used for the correct annual cycle
    ! remember the mismatch between calendar and
    ! astronomical year accumulated until 1582
    !
    ! this version gives adjustment until 1583
    !
    TYPE (t_ly360_date), INTENT(IN) :: my_day
    INTEGER, INTENT(OUT)            :: ky
    INTEGER, INTENT(OUT)            :: kd
    INTEGER, INTENT(out)            :: ksec

    ky = INT(my_day%day/360_i8)
    kd  = INT(MOD(my_day%day, 360_i8))

    ksec = frac2sec(my_day%fraction)


  END SUBROUTINE get_ly360_yearday

  INTEGER FUNCTION get_ly360_yearlen()

    ! get_ly360_yearlen  [function, integer]
    !    get the length of a Ly360 year in days
    !    (
    !    )
    !

    get_ly360_yearlen = 360

  END FUNCTION get_ly360_yearlen

  INTEGER FUNCTION get_ly360_monlen()

    ! get_ly360_monlen [function, integer]
    !    get the length of a months in a Ly360 year
    !    (
    !    )
    !

    get_ly360_monlen = 30

  END FUNCTION get_ly360_monlen

  SUBROUTINE print_ly360_day(my_day)

    ! print_ly360_day [subroutine] interface [print]
    !   print Ly360/model day on standard output
    !   (
    !   date [t_ly360_date] input (model/Ly360 calendar day)
    !   )
    !

    TYPE (t_ly360_date), INTENT(in) :: my_day

    WRITE(message_text,*) REAL(my_day%day,wp)+my_day%fraction
    CALL message('print_ly360_day', message_text)

  END SUBROUTINE print_ly360_day

  FUNCTION sec2frac(isec) RESULT (zfrac)

    ! sec2frac [function, real_wp]
    !    convert seconds of day into fraction
    !    (
    !    second [integer] input (seconds of the day)
    !    )
    !

    INTEGER , INTENT(in) :: isec
    REAL(wp) :: zfrac

    SELECT CASE (get_calendar_type())
    CASE (julian)
      zfrac = (REAL(isec,wp)+0.000001_wp)/rdaylen
    CASE (cyl360)
      zfrac = REAL(isec,wp)/rdaylen
    CASE DEFAULT
      zfrac = 0.0_wp
      CALL finish('mo_time_base:sec2frac','calendar_type must be julian or cyl360')
    END SELECT

  END FUNCTION sec2frac

  FUNCTION frac2sec(zfrac) RESULT (isec)

    ! frac2sec [function, integer]
    !    convert fraction of day into seconds
    !    (
    !    fraction [real_wp] input (fraction of the day)
    !    )
    !

    REAL(wp), INTENT(in) :: zfrac
    INTEGER :: isec

    isec = INT(zfrac*REAL(IDAYLEN,wp))

  END FUNCTION frac2sec

  SUBROUTINE add_time (days, seconds, my_day)
    INTEGER,              INTENT(in)    :: days, seconds
    TYPE (t_julian_date), INTENT(inout) :: my_day

    REAL(wp)  :: zdayl
    INTEGER   :: idays, isecs

!write (0,*) '1: ', days, seconds
!write (0,*) '1: ', my_day%day, my_day%fraction

    zdayl = rdaylen
    isecs = frac2sec(my_day%fraction) + seconds
    IF (isecs < 0) THEN
       idays = INT((REAL(isecs,wp)-0.001_wp)/zdayl)
    ELSE
       idays = INT((REAL(isecs,wp)+0.001_wp)/zdayl)
    END IF
    isecs = isecs - idays*idaylen
    idays = INT(my_day%day) + days + idays

    IF (isecs < 0) THEN
       isecs = idaylen + isecs
       idays = idays - 1
    END IF

!write (0,*) '1: ', idays, isecs

    my_day%day      = REAL(idays,wp)
    my_day%fraction = sec2frac(isecs)

!write (0,*) '1: ', my_day%day, my_day%fraction

  END SUBROUTINE add_time

END MODULE mo_time_base
