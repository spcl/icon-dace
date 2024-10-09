! mo_time_base.f90 - Routines to handle Julian dates and times
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

#ifdef __xlC__
@PROCESS STRICT
#endif
MODULE mo_time_base
  !+
  !
  ! mo_time_base [module]
  !    routines to handle Julian dates and times
  !
  ! Version:    E5/R1.07+  04-February-2003
  !
  ! Authors:
  !   L. Kornblueh, MPI                prepared basic version
  !   I. Kirchner,  MPI July 2000      large revision
  !   L. Kornblueh, MPI August 2000    adapting for right real kind (dp)
  !                                    and other minor corrections, removed
  !                                    30 day support, removed print overload.
  !   L. Kornblueh, MPI February 2003  adaption to PRISM calendar
  !                                    dates before the change to the gregorian
  !                                    calendar are treated as gregorian as 
  !                                    well to have a consitent view on the 
  !                                    seasons in climte sense - introduce
  !                                    preprocessor flag IAU for the original
  !                                    code (IAU - International Astrophysical
  !                                    Union).
  !   I. Kirchner,  FUB February 2003  code review
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.
  !
  !-
#undef IAU
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
 
  USE mo_kind,      ONLY: dp, i8
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
  PUBLIC :: Set_JulianDay      ! convert calendar date/time into Julian day
  PUBLIC :: Set_JulianCalendar ! convert Julian day into calendar date/time
  PUBLIC :: Get_JulianYearDay  ! select the day number of the year
  PUBLIC :: Get_JulianYearLen  ! get the length of the year
  PUBLIC :: Get_JulianMonLen   ! get the length of the months
  PUBLIC :: Print_JulianDay    ! print out the date information
  PUBLIC :: Set_Ly360Day       ! convert calendar date/time into Ly360 day
  PUBLIC :: Set_Ly360Calendar  ! convert Ly360 day into calendar date/time
  PUBLIC :: Get_Ly360YearDay   ! select the day number of the year
  PUBLIC :: Get_Ly360YearLen   ! get the length of the year
  PUBLIC :: Get_Ly360MonLen    ! get the length of the months
  PUBLIC :: Print_Ly360Day     ! print out the date information
  PUBLIC :: sec2frac, frac2sec ! convert sec to fraction and vice versa

  !----------------------------------------------------------------------
  ! Structure declarations, variables
  !
  TYPE, PUBLIC :: julian_date
    !+
    !
    ! julian_date [structure]
    !     day      [real_dp]  (day in calendar)
    !     fraction [real_dp]  (fraction of the day)
    !
    !-
    REAL(dp) :: day       ! day since 1. January 4713 BC, 12UTC
    REAL(dp) :: fraction  ! fraction of day
  END TYPE julian_date

  TYPE, PUBLIC :: Ly360_date
    !+
    !
    ! Ly360_date [structure]
    !     day      [real_dp]  (day in calendar)
    !     fraction [real_dp]  (fraction of the day)
    !
    !-
    INTEGER(i8) :: day       ! days since start of experiment
    REAL(dp)    :: fraction  ! fraction of day
  END TYPE Ly360_date

  !+
  !
  ! idaylen [integer, constant]
  !    length of a day in seconds
  !
  !-
  INTEGER, PUBLIC, PARAMETER :: IDAYLEN = 86400

  CHARACTER(len=256) :: message_text

  INTEGER, PUBLIC, PARAMETER :: JULIAN = 0   ! Julian day based full fledged 
  INTEGER, PUBLIC, PARAMETER :: CYL365 = 1   ! constant year length of 365 days
  INTEGER, PUBLIC, PARAMETER :: CYL360 = 2   ! constant year length of 360 days

  ! used calendar predefined as JULIAN

  INTEGER, SAVE :: used_calendar = JULIAN  

CONTAINS

  SUBROUTINE set_calendar_type (ktype)

    INTEGER, INTENT(in) :: ktype

    used_calendar = ktype

  END SUBROUTINE set_calendar_type
  
  INTEGER FUNCTION get_calendar_type ()

    get_calendar_type = used_calendar

  END FUNCTION get_calendar_type

  !----------------------------------------------------------------------
  !

  SUBROUTINE Set_JulianDay(ky, km, kd, ksec, my_day, lperpetual_year)
    !+
    !
    ! Set_JulianDay  [subroutine]
    !    convert year, month, day, seconds into Julian calendar day
    !    (
    !    year   [integer] input (calendar year)
    !    month  [integer] input (month of the year)
    !    day    [integer] input (day of the month)
    !    second [integer] input (seconds of the day)
    !    date   [julian_date] output (Julian day)
    !    )
    !
    !-
    !
    INTEGER, INTENT(in) :: ky
    INTEGER, INTENT(in) :: km
    INTEGER, INTENT(in) :: kd
    INTEGER, INTENT(in) :: ksec
    LOGICAL, INTENT(in), OPTIONAL :: lperpetual_year
    !
    ! for reference: 1. January 1998 00 UTC === Julian Day 2450814.5
    !
    TYPE (julian_date), INTENT(out) :: my_day

    LOGICAL :: lperp
    INTEGER :: ib, iy, im, idmax
    REAL(dp) :: zd, zsec

    IF ( ksec < 0 .OR. ksec > 86399) THEN
      CALL finish ('mo_time_base:Set_JulianDay', 'invalid number of seconds')
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

#ifdef IAU
    IF (ky > 1582 .OR. (ky == 1582 .AND. km > 10 &
         .OR. (km == 10 .AND. kd >= 15))) THEN
      ! 15th October 1582 AD or later no changes

    ELSE IF (ky == 1582 .AND. km == 10 .AND. (4 < kd .AND. kd < 15) ) THEN
      ! a small pice from the history:
      !     Papst Gregor XIII corrected the calendar,
      !     he skipped 10 days between the 4th and the 15th of October 1582
      !
      WRITE (message_text, '(a,i2.2,a1,i2.2,a1,i4.4,a)') &
           'date ', km, '/', kd, '/', ky, ' not valid'
      CALL finish ('mo_time_base:Set_JulianDay', message_text)
    ELSE 

      ! 4th October 1582 AD or earlier
      ib = -2

    ENDIF
#endif

    lperp = .false.
    IF (PRESENT(lperpetual_year)) lperp = lperpetual_year
    IF (.NOT. lperp) THEN
      ! check the length of the month
      idmax = Get_JulianMonLen (ky, km)

      IF (kd < 1 .OR. idmax < kd) &
           CALL finish ('mo_time_base:Set_JulianDay' ,'day in months invalid')
    ENDIF

    zd = FLOOR(365.25_dp*iy)+INT(30.6001_dp*(im+1)) &
        +REAL(ib,dp)+1720996.5_dp+REAL(kd,dp)+zsec

    my_day %day      = AINT(zd,dp)
    my_day %fraction = zd-AINT(zd,dp)

  END SUBROUTINE Set_JulianDay

  SUBROUTINE Set_JulianCalendar(my_day, ky, km, kd, ksec)
    !+
    !
    ! Set_JulianCalendar [subroutine]
    !    convert Julian date into year, months, day, seconds of day
    !    (
    !    date   [julian_date]  input (Julian day)
    !    year   [integer] output (calendar year)
    !    month  [integer] output (month of the year)
    !    day    [integer] output (day of the month)
    !    second [integer] output (seconds of the day)
    !    )
    !
    !-
    TYPE (julian_date), INTENT(IN) :: my_day
    INTEGER, INTENT(OUT)           :: ky
    INTEGER, INTENT(OUT)           :: km
    INTEGER, INTENT(OUT)           :: kd
    INTEGER, INTENT(OUT)           :: ksec
    
    REAL(dp) :: za, zb, zc, zd, ze, zf, zday, zfrac, zsec

    zday  = my_day%day
    zfrac = my_day%fraction

    za   = FLOOR(zday+zfrac+0.5_dp)

    zb = FLOOR((za-1867216.25_dp)/36524.25_dp)
    zc = za+zb-FLOOR(zb/4)+1525

#ifdef IAU
    IF (za < 2299161.0_dp) zc = za+1524.0_dp
#endif

    zd = FLOOR((zc-122.1_dp)/365.25_dp)
    ze = FLOOR(365.25_dp*zd)
    zf = FLOOR((zc-ze)/30.6001_dp)

    kd = INT(zc-ze- FLOOR(30.6001_dp*zf))
    km = INT(zf-1-12*FLOOR(zf/14))
    ky = INT(zd-4715-((7+km)/10))

    ! convert fraction in seconds of the day
    IF (zfrac < -0.5_dp) THEN
       zsec = 1.5_dp + zfrac
    ELSE IF (zfrac < 0.5_dp) THEN
       zsec = 0.5_dp + zfrac
    ELSE
       zsec = zfrac - 0.5_dp
    END IF
    ksec = frac2sec(zsec)

  END SUBROUTINE Set_JulianCalendar

  SUBROUTINE Get_JulianYearDay(my_day, kjy, kd, ksec)
    !+
    !
    ! Get_JulianYearDay [subroutine]
    !    get Julian year, day of year and seconds of day from Julian date
    !    (
    !    date   [julian_date]  input (Julian day)
    !    year   [integer] output (Julian calendar year)
    !    day    [integer] output (day of the year)
    !    second [integer] output (seconds of the day)
    !    )
    !
    !-
    ! Day of Astronomical Year (approximated)
    !
    ! it is used for the correct annual cycle
    ! remember the mismatch between calendar and
    ! astronomical year accumulated until 1582
    !
    ! this version gives adjustment until 1583
    !
    TYPE (julian_date), INTENT(IN) :: my_day
    INTEGER, INTENT(OUT)           :: kjy
    INTEGER, INTENT(OUT)           :: kd
    INTEGER, INTENT(out)           :: ksec

    TYPE (julian_date) :: first_jan, present_day
    INTEGER :: iy, id, im

    ! AD 1 is year 4713 in Julian calendar, 31st of December 1997 the Julian 
    ! day 2450814.0 is elapsed at 12 UTC and the Julian year is 6711, so:
    ! year_len = 2450814. / 6710. = 365.25

    ! get the date of the present julian day
    CALL Set_JulianCalendar(my_day, iy, im, id, ksec)
    
    ! find the first of January 00UTC
    CALL Set_JulianDay(iy,  1,  1, 0, first_jan)

    ! find the present Julian Day at 00UTC
    CALL Set_JulianDay(iy, im, id, 0, present_day)

    ! the day in the year is given
    kd = present_day%day-first_jan%day+1

    ! get the Julian year
    kjy = iy+4712

  END SUBROUTINE Get_JulianYearDay

  INTEGER FUNCTION Get_JulianYearLen(ky)
    !+
    !
    ! Get_JulianYearLen  [function, integer]
    !    get the length of a Julian year in days
    !    (
    !    year [integer] input (Calendar year)
    !    )
    !
    !-
    INTEGER, INTENT(in) :: ky

    INTEGER :: ylen

#ifdef IAU
    IF (ky == 1582) THEN
       ylen = 355
    ELSE IF ( (MOD(ky,4)==0 .AND. MOD(ky,100)/=0) .OR. MOD(ky,400)==0 ) THEN
       ylen = 366
    ELSE
       ylen = 365
    END IF
#else
    IF ( (MOD(ky,4)==0 .AND. MOD(ky,100)/=0) .OR. MOD(ky,400)==0 ) THEN
       ylen = 366
    ELSE
       ylen = 365
    END IF
#endif
    Get_JulianYearLen = ylen

  END FUNCTION Get_JulianYearLen

  INTEGER FUNCTION Get_JulianMonLen(ky, km)
    !+
    !
    ! Get_JulianMonLen [function, integer]
    !    get the length of a months in a Julian year
    !    (
    !    year  [integer] input (Calendar year)
    !    month [integer] input (month of the year)
    !    )
    !
    !-
    INTEGER, INTENT(in) :: km, ky

    INTEGER :: idmax

    SELECT CASE(km)
    CASE(1,3,5,7,8,10,12);  idmax = 31
    CASE(4,6,9,11);                  idmax = 30
    CASE(2)
      IF ( (MOD(ky,4)==0 .AND. MOD(ky,100)/=0) .OR. MOD(ky,400)==0 ) THEN
        ! leap year found
        idmax = 29
      ELSE
        idmax = 28
      END IF

    CASE default
      CALL finish('mo_time_base:Get_JulianMonLen', 'month invalid')

    END SELECT
    Get_JulianMonLen = idmax

  END FUNCTION Get_JulianMonLen

  SUBROUTINE Print_JulianDay(my_day)
    !+
    !
    ! Print_JulianDay [subroutine] interface [print]
    !   print Julian/model day on standard output
    !   (
    !   date [julian_date] input (model/Julian calendar day)
    !   )
    !
    !-
    TYPE (julian_date), INTENT(in) :: my_day

    WRITE(message_text,*) my_day%day+my_day%fraction
    CALL message('Print_JulianDay', message_text)

  END SUBROUTINE Print_JulianDay        

  SUBROUTINE Set_Ly360Day(ky, km, kd, ksec, my_day)
    !+
    !
    ! Set_Ly360Day  [subroutine]
    !    convert year, month, day, seconds into Ly360 calendar day
    !    (
    !    year   [integer] input (calendar year)
    !    month  [integer] input (month of the year)
    !    day    [integer] input (day of the month)
    !    second [integer] input (seconds of the day)
    !    date   [Ly360_date] output (Ly360 day)
    !    )
    !
    !-
    !
    INTEGER, INTENT(IN) :: ky
    INTEGER, INTENT(IN) :: km
    INTEGER, INTENT(IN) :: kd
    INTEGER, INTENT(IN) :: ksec
    !
    TYPE (Ly360_date), INTENT(out) :: my_day

    IF ( ksec < 0 .OR. ksec > 86399) THEN
      CALL finish ('mo_time_base:Set_Ly360Day', 'invalid number of seconds')
    ENDIF
    
    my_day %day      = 360*ky+(km-1)*30+(kd-1)
    my_day %fraction = sec2frac(ksec)

#ifdef DEBUG
    write (0,*) 'date2internal: ', ky, km, kd, ksec, ' -> ', my_day%day, my_day%fraction
#endif

  END SUBROUTINE Set_Ly360Day

  SUBROUTINE Set_Ly360Calendar(my_day, ky, km, kd, ksec)
    !+
    !
    ! Set_Ly360Calendar [subroutine]
    !    convert Ly360 date into year, months, day, seconds of day
    !    (
    !    date   [Ly360_date]  input (Ly360 day)
    !    year   [integer] output (calendar year)
    !    month  [integer] output (month of the year)
    !    day    [integer] output (day of the month)
    !    second [integer] output (seconds of the day)
    !    )
    !
    !-
    TYPE (Ly360_date), INTENT(IN) :: my_day
    INTEGER, INTENT(OUT)           :: ky
    INTEGER, INTENT(OUT)           :: km
    INTEGER, INTENT(OUT)           :: kd
    INTEGER, INTENT(OUT)           :: ksec
    
    INTEGER :: ir

    IF (my_day%day < 0) THEN
      ky = my_day%day/360-1
      ir  = MOD(my_day%day, 360_i8)
      km = 12-ir/30
      kd = 30-(MOD(ir,30)+1)
    ELSE
      ky = my_day%day/360
      ir  = MOD(my_day%day, 360_i8)
      km = ir/30+1
      kd = MOD(ir,30)+1
    ENDIF

    ksec = frac2sec(my_day%fraction)

#ifdef DEBUG
    write (0,*) 'internal2date: ', my_day%day, my_day%fraction,  ' -> ', ky, km, kd, ksec
#endif

  END SUBROUTINE Set_Ly360Calendar

  SUBROUTINE Get_Ly360YearDay(my_day, ky, kd, ksec)
    !+
    !
    ! Get_Ly360YearDay [subroutine]
    !    get Ly360 year, day of year and seconds of day from Ly360 date
    !    (
    !    date   [Ly360_date]  input (Ly360 day)
    !    year   [integer] output (Ly360 calendar year)
    !    day    [integer] output (day of the year)
    !    second [integer] output (seconds of the day)
    !    )
    !
    !-
    ! Day of Astronomical Year (approximated)
    !
    ! it is used for the correct annual cycle
    ! remember the mismatch between calendar and
    ! astronomical year accumulated until 1582
    !
    ! this version gives adjustment until 1583
    !
    TYPE (Ly360_date), INTENT(IN) :: my_day
    INTEGER, INTENT(OUT)           :: ky
    INTEGER, INTENT(OUT)           :: kd
    INTEGER, INTENT(out)           :: ksec

    ky = my_day%day/360
    kd  = MOD(my_day%day, 360_i8)

    ksec = frac2sec(my_day%fraction)

#ifdef DEBUG
    write (0,*) 'internal2date: ', my_day%day, my_day%fraction,  ' -> (2) ', ky, kd, ksec 
#endif

  END SUBROUTINE Get_Ly360YearDay

  INTEGER FUNCTION Get_Ly360YearLen()
    !+
    !
    ! Get_Ly360YearLen  [function, integer]
    !    get the length of a Ly360 year in days
    !    (
    !    )
    !
    !-

    Get_Ly360YearLen = 360

  END FUNCTION Get_Ly360YearLen

  INTEGER FUNCTION Get_Ly360MonLen()
    !+
    !
    ! Get_Ly360MonLen [function, integer]
    !    get the length of a months in a Ly360 year
    !    (
    !    )
    !
    !-

    Get_Ly360MonLen = 30

  END FUNCTION Get_Ly360MonLen

  SUBROUTINE Print_Ly360Day(my_day)
    !+
    !
    ! Print_Ly360Day [subroutine] interface [print]
    !   print Ly360/model day on standard output
    !   (
    !   date [Ly360_date] input (model/Ly360 calendar day)
    !   )
    !
    !-
    TYPE (Ly360_date), INTENT(in) :: my_day

    WRITE(message_text,*) my_day%day+my_day%fraction
    CALL message('Print_Ly360Day', message_text)

  END SUBROUTINE Print_Ly360Day 

  FUNCTION sec2frac(isec) RESULT (zfrac)
    !+
    !
    ! sec2frac [function, real_dp]
    !    convert seconds of day into fraction
    !    (
    !    second [integer] input (seconds of the day)
    !    )
    !
    !-
    INTEGER , INTENT(in) :: isec
    REAL(dp) :: zfrac

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      zfrac = (REAL(isec,dp)+0.000001_dp)/REAL(IDAYLEN,dp)
    CASE (CYL360)
      zfrac = REAL(isec,dp)/REAL(IDAYLEN,dp)
    END SELECT

  END FUNCTION Sec2frac

  FUNCTION frac2sec(zfrac) RESULT (isec)
    !+
    !
    ! frac2sec [function, integer]
    !    convert fraction of day into seconds
    !    (
    !    fraction [real_dp] input (fraction of the day)
    !    )
    !
    !-
    REAL(dp), INTENT(in) :: zfrac
    INTEGER :: isec

    isec = INT(zfrac*REAL(IDAYLEN,dp))

  END FUNCTION frac2sec

END MODULE mo_time_base
