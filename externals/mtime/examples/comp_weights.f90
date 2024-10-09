!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM comp_weights

  USE mtime

  IMPLICIT NONE

  INTEGER, PARAMETER :: pd = 12
  INTEGER, PARAMETER :: rd = 307
  !
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(pd, rd) !< double precision

  INTEGER  :: m1, m2
  REAL(wp) :: pw1, pw2

  TYPE(datetime), POINTER :: test_date => NULL()

  CALL month2hour(2000, 3, 5, 13, m1, m2, pw1, pw2)
  WRITE (0, *) 'month2hour 1: ', m1, m2, pw1, pw2

  CALL time_weights_limm(2000, 3, 5, 13, m1, m2, pw1, pw2)
  WRITE (0, *) 'time_intpl 1: ', m1, m2, pw1, pw2

  CALL setCalendar(proleptic_gregorian)

  test_date => newDatetime('2000-03-05T13:00:00')
  CALL calculate_time_interpolation_weights(test_date, m1, m2, pw1, pw2)
  WRITE (0, *) 'mtime      1: ', m1, m2, pw1, pw2
  CALL deallocateDatetime(test_date)

  CALL month2hour(2000, 7, 23, 7, m1, m2, pw1, pw2)
  WRITE (0, *) 'month2hour 1: ', m1, m2, pw1, pw2

  CALL time_weights_limm(2000, 7, 23, 7, m1, m2, pw1, pw2)
  WRITE (0, *) 'time_intpl 1: ', m1, m2, pw1, pw2

  CALL setCalendar(proleptic_gregorian)

  test_date => newDatetime('2000-07-23T07:00:00')
  CALL calculate_time_interpolation_weights(test_date, m1, m2, pw1, pw2)
  WRITE (0, *) 'mtime      1: ', m1, m2, pw1, pw2
  CALL deallocateDatetime(test_date)

CONTAINS

  INTEGER FUNCTION mleapy(myy)
    INTEGER, INTENT(in) :: myy
    mleapy = MAX(1 - MODULO(myy, 4), 0) - MAX(1 - MODULO(myy, 100), 0) + MAX(1 - MODULO(myy, 400), 0)
  END FUNCTION mleapy

  SUBROUTINE month2hour(year, month, day, hour, m1, m2, pw1, pw2)
    INTEGER, INTENT(in)  :: year, month, day, hour
    INTEGER, INTENT(out) :: m1, m2     ! indices of nearest months
    REAL(wp), INTENT(out) :: pw1, pw2   ! weights of nearest months

    INTEGER :: dayhour,                & ! actual date (in hours of month)
         &     midthhours,             & ! midth of month in hours
         &     i, ip1                    ! month indices (ip1=i+1)

    REAL(wp) :: zdiff(12),             & ! difference between midth of following months in days
         &      zhalf(12),             & ! number of days for half month
         &      zact                     ! actual time in hours of month

    ! number of days for each month
    INTEGER :: month_days(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    ! compute half of each month (in days)

    ! leap year ??
    month_days(2) = 28 + mleapy(year)
    zhalf(:) = 0.5_wp*REAL(month_days(:), wp)

    ! compute difference between the midth of actual month and the
    ! following one (in days)

    DO i = 1, 12
      ip1 = MOD(i, 12) + 1
      zdiff(i) = zhalf(i) + zhalf(ip1)
    END DO

    ! compute actual date (day and hours) and midth of actual month in hours
    dayhour = (day - 1)*24 + hour
    midthhours = NINT(zhalf(month)*24._wp)

    ! determine the months needed for interpolation of current values
    ! search for the position of date in relation to first of month.
    ! the original data are valid for the mid-month.
    !
    ! example 1
    !        march    !  april     !   may             x : aerosol data
    !       ----x-----!-----x----o-!-----x-----        ! : first of month
    !                       !    ^       !             o : current date
    !                       !    ^ interpolation for that point in time
    !                       !  zdiff(4)  !
    !                       !zact!
    !
    ! example 2
    !        march    !  april     !   may             x : ndvi_ratio
    !       ----x-----!-----x------!----ox-----        ! : first of month
    !                       !           ^              o : current date
    !                       !      interpolation for that point in time
    !                       !zhalf !
    !                       !  zdiff(4)  !
    !                       !   zact    !
    !
    !

    IF (dayhour < midthhours) THEN
      ! point is in first half of month (example 2)
      m1 = month - 1
      IF (month == 1) m1 = 12
      zact = zhalf(m1) + REAL(dayhour, wp)/24._wp
    ELSE
      ! point is in second half of month (example 1)
      m1 = month
      zact = REAL(dayhour - midthhours, wp)/24._wp
    END IF
    m2 = MOD(m1, 12) + 1
    pw2 = zact/zdiff(m1)
    pw1 = 1.0_wp - pw2

  END SUBROUTINE month2hour

  SUBROUTINE time_weights_limm(year, month, day, hour, m1, m2, pw1, pw2)
    INTEGER, INTENT(in)  :: year, month, day, hour
    INTEGER, INTENT(out) :: m1, m2     ! indices of nearest months
    REAL(wp), INTENT(out) :: pw1, pw2   ! weights of nearest months

    REAL(wp) :: zcmonfrc ! current month fraction
    REAL(wp) :: zevent_tim !time in month
    REAL(wp) :: zcmlen2, znmlen2 !half of current/nearest month length

    INTEGER :: my_year, my_month, my_day, my_hour, my_monlen

    INTEGER :: month_days(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    month_days(2) = 28 + mleapy(year)

    ! compute half of each month (in days)

    zcmlen2 = 0.5_wp*REAL(month_days(month), wp)
    zevent_tim = (24.0_wp*REAL(day - 1, wp) + hour)/24.0_wp ! event time since first of month in frac. days
    zcmonfrc = zevent_tim/REAL(month_days(month), wp)

    my_year = year
    my_month = month
    my_day = 1
    my_hour = hour
    IF (zcmonfrc <= 0.5_wp) THEN
      ! interpolate between value of previous and current month,
      m1 = month - 1
      m2 = month
      IF (m1 == 0) THEN
        m1 = 12
      END IF
      my_monlen = month_days(m1)
      znmlen2 = my_monlen*0.5_wp
      pw1 = (zcmlen2 - zevent_tim)/(zcmlen2 + znmlen2)
      pw2 = 1._wp - pw1
    ELSE
      ! interpolate between value of currrent and next month,
      m1 = month
      m2 = month + 1
      IF (m2 == 13) THEN
        m2 = 1
      END IF
      my_monlen = month_days(m1)
      znmlen2 = my_monlen*0.5_wp
      pw2 = (zevent_tim - zcmlen2)/(zcmlen2 + znmlen2)
      pw1 = 1._wp - pw2
    END IF

  END SUBROUTINE time_weights_limm

  SUBROUTINE calculate_time_interpolation_weights(current_date, month1, month2, weight1, weight2)
    TYPE(datetime), POINTER, INTENT(in)  :: current_date
    INTEGER, INTENT(out) :: month1, month2
    REAL(wp), INTENT(out) :: weight1, weight2

    TYPE(datetime), POINTER :: next_month => NULL(), previous_month => NULL()
    TYPE(timedelta), POINTER :: one_month => NULL()

    INTEGER :: days_in_previous_month, days_in_month, days_in_next_month
    INTEGER :: seconds_in_month, seconds_in_middle_of_previous_month, seconds_in_middle_of_month, seconds_in_middle_of_next_month

    days_in_month = getNoOfDaysInMonthDateTime(current_date)
    seconds_in_middle_of_month = 43200*days_in_month ! 86400 * my_month_len / 2

    seconds_in_month = INT(getNoOfSecondsElapsedInMonthDateTime(current_date))

    IF (seconds_in_month <= seconds_in_middle_of_month) THEN
      ! first half of month
      one_month => newTimedelta('-P1M')
      previous_month => newDatetime(current_date)
      previous_month = current_date + one_month
      days_in_previous_month = getNoOfDaysInMonthDateTime(previous_month)
      seconds_in_middle_of_previous_month = 43200*days_in_previous_month ! 86400 * my_month_len / 2
      ! simple linear interpolation
      weight1 = REAL(seconds_in_middle_of_month - seconds_in_month, wp) &
           &   /REAL(seconds_in_middle_of_month + seconds_in_middle_of_previous_month, wp)
      weight2 = 1.0_wp - weight1
      month1 = previous_month%date%month
      month2 = current_date%date%month
    ELSE
      ! second half of month

      one_month => newTimedelta('P1M')
      next_month => newDatetime(current_date)
      next_month = current_date + one_month
      days_in_next_month = getNoOfDaysInMonthDateTime(next_month)
      seconds_in_middle_of_next_month = 43200*days_in_next_month ! 86400 * my_month_len / 2
      ! simple linear interpolation
      weight2 = REAL(seconds_in_month - seconds_in_middle_of_month, wp) &
           &   /REAL(seconds_in_middle_of_month + seconds_in_middle_of_next_month, wp)
      weight1 = 1.0_wp - weight2
      month1 = current_date%date%month
      month2 = next_month%date%month
    END IF

    CALL deallocateTimedelta(one_month)
    IF (ASSOCIATED(previous_month)) CALL deallocateDatetime(previous_month)
    IF (ASSOCIATED(next_month)) CALL deallocateDatetime(next_month)

  END SUBROUTINE calculate_time_interpolation_weights

END PROGRAM comp_weights
