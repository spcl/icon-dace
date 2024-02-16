!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
program comp_weights

  use mtime
  
  implicit none

  integer, parameter :: pd =  12
  integer, parameter :: rd = 307
  !
  integer, parameter :: wp = selected_real_kind(pd,rd) !< double precision

  integer  :: m1, m2
  real(wp) :: pw1, pw2

  type(datetime), pointer :: test_date => null()
  
  call month2hour(2000,  3,  5, 13, m1, m2, pw1, pw2)
  write(0,*) 'month2hour 1: ', m1, m2, pw1, pw2

  call time_weights_limm(2000, 3, 5, 13, m1, m2, pw1, pw2)
  write(0,*) 'time_intpl 1: ', m1, m2, pw1, pw2

  call setCalendar(proleptic_gregorian)
  
  test_date => newDatetime('2000-03-05T13:00:00')
  call calculate_time_interpolation_weights(test_date, m1, m2, pw1, pw2)
  write(0,*) 'mtime      1: ', m1, m2, pw1, pw2
  call deallocateDatetime(test_date)
  
  call month2hour(2000,  7,  23, 7, m1, m2, pw1, pw2)
  write(0,*) 'month2hour 1: ', m1, m2, pw1, pw2

  call time_weights_limm(2000, 7, 23, 7, m1, m2, pw1, pw2)
  write(0,*) 'time_intpl 1: ', m1, m2, pw1, pw2

  call setCalendar(proleptic_gregorian)
  
  test_date => newDatetime('2000-07-23T07:00:00')
  call calculate_time_interpolation_weights(test_date, m1, m2, pw1, pw2)
  write(0,*) 'mtime      1: ', m1, m2, pw1, pw2
  call deallocateDatetime(test_date)
  
contains
  
  integer function mleapy(myy)
    integer, intent(in) :: myy
    mleapy = max(1-modulo(myy,4),0)-max(1-modulo(myy,100),0)+max(1-modulo(myy,400),0)
  end function mleapy

  subroutine month2hour( year, month, day, hour, m1, m2, pw1, pw2 )
    integer,  intent(in)  :: year, month, day, hour
    integer,  intent(out) :: m1, m2     ! indices of nearest months
    real(wp), intent(out) :: pw1, pw2   ! weights of nearest months
    
    integer :: dayhour,                & ! actual date (in hours of month)
         &     midthhours,             & ! midth of month in hours
         &     i, ip1                    ! month indices (ip1=i+1)
    
    real(wp) :: zdiff(12),             & ! difference between midth of following months in days
         &      zhalf(12),             & ! number of days for half month
         &      zact                     ! actual time in hours of month
    
    ! number of days for each month
    integer :: month_days(12) = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
    
    ! compute half of each month (in days)

    ! leap year ??
    month_days (2) = 28 + mleapy(year)
    zhalf(:) = 0.5_wp * real(month_days(:),wp)
    
    ! compute difference between the midth of actual month and the
    ! following one (in days)

    do i = 1,12
      ip1 = mod(i,12)+1
      zdiff(i) = zhalf(i) + zhalf(ip1)
    enddo
    
    ! compute actual date (day and hours) and midth of actual month in hours
    dayhour = (day-1)*24 + hour
    midthhours = nint(zhalf(month)*24._wp)
    
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
    
    if( dayhour < midthhours) then
      ! point is in first half of month (example 2)
      m1 = month - 1
      if(month == 1) m1 = 12
      zact = zhalf(m1) + real(dayhour,wp)/24._wp
    else
      ! point is in second half of month (example 1)
      m1 = month
      zact = real(dayhour-midthhours,wp)/24._wp
    endif
    m2 = mod(m1,12) + 1
    pw2 =  zact/zdiff(m1)
    pw1 = 1.0_wp-pw2
    
  end subroutine month2hour

  subroutine time_weights_limm(year, month, day, hour, m1, m2, pw1, pw2)
    integer,  intent(in)  :: year, month, day, hour
    integer,  intent(out) :: m1, m2     ! indices of nearest months
    real(wp), intent(out) :: pw1, pw2   ! weights of nearest months

    real(wp) :: zcmonfrc ! current month fraction
    real(wp) :: zevent_tim !time in month
    real(wp) :: zcmlen2, znmlen2 !half of current/nearest month length

    integer :: my_year, my_month, my_day, my_hour, my_monlen
    
    integer :: month_days(12) = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
    month_days (2) = 28 + mleapy(year)

    ! compute half of each month (in days)

    zcmlen2 = 0.5_wp * real(month_days(month),wp)
    zevent_tim = (24.0_wp * real(day-1,wp) + hour) / 24.0_wp ! event time since first of month in frac. days
    zcmonfrc = zevent_tim / real(month_days(month),wp)

    my_year = year
    my_month = month
    my_day = 1
    my_hour = hour
    if (zcmonfrc <= 0.5_wp) then
      ! interpolate between value of previous and current month, 
      m1 = month-1  
      m2 = month    
      if (m1 == 0) then
        m1 = 12
      end if
      my_monlen = real(month_days(m1),wp)
      znmlen2 = my_monlen*0.5_wp
      pw1 = (zcmlen2-zevent_tim)/(zcmlen2+znmlen2)
      pw2 = 1._wp-pw1
    else
      ! interpolate between value of currrent and next month, 
      m1 = month
      m2 = month+1
      if (m2 == 13) then
        m2 = 1
      end if
      my_monlen = real(month_days(m1),wp)      
      znmlen2 = my_monlen*0.5_wp
      pw2 = (zevent_tim-zcmlen2)/(zcmlen2+znmlen2)
      pw1 = 1._wp-pw2
    end if
    
  end subroutine time_weights_limm

  subroutine calculate_time_interpolation_weights(current_date, month1, month2, weight1, weight2)
    type(datetime), pointer, intent(in)  :: current_date
    integer,                 intent(out) :: month1, month2
    real(wp),                intent(out) :: weight1, weight2

    type(datetime), pointer :: next_month => null(), previous_month => null()
    type(timedelta), pointer :: one_month => null()

    integer :: days_in_previous_month, days_in_month, days_in_next_month
    integer :: seconds_in_month, seconds_in_middle_of_previous_month, seconds_in_middle_of_month, seconds_in_middle_of_next_month

    integer :: errno
    
    days_in_month = getNoOfDaysInMonthDateTime(current_date) 
    seconds_in_middle_of_month = 43200 * days_in_month ! 86400 * my_month_len / 2
    
    seconds_in_month = getNoOfSecondsElapsedInMonthDateTime(current_date)

    if (seconds_in_month <= seconds_in_middle_of_month) then
      ! first half of month 
      one_month => newTimedelta('-P1M')
      previous_month => newDatetime(current_date)
      previous_month = current_date + one_month
      days_in_previous_month = getNoOfDaysInMonthDateTime(previous_month)
      seconds_in_middle_of_previous_month = 43200 * days_in_previous_month ! 86400 * my_month_len / 2
      ! simple linear interpolation
      weight1 = real(seconds_in_middle_of_month - seconds_in_month,wp) &
           &   /real(seconds_in_middle_of_month + seconds_in_middle_of_previous_month,wp)
      weight2 = 1.0_wp - weight1
      month1 = previous_month%date%month
      month2 = current_date%date%month
    else
      ! second half of month

      one_month => newTimedelta('P1M')
      next_month => newDatetime(current_date)
      next_month = current_date + one_month      
      days_in_next_month = getNoOfDaysInMonthDateTime(next_month)
      seconds_in_middle_of_next_month = 43200 * days_in_next_month ! 86400 * my_month_len / 2
      ! simple linear interpolation
      weight2 = real(seconds_in_month - seconds_in_middle_of_month,wp) &
           &   /real(seconds_in_middle_of_month + seconds_in_middle_of_next_month,wp)
      weight1 = 1.0_wp - weight2
      month1 = current_date%date%month
      month2 = next_month%date%month
    endif

    call deallocateTimedelta(one_month)
    if (associated(previous_month)) call deallocateDatetime(previous_month)
    if (associated(next_month)) call deallocateDatetime(next_month)
    
  end subroutine calculate_time_interpolation_weights
    
end program comp_weights
