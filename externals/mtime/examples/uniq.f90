!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
program output_test

  use mtime
  
  implicit none

  integer, parameter :: i8 = selected_int_kind(14)
  
  integer, parameter :: nintvls = 8
  
  character(len=24) :: begin_str(nintvls)
  character(len=24) :: end_str(nintvls)
  character(len=24) :: intvl_str(nintvls)

  type(datetime),  pointer :: mtime_begin => null()
  type(datetime),  pointer :: mtime_end => null()
  type(timedelta), pointer :: delta => null()
  type(datetime),  pointer :: mtime_date => null()

  type(julianday), pointer :: tmp_jd => null()
  
  type tmp_container
    integer(i8) :: day
    integer(i8) :: ms
  end type tmp_container

  type(tmp_container), allocatable, target :: mtime_date_container_a(:)
  type(tmp_container), allocatable, target :: mtime_date_container_b(:)  
  type(tmp_container), allocatable, target :: tmp(:)
  type(tmp_container), allocatable :: mtime_date_uniq(:)
  
  type(tmp_container), pointer :: mtime_date_container(:) => null()
  
  integer :: indices_to_use(nintvls)
  integer :: remaining_intvls, iselected_intvl
  
  integer :: iintvl
  integer :: n_event_steps, n_event_steps_a, n_event_steps_b, remaining_event_steps, ierrstat   

  logical :: l_active

  call setCalendar(proleptic_gregorian)
  
  begin_str(1) = "1979-01-01T00:00:00Z"
  end_str(1)   = "2009-01-01T00:00:00Z"
!  end_str(1)   = "1979-01-02T00:00:00Z"  
  intvl_str(1) = "PT6H"

  begin_str(2) = "1979-01-01T00:00:00Z"
  end_str(2)   = "2009-01-01T00:00:00Z"
!  end_str(2)   = "1979-01-02T00:00:00Z"  
  intvl_str(2) = "PT2H"

  begin_str(3) = "1979-01-01T00:00:00Z"
  end_str(3)   = "2009-01-01T00:00:00Z"
!  end_str(3)   = "1979-01-02T00:00:00Z"  
  intvl_str(3) = "PT6H"

  begin_str(4) = "1979-01-01T00:00:00Z"
  end_str(4)   = "2009-01-01T00:00:00Z"
!  end_str(4)   = "1979-01-02T00:00:00Z"  
  intvl_str(4) = "PT2H"

  begin_str(5) = "1979-01-01T00:00:00Z"
  end_str(5)   = "2009-01-01T00:00:00Z"
!  end_str(5)   = "1979-01-02T00:00:00Z"  
  intvl_str(5) = "PT2H"

  begin_str(6) = "1979-01-01T00:00:00Z"
  end_str(6)   = "2009-01-01T00:00:00Z"
!  end_str(6)   = "1979-01-02T00:00:00Z"  
  intvl_str(6) = "PT6H"

  begin_str(7) = "1979-01-01T00:00:00Z"
  end_str(7)   = "2009-01-01T00:00:00Z"
!  end_str(7)   = "1979-01-02T00:00:00Z"  
  intvl_str(7) = "PT1H"

  begin_str(8) = "1979-01-01T00:00:00Z"
  end_str(8)   = "2009-01-01T00:00:00Z"
!  end_str(8)   = "1979-01-02T00:00:00Z"  
  intvl_str(8) = "PT6H"

  call remove_duplicate_intervals(begin_str, end_str, intvl_str, nintvls, indices_to_use, remaining_intvls)

  n_event_steps = 0

  tmp_jd => newJulianday(0_i8, 0_i8);

  allocate(mtime_date_container_a(256), stat=ierrstat)
  allocate(mtime_date_container_b(256), stat=ierrstat)      

  mtime_date_container => mtime_date_container_a
  
  interval_loop:  do iselected_intvl = 1, remaining_intvls
    iintvl = indices_to_use(iselected_intvl)
    
    write (0,*) 'Handling set ', iintvl
    
    mtime_begin => newDatetime(trim(begin_str(iintvl)))
    mtime_end   => newDatetime(trim(end_str(iintvl)))
    delta       => newTimedelta(trim(intvl_str(iintvl)))
    
    mtime_date  => mtime_begin
    
    event_loop: do
      
      n_event_steps = n_event_steps + 1
      
      if (n_event_steps > size(mtime_date_container)) then
        allocate(tmp(2*size(mtime_date_container)), stat=ierrstat)
        if (ierrstat /= 0) stop 'allocate failed'
        tmp(1:size(mtime_date_container)) = mtime_date_container(:)
        if (associated(mtime_date_container, mtime_date_container_a)) then
          call move_alloc(tmp, mtime_date_container_a)
          mtime_date_container => mtime_date_container_a
        else
          call move_alloc(tmp, mtime_date_container_b)
          mtime_date_container => mtime_date_container_b
        endif
      endif

      call getJulianDayFromDatetime(mtime_date, tmp_jd)
      mtime_date_container(n_event_steps)%day = tmp_jd%day
      mtime_date_container(n_event_steps)%ms = tmp_jd%ms

      mtime_date = mtime_date + delta
      
      l_active = .not. (mtime_date > mtime_end)
      
      if (.not. l_active) exit event_loop
      
    enddo event_loop

    write (0,*) '   events: ', n_event_steps

    ! for first event set we do not remove douplicates

    if (iintvl == 1) then
      n_event_steps_a = n_event_steps
      mtime_date_container => mtime_date_container_b
      n_event_steps = 0
      cycle interval_loop
    else
      n_event_steps_b = n_event_steps
      n_event_steps = 0      

      call merge2SortedAndRemoveDublicates(mtime_date_container_a, n_event_steps_a, &
           &                               mtime_date_container_b, n_event_steps_b, &
           &                               mtime_date_uniq, remaining_event_steps)
    endif

    write (0,*) '     remaining events: ', remaining_event_steps

    if (remaining_event_steps > size(mtime_date_container_a)) then
      allocate(tmp(size(mtime_date_container_a)), stat=ierrstat)
      if (ierrstat /= 0) stop 'allocate failed'
      tmp(1:remaining_event_steps) = mtime_date_uniq(1:remaining_event_steps)
      call move_alloc(tmp, mtime_date_container_a)
    endif
    
    n_event_steps_a = remaining_event_steps

    deallocate(mtime_date_uniq)
    
  enddo interval_loop

  deallocate(mtime_date_container_a)
  deallocate(mtime_date_container_b)  
  deallocate(tmp_jd)
  
contains

  subroutine merge2SortedAndRemoveDublicates(InputArray1, nsize_IA1, &
       &                                     InputArray2, nsize_IA2, &
       &                                     OutputArray, nsize_OA)
    type(tmp_container), intent(in) :: InputArray1(:)
    type(tmp_container), intent(in) :: InputArray2(:)
    type(tmp_container), allocatable, intent(out) :: OutputArray(:)
    integer, intent(in) :: nsize_IA1
    integer, intent(in) :: nsize_IA2
    integer, intent(out) :: nsize_OA

    integer(i8) :: diff
    
    integer :: n, na, nb
    integer :: i, j, k
    
    na = nsize_IA1
    nb = nsize_IA2 
    n = na + nb

    allocate(OutputArray(n))
    
    i = 1
    j = 1
    k = 1
    
    do while(i <= na .and. j <= nb)
      diff = 86400000_i8 * (InputArray1(i)%day - InputArray2(j)%day) + InputArray1(i)%ms - InputArray2(j)%ms
      if (diff < 0_i8) then
        if (k == 1 .or. ((InputArray1(i)%day /= OutputArray(k-1)%day) .or. (InputArray1(i)%ms /= OutputArray(k-1)%ms))) then
          OutputArray(k) = InputArray1(i)
          k = k+1
        endif
        i=i+1
      else if (diff > 0_i8) then
        if (k == 1 .or. ((InputArray2(j)%day /= OutputArray(k-1)%day) .or. (InputArray2(j)%ms /= OutputArray(k-1)%ms))) then
          OutputArray(k) = InputArray2(j)
          k = k+1
        endif
        j = j+1
      else
        if (k == 1 .or. ((InputArray1(i)%day /= OutputArray(k-1)%day) .or. (InputArray1(i)%ms /= OutputArray(k-1)%ms))) then
          OutputArray(k) = InputArray1(i)
          k = k+1
        endif
        i = i+1
        j = j+1
      endif
    enddo
    
    do while (i <= na)
      if ((InputArray1(i)%day /= OutputArray(k-1)%day) .or. (InputArray1(i)%ms /= OutputArray(k-1)%ms)) then
        OutputArray(k) = InputArray1(i)
        k = k+1
        i = i+1
      else
        i = i+1
      endif
    enddo
    
    do while (j <= nb)
      if ((InputArray2(j)%day /= OutputArray(k-1)%day) .or. (InputArray2(j)%ms /= OutputArray(k-1)%ms)) then
        OutputArray(k) = InputArray2(j)
        k = k+1
        j = j+1
      else
        j = j+1
      endif
    enddo
    
    ! do i = 1, k-1
    !   write (0,*) OutputArray(i)
    ! enddo

    nsize_OA = k-1
    
  end subroutine merge2SortedAndRemoveDublicates

  subroutine remove_duplicate_intervals(starts, ends, intvls, n, indices_to_use, remaining)
    character(len=*), intent(in) :: starts(:)
    character(len=*), intent(in) :: ends(:)
    character(len=*), intent(in) :: intvls(:)
    integer, intent(in) :: n         
    integer, intent(out) :: indices_to_use(n)
    integer, intent(out) :: remaining
    
    integer :: i, j
    
    remaining = 1
    indices_to_use(1) = 1
    
    outer: do i = 2, n
      do j = 1, remaining
        if (trim(starts(j)) == trim(starts(i))) then
          if (trim(ends(j)) == trim(ends(i))) then
            if (trim(intvls(j)) == trim(intvls(i))) then
              ! found a match so start looking again
              cycle outer
            endif
          endif
        endif
      enddo
      ! no match found so add it to the output
      remaining = remaining + 1
      indices_to_use(remaining) = i
    end do outer
    
  end subroutine remove_duplicate_intervals

end program output_test
