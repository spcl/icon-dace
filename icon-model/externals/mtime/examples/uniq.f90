!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM output_test

  USE mtime

  IMPLICIT NONE

  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)

  INTEGER, PARAMETER :: nintvls = 8

  CHARACTER(len=24) :: begin_str(nintvls)
  CHARACTER(len=24) :: end_str(nintvls)
  CHARACTER(len=24) :: intvl_str(nintvls)

  TYPE(datetime), POINTER :: mtime_begin => NULL()
  TYPE(datetime), POINTER :: mtime_end => NULL()
  TYPE(timedelta), POINTER :: delta => NULL()
  TYPE(datetime), POINTER :: mtime_date => NULL()

  TYPE(julianday), POINTER :: tmp_jd => NULL()

  TYPE tmp_container
    INTEGER(i8) :: day
    INTEGER(i8) :: ms
  END TYPE tmp_container

  TYPE(tmp_container), ALLOCATABLE, TARGET :: mtime_date_container_a(:)
  TYPE(tmp_container), ALLOCATABLE, TARGET :: mtime_date_container_b(:)
  TYPE(tmp_container), ALLOCATABLE, TARGET :: tmp(:)
  TYPE(tmp_container), ALLOCATABLE :: mtime_date_uniq(:)

  TYPE(tmp_container), POINTER :: mtime_date_container(:) => NULL()

  INTEGER :: indices_to_use(nintvls)
  INTEGER :: remaining_intvls, iselected_intvl

  INTEGER :: iintvl
  INTEGER :: n_event_steps, n_event_steps_a, n_event_steps_b, remaining_event_steps, ierrstat

  LOGICAL :: l_active

  CALL setCalendar(proleptic_gregorian)

  begin_str(1) = "1979-01-01T00:00:00Z"
  end_str(1) = "2009-01-01T00:00:00Z"
!  end_str(1)   = "1979-01-02T00:00:00Z"
  intvl_str(1) = "PT6H"

  begin_str(2) = "1979-01-01T00:00:00Z"
  end_str(2) = "2009-01-01T00:00:00Z"
!  end_str(2)   = "1979-01-02T00:00:00Z"
  intvl_str(2) = "PT2H"

  begin_str(3) = "1979-01-01T00:00:00Z"
  end_str(3) = "2009-01-01T00:00:00Z"
!  end_str(3)   = "1979-01-02T00:00:00Z"
  intvl_str(3) = "PT6H"

  begin_str(4) = "1979-01-01T00:00:00Z"
  end_str(4) = "2009-01-01T00:00:00Z"
!  end_str(4)   = "1979-01-02T00:00:00Z"
  intvl_str(4) = "PT2H"

  begin_str(5) = "1979-01-01T00:00:00Z"
  end_str(5) = "2009-01-01T00:00:00Z"
!  end_str(5)   = "1979-01-02T00:00:00Z"
  intvl_str(5) = "PT2H"

  begin_str(6) = "1979-01-01T00:00:00Z"
  end_str(6) = "2009-01-01T00:00:00Z"
!  end_str(6)   = "1979-01-02T00:00:00Z"
  intvl_str(6) = "PT6H"

  begin_str(7) = "1979-01-01T00:00:00Z"
  end_str(7) = "2009-01-01T00:00:00Z"
!  end_str(7)   = "1979-01-02T00:00:00Z"
  intvl_str(7) = "PT1H"

  begin_str(8) = "1979-01-01T00:00:00Z"
  end_str(8) = "2009-01-01T00:00:00Z"
!  end_str(8)   = "1979-01-02T00:00:00Z"
  intvl_str(8) = "PT6H"

  CALL remove_duplicate_intervals(begin_str, end_str, intvl_str, nintvls, indices_to_use, remaining_intvls)

  n_event_steps = 0

  tmp_jd => newJulianday(0_i8, 0_i8); 
  ALLOCATE (mtime_date_container_a(256), stat=ierrstat)
  ALLOCATE (mtime_date_container_b(256), stat=ierrstat)

  mtime_date_container => mtime_date_container_a

  interval_loop: DO iselected_intvl = 1, remaining_intvls
    iintvl = indices_to_use(iselected_intvl)

    WRITE (0, *) 'Handling set ', iintvl

    mtime_begin => newDatetime(TRIM(begin_str(iintvl)))
    mtime_end => newDatetime(TRIM(end_str(iintvl)))
    delta => newTimedelta(TRIM(intvl_str(iintvl)))

    mtime_date => mtime_begin

    event_loop: DO

      n_event_steps = n_event_steps + 1

      IF (n_event_steps > SIZE(mtime_date_container)) THEN
        ALLOCATE (tmp(2*SIZE(mtime_date_container)), stat=ierrstat)
        IF (ierrstat /= 0) STOP 'allocate failed'
        tmp(1:SIZE(mtime_date_container)) = mtime_date_container(:)
        IF (ASSOCIATED(mtime_date_container, mtime_date_container_a)) THEN
          CALL MOVE_ALLOC(tmp, mtime_date_container_a)
          mtime_date_container => mtime_date_container_a
        ELSE
          CALL MOVE_ALLOC(tmp, mtime_date_container_b)
          mtime_date_container => mtime_date_container_b
        END IF
      END IF

      CALL getJulianDayFromDatetime(mtime_date, tmp_jd)
      mtime_date_container(n_event_steps)%day = tmp_jd%day
      mtime_date_container(n_event_steps)%ms = tmp_jd%ms

      mtime_date = mtime_date + delta

      l_active = .NOT. (mtime_date > mtime_end)

      IF (.NOT. l_active) EXIT event_loop

    END DO event_loop

    WRITE (0, *) '   events: ', n_event_steps

    ! for first event set we do not remove douplicates

    IF (iintvl == 1) THEN
      n_event_steps_a = n_event_steps
      mtime_date_container => mtime_date_container_b
      n_event_steps = 0
      CYCLE interval_loop
    ELSE
      n_event_steps_b = n_event_steps
      n_event_steps = 0

      CALL merge2SortedAndRemoveDublicates(mtime_date_container_a, n_event_steps_a, &
           &                               mtime_date_container_b, n_event_steps_b, &
           &                               mtime_date_uniq, remaining_event_steps)
    END IF

    WRITE (0, *) '     remaining events: ', remaining_event_steps

    IF (remaining_event_steps > SIZE(mtime_date_container_a)) THEN
      ALLOCATE (tmp(SIZE(mtime_date_container_a)), stat=ierrstat)
      IF (ierrstat /= 0) STOP 'allocate failed'
      tmp(1:remaining_event_steps) = mtime_date_uniq(1:remaining_event_steps)
      CALL MOVE_ALLOC(tmp, mtime_date_container_a)
    END IF

    n_event_steps_a = remaining_event_steps

    DEALLOCATE (mtime_date_uniq)

  END DO interval_loop

  DEALLOCATE (mtime_date_container_a)
  DEALLOCATE (mtime_date_container_b)
  DEALLOCATE (tmp_jd)

CONTAINS

  SUBROUTINE merge2SortedAndRemoveDublicates(InputArray1, nsize_IA1, &
       &                                     InputArray2, nsize_IA2, &
       &                                     OutputArray, nsize_OA)
    TYPE(tmp_container), INTENT(in) :: InputArray1(:)
    TYPE(tmp_container), INTENT(in) :: InputArray2(:)
    TYPE(tmp_container), ALLOCATABLE, INTENT(out) :: OutputArray(:)
    INTEGER, INTENT(in) :: nsize_IA1
    INTEGER, INTENT(in) :: nsize_IA2
    INTEGER, INTENT(out) :: nsize_OA

    INTEGER(i8) :: diff

    INTEGER :: n, na, nb
    INTEGER :: i, j, k

    na = nsize_IA1
    nb = nsize_IA2
    n = na + nb

    ALLOCATE (OutputArray(n))

    i = 1
    j = 1
    k = 1

    DO WHILE (i <= na .AND. j <= nb)
      diff = 86400000_i8*(InputArray1(i)%day - InputArray2(j)%day) + InputArray1(i)%ms - InputArray2(j)%ms
      IF (diff < 0_i8) THEN
        IF (k == 1 .OR. ((InputArray1(i)%day /= OutputArray(k - 1)%day) .OR. (InputArray1(i)%ms /= OutputArray(k - 1)%ms))) THEN
          OutputArray(k) = InputArray1(i)
          k = k + 1
        END IF
        i = i + 1
      ELSE IF (diff > 0_i8) THEN
        IF (k == 1 .OR. ((InputArray2(j)%day /= OutputArray(k - 1)%day) .OR. (InputArray2(j)%ms /= OutputArray(k - 1)%ms))) THEN
          OutputArray(k) = InputArray2(j)
          k = k + 1
        END IF
        j = j + 1
      ELSE
        IF (k == 1 .OR. ((InputArray1(i)%day /= OutputArray(k - 1)%day) .OR. (InputArray1(i)%ms /= OutputArray(k - 1)%ms))) THEN
          OutputArray(k) = InputArray1(i)
          k = k + 1
        END IF
        i = i + 1
        j = j + 1
      END IF
    END DO

    DO WHILE (i <= na)
      IF ((InputArray1(i)%day /= OutputArray(k - 1)%day) .OR. (InputArray1(i)%ms /= OutputArray(k - 1)%ms)) THEN
        OutputArray(k) = InputArray1(i)
        k = k + 1
        i = i + 1
      ELSE
        i = i + 1
      END IF
    END DO

    DO WHILE (j <= nb)
      IF ((InputArray2(j)%day /= OutputArray(k - 1)%day) .OR. (InputArray2(j)%ms /= OutputArray(k - 1)%ms)) THEN
        OutputArray(k) = InputArray2(j)
        k = k + 1
        j = j + 1
      ELSE
        j = j + 1
      END IF
    END DO

    ! do i = 1, k-1
    !   write (0,*) OutputArray(i)
    ! enddo

    nsize_OA = k - 1

  END SUBROUTINE merge2SortedAndRemoveDublicates

  SUBROUTINE remove_duplicate_intervals(starts, ends, intvls, n, indices_to_use, remaining)
    CHARACTER(len=*), INTENT(in) :: starts(:)
    CHARACTER(len=*), INTENT(in) :: ends(:)
    CHARACTER(len=*), INTENT(in) :: intvls(:)
    INTEGER, INTENT(in) :: n
    INTEGER, INTENT(out) :: indices_to_use(n)
    INTEGER, INTENT(out) :: remaining

    INTEGER :: i, j

    remaining = 1
    indices_to_use(1) = 1

    outer: DO i = 2, n
      DO j = 1, remaining
        IF (TRIM(starts(j)) == TRIM(starts(i))) THEN
          IF (TRIM(ends(j)) == TRIM(ends(i))) THEN
            IF (TRIM(intvls(j)) == TRIM(intvls(i))) THEN
              ! found a match so start looking again
              CYCLE outer
            END IF
          END IF
        END IF
      END DO
      ! no match found so add it to the output
      remaining = remaining + 1
      indices_to_use(remaining) = i
    END DO outer

  END SUBROUTINE remove_duplicate_intervals

END PROGRAM output_test
