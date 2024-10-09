!! Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
!!
!! SPDX-License-Identifier: BSD-3-Clause
!!
PROGRAM repetitor

  USE mtime

  IMPLICIT NONE

  CHARACTER(len=64) :: case1 = 'R2/20130304T000000/20130504T030405'
  CHARACTER(len=64) :: case2 = 'R3/8000101T230516/P10D'
  CHARACTER(len=64) :: case3 = 'R/PT2M/-34560101T040404'
  CHARACTER(len=64) :: case4 = 'R4/P1Y'
  CHARACTER(len=64) :: case5 = 'R012/2013-03-04T00:00:00.546/2013-05-04T03:04:05'
  CHARACTER(len=64) :: case6 = 'R3/800-01-01T23:05:16/P10D'
  CHARACTER(len=64) :: case7 = 'R/PT2M/-347856-01-01T04:04:04'
  CHARACTER(len=64) :: case8 = 'R4/P01Y'

  CALL setCalendar(proleptic_gregorian); 
  CALL testCase(case1)
  CALL testCase(case2)
  CALL testCase(case3)
  CALL testCase(case4)
  CALL testCase(case5)
  CALL testCase(case6)
  CALL testCase(case7)
  CALL testCase(case8)

CONTAINS

  SUBROUTINE testCase(caseX)
    CHARACTER(len=*), INTENT(in) :: caseX

    CHARACTER(len=MAX_REPETITION_STR_LEN) :: repetitor
    CHARACTER(len=MAX_REPETITION_STR_LEN) :: start
    CHARACTER(len=MAX_REPETITION_STR_LEN) :: END
    CHARACTER(len=MAX_REPETITION_STR_LEN) :: duration

    TYPE(timedelta), POINTER :: d
    TYPE(datetime), POINTER :: s, e

    CHARACTER(len=MAX_DATETIME_STR_LEN) :: dstring
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: tdstring

    LOGICAL :: lrepetitor, lstart, lend, lduration

    PRINT *, 'Input ', caseX

    CALL splitRepetitionString(caseX, repetitor, start, END, duration, &
         &                            lrepetitor, lstart, lend, lduration); 
    IF (lrepetitor) THEN
      PRINT *, 'Repetitor: ', TRIM(repetitor)
      PRINT *, 'Repetitions ', getRepetitions(repetitor)
    END IF
    IF (lstart) THEN
      PRINT *, 'Start: ', TRIM(start)
      s => newDateTime(start)
      CALL datetimeToString(s, dstring)
      PRINT *, 'Start rev: ', TRIM(dstring)
      CALL deallocateDateTime(s)
    END IF
    IF (lend) THEN
      PRINT *, 'End: ', TRIM(END)
      e => newDateTime(END)
      CALL datetimeToString(e, dstring)
      PRINT *, 'End rev: ', TRIM(dstring)
      CALL deallocateDateTime(e)
    END IF
    IF (lduration) THEN
      PRINT *, 'Duration: ', TRIM(duration)
      d => newTimeDelta(duration)
      CALL timedeltaToString(d, tdstring)
      PRINT *, 'Duration rev: ', TRIM(tdstring)
      CALL deallocateTimeDelta(d)
    END IF

    PRINT *, ''

  END SUBROUTINE testCase

END PROGRAM repetitor
