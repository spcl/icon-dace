! mo_util_string.f90 - String conversion utilities
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_util_string
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.
  !----------------------------------------------
  ! This module holds string conversion utilities
  !----------------------------------------------
  USE mo_kind, ONLY: dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tolower        ! Conversion   : 'ABCXYZ' -> 'abcxyz'   
  PUBLIC :: toupper        ! Conversion   : 'abcxyz' -> 'ABCXYZ'
  PUBLIC :: char2          ! Conversion   : INTEGER  -> CHAR (LEN=2)
  PUBLIC :: separator      ! Format string: (/"-----...-----"/)
  PUBLIC :: int2string     ! returns integer n as a string
  PUBLIC :: real2string    ! returns real n as a string
  PUBLIC :: logical2string ! returns logical n as a string
  !-----------------
  ! module variables
  !-----------------
  CHARACTER(len=*), PARAMETER :: separator = REPEAT('-',78)
!==============================================================================
CONTAINS
!==============================================================================
  FUNCTION tolower (upper)
    !-----------------------------------
    ! Conversion: Uppercase -> Lowercase
    !-----------------------------------
    CHARACTER(LEN=*)              ,INTENT(in) :: upper
    CHARACTER(LEN=LEN_TRIM(upper))            :: tolower

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')

    DO i=1,LEN_TRIM(upper)
      IF (ICHAR(upper(i:i)) >= ICHAR('A') .AND. &
          ICHAR(upper(i:i)) <= ICHAR('Z')) THEN
        tolower(i:i) = CHAR( ICHAR(upper(i:i)) + idel )
      ELSE
        tolower(i:i) = upper(i:i)
      END IF
    END DO

  END FUNCTION tolower
!------------------------------------------------------------------------------
  FUNCTION toupper (lower)
    !-----------------------------------
    ! Conversion: Lowercase -> Uppercase
    !-----------------------------------
    CHARACTER(LEN=*)              ,INTENT(in) :: lower
    CHARACTER(LEN=LEN_TRIM(lower))            :: toupper

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('A')-ICHAR('a')

    DO i=1,LEN_TRIM(lower)
      IF (ICHAR(lower(i:i)) >= ICHAR('a') .AND. &
          ICHAR(lower(i:i)) <= ICHAR('z')) THEN
        toupper(i:i) = CHAR( ICHAR(lower(i:i)) + idel )
      ELSE
        toupper(i:i) = lower(i:i)
      END IF
    END DO

  END FUNCTION toupper
!------------------------------------------------------------------------------
  FUNCTION char2 (i, zero)
    !----------------------------------------
    ! Conversion: INTEGER -> CHARACTER(LEN=2)
    !----------------------------------------
    CHARACTER(LEN=2)                       :: char2 ! result
    INTEGER          ,INTENT(in)           :: i     ! argument
    CHARACTER        ,INTENT(in) ,OPTIONAL :: zero  ! padding instead of '0'

    INTEGER ,PARAMETER :: i0 = ICHAR ('0')

    IF (i>99 .OR. i<0) THEN
      char2 = '**'
    ELSE
      char2(1:1) = CHAR(    i/10  + i0)
      char2(2:2) = CHAR(MOD(i,10) + i0)
    ENDIF

    IF(PRESENT(zero)) THEN
      IF(char2(1:1) == '0') char2(1:1) = zero
      IF(char2(2:2) == '0') char2(2:2) = zero
    ENDIF
  END FUNCTION char2
!------------------------------------------------------------------------------
  FUNCTION int2string(n) ! returns integer n as a string (often needed in printing messages)

    INTEGER           :: n
    CHARACTER(len=10) :: int2string

    WRITE(int2string,'(I10)') n
    int2string = ADJUSTL(int2string)

  END FUNCTION int2string
!------------------------------------------------------------------------------
  CHARACTER(len=32) FUNCTION real2string(n) ! returns real n as a string (often needed in printing messages)

    REAL(dp), INTENT(in) :: n

    WRITE(real2string,'(G32.5)') n
    real2string = ADJUSTL(real2string)

  END FUNCTION real2string
!------------------------------------------------------------------------------
  CHARACTER(len=10) FUNCTION logical2string(n) ! returns integer n as a string (often needed in printing messages)

    LOGICAL, INTENT(in)    :: n

    WRITE(logical2string,'(L10)') n
    logical2string = ADJUSTL(logical2string)

  END FUNCTION logical2string

END MODULE mo_util_string
