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

!>
!!                 Sets the standard I/O units.
!!

MODULE mo_io_units
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
  IMPLICIT NONE

  PUBLIC

! This paramter is taken from /usr/include/stdio.h (ANSI C standard). If problems
! with filename length appear, check the before mentioned file.

  INTEGER, PARAMETER :: filename_max = 1024

! Standard I/O-units

#ifdef hpux
  INTEGER, PARAMETER :: nerr = 7 ! error output
#else
  INTEGER, PARAMETER :: nerr = 0 ! error output
#endif
  INTEGER, PARAMETER :: nlog = 1 ! standard log file unit
  INTEGER, PARAMETER :: nnml = 2 ! standard namelist file unit
  INTEGER, PARAMETER :: nstat = 3 ! standard statistics file unit
  INTEGER, PARAMETER :: ngmt = 4 ! standard GMT output file unit
  INTEGER, PARAMETER :: nin = 5 ! standard input
  INTEGER, PARAMETER :: nout = 6 ! standard output

  INTEGER, PARAMETER, PRIVATE :: NONE = -1 ! unit given back, when nothing
  ! in the allowed range is available

  INTEGER :: nnml_output ! unit of the ASCII output that contains the
  ! namelist variables and their actual values.
!-------------------------------------------------------------------------

CONTAINS
  !
  FUNCTION find_next_free_unit(istart, istop) RESULT(iunit)
    INTEGER :: iunit
    INTEGER, INTENT(in) :: istart, istop
    !
    INTEGER :: kstart, kstop
    LOGICAL :: lfound, lopened
    INTEGER :: i
    !
    lfound = .FALSE.
    !
    kstart = istart
    kstop = istop
    IF (kstart < 10) kstart = 10
    IF (kstop <= kstart) kstop = kstart + 10
    !
    DO i = kstart, kstop
      INQUIRE (unit=i, opened=lopened)
      IF (.NOT. lopened) THEN
        iunit = i
        lfound = .TRUE.
        EXIT
      END IF
    END DO
    !
    IF (.NOT. lfound) THEN
      iunit = NONE
    END IF
    !
  END FUNCTION find_next_free_unit
  !
END MODULE mo_io_units
