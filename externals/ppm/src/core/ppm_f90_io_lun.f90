!>
!! @file ppm_f90_io_lun.f90
!! @brief query next free unit number
!!
!! @copyright Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Maintainer: Thomas Jahns <jahns@dkrz.de>
! URL: https://www.dkrz.de/redmine/projects/scales-ppm
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!> module ppm_f90_io_lun
!! provides an abstraction layer so that unique unit
!! numbers can be retrieved without filling a global table
!! cross-referenced by humans
!!
!! use one of the next_free_unit* functions to find a number suitable for
!! use as logical unit number (lun) and not yet opened or reserved
!! this is further referred to as available
!!
!! @author Thomas Jahns <jahns@dkrz.de>
!! @date 2009-04-09

MODULE ppm_f90_io_lun
  USE ppm_base, ONLY: abort_ppm
  IMPLICIT NONE
  PRIVATE
  INTEGER :: default_min_lun=1, default_max_lun=100
  INTEGER, PARAMETER :: default_min_reserved_lun=1, &
       default_max_reserved_lun=10
  INTEGER :: min_reserved_lun=default_min_reserved_lun, &
       max_reserved_lun=default_max_reserved_lun
  LOGICAL, ALLOCATABLE, PRIVATE :: reserved_luns(:)
  PUBLIC add_lun_reservation, remove_lun_reservation, add_lun_reservations, &
       remove_lun_reservations, set_reserved_luns, next_free_unit_in_range, &
       next_free_unit, default_min_lun, default_max_lun
  PUBLIC :: setup_lun_table, take_down_lun_table
  CHARACTER(len=*), PARAMETER :: filename = 'ppm_f90_io_lun.f90'
CONTAINS

  !> allocate array of reserved luns, reserves unit from 1 to 10 by default
  SUBROUTINE setup_lun_table
    ALLOCATE(reserved_luns(default_min_reserved_lun:default_max_reserved_lun))
    reserved_luns=.TRUE.
  END SUBROUTINE setup_lun_table

  SUBROUTINE take_down_lun_table
    DEALLOCATE(reserved_luns)
  END SUBROUTINE take_down_lun_table

  !> excludes lun from list of automatically probed luns
  !> @param lun number of lun to be marked as unavailable
  SUBROUTINE add_lun_reservation(lun)
    INTEGER, INTENT(in) :: lun
    CALL adjust_reserved_map(lun)
    reserved_luns(lun) = .TRUE.
  END SUBROUTINE add_lun_reservation

  !> includes lun in list of automatically probed luns
  !> (i.e. cancels previous reservation)
  !> @param lun number of lun to be marked as available
  SUBROUTINE remove_lun_reservation(lun)
    INTEGER, INTENT(in) :: lun
    IF (lun .GE. min_reserved_lun .AND. lun .LE. max_reserved_lun) THEN
      reserved_luns(lun) = .TRUE.
    END IF
  END SUBROUTINE remove_lun_reservation

  !> excludes luns from list of automatically probed luns
  !> @param luns array holding numbers of luns to be marked as unavailable
  !> @param nluns size of array luns
  SUBROUTINE add_lun_reservations(luns, nluns)
    INTEGER, INTENT(in) :: nluns
    INTEGER, DIMENSION(nluns), INTENT(in) :: luns
    INTEGER :: i
    DO i=1,nluns
      CALL add_lun_reservation(luns(i))
    END DO
  END SUBROUTINE add_lun_reservations

  !> includes luns in list of automatically probed luns
  !> @param luns array holding numbers of luns to be marked as available
  !> @param nluns size of array luns
  SUBROUTINE remove_lun_reservations(luns, nluns)
    INTEGER, INTENT(in) :: nluns
    INTEGER, DIMENSION(nluns), INTENT(in) :: luns
    INTEGER :: i
    DO i=1,nluns
      CALL remove_lun_reservation(luns(i))
    END DO
  END SUBROUTINE remove_lun_reservations

  !> set list of luns excluded from automatically probed luns
  !> @param luns array holding numbers of luns to be marked as unavailable
  !> @param nluns size of array luns
  SUBROUTINE set_reserved_luns(luns, nluns)
    INTEGER, INTENT(in) :: nluns
    INTEGER, DIMENSION(nluns), INTENT(in) :: luns
    INTEGER :: i
    max_reserved_lun=MAXVAL(luns)
    min_reserved_lun=MINVAL(luns)
    DEALLOCATE(reserved_luns)
    ALLOCATE(reserved_luns(min_reserved_lun:max_reserved_lun))
    DO i=1,nluns
      CALL add_lun_reservation(luns(i))
    END DO
  END SUBROUTINE set_reserved_luns

  !> returns next available unit in range [min_lun,max_lun] in param unit
  !> @param min_lun do not search for luns less than this
  !> @param max_lun do not search for luns greater than this
  !> @param found set to .false. if no available unit in
  !>   given range could be found
  !> @param unit set to number of available unit iff found is .true.
  SUBROUTINE next_free_unit_in_range(min_lun, max_lun, found, unit)
    INTEGER, INTENT(in) :: min_lun, max_lun
    INTEGER, INTENT(out) :: unit
    LOGICAL, INTENT(out) :: found
    LOGICAL :: opened
    INTEGER :: i
    found = .FALSE.
    DO i=min_lun,max_lun
      IF ((i .LE. max_reserved_lun .AND. i .GE. min_reserved_lun)) THEN
        IF (reserved_luns(i)) CYCLE
      END IF
      INQUIRE(unit=i,opened=opened)
      IF (.NOT.opened) THEN
        unit = i
        found = .TRUE.
        EXIT
      END IF
    END DO
  END SUBROUTINE next_free_unit_in_range

  !> find currently available unit number, the range is determined
  !> by default_min_lun and default_max_lun
  !> @return next available unit
  FUNCTION next_free_unit()
    INTEGER :: next_free_unit
    LOGICAL :: found
    CALL next_free_unit_in_range(default_min_lun, &
         default_max_lun, found, next_free_unit)
    IF (.NOT. found) THEN
      CALL abort_ppm('No free logical unit available for open.', &
           source=filename, line=173)
    END IF
  END FUNCTION next_free_unit

  !> find currently available unit number,  the range is determined
  !> by default_min_lun and default_max_lun and the returned unit is also
  !> immediately reserved
  !> @return next available unit
  FUNCTION reserve_and_get_next_free_unit()
    INTEGER :: reserve_and_get_next_free_unit
    reserve_and_get_next_free_unit = next_free_unit()
    CALL add_lun_reservation(reserve_and_get_next_free_unit)
  END FUNCTION reserve_and_get_next_free_unit

  SUBROUTINE adjust_reserved_map(lun)
    INTEGER, INTENT(in) :: lun
    LOGICAL, ALLOCATABLE :: reserved_luns_copy(:)
    INTEGER :: newmin, newmax
    IF (ALLOCATED(reserved_luns)) THEN
      IF (lun < min_reserved_lun .OR. lun > max_reserved_lun) THEN
        ALLOCATE(reserved_luns_copy(min_reserved_lun:max_reserved_lun))
        reserved_luns_copy = reserved_luns
        newmin=MIN(min_reserved_lun, lun)
        newmax=MAX(max_reserved_lun, lun)
        DEALLOCATE(reserved_luns)
        ALLOCATE(reserved_luns(newmin:newmax))
        reserved_luns=.false.
        reserved_luns(min_reserved_lun:max_reserved_lun)=reserved_luns_copy
        min_reserved_lun=newmin
        max_reserved_lun=newmax
      END IF
    ELSE
      ALLOCATE(reserved_luns(lun:lun))
    END IF
  END SUBROUTINE adjust_reserved_map

! the next function is a port from ECHAM but needs the backtrace
! generating function 'finish' not currently implemented in MPIOM

END MODULE ppm_f90_io_lun
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
