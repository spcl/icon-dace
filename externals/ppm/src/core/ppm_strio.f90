!> @file ppm_strio.f90
!! @brief routines for quick, C-like string-parsing
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
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

MODULE ppm_strio_internal
  USE ppm_std_type_kinds, ONLY: i4
  IMPLICIT NONE
  PRIVATE
  INCLUDE 'ftype_size.inc'
  EXTERNAL :: ppm_get_address

  TYPE fmt_elem
    SEQUENCE
    INTEGER(i4) :: argtype, flags
    INTEGER(ppm_address_kind) :: addr
  END TYPE fmt_elem

  INTERFACE
    FUNCTION ppm_sscana(str, out, out_size, ierror) RESULT(count)
      USE ppm_std_type_kinds, ONLY: i4
      IMPORT :: fmt_elem
      CHARACTER(len=*), INTENT(in) :: str
      TYPE(fmt_elem), INTENT(inout) :: out(*)
      INTEGER(i4), INTENT(in) :: out_size
      INTEGER(i4), INTENT(out) :: ierror
      INTEGER(i4) :: count
    END FUNCTION ppm_sscana
  END INTERFACE

  PUBLIC :: ppm_sscana, ppm_address_kind, ppm_get_address, fmt_elem
END MODULE ppm_strio_internal
MODULE ppm_strio
  USE ppm_base, ONLY: abort_ppm
  USE ppm_strio_internal, ONLY: get_address => ppm_get_address, &
       ppm_sscana, fmt_elem, ppm_address_kind
  USE ppm_std_type_kinds, ONLY: i4
  IMPLICIT NONE
  PRIVATE


  INTEGER(i4), PUBLIC, PARAMETER :: arg_i4=1, arg_i8=2, arg_sp=3, arg_dp=4

  PUBLIC :: sscana, ppm_address_kind, get_address, fmt_elem
  CHARACTER(len=*), PARAMETER :: filename = 'ppm_strio.f90'
CONTAINS
  SUBROUTINE sscana(str, out, count, ierror)
    CHARACTER(len=*), INTENT(in) :: str
    TYPE(fmt_elem), INTENT(inout) :: out(:)
    INTEGER, OPTIONAL, INTENT(out) :: count, ierror

    INTEGER(i4) :: cnt, ierr

    cnt = ppm_sscana(str, out, INT(SIZE(out), i4), ierr)
    IF (PRESENT(count)) count = cnt
    IF (PRESENT(ierror)) THEN
      ierror = ierr
    ELSE
      IF (ierr /= 0) CALL abort_ppm('error in sscana', filename, 92)
    END IF
  END SUBROUTINE sscana

END MODULE ppm_strio
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! license-markup: "doxygen"
! End:
!
