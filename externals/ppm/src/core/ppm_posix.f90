!>
!! @file ppm_posix.f90
!! @brief Fortran wrapper for POSIX C interface functions
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
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

MODULE ppm_posix_types
  USE ppm_std_type_kinds, ONLY: i4, i8
  IMPLICIT NONE
  PUBLIC
  TYPE ppm_stat
    SEQUENCE
    INTEGER(i4) :: st_dev, st_mode
    INTEGER(i8) :: st_ino
    INTEGER(i4) :: st_nlink, st_uid, st_gid, st_rdev
    INTEGER(i8) :: st_size
    INTEGER(i4) :: st_blksize, fill
    INTEGER(i8) :: st_blocks, st_atime, st_mtime, st_ctime
  END TYPE ppm_stat
END MODULE ppm_posix_types

MODULE ppm_posix
  USE iso_c_binding, ONLY: c_int
  USE ppm_std_type_kinds, ONLY: i4
  USE ppm_posix_types
  IMPLICIT NONE
  PRIVATE
  INCLUDE 'ftype_size.inc'
  CHARACTER(1), PARAMETER :: dir_sep='/'
  INTEGER(i4), PARAMETER :: ppm_posix_success=0
  INTEGER(i4), PARAMETER :: max_strerror_len=255
  INTEGER(i4), PARAMETER :: path_max=c_limit_path_max
  PUBLIC :: mkdir, rmdir, stat, strerror
  PUBLIC :: dir_sep, ppm_posix_success
  PUBLIC :: ppm_stat, path_max
  INTERFACE
    FUNCTION setjmp(env) BIND(c, name='setjmp') RESULT(longjmp_return)
      IMPORT :: c_int, jmp_buf_isize
      INTEGER(c_int) :: longjmp_return
      INTEGER(c_int), INTENT(out) :: env(jmp_buf_isize)
    END FUNCTION setjmp
    SUBROUTINE longjmp(env, retval) BIND(c, name='longjmp')
      IMPORT :: c_int, jmp_buf_isize
      INTEGER(c_int), INTENT(in) :: env(jmp_buf_isize)
      INTEGER(c_int), VALUE :: retval
    END SUBROUTINE longjmp
  END INTERFACE
  PUBLIC :: setjmp, longjmp
  PUBLIC :: jmp_buf_isize
CONTAINS
  SUBROUTINE mkdir(path, ierr, mode)
    CHARACTER(len=*), INTENT(in) :: path
    INTEGER(i4), INTENT(out) :: ierr
    INTEGER(i4), OPTIONAL, INTENT(in) :: mode
    INTEGER(i4) :: emode
    INTERFACE
      SUBROUTINE ppm_mkdir_f(path, mode, ierr)
        USE ppm_std_type_kinds, ONLY: i4
        CHARACTER(len=*), INTENT(in) :: path
        INTEGER(i4), INTENT(in) :: mode
        INTEGER(i4), INTENT(out) :: ierr
      END SUBROUTINE ppm_mkdir_f
    END INTERFACE
    IF (PRESENT(mode)) THEN
      emode = mode
    ELSE
      ! 511 equals 777 octal
      emode = 511_i4
    END IF
    CALL ppm_mkdir_f(path, emode, ierr)
  END SUBROUTINE mkdir

  SUBROUTINE rmdir(path, ierr)
    CHARACTER(len=*), INTENT(in) :: path
    INTEGER(i4), INTENT(out) :: ierr
    INTERFACE
      SUBROUTINE ppm_rmdir_f(path, ierr)
        USE ppm_std_type_kinds, ONLY: i4
        CHARACTER(len=*), INTENT(in) :: path
        INTEGER(i4), INTENT(out) :: ierr
      END SUBROUTINE ppm_rmdir_f
    END INTERFACE
    CALL ppm_rmdir_f(path, ierr)
  END SUBROUTINE rmdir

  SUBROUTINE stat(path, buf, ierr)
    CHARACTER(len=*), INTENT(in) :: path
    TYPE(ppm_stat), INTENT(out) :: buf
    INTEGER(i4), INTENT(out) :: ierr
    INTERFACE
      SUBROUTINE ppm_stat_f(path, buf, ierr)
        USE ppm_std_type_kinds, ONLY: i4
        USE ppm_posix_types, ONLY: ppm_stat
        CHARACTER(len=*), INTENT(in) :: path
        TYPE(ppm_stat), INTENT(out) :: buf
        INTEGER(i4), INTENT(out) :: ierr
      END SUBROUTINE ppm_stat_f
    END INTERFACE
    CALL ppm_stat_f(path, buf, ierr)
  END SUBROUTINE stat

  FUNCTION is_dir(stats)
    LOGICAL :: is_dir
    TYPE(ppm_stat), INTENT(in) :: stats
    INTERFACE
      FUNCTION ppm_is_dir_f(stats)
        USE ppm_posix_types, ONLY: ppm_stat
        TYPE(ppm_stat), INTENT(in) :: stats
        LOGICAL :: ppm_is_dir_f
      END FUNCTION ppm_is_dir_f
    END INTERFACE
    is_dir = ppm_is_dir_f(stats)
  END FUNCTION is_dir

  FUNCTION strerror(ierr, full_len)
    CHARACTER(len=max_strerror_len) :: strerror
    INTEGER(i4), INTENT(in) :: ierr
    INTEGER(i4), OPTIONAL, INTENT(out) :: full_len
    INTERFACE
      SUBROUTINE ppm_strerror_f(buf, buf_len, ierr, result_len)
        USE ppm_std_type_kinds, ONLY: i4
        CHARACTER, INTENT(inout) :: buf
        INTEGER(i4), INTENT(in) :: ierr, buf_len
        INTEGER(i4), INTENT(out) :: result_len
      END SUBROUTINE ppm_strerror_f
    END INTERFACE
    INTEGER(i4) :: n

    strerror = " "
    CALL ppm_strerror_f(strerror, max_strerror_len, ierr, n)
    IF (PRESENT(full_len)) THEN
      full_len = n
    END IF
  END FUNCTION strerror

END MODULE ppm_posix
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
