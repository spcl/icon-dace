!> @file ppm_checksum.f90
!! @brief support checksumming of data in Fortran
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>

!
! Keywords: checksums, hash digests
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
#include "fc_feature_defs.inc"
MODULE ppm_checksum
  USE ppm_std_type_kinds, ONLY: i4, i8, dp, sp
  IMPLICIT NONE
  PRIVATE
  INCLUDE 'ftype_size.inc'

  TYPE digest_description
    SEQUENCE
    INTEGER(ppm_address_kind) :: helper
    INTEGER(i4) :: size, family
  END TYPE digest_description

  INTEGER, PARAMETER :: ppm_md5=1, ppm_sha1=2

  TYPE(digest_description) :: hashes(2)

  !> Perform full cryptographic checksum computation
  INTERFACE hex_checksum
    MODULE PROCEDURE hex_checksum_dp_1d
    MODULE PROCEDURE hex_checksum_dp_2d
    MODULE PROCEDURE hex_checksum_dp_3d
    MODULE PROCEDURE hex_checksum_sp_1d
    MODULE PROCEDURE hex_checksum_sp_2d
    MODULE PROCEDURE hex_checksum_sp_3d
    MODULE PROCEDURE hex_checksum_i4_1d
    MODULE PROCEDURE hex_checksum_i4_2d
    MODULE PROCEDURE hex_checksum_i4_3d
    MODULE PROCEDURE hex_checksum_i8_1d
    MODULE PROCEDURE hex_checksum_i8_2d
    MODULE PROCEDURE hex_checksum_i8_3d
    MODULE PROCEDURE hex_checksum_char
  END INTERFACE hex_checksum

  !> Compute
  !! SUM(a)/SIZE(a) + (SUM(a - SUM(a)/SIZE(a))
  !!                   + SUM((a - SUM(a)/SIZE(a))**2))/SIZE(a)
  !! REAL type array argument a.
  !!
  !! This summation is (almost)
  !! permutation-invariant and will give similar results for arrays
  !! with similar composition
  INTERFACE deviation_controlsum
    MODULE PROCEDURE deviation_controlsum_dp_1d
    MODULE PROCEDURE deviation_controlsum_dp_2d
    MODULE PROCEDURE deviation_controlsum_dp_3d
    MODULE PROCEDURE deviation_controlsum_sp_1d
    MODULE PROCEDURE deviation_controlsum_sp_2d
    MODULE PROCEDURE deviation_controlsum_sp_3d
  END INTERFACE deviation_controlsum
  PUBLIC :: hex_checksum, init_digests
  PUBLIC :: ppm_md5, ppm_sha1, digest_description
  PUBLIC :: hashes

#if defined __GNUC__ && __GNUC__ > 4
  INTERFACE
    SUBROUTINE ppm_hex_checksum_f(a,n,es,hd,hex)
      IMPORT :: digest_description
      INTEGER, INTENT(in) :: n, es
      INTEGER, INTENT(in) :: a(*)
!DEC$ ATTRIBUTES NO_ARG_CHECK :: a
!GCC$ ATTRIBUTES NO_ARG_CHECK :: a
!$PRAGMA IGNORE_TKR a
!DIR$ IGNORE_TKR a
!IBM* IGNORE_TKR a
      TYPE(digest_description), INTENT(in) :: hd
      CHARACTER(len=*) :: hex
    END SUBROUTINE ppm_hex_checksum_f
  END INTERFACE
#endif
CONTAINS
  SUBROUTINE init_digests
    INTERFACE
      SUBROUTINE ppm_describe_digest(hash_type, hd)
        IMPORT :: digest_description
        INTEGER, INTENT(in) :: hash_type
        TYPE(digest_description), INTENT(out) :: hd
      END SUBROUTINE ppm_describe_digest
    END INTERFACE
    INTEGER :: i
    DO i = ppm_md5, ppm_sha1
      CALL ppm_describe_digest(i, hashes(i))
    END DO
  END SUBROUTINE init_digests

  FUNCTION hex_checksum_dp_1d(a, hash_type) RESULT(hex)
    REAL(dp), INTENT(in) :: a(:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 8, hashes(hash_type), hex)
  END FUNCTION hex_checksum_dp_1d

  FUNCTION hex_checksum_dp_2d(a, hash_type) RESULT(hex)
    REAL(dp), INTENT(in) :: a(:,:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 8, hashes(hash_type), hex)
  END FUNCTION hex_checksum_dp_2d

  FUNCTION hex_checksum_dp_3d(a, hash_type) RESULT(hex)
    REAL(dp), INTENT(in) :: a(:,:,:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 8, hashes(hash_type), hex)
  END FUNCTION hex_checksum_dp_3d

  FUNCTION hex_checksum_sp_1d(a, hash_type) RESULT(hex)
    REAL(sp), INTENT(in) :: a(:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 4, hashes(hash_type), hex)
  END FUNCTION hex_checksum_sp_1d

  FUNCTION hex_checksum_sp_2d(a, hash_type) RESULT(hex)
    REAL(sp), INTENT(in) :: a(:,:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 4, hashes(hash_type), hex)
  END FUNCTION hex_checksum_sp_2d

  FUNCTION hex_checksum_sp_3d(a, hash_type) RESULT(hex)
    REAL(sp), INTENT(in) :: a(:,:,:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 4, hashes(hash_type), hex)
  END FUNCTION hex_checksum_sp_3d

  FUNCTION hex_checksum_i4_1d(a, hash_type) RESULT(hex)
    INTEGER(i4), INTENT(in) :: a(:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 4, hashes(hash_type), hex)
  END FUNCTION hex_checksum_i4_1d

  FUNCTION hex_checksum_i4_2d(a, hash_type) RESULT(hex)
    INTEGER(i4), INTENT(in) :: a(:,:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 4, hashes(hash_type), hex)
  END FUNCTION hex_checksum_i4_2d

  FUNCTION hex_checksum_i4_3d(a, hash_type) RESULT(hex)
    INTEGER(i4), INTENT(in) :: a(:,:,:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 4, hashes(hash_type), hex)
  END FUNCTION hex_checksum_i4_3d

  FUNCTION hex_checksum_i8_1d(a, hash_type) RESULT(hex)
    INTEGER(i8), INTENT(in) :: a(:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 8, hashes(hash_type), hex)
  END FUNCTION hex_checksum_i8_1d

  FUNCTION hex_checksum_i8_2d(a, hash_type) RESULT(hex)
    INTEGER(i8), INTENT(in) :: a(:,:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 8, hashes(hash_type), hex)
  END FUNCTION hex_checksum_i8_2d

  FUNCTION hex_checksum_i8_3d(a, hash_type) RESULT(hex)
    INTEGER(i8), INTENT(in) :: a(:,:,:)
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_f(a, SIZE(a), 8, hashes(hash_type), hex)
  END FUNCTION hex_checksum_i8_3d

  FUNCTION hex_checksum_char(a, hash_type) RESULT(hex)
    CHARACTER(*), INTENT(in) :: a
    INTEGER, INTENT(in) :: hash_type
    CHARACTER(hashes(hash_type)%size * 2) :: hex
    hex = 'h'
    CALL ppm_hex_checksum_str_f(a, LEN(a), hashes(hash_type), hex)
  END FUNCTION hex_checksum_char

  FUNCTION deviation_controlsum_dp_1d(a) RESULT(checksum)
    REAL(dp), INTENT(in) :: a(:)
    REAL(dp) :: avg, checksum, deviation, dev_sum, dev_sqr_sum
    INTEGER :: i

    avg = SUM(a) / REAL(SIZE(a), dp)
    dev_sum = 0.0_dp
    dev_sqr_sum = 0.0_dp
    DO i = 1, SIZE(a)
      deviation = a(i) - avg
      dev_sum = dev_sum + deviation
      dev_sqr_sum = dev_sqr_sum + deviation**2
    END DO
    checksum = avg + (dev_sum + dev_sqr_sum) / REAL(SIZE(a), dp)
  END FUNCTION deviation_controlsum_dp_1d

  FUNCTION deviation_controlsum_dp_2d(a) RESULT(checksum)
    REAL(dp), INTENT(in) :: a(:,:)
    REAL(dp) :: avg, checksum, deviation, dev_sum, dev_sqr_sum
    INTEGER :: i, j

    avg = SUM(a) / REAL(SIZE(a), dp)
    dev_sum = 0.0_dp
    dev_sqr_sum = 0.0_dp
    DO j = 1, SIZE(a, 2)
      DO i = 1, SIZE(a, 1)
        deviation = a(i, j) - avg
        dev_sum = dev_sum + deviation
        dev_sqr_sum = dev_sqr_sum + deviation**2
      END DO
    END DO
    checksum = avg + (dev_sum + dev_sqr_sum) / REAL(SIZE(a), dp)
  END FUNCTION deviation_controlsum_dp_2d

  FUNCTION deviation_controlsum_dp_3d(a) RESULT(checksum)
    REAL(dp), INTENT(in) :: a(:,:,:)
    REAL(dp) :: avg, checksum, deviation, dev_sum, dev_sqr_sum
    INTEGER :: i, j, k

    avg = SUM(a) / REAL(SIZE(a), dp)
    dev_sum = 0.0_dp
    dev_sqr_sum = 0.0_dp
    DO k = 1, SIZE(a, 3)
      DO j = 1, SIZE(a, 2)
        DO i = 1, SIZE(a, 1)
          deviation = a(i, j, k) - avg
          dev_sum = dev_sum + deviation
          dev_sqr_sum = dev_sqr_sum + deviation**2
        END DO
      END DO
    END DO
    checksum = avg + (dev_sum + dev_sqr_sum) / REAL(SIZE(a), dp)
  END FUNCTION deviation_controlsum_dp_3d

  FUNCTION deviation_controlsum_sp_1d(a) RESULT(checksum)
    REAL(sp), INTENT(in) :: a(:)
    REAL(sp) :: avg, checksum, deviation, dev_sum, dev_sqr_sum
    INTEGER :: i

    avg = SUM(a) / REAL(SIZE(a), sp)
    dev_sum = 0.0_sp
    dev_sqr_sum = 0.0_sp
    DO i = 1, SIZE(a)
      deviation = a(i) - avg
      dev_sum = dev_sum + deviation
      dev_sqr_sum = dev_sqr_sum + deviation**2
    END DO
    checksum = avg + (dev_sum + dev_sqr_sum) / REAL(SIZE(a), sp)
  END FUNCTION deviation_controlsum_sp_1d

  FUNCTION deviation_controlsum_sp_2d(a) RESULT(checksum)
    REAL(sp), INTENT(in) :: a(:,:)
    REAL(sp) :: avg, checksum, deviation, dev_sum, dev_sqr_sum
    INTEGER :: i, j

    avg = SUM(a) / REAL(SIZE(a), sp)
    dev_sum = 0.0_sp
    dev_sqr_sum = 0.0_sp
    DO j = 1, SIZE(a, 2)
      DO i = 1, SIZE(a, 1)
        deviation = a(i, j) - avg
        dev_sum = dev_sum + deviation
        dev_sqr_sum = dev_sqr_sum + deviation**2
      END DO
    END DO
    checksum = avg + (dev_sum + dev_sqr_sum) / REAL(SIZE(a), sp)
  END FUNCTION deviation_controlsum_sp_2d

  FUNCTION deviation_controlsum_sp_3d(a) RESULT(checksum)
    REAL(sp), INTENT(in) :: a(:,:,:)
    REAL(sp) :: avg, checksum, deviation, dev_sum, dev_sqr_sum
    INTEGER :: i, j, k

    avg = SUM(a) / REAL(SIZE(a), sp)
    dev_sum = 0.0_sp
    dev_sqr_sum = 0.0_sp
    DO k = 1, SIZE(a, 3)
      DO j = 1, SIZE(a, 2)
        DO i = 1, SIZE(a, 1)
          deviation = a(i, j, k) - avg
          dev_sum = dev_sum + deviation
          dev_sqr_sum = dev_sqr_sum + deviation**2
        END DO
      END DO
    END DO
    checksum = avg + (dev_sum + dev_sqr_sum) / REAL(SIZE(a), sp)
  END FUNCTION deviation_controlsum_sp_3d

END MODULE ppm_checksum
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
