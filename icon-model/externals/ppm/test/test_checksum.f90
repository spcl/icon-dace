!> @file test_checksum.f90
!! @brief Test checksumming algorithms
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>

!
! Keywords:
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
PROGRAM test_checksum
  USE ppm_base, ONLY: abort_ppm, assertion
  USE ppm_checksum, ONLY: init_digests, hex_checksum, ppm_md5, &
       ppm_sha1, hashes
  USE ppm_std_type_kinds, ONLY: i4
  IMPLICIT NONE

  INTEGER(i4), PARAMETER :: test_val(2) = (/ 4711, 42 /)
  INTEGER, PARAMETER :: max_checksum_len=512/8
  CHARACTER(max_checksum_len * 2) :: checksum
  INTEGER, PARAMETER :: be = 1, le = 2
  CHARACTER(max_checksum_len * 2) :: &
       checksum_ref(2, 2)
  INTEGER :: endian
  CHARACTER(len=*), PARAMETER :: filename = 'test_checksum.f90'

  checksum_ref(ppm_md5, be) = '49f7c48d7b07242588e58e90a8d7360e'
  checksum_ref(ppm_md5, le) = '7bd493910ea0d0b160e844c1be457da7'
  checksum_ref(ppm_sha1, be) = 'b90fc3f44d2e0b5b1046664027d84ebbee51edca'
  checksum_ref(ppm_sha1, le) = '66b6f33a5192a93f6fa1c141f3270c20ea71196d'

  CALL init_digests
  checksum = hex_checksum("The quick brown fox jumps over the lazy dog", &
       ppm_md5)
  CALL assertion(checksum(1:hashes(ppm_md5)%size * 2) &
       == '9e107d9d372bb6826bd81d3542a419d6', &
       filename, __LINE__, 'checksum reference mismatch')
  checksum = hex_checksum(test_val, ppm_md5)
  IF (checksum &
       == checksum_ref(ppm_md5, be)(1:hashes(ppm_md5)%size * 2)) THEN
    endian = be
  ELSE IF (checksum &
       == checksum_ref(ppm_md5, le)(1:hashes(ppm_md5)%size * 2)) THEN
    endian = le
  ELSE
    CALL abort_ppm('unexpected checksum result', filename, __LINE__)
    STOP 1
  END IF
  IF (hashes(ppm_sha1)%size > 0) THEN
    checksum = hex_checksum(test_val, ppm_sha1)
    CALL assertion(checksum_ref(ppm_sha1, endian)(1:hashes(ppm_sha1)%size * 2) &
         == checksum(1:hashes(ppm_sha1)%size * 2), &
         filename, __LINE__, 'checksum reference mismatch')
  END IF
END PROGRAM test_checksum
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
