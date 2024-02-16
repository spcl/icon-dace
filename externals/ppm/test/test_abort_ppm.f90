!>
!! @file test_abort_ppm.f90
!! @brief test functions to exchange abort handler from Fortran
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
#include "fc_feature_defs.inc"
MODULE test_abort_ppm_static_data
  USE iso_c_binding, ONLY: c_int
  USE ppm_posix, ONLY: jmp_buf_isize
  IMPLICIT NONE
  PRIVATE
  INTEGER, SAVE :: testi = 0, testcomm = 1
  CHARACTER(len=255), SAVE :: testmsg="foo", testsource="test_abort_ppm.f90"
  INTEGER(c_int) :: env(jmp_buf_isize)
  PUBLIC :: testi, testcomm, testmsg, testsource, env
END MODULE test_abort_ppm_static_data

PROGRAM test_abort_ppm
  USE iso_c_binding, ONLY: c_int
  USE ppm_base, ONLY: abort_ppm, set_abort_handler, &
       restore_default_abort_handler, assertion
  USE ppm_posix, ONLY: setjmp
  USE test_abort_ppm_static_data, ONLY: env, testi, testmsg
  IMPLICIT NONE
  INTEGER(c_int) :: rv
  CHARACTER(len=*), PARAMETER :: filename = 'test_abort_ppm.f90'
  EXTERNAL :: my_handler
  CALL assertion(testi == 0 .AND. testmsg == 'foo', filename, __LINE__, &
       'initial values mangled')
  CALL set_abort_handler(my_handler)
  env = 0_c_int
  rv = setjmp(env)
  IF (rv == 0_c_int) THEN
    CALL abort_ppm("bar", "baz", 42)
  END IF
  CALL restore_default_abort_handler
  CALL assertion(testi == 42 .AND. testmsg == 'bar', filename, __LINE__, &
       'value copy incomplete')
END PROGRAM test_abort_ppm

SUBROUTINE my_handler(comm, msg, source, line)
  USE ppm_posix, ONLY: longjmp
  USE iso_c_binding, ONLY: c_int
  USE test_abort_ppm_static_data, ONLY: env, testi, testcomm, testmsg, &
       testsource
  INTEGER, INTENT(in) :: comm, line
  CHARACTER(len=*), INTENT(in) :: msg, source
  testi = line
  testmsg = msg
  testcomm = comm
  testsource = source
  CALL longjmp(env, 1_c_int)
END SUBROUTINE my_handler
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
