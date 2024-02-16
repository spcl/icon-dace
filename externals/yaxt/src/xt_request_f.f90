!>
!! @file xt_request_f.f90
!! @brief xt_request-related procedures of Fortran interface
!!
!! @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
MODULE xt_requests
  USE iso_c_binding, ONLY: c_null_ptr, c_ptr, c_associated, c_int
  IMPLICIT NONE
  PRIVATE

  TYPE, BIND(C), PUBLIC :: xt_request
#ifndef __G95__
    PRIVATE
#endif
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_request

  PUBLIC :: xt_request_init, xt_request_wait, xt_request_test, xt_request_f2c, &
            xt_is_null

  INTERFACE
    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_request_f2c(request) BIND(c, name='xt_request_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_request
      IMPLICIT NONE
      TYPE(xt_request), INTENT(in) :: request
      TYPE(c_ptr) :: p
    END FUNCTION xt_request_f2c

    SUBROUTINE xt_request_wait(request) BIND(C, name='xt_request_wait')
      IMPORT :: xt_request
      TYPE(xt_request), INTENT(inout) :: request
    END SUBROUTINE xt_request_wait

  END INTERFACE

  TYPE(xt_request), PARAMETER, PUBLIC :: xt_request_null = xt_request(c_null_ptr)

  INTERFACE xt_is_null
    MODULE PROCEDURE xt_request_is_null
  END INTERFACE xt_is_null

CONTAINS

  SUBROUTINE xt_request_init(request, cptr)
    TYPE(xt_request),INTENT(out) :: request
    TYPE(c_ptr), INTENT(in) :: cptr
    request%cptr = cptr
  END SUBROUTINE xt_request_init


  SUBROUTINE xt_request_test(request, flag)
    TYPE(xt_request), INTENT(inout) :: request
    LOGICAL, INTENT(out) :: flag
    INTEGER(c_int) :: flag_c
    INTERFACE
      SUBROUTINE xt_request_test_c(request_c, flag_c) &
          BIND(C, name='xt_request_test')
        IMPORT:: c_ptr, c_int
        TYPE(c_ptr), INTENT(inout) :: request_c
        INTEGER(c_int), INTENT(out) :: flag_c
      END SUBROUTINE xt_request_test_c
    END INTERFACE
    CALL xt_request_test_c(request%cptr, flag_c)
    flag = flag_c /= 0
  END SUBROUTINE xt_request_test

  FUNCTION xt_request_is_null(request) RESULT(p)
    TYPE(xt_request), INTENT(in) :: request
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(request%cptr)
  END FUNCTION xt_request_is_null

END MODULE xt_requests
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
