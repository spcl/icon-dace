!>
!! @file xt_idxlist_collection_f.f90
!! @brief Fortran interface to xt_idxlist_collection constructors
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
MODULE xt_idxlist_collection
  USE xt_core, ONLY: xt_abort, xt_get_default_comm
  USE xt_idxlist_abstract, ONLY: xt_idxlist, xt_idxlist_c2f
  USE iso_c_binding, ONLY: c_int, c_ptr
  IMPLICIT NONE
  PRIVATE
  INTERFACE
    FUNCTION xt_idxlist_collection_new_c(idxlists, num_lists) &
         BIND(C, name='xt_idxlist_collection_new') RESULT(res_ptr)
      IMPORT:: c_ptr, c_int, xt_idxlist
      IMPLICIT NONE
      TYPE(xt_idxlist), INTENT(in) :: idxlists(*)
      INTEGER(c_int), VALUE, INTENT(in) :: num_lists
      TYPE(c_ptr) :: res_ptr
    END FUNCTION xt_idxlist_collection_new_c
  END INTERFACE

  INTERFACE xt_idxlist_collection_new
    MODULE PROCEDURE xt_idxlist_collection_new_a1d
    MODULE PROCEDURE xt_idxlist_collection_new_a2d
  END INTERFACE xt_idxlist_collection_new

  PUBLIC :: xt_idxlist_collection_new
  CHARACTER(len=*), PARAMETER :: filename = 'xt_idxlist_collection_f.f90'
CONTAINS

  FUNCTION xt_idxlist_collection_new_a1d(idxlists) RESULT(res)
    TYPE(xt_idxlist), INTENT(in) :: idxlists(:)
    TYPE(xt_idxlist) :: res
    INTEGER(c_int) :: num_idxlists_c

    IF (SIZE(idxlists) > HUGE(1_c_int)) &
         CALL xt_abort(xt_get_default_comm(), "idxlists array too large", &
         filename, __LINE__)
    num_idxlists_c = INT(SIZE(idxlists), c_int)
    res = xt_idxlist_c2f(xt_idxlist_collection_new_c(idxlists, num_idxlists_c))
  END FUNCTION xt_idxlist_collection_new_a1d

  FUNCTION xt_idxlist_collection_new_a2d(idxlists) RESULT(res)
    TYPE(xt_idxlist), INTENT(in) :: idxlists(:,:)
    TYPE(xt_idxlist) :: res
    INTEGER(c_int) :: num_idxlists_c

    IF (SIZE(idxlists) > HUGE(1_c_int)) &
         CALL xt_abort(xt_get_default_comm(), "idxlists array too large", &
         filename, __LINE__)
    num_idxlists_c = INT(SIZE(idxlists), c_int)
    res = xt_idxlist_c2f(xt_idxlist_collection_new_c(idxlists, num_idxlists_c))
  END FUNCTION xt_idxlist_collection_new_a2d

END MODULE xt_idxlist_collection
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
