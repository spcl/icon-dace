!> @file xt_idxsection_f.f90
!! @brief Fortran interface to yaxt implementation
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
MODULE xt_idxsection
  USE iso_c_binding, ONLY: c_int, c_ptr
  USE xt_core, ONLY: xt_int_kind, xt_abort, i2, i4, i8
  USE xt_idxlist_abstract, ONLY: xt_idxlist, xt_idxlist_c2f
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xt_idxsection_new, xt_idxfsection_new

  INTERFACE xt_idxsection_new
    MODULE PROCEDURE xt_idxsection_new_a
    MODULE PROCEDURE xt_idxsection_new_i2
    MODULE PROCEDURE xt_idxsection_new_i4
    MODULE PROCEDURE xt_idxsection_new_i8
  END INTERFACE xt_idxsection_new

  INTERFACE
    FUNCTION xt_idxsection_new_c(start, num_dimensions, global_size, &
         local_size, local_start) BIND(c, name='xt_idxsection_new') &
         RESULT(idxsection)
      IMPORT :: c_int, c_ptr, xt_idxlist, xt_int_kind
      INTEGER(xt_int_kind), VALUE, INTENT(in) :: start
      INTEGER(c_int), VALUE, INTENT(in) :: num_dimensions
      INTEGER(xt_int_kind), INTENT(in) :: global_size(num_dimensions), &
           local_start(num_dimensions)
      INTEGER(c_int), INTENT(in) :: local_size(num_dimensions)
      TYPE(c_ptr) :: idxsection
    END FUNCTION xt_idxsection_new_c
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: filename = 'xt_idxsection_f.f90'
CONTAINS

  FUNCTION xt_idxsection_new_a(start, global_size, local_size, local_start) &
       RESULT(idxsection)
    INTEGER(xt_int_kind), INTENT(in) :: start, local_start(:), global_size(:)
    INTEGER, INTENT(in) :: local_size(:)
    TYPE(xt_idxlist) :: idxsection
    INTEGER :: num_dimensions
    INTEGER(c_int) :: num_dimensions_c
    num_dimensions = SIZE(global_size)
    IF (SIZE(local_size) /= num_dimensions &
         .OR. SIZE(local_start) /= num_dimensions) &
         CALL xt_abort("non-matching array sizes", filename, __LINE__)
    num_dimensions_c = INT(num_dimensions, c_int)
    idxsection = xt_idxlist_c2f(&
         xt_idxsection_new_c(start, num_dimensions_c, global_size, &
         &                   INT(local_size, c_int), local_start))
  END FUNCTION xt_idxsection_new_a

  FUNCTION xt_idxsection_new_i2(start, num_dimensions, global_size, &
       local_size, local_start) RESULT(idxsection)
    INTEGER(i2), INTENT(in) :: num_dimensions
    INTEGER(xt_int_kind), INTENT(in) :: start, global_size(num_dimensions), &
         local_start(num_dimensions)
    INTEGER, INTENT(in) :: local_size(num_dimensions)
    TYPE(xt_idxlist) :: idxsection
    INTEGER(c_int) :: num_dimensions_c

    num_dimensions_c = INT(num_dimensions, c_int)
    idxsection = xt_idxlist_c2f(&
         xt_idxsection_new_c(start, num_dimensions_c, global_size, &
         &                   INT(local_size, c_int), local_start))
  END FUNCTION xt_idxsection_new_i2

  FUNCTION xt_idxsection_new_i4(start, num_dimensions, global_size, &
       local_size, local_start) RESULT(idxsection)
    INTEGER(i4), INTENT(in) :: num_dimensions
    INTEGER(xt_int_kind), INTENT(in) :: start, global_size(num_dimensions), &
         local_start(num_dimensions)
    INTEGER, INTENT(in) :: local_size(num_dimensions)
    TYPE(xt_idxlist) :: idxsection
    INTEGER(c_int), PARAMETER :: dummy = 1
    INTEGER(c_int) :: num_dimensions_c

    IF (num_dimensions > HUGE(dummy)) &
         CALL xt_abort("num_dimensions too large", filename, __LINE__)
    num_dimensions_c = INT(num_dimensions, c_int)
    idxsection = xt_idxlist_c2f(&
         xt_idxsection_new_c(start, num_dimensions_c, global_size, &
         &                   INT(local_size, c_int), local_start))
  END FUNCTION xt_idxsection_new_i4

  FUNCTION xt_idxsection_new_i8(start, num_dimensions, global_size, &
       local_size, local_start) RESULT(idxsection)
    INTEGER(i8), INTENT(in) :: num_dimensions
    INTEGER(xt_int_kind), INTENT(in) :: start, global_size(num_dimensions), &
         local_start(num_dimensions)
    INTEGER, INTENT(in) :: local_size(num_dimensions)
    TYPE(xt_idxlist) :: idxsection
    INTEGER(c_int), PARAMETER :: dummy = 1
    INTEGER(c_int) :: num_dimensions_c

    IF (num_dimensions > HUGE(dummy)) &
         CALL xt_abort("num_dimensions too large", filename, __LINE__)
    num_dimensions_c = INT(num_dimensions, c_int)
    idxsection = xt_idxlist_c2f(&
         xt_idxsection_new_c(start, num_dimensions_c, global_size, &
         &                   INT(local_size, c_int), local_start))
  END FUNCTION xt_idxsection_new_i8

  !> Fortran style version of \ref xt_idxsection_new. Compared to xt_idxsection_new, here
  !! the elements of the vector arguments are used in reversed order and the values of the elements
  !! of local_start are shifted by one. This means that, e.g.,  to start your local section with
  !! the global start index you have to set all coords in local_start to ONE instead of ZERO (as it would be required
  !! in xt_idxsection_new). The local section must be contained within the global index space.
  !! @param[in] start       start index of the global index space
  !! @param[in] global_size vector holding the global size for each dimension
  !! @param[in] local_size  vector holding the local section size for each dimension
  !! @param[in] local_start vector holding the coordinates of the section start; lowest coodinate is ONE for each dimension
  FUNCTION xt_idxfsection_new(start, global_size, local_size, local_start) &
       RESULT(idxfsection)
    INTEGER(xt_int_kind), INTENT(in) :: start, global_size(:), local_start(:)
    INTEGER, INTENT(in) :: local_size(:)
    TYPE(xt_idxlist) :: idxfsection

    INTEGER :: idim, ndim
    LOGICAL :: err_state

    ndim = SIZE(global_size)
    IF (SIZE(local_size) /= ndim .OR. SIZE(local_start) /= ndim) &
         CALL xt_abort("non-matching array sizes", filename, __LINE__)

    ! check if local indices are a subset of global indices:
    err_state = .FALSE.
    DO idim = 1, ndim
      err_state = err_state .OR. (local_start(idim) < 1) .OR. &
           (local_start(idim) + local_size(idim) - 1 >  global_size(idim))
    ENDDO
    IF (err_state) CALL xt_abort("local indices out of global index space", &
         filename, __LINE__)

    ! Fortran style map of mult-dim coords to indices:
    ! => reverse order of dimensions and coords starting at 1 (instead of 0 as in c)
    idxfsection = xt_idxsection_new(start, &
         global_size(ndim:1:-1), &
         local_size(ndim:1:-1), &
         local_start(ndim:1:-1) - 1_xt_int_kind )

  END FUNCTION xt_idxfsection_new

END MODULE xt_idxsection
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
