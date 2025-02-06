!>
!! @file ppm_strided_extents.f90
!! @brief extend extent type to handle strided access
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
MODULE ppm_strided_extents
#ifdef USE_MPI
  USE ppm_base, ONLY: assertion
#endif
  USE ppm_extents, ONLY: extent, extent_size, extent_start, &
       extent_end, OPERATOR(==)
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif
  ! succinct description of subarrays
  TYPE strided_extent
    TYPE(extent) :: ext
    INTEGER :: stride
  END TYPE strided_extent
  !> string representation of extent size/position takes
  !! 11 decimal places (10 + sign)
  INTEGER, PARAMETER :: sext_i2s_len=11

  PUBLIC :: strided_extent, extent_size, char, extent_start, extent_end, &
       OPERATOR(==)

  INTERFACE extent_size
    MODULE PROCEDURE strided_extent_size_1d
    MODULE PROCEDURE strided_extent_size_nd
  END INTERFACE

  INTERFACE extent_start
    MODULE PROCEDURE strided_extent_start_1d
  END INTERFACE

  INTERFACE extent_end
    MODULE PROCEDURE strided_extent_end_1d
  END INTERFACE

  INTERFACE char
    MODULE PROCEDURE char_auto
  END INTERFACE

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE strided_extent_equality
  END INTERFACE

#ifdef USE_MPI
  PUBLIC :: subarray_mpi_datatype
#endif /* USE_MPI */
  CHARACTER(len=*), PARAMETER :: filename = 'ppm_strided_extents.f90'
CONTAINS

  FUNCTION strided_extent_size_1d(sext) RESULT(sext_size)
    INTEGER :: sext_size
    TYPE(strided_extent), INTENT(in) :: sext

    sext_size = (extent_size(sext%ext) + sext%stride - SIGN(1, sext%stride)) &
         / sext%stride
  END FUNCTION strided_extent_size_1d

  FUNCTION strided_extent_size_nd(sext) RESULT(sext_size)
    INTEGER :: sext_size
    TYPE(strided_extent), INTENT(in) :: sext(:)

    INTEGER :: i
    sext_size = 1
    DO i = 1, SIZE(sext)
      sext_size = sext_size * extent_size(sext(i)%ext)/sext(i)%stride
    END DO
  END FUNCTION strided_extent_size_nd

  ELEMENTAL FUNCTION strided_extent_start_1d(sext) RESULT(sext_start)
    INTEGER :: sext_start
    TYPE(strided_extent), INTENT(in) :: sext
    sext_start = extent_start(sext%ext)
  END FUNCTION strided_extent_start_1d

  ELEMENTAL FUNCTION strided_extent_end_1d(sext) RESULT(sext_end)
    INTEGER :: sext_end
    TYPE(strided_extent), INTENT(in) :: sext
    sext_end = extent_end(sext%ext) - MOD(extent_end(sext%ext), sext%stride)
  END FUNCTION strided_extent_end_1d

  ELEMENTAL FUNCTION char_auto(sext) RESULT(str)
    CHARACTER(len=3*sext_i2s_len+5) :: str
    TYPE(strided_extent), INTENT(in) :: sext
    IF (sext%ext%size /= 0) THEN
      WRITE (str, '(3(a,i0),a)') '[', sext%ext%first, ',', extent_end(sext), &
           ',', sext%stride, ']'
    ELSE
      str = '{}'
    END IF
  END FUNCTION char_auto

  ELEMENTAL FUNCTION strided_extent_equality(a, b) RESULT(l)
    TYPE(strided_extent), INTENT(in) :: a, b
    LOGICAL :: l
    l = a%ext == b%ext .AND. a%stride == b%stride
  END FUNCTION strided_extent_equality

#ifdef USE_MPI
  FUNCTION subarray_mpi_datatype(base_type, outer_shape, &
       inner_shape, rep_offset) RESULT(subarray_dt)
    INTEGER :: subarray_dt

    INTEGER, INTENT(in) :: base_type
    TYPE(extent), INTENT(in) :: outer_shape(:)
    TYPE(strided_extent), INTENT(in) :: inner_shape(:)
    INTEGER, OPTIONAL, INTENT(in) :: rep_offset(:)

    !      use MPI_handles
    INTEGER :: i, ranks, outer_ranks, offset_ranks
    INTEGER :: intermediate_handles(SIZE(inner_shape)), resize_handle
    INTEGER :: external_sizes(SIZE(outer_shape))
    INTEGER(mpi_address_kind) :: base_type_lb, base_type_extent, new_extent
    INTEGER :: ierror

    ranks = SIZE(inner_shape)
    outer_ranks = SIZE(outer_shape)
    CALL assertion(ranks <= SIZE(outer_shape) .AND. ranks >= 0, filename, &
         __LINE__, &
         msg='inner shape rank number must not exceed outer shape ranks')

    IF (ranks == 0) THEN
      ! the user inexplicably requested a single scalar?
      subarray_dt = base_type
    ELSE
      CALL mpi_type_vector((extent_size(inner_shape(1)%ext) &
           + inner_shape(1)%stride - SIGN(1, inner_shape(1)%stride)) &
           / inner_shape(1)%stride, &
           1, inner_shape(1)%stride, base_type, intermediate_handles(1), &
           ierror)
      CALL assertion(ierror == mpi_success, filename, __LINE__, &
           msg="subarray datatype creation failed")
      external_sizes(1) = extent_size(outer_shape(1))
      CALL mpi_type_get_extent(base_type, base_type_lb, base_type_extent, &
           ierror)

      DO i = 2, ranks
        CALL mpi_type_create_hvector((extent_size(inner_shape(i)%ext) &
             + inner_shape(1)%stride - SIGN(1, inner_shape(1)%stride)) &
             / inner_shape(i)%stride, 1, &
             INT(external_sizes(i - 1), mpi_address_kind) &
             * base_type_extent * INT(inner_shape(i)%stride, mpi_address_kind),&
             intermediate_handles(i - 1), intermediate_handles(i), &
             ierror)
        CALL assertion(ierror == mpi_success, filename, __LINE__, &
             msg="subarray datatype creation failed")
        external_sizes(i) = external_sizes(i - 1) * extent_size(outer_shape(i))
      END DO
      DO i = ranks + 1, outer_ranks
        external_sizes(i) = external_sizes(i - 1) * extent_size(outer_shape(i))
      END DO
      subarray_dt = intermediate_handles(ranks)

    END IF
    ! assuming the first index tuple used in the transfer is i,j,k,...
    ! let the next repeated element begin at
    ! i+offset_ranks(1),j+offset_ranks(2),k+offset_ranks(3),...
    ! when using the type in a transfer with count > 1
    offset_ranks = 0
    IF (PRESENT(rep_offset)) THEN
      offset_ranks = SIZE(rep_offset)
      CALL assertion(offset_ranks <= outer_ranks, line=__LINE__, &
           source=filename, &
           msg='input error: SIZE(rep_offset) not le SIZE(outer_shape)')
      IF (offset_ranks > 0) THEN
        new_extent = base_type_extent &
             * (INT(rep_offset(1), mpi_address_kind) &
             &  + SUM(INT(external_sizes(1:offset_ranks-1), mpi_address_kind) &
             &        * INT(rep_offset(2:offset_ranks), mpi_address_kind)))
        CALL mpi_type_create_resized(subarray_dt, base_type_lb, new_extent,  &
             resize_handle, ierror)
        CALL assertion(ierror == mpi_success, filename, __LINE__, &
             msg="subarray datatype creation failed")
        subarray_dt = resize_handle
      END IF
    END IF

    CALL mpi_type_commit(subarray_dt, ierror)
    CALL assertion(ierror == mpi_success, filename, __LINE__, &
         msg="subarray datatype creation failed")

    ! clean up intermediate types
    DO i = 1, ranks - 1
      CALL mpi_type_free(intermediate_handles(i), ierror)
      CALL assertion(ierror == mpi_success, filename, __LINE__, &
           msg="subarray datatype creation failed")
    END DO
    IF (PRESENT(rep_offset) .AND. offset_ranks > 0 .AND. ranks > 0) THEN
      CALL mpi_type_free(intermediate_handles(ranks), ierror)
      CALL assertion(ierror == mpi_success, filename, __LINE__, &
           msg="subarray datatype creation failed")
    END IF
  END FUNCTION subarray_mpi_datatype
#endif /* USE_MPI */

END MODULE ppm_strided_extents
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
