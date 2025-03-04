!> @file ppm_compare.f90
!! @brief comparison utility routines e.g. as used in sorting
!!
!! @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
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
!> collection of simple comparison functions

MODULE ppm_compare
  USE ppm_std_type_kinds, ONLY: i4, i8, dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cmp_i4, cmp_i4_indirect_i4, cmp_i4_indirect_i8, &
       cmp_i4_indirect_dp
  PUBLIC :: cmp_i8
  PUBLIC :: rcmp_i4, rcmp_i4_indirect_dp, rcmp_i4_indirect_i8
  PUBLIC :: cmp_dp, rcmp_dp
CONTAINS
  !> compare two integers
  !! @param a first number to compare
  !! @param b second number to compare
  !! @return  0 if a == b
  !!          1 if a > b
  !!         -1 if a < b
  PURE FUNCTION cmp_i4(a, b) RESULT(cmp)
    INTEGER(i4), INTENT(in) :: a, b
    INTEGER :: cmp
    cmp = MERGE(1, MERGE(-1, 0, a < b), a > b)
  END FUNCTION cmp_i4

  !> compare two integers, reverse result
  !! @param a first number to compare
  !! @param b second number to compare
  !! @return  0 if a == b
  !!         -1 if a > b
  !!          1 if a < b
  PURE FUNCTION rcmp_i4(a, b) RESULT(cmp)
    INTEGER(i4), INTENT(in) :: a, b
    INTEGER :: cmp
    cmp = MERGE(-1, MERGE(1, 0, a < b), a > b)
  END FUNCTION rcmp_i4

  !> compare two integers indirectly
  !! @param a first index to compare
  !! @param b second index to compare
  !! @param indirection table to reference for values
  !! @return  0 if indirection(a) == indirection(b)
  !!          1 if indirection(a) > indirection(b)
  !!         -1 if indirection(a) < indirection(b)
  PURE FUNCTION cmp_i4_indirect_i4(a, b, indirection) RESULT(cmp)
    INTEGER(i4), INTENT(in) :: a, b, indirection(*)
    INTEGER :: cmp

    cmp = MERGE(1, MERGE(-1, 0, indirection(a) < indirection(b)), &
         &      indirection(a) > indirection(b))
  END FUNCTION cmp_i4_indirect_i4

  !> compare two integers indirectly
  !! @param a first index to compare
  !! @param b second index to compare
  !! @param indirection table to reference for values
  !! @return  0 if indirection(a) == indirection(b)
  !!          1 if indirection(a) > indirection(b)
  !!         -1 if indirection(a) < indirection(b)
  PURE FUNCTION cmp_i4_indirect_i8(a, b, indirection) RESULT(cmp)
    INTEGER(i4), INTENT(in) :: a, b
    INTEGER(i8), INTENT(in) :: indirection(*)
    INTEGER :: cmp

    cmp = MERGE(1, MERGE(-1, 0, indirection(a) < indirection(b)), &
         &      indirection(a) > indirection(b))
  END FUNCTION cmp_i4_indirect_i8

  !> compare two integers indirectly, reverse result
  !! @param a first index to compare
  !! @param b second index to compare
  !! @param indirection table to reference for values
  !! @return  0 if indirection(a) == indirection(b)
  !!         -1 if indirection(a) > indirection(b)
  !!          1 if indirection(a) < indirection(b)
  PURE FUNCTION rcmp_i4_indirect_i8(a, b, indirection) RESULT(cmp)
    INTEGER(i4), INTENT(in) :: a, b
    INTEGER(i8), INTENT(in) :: indirection(*)
    INTEGER :: cmp

    cmp = MERGE(-1, MERGE(1, 0, indirection(a) < indirection(b)), &
         &      indirection(a) > indirection(b))
  END FUNCTION rcmp_i4_indirect_i8

  !> compare two integers
  !! @param a first number to compare
  !! @param b second number to compare
  !! @return  0 if a == b
  !!          1 if a > b
  !!         -1 if a < b
  PURE FUNCTION cmp_i8(a, b) RESULT(cmp)
    INTEGER(i8), INTENT(in) :: a, b
    INTEGER :: cmp
    cmp = MERGE(1, MERGE(-1, 0, a < b), a > b)
  END FUNCTION cmp_i8

  !> compare two double precision reals
  !! @param a first number to compare
  !! @param b second number to compare
  !! @return  0 if a == b
  !!          1 if a > b
  !!         -1 if a < b
  PURE FUNCTION cmp_dp(a, b) RESULT(cmp)
    REAL(dp), INTENT(in) :: a, b
    INTEGER :: cmp
    cmp = MERGE(1, MERGE(-1, 0, a < b), a > b)
  END FUNCTION cmp_dp

  !> compare two double precision reals, reverse result
  !! @param a first number to compare
  !! @param b second number to compare
  !! @return  0 if a == b
  !!         -1 if a > b
  !!          1 if a < b
  PURE FUNCTION rcmp_dp(a, b) RESULT(cmp)
    REAL(dp), INTENT(in) :: a, b
    INTEGER :: cmp
    cmp = MERGE(-1, MERGE(1, 0, a < b), a > b)
  END FUNCTION rcmp_dp

  !> compare two integers indirectly
  !! @param a first index to compare
  !! @param b second index to compare
  !! @param indirection table to reference for values
  !! @return  0 if indirection(a) == indirection(b)
  !!          1 if indirection(a) > indirection(b)
  !!         -1 if indirection(a) < indirection(b)
  PURE FUNCTION cmp_i4_indirect_dp(a, b, indirection) RESULT(cmp)
    INTEGER(i4), INTENT(in) :: a, b
    REAL(dp), INTENT(in) :: indirection(*)
    INTEGER :: cmp

    cmp = MERGE(1, MERGE(-1, 0, indirection(a) < indirection(b)), &
         &      indirection(a) > indirection(b))
  END FUNCTION cmp_i4_indirect_dp

  !> compare two integers indirectly, reverse result
  !! @param a first index to compare
  !! @param b second index to compare
  !! @param indirection table to reference for values
  !! @return  0 if indirection(a) == indirection(b)
  !!         -1 if indirection(a) > indirection(b)
  !!          1 if indirection(a) < indirection(b)
  PURE FUNCTION rcmp_i4_indirect_dp(a, b, indirection) RESULT(cmp)
    INTEGER(i4), INTENT(in) :: a, b
    REAL(dp), INTENT(in) :: indirection(*)
    INTEGER :: cmp

    cmp = MERGE(-1, MERGE(1, 0, indirection(a) < indirection(b)), &
         &      indirection(a) > indirection(b))
  END FUNCTION rcmp_i4_indirect_dp
END MODULE ppm_compare
!
! Local Variables:
! license-markup: "doxygen"
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
!
