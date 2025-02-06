! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_util_stride

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: util_stride_1d, util_stride_2d, util_get_ptrdiff

  INTERFACE
    SUBROUTINE util_stride_1d(f_out, elemsize, p1, p2) BIND(C)
      USE ISO_C_BINDING, ONLY: c_int, c_ptr
      INTEGER(c_int), INTENT(OUT) :: f_out
      INTEGER(c_int), VALUE, INTENT(IN) :: elemsize
      TYPE(c_ptr), VALUE, INTENT(IN) :: p1, p2
    END SUBROUTINE util_stride_1d

    SUBROUTINE util_stride_2d(f_out, elemsize, p1, p2, p3) BIND(C)
      USE ISO_C_BINDING, ONLY: c_int, c_ptr
      INTEGER(c_int), INTENT(OUT) :: f_out(2)
      INTEGER(c_int), VALUE, INTENT(IN) :: elemsize
      TYPE(c_ptr), VALUE, INTENT(IN) :: p1, p2, p3
    END SUBROUTINE util_stride_2d

    FUNCTION util_get_ptrdiff(a, b) RESULT(s) BIND(C, NAME='util_get_ptrdiff')
      USE ISO_C_BINDING, ONLY: c_size_t
      INTEGER(c_size_t) :: s, a, b
    END FUNCTION util_get_ptrdiff
  END INTERFACE

END MODULE mo_util_stride

