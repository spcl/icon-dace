dnl acx_fc_have_c_loc_of_slice.m4 --- test whether the compiler supports
dnl c_loc for arrays slices (a Fortran 2008 feature)
dnl
dnl Copyright  (C)  2019  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords:
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://www.dkrz.de/redmine/projects/scales-ppm
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are  permitted provided that the following conditions are
dnl met:
dnl
dnl Redistributions of source code must retain the above copyright notice,
dnl this list of conditions and the following disclaimer.
dnl
dnl Redistributions in binary form must reproduce the above copyright
dnl notice, this list of conditions and the following disclaimer in the
dnl documentation and/or other materials provided with the distribution.
dnl
dnl Neither the name of the DKRZ GmbH nor the names of its contributors
dnl may be used to endorse or promote products derived from this software
dnl without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
dnl IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
dnl TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
dnl PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
dnl OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
dnl EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
dnl PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
dnl PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
dnl LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
dnl NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
dnl
dnl ACX_FC_HAVE_C_LOC_OF_SLICE([ACTION-IF-WORKING],[ACTION-IF-NOT-WORKING])
AC_DEFUN([ACX_FC_HAVE_C_LOC_OF_SLICE],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_CACHE_CHECK([whether $FC supports C_LOC for a slice],
     [acx_cv_fc_have_c_loc_of_slice],
     [AC_LANG_PUSH([Fortran])
      AC_LANG_CONFTEST([AC_LANG_PROGRAM(
      [
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_PTR, C_LOC, &
        C_INTPTR_T
      IMPLICIT NONE
      ! this type is interoperable but gfortran 4.4 fails to handle that
      TYPE, BIND(C) :: t
      INTEGER(C_INT) :: bar
      END TYPE t
      INTERFACE
      SUBROUTINE POSIX_EXIT(STATUS) BIND(C, NAME='exit')
      IMPORT :: C_INT
      INTEGER(C_INT), VALUE :: STATUS
      END SUBROUTINE POSIX_EXIT
      END INTERFACE
      INTEGER, TARGET :: II(10,2)
      TYPE(C_PTR) :: P_II, P_II1, P_II2
      INTEGER(C_INTPTR_T) :: I_P_II, I_P_II1, I_P_II2
      TYPE(T), TARGET :: VT(3)
      P_II = C_LOC(II)
      P_II1 = C_LOC(II(:,1))
      P_II2 = C_LOC(II(:,2))
      I_P_II = TRANSFER(P_II, I_P_II)
      I_P_II1 = TRANSFER(P_II1, I_P_II1)
      I_P_II2 = TRANSFER(P_II2, I_P_II2)
      IF (I_P_II /= I_P_II1 .OR. I_P_II1 == I_P_II2) &
      CALL POSIX_EXIT(1_C_INT)
      CALL F(VT, P_II)
      CONTAINS
      SUBROUTINE F(A, P)
      TYPE(T), TARGET, INTENT(IN) :: A(:)
      TYPE(C_PTR), INTENT(OUT) :: P
      P = C_LOC(A(1))
      END SUBROUTINE F
])])
      AC_LINK_IFELSE(,
        [AC_RUN_IFELSE(,
          [acx_cv_fc_have_c_loc_of_slice=yes],
          [acx_cv_fc_have_c_loc_of_slice=no],
          [acx_cv_fc_have_c_loc_of_slice=yes])],
        [acx_cv_fc_have_c_loc_of_slice=no])])
   m4_ifval([$1$2],
     [AS_IF([test $acx_cv_fc_have_c_loc_of_slice = yes],
        [$1],[$2])])
   AC_LANG_POP([Fortran])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
