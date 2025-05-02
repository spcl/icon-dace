dnl acx_fc_have_is_contiguous.m4 --- test whether the compiler supports
dnl is_contiguous for arrays and pointers
dnl
dnl Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
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
dnl ACX_FC_HAVE_IS_CONTIGUOUS([ACTION-IF-WORKING],[ACTION-IF-NOT-WORKING])
AC_DEFUN([ACX_FC_HAVE_IS_CONTIGUOUS],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_CACHE_CHECK([whether $FC supports IS_CONTIGUOUS query],
     [acx_cv_fc_have_is_contiguous],
     [AC_LANG_PUSH([Fortran])
      AC_LANG_CONFTEST([AC_LANG_PROGRAM(
      [
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
      IMPLICIT NONE
      INTERFACE
      SUBROUTINE POSIX_EXIT(STATUS) BIND(C, NAME='exit')
      IMPORT :: C_INT
      INTEGER(C_INT), VALUE :: STATUS
      END SUBROUTINE POSIX_EXIT
      END INTERFACE
      INTEGER, TARGET :: II(10)
      INTEGER, POINTER :: P(:)
      DOUBLE PRECISION :: D(2)
      LOGICAL :: IS_II_CONT, IS_P_CONT, IS_D_CONT
      IS_D_CONT = TEST_ARG_CONT(D)
      IS_II_CONT = IS_CONTIGUOUS(II)
      P => II(::2)
      IS_P_CONT = IS_CONTIGUOUS(P)
      IF (.NOT. IS_D_CONT .OR. .NOT. IS_II_CONT .OR. IS_P_CONT) &
      CALL POSIX_EXIT(1_C_INT)
      CONTAINS
      FUNCTION TEST_ARG_CONT(A) RESULT(IS_CONT)
      DOUBLE PRECISION, INTENT(IN) :: A(:)
      LOGICAL :: IS_CONT
      IS_CONT = IS_CONTIGUOUS(A)
      END FUNCTION TEST_ARG_CONT
])])
      AC_LINK_IFELSE([],
        [AC_RUN_IFELSE([],
          [acx_cv_fc_have_is_contiguous=yes],
          [acx_cv_fc_have_is_contiguous=no],
          [acx_cv_fc_have_is_contiguous=yes])],
        [acx_cv_fc_have_is_contiguous=no])])
   m4_ifval([$1$2],
     [AS_IF([test $acx_cv_fc_have_is_contiguous = yes],
        [$1],[$2])])
   AC_LANG_POP([Fortran])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
