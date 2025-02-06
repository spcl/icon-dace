dnl acx_fc_attribute_contiguous.m4 --- test wether the compiler supports
dnl                              pointer and assumed shape declarations
dnl                              that include the CONTIGUOUS attribute
dnl
dnl Copyright  (C)  2015  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1
dnl Keywords:
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://www.dkrz.de/redmine/projects/show/scales-ppm
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
dnl SPDX-License-Identifier: BSD-3-Clause
dnl
dnl
AC_DEFUN([ACX_FC_ATTRIBUTE_CONTIGUOUS],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_CACHE_CHECK([whether $FC supports CONTIGUOUS attribute],
     [acx_cv_fc_attribute_contiguous],
     [AC_LANG_PUSH([Fortran])
      AC_COMPILE_IFELSE([      SUBROUTINE f(x)
      REAL, CONTIGUOUS, TARGET, INTENT(IN) :: x(:, :)
      REAL, CONTIGUOUS, POINTER :: p(:, :)
      INTEGER :: i, j
      p => x
      DO j = 1, SIZE(p, 2)
      DO i = 1, SIZE(p, 1)
      p(i, j) = p(i, j) * 0.1
      END DO
      END DO
      END SUBROUTINE f

dnl some compilers have problems when dereferencing contiguous arguments in a
dnl old-style call (e.g. MPI functions)
      SUBROUTINE bar(b)
      REAL, POINTER, CONTIGUOUS :: b(:)
      CALL foo(b(1))
      END SUBROUTINE bar

      SUBROUTINE foo(x)
      REAL, INTENT(in) :: x
      WRITE (0, *) x
      END SUBROUTINE foo
],
        [acx_cv_fc_attribute_contiguous=yes],
        [acx_cv_fc_attribute_contiguous=no])
      AC_LANG_POP([Fortran])])
    AS_IF([test x"$acx_cv_fc_attribute_contiguous" = xyes],
      [$1],
      [$2])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
