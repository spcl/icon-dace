dnl
dnl Copyright  (C)  2019  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: configure configure.ac autotools
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
dnl SPDX-License-Identifier: BSD-3-Clause
dnl
dnl _ACX_CHECK_INTEGER_EXPRESSION_SOURCE(Fortran)([DECLARATION],
dnl    [TEST_EXPRESSION],[VALUE_TO_CHECK_FOR_EQUALITY],[EXTRA_DECLARATION])
dnl
dnl
dnl                           [ACTION-IF-TRUE], [ACTION-IF-FALSE])
dnl
dnl  produce source that fails to compile if TEST_EXPRESSION equals
dnl  VALUE_TO_CHECK_FOR_EQUALITY
dnl
m4_define([_ACX_LANG_CHECK_INTEGER_EXPRESSION_SOURCE(Fortran)],
[AC_LANG_SOURCE([      module conftest]m4_ifval([$4],[
      $4
])[      contains
      integer function callee(a)
      integer a(:)
      integer i
      $1
      do i=1,1024,($2 - $3)
      a(i) = a(i)
      end do
      callee = a(1)
      end function callee
      subroutine caller
      integer b(1024)
      $1
      b(1) = callee(b(::($2 - $3)))
      end subroutine
      end module])])
m4_define([_ACX_LANG_CHECK_INTEGER_EXPRESSION_SOURCE],
  [_AC_LANG_DISPATCH([$0], _AC_LANG, $@)])
dnl
dnl _ACX_FORTRAN_CHECK_INTEGER_EXPRESSION_COMPILE(
dnl                           [DECL],[TEST_EXPRESSION],[EXTRA-DECL],
dnl                           [LOW-BOUND],[HIGH-BOUND],[INCREMENT]
dnl                           [VAR-TO-SET-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl
dnl
m4_define([_ACX_FORTRAN_CHECK_INTEGER_EXPRESSION_COMPILE],dnl
  [ac_lo=m4_ifval([$4],[$4],) ac_hi=m4_ifval([$5],[$5],)
   while :; do
   AC_COMPILE_IFELSE([_ACX_LANG_CHECK_INTEGER_EXPRESSION_SOURCE([$1], [$2], [$ac_lo], [$3])],
     [_AC_MSG_LOG_CONFTEST
      ac_lo=`expr $ac_lo + $6`
      AS_IF([test $ac_lo -gt $ac_hi],
        [ac_lo= ac_hi= ; break])],
     [break])
   done
   AS_IF([test -z "$ac_lo"],
     m4_ifval([$8],[$8],[AC_MSG_FAILURE([cannot determine integer expression at compile time], 77)]),
     m4_ifval([$7],[AS_VAR_SET([$7], [$ac_lo])]))])
dnl
dnl _ACX_FORTRAN_CHECK_INTEGER_EXPRESSION_RUN(
dnl             [DECL],[TEST_EXPRESSION],[EXTRA-DECL],
dnl             [VAR-TO-SET-IF-FOUND], [ACTION-IF-NOT-FOUND])
m4_define([_ACX_FORTRAN_CHECK_INTEGER_EXPRESSION_RUN],dnl
  [AC_RUN_IFELSE([AC_LANG_PROGRAM([],m4_ifval([$3],[      $3])
[      $1
      open(10,file="conftestval")
      write(10,*)$2
      close(10)])],
     [ac_lo=`sed -e 's/^ *//' conftestval`
      AS_VAR_SET([$4],[$ac_lo])],
     m4_ifval([$5],[$5],[AC_MSG_FAILURE([failed to determine integer expression at run-time])]))
   /bin/rm -f conftestval])
dnl
dnl ACX_FORTRAN_CHECK_INTEGER_EXPRESSION(
dnl             [DECL],[TEST_EXPRESSION],[EXTRA-DECL],
dnl             [LOW-BOUND],[HIGH-BOUND],[INCREMENT]
dnl             [VAR-TO-SET-IF-FOUND], [ACTION-IF-NOT-FOUND])
AC_DEFUN([ACX_FORTRAN_CHECK_INTEGER_EXPRESSION],
  [AS_IF([test x"$cross_compiling" = xyes],
     [_ACX_FORTRAN_CHECK_INTEGER_EXPRESSION_COMPILE(
        [$1],[$2],[$3],[$4],[$5],[$6],
        [$7],[$8])],
     [_ACX_FORTRAN_CHECK_INTEGER_EXPRESSION_RUN([$1],[$2],[$3],
        [$7],[$8])])])

dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
