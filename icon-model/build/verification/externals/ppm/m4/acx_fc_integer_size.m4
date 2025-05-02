dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
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
dnl Commentary:
dnl
dnl
dnl
dnl Code:
dnl
dnl
dnl
dnl _ACX_CHECK_SIZEOF_PROGRAM(Fortran)(type, size_to_test,
dnl                           [ACTION-IF-TRUE], [ACTION-IF-FALSE])
dnl
dnl  produce program that fails to compile if size_to_test equals
dnl  the size of type
dnl
m4_define([_ACX_LANG_CHECK_SIZEOF_INTEGRAL_TYPE_PROGRAM(Fortran)],
[AC_LANG_PROGRAM([],
[      contains]m4_ifval([$3],[
      $3])[
      integer function callee(a)
      integer a(:)
      integer i
      $1 c
      do i=1,1024,(bit_size(c) - $2)
      a(i) = a(i)
      end do
      callee = a(1)
      end function callee
      subroutine caller
      integer b(1024)
      $1 c
      b(1) = callee(b(::(bit_size(c) - $2)))
      end subroutine])])
dnl
dnl _ACX_
m4_define([_ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE_COMPILE],dnl
  [ac_lo=$acx_cv_c_char_bits ac_hi=`expr 32 \* $acx_cv_c_char_bits`
   while :; do
   AC_COMPILE_IFELSE([_AC_LANG_DISPATCH(_ACX_LANG_CHECK_SIZEOF_INTEGRAL_TYPE_PROGRAM,
      _AC_LANG, [$1], $ac_lo, $2)],
     [_AC_MSG_LOG_CONFTEST
      ac_lo=`expr $ac_lo + $acx_cv_c_char_bits`
      AS_IF([test $ac_lo -gt $ac_hi],
        [ac_lo= ac_hi= ; break])],
     [ac_lo=`expr $ac_lo / $acx_cv_c_char_bits` ; break])
   done
   AS_IF([test -z "$ac_lo"],
     [AC_MSG_FAILURE([cannot compute sizeof ($1)], 77)],
     [AS_VAR_SET([acx_fortran_Sizeof], [$ac_lo])])])
dnl
dnl _ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE_RUN(TYPE, [PROLOGUE])
m4_define([_ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE_RUN],dnl
  [AC_RUN_IFELSE([AC_LANG_PROGRAM([$2],
[      $1 :: itest
      integer :: integer_bits
      integer_bits=bit_size(itest)
      open(10,file="conftestval")
      write(10,*)integer_bits
      close(10)])],
     [ac_lo=`cat conftestval`
      AS_VAR_SET([acx_fortran_Sizeof],[`expr $ac_lo / $acx_cv_c_char_bits`])],
     [AS_VAR_SET([acx_fortran_Sizeof],[unavailable])])
   /bin/rm -f conftestval])
dnl
dnl
dnl ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE(type, [optinal-prologue],[subst?])
dnl
dnl sets shell variable acx_cv_fortran_sizeof_mangled_type
dnl to size of type in bytes and defines cpp macro named like
dnl FORTRAN_SIZEOF_type
dnl where type is a mangled version of the fortran type declaration
dnl also defines an AC_SUBST if subst is not the empty string
dnl
AC_DEFUN([ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE],
[
  AC_REQUIRE([ACX_C_CHAR_BITS])
  AS_VAR_PUSHDEF([acx_fortran_Sizeof], [acx_cv_fortran_sizeof_$1])dnl
  AS_VAR_PUSHDEF([acx_fortran_Type], [acx_cv_fortran_type_$1])dnl
  ACX_FORTRAN_CHECK_TYPE([$1])
  AS_IF([test x]AS_VAR_GET([acx_fortran_Type])[ = xyes],[
     AC_CACHE_CHECK([size of Fortran type $1],[acx_fortran_Sizeof],[
        AC_LANG_PUSH([Fortran])
        AS_IF([test "$cross_compiling" = yes],dnl
          [_ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE_COMPILE([$1],[$2])],dnl
          [_ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE_RUN([$1],[$2])])
        AC_LANG_POP([Fortran])
       ])
    m4_ifval([$3],dnl
      [AS_TR_CPP([FORTRAN_SIZEOF_$1])=AS_VAR_GET([acx_fortran_Sizeof])
       AC_SUBST(AS_TR_CPP([FORTRAN_SIZEOF_$1]))])
    AC_DEFINE_UNQUOTED(AS_TR_CPP(fortran_sizeof_$1),dnl
       AS_VAR_GET([acx_fortran_Sizeof]),dnl
       [The size of `$1', as computed by sizeof.])])
  AS_VAR_POPDEF([acx_fortran_Type])dnl
  AS_VAR_POPDEF([acx_fortran_Sizeof])dnl
])
dnl Local Variables:
dnl mode: autoconf
dnl End:
