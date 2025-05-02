dnl acx_fc_integer_exp_is_precise.m4 --- test whether the Fortran processor executes
dnl integer power of integer via imprecise pow() or precisely via multiplication
dnl
dnl Copyright  (C)  2022  Thomas Jahns <jahns@dkrz.de>
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
dnl ACX_FC_INTEGER_EXP_IS_PRECISE([ACTION-IF-WORKING],[ACTION-IF-NOT-WORKING])
AC_DEFUN([ACX_FC_INTEGER_EXP_IS_PRECISE],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_CACHE_CHECK([whether $FC evaluates small integer power of large integers accurately],
     [acx_cv_fc_integer_exp_is_precise],
     [AC_LANG_PUSH([Fortran])
      AC_LANG_CONFTEST([AC_LANG_PROGRAM(
      [      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
      INTEGER, PARAMETER :: i8 =  SELECTED_INT_KIND(14)
      INTEGER(i8), PARAMETER :: llim = 2097143_i8, ulim=2097151_i8
      INTEGER(i8) :: cubed, i
      INTERFACE
        SUBROUTINE POSIX_EXIT(STATUS) BIND(C, NAME='exit')
          IMPORT :: C_INT
          INTEGER(C_INT), VALUE :: STATUS
        END SUBROUTINE POSIX_EXIT
      END INTERFACE
      DO i = llim, ulim
        cubed = i**3
        IF (cubed /= i * i * i) THEN
          WRITE (0, '(a)') 'incorrect cube operation detected!'
          CALL posix_exit(1_c_int)
        END IF
      END DO
])])
      AC_LINK_IFELSE(,
        [AC_RUN_IFELSE(,
          [acx_cv_fc_integer_exp_is_precise=yes],
          [acx_cv_fc_integer_exp_is_precise=no],
          [acx_cv_fc_integer_exp_is_precise=yes])],
        [acx_cv_fc_integer_exp_is_precise=no])])
   m4_ifval([$1$2],
     [AS_IF([test $acx_cv_fc_integer_exp_is_precise = yes],
        [$1],[$2])])
   AC_LANG_POP([Fortran])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
