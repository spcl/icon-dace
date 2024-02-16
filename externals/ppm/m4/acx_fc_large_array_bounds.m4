dnl acx_fc_large_array_bounds.m4 --- test whether the compiler arrays with large
dnl                                  values for the bounds
dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
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
dnl  ACX_FC_LARGE_ARRAY_BOUNDS([ACTION-IF-UNCONDITIONALLY-SUPPORTED],
dnl                            [ACTION-IF-SUPPORTED-WITH-EXTRA-OPTION],
dnl                            [ACTION-IF-UNSUPPORTED],
dnl                            [VARIABLE-TO-SET-TO-EXTRA-OPTION],
dnl                            [WHAT-TO-ASSUME-IN-CROSS-COMPILATION-MODE])
AC_DEFUN([ACX_FC_LARGE_ARRAY_BOUNDS],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_CACHE_CHECK([whether $FC supports small arrays with large bounds],
     [acx_cv_fc_large_array_bounds],
     [AS_IF([test $cross_compiling = yes],
        [acx_cv_fc_large_array_bounds=m4_default([$5],[no])],
        [AC_LANG_PUSH([Fortran])
         _AC_RUN_IFELSE([AC_LANG_PROGRAM([],[      IMPLICIT NONE
      INTEGER :: strt(3), ends(3)
      REAL, ALLOCATABLE :: a(:, :, :), b(:, :, :)
      LOGICAL, ALLOCATABLE :: m(:, :, :)
      strt(1) = 1456249828
      strt(2) = -1955192049
      strt(3) = 1451284827
      ends(1) = 1456250114
      ends(2) = -1955191456
      ends(3) = 1451284985
      ALLOCATE(a(strt(1):ends(1), strt(2):ends(2), strt(3):ends(3)), &
           b(strt(1):ends(1), strt(2):ends(2), strt(3):ends(3)), &
           m(strt(1):ends(1), strt(2):ends(2), strt(3):ends(3)))
      CALL random_number(a)
      b = a
      m = a > 0.5])],
           [acx_cv_fc_large_array_bounds=yes],
           [FCFLAGS_save=$FCFLAGS
            FCFLAGS="${FCFLAGS+$FCFLAGS }-Mlarge_arrays"
            _AC_RUN_IFELSE(,[acx_cv_fc_large_array_bounds=-Mlarge_arrays],
              [acx_cv_fc_large_array_bounds=no])
            FCFLAGS=$FCFLAGS_save])
      AC_LANG_POP([Fortran])])])
    AS_IF([test x"$acx_cv_fc_large_array_bounds" = xyes],
      [$1],
      [test x"$acx_cv_fc_large_array_bounds" != xno],
      [m4_ifval([$4],[AS_VAR_SET([$4],[$acx_cv_fc_large_array_bounds]) ])m4_ifval([$2],[; $2])],[$3])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
