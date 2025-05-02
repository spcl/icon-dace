dnl acx_fc_pointer_intent.m4 --- test wether the compiler supports giving both,
dnl                              intent and pointer statements on argument
dnl                              declarations
dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords:
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://swprojects.dkrz.de/redmine/projects/scales-ppm
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
dnl ACX_FC_POINTER_BOUNDS_REMAP([ACTION-IF-SUPPORTED],[ACTION-IF-UNSUPPORTED])
dnl
AC_DEFUN([ACX_FC_POINTER_BOUNDS_REMAP],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_CACHE_CHECK([whether $FC supports pointer bounds remapping],
     [acx_cv_fc_pointer_bounds_remap],
     [AC_LANG_PUSH([Fortran])
      AC_COMPILE_IFELSE([_ACX_FC_POINTER_BOUNDS_REMAP_TEST],
        [acx_cv_fc_pointer_bounds_remap=yes],
        [acx_cv_fc_pointer_bounds_remap=no])
      AC_LANG_POP([Fortran])])
    AS_IF([test x"$acx_cv_fc_pointer_bounds_remap" = xyes],
      [$1],
      [$2])])

AC_DEFUN([_ACX_FC_POINTER_BOUNDS_REMAP_TEST],
  [AC_LANG_PROGRAM(,
[      IMPLICIT NONE
      INTEGER, PARAMETER :: nlb=-2147483471, sz=5
      INTEGER, PARAMETER :: nub=nlb+sz-1
      INTEGER, TARGET :: a(sz, sz)
      REAL :: rnd(sz, sz)
      INTEGER, ALLOCATABLE, TARGET :: b(:, :, :)
      INTEGER, POINTER :: p(:, :)
      INTEGER :: i
      p(nlb:, nlb:) => a
      a = 7
      DO i = nlb, nub
        PRINT *, p(:, i)
      END DO
      ALLOCATE(b(nlb:nub, nlb:nub, nlb:nub))
      b = 1
      p(nlb:, nlb:) => b(:, :, nlb)
      DO i = nlb, nub
        PRINT *, p(:, i)
      END DO
      CALL RANDOM_NUMBER(rnd)
      WHERE (rnd > 0.5) p = p + 1
      DO i = nlb, nub
        PRINT *, p(:, i)
      END DO])])

AC_DEFUN([_ACX_FC_POINTER_HUGE_BOUNDS_REMAP_VERSION_TEST],
  [AS_IF([test x${ac_compiler_gnu} = xyes],
     [acx_temp=`$FC --version 2>&1 | grep '^GNU Fortran'`
      acx_fc_major=`echo "$acx_temp " | sed -e 's/.* *\([0-9]*\)[0-9.]* .*/\1/'`
      acx_fc_minor=`echo "$acx_temp " | sed -e 's/.* *[0-9]*\.\([0-9]*\)[0-9.]* .*/\1/'`
      AS_IF([test $acx_fc_major -eq 4 && test $acx_fc_minor -eq 9],
        [acx_cv_fc_pointer_huge_bounds_remap=no])],
     [acx_temp=`$FC -V 2>&1` \
      && echo "$acx_temp" \
      | grep '^Copyright.*\(The Portland Group\|NVIDIA CORPORATION\)' \
      >/dev/null],
     [AS_CASE([" $FC $FCFLAGS "],
        [*\ -Mlarge_arrays\ *],
        [],
        [AC_MSG_NOTICE([$FC will need option -Mlarge_arrays for large values of lower bounds!])])])])
dnl
dnl
AC_DEFUN([ACX_FC_POINTER_HUGE_BOUNDS_REMAP],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_CACHE_CHECK([whether $FC supports pointer bounds remapping for huge bounds values],
     [acx_cv_fc_pointer_huge_bounds_remap],
     [acx_cv_fc_pointer_huge_bounds_remap=no
      AC_LANG_PUSH([Fortran])dnl
      AC_LINK_IFELSE([_ACX_FC_POINTER_BOUNDS_REMAP_TEST],
        [AS_IF([test "$cross_compiling" = yes],
           [_ACX_FC_POINTER_HUGE_BOUNDS_REMAP_VERSION_TEST],
           [AS_IF([expr "$ac_link" : '.*/libtool --mode=link' >/dev/null],
              [acx_temp=`echo "$ac_link" | sed -e 's@\(.*/libtool --mode=\)link.*@\1@'`"execute ./conftest$ac_exeext"])
            AS_IF([_AC_RUN_LOG([LIBC_FATAL_STDERR_=1 $acx_temp >&2 || exit 1],[echo "running $acx_temp"])],
              [acx_cv_fc_pointer_huge_bounds_remap=yes],
              [_AC_MSG_LOG_CONFTEST])])])
      AC_LANG_POP([Fortran])])
    AS_IF([test x"$acx_cv_fc_pointer_huge_bounds_remap" = xyes],
      [$1],
      [$2])])
dnl
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
