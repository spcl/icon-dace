dnl acx_c_strtod_denormal.m4 --- test strtod works reliably for denormals
dnl
dnl Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: http://
dnl
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
dnl
m4_define([_ACX_C_STRTOLDF_DENORMAL],
  [AC_CACHE_CHECK([if $1 can read denormals],
     [acx_cv_c_$1_denormal],
     [AC_LANG_PUSH([C])
      saved_LIBS=$LIBS
      LIBS="$LIBS ${LIBM--lm}"
      AC_RUN_IFELSE([AC_LANG_PROGRAM([@%:@include <errno.h>
@%:@include <float.h>
@%:@include <math.h>
@%:@include <stdlib.h>
@%:@include <stdio.h>
@%:@include <string.h>

@%:@pragma float_control (precise,on)
@%:@pragma fenv_access (on)
],[$5
  char *endptr;
  errno = 0;
  $2 val = $1(str, &endptr);
  if ((errno == ERANGE && (val == $3 || val == -$3))
      || (errno != 0 && val == 0.0) || (endptr == str)) {
    if (errno)
      perror("parse error");
    else
      fputs("parse error\n", stderr);
    return EXIT_FAILURE;
  }
  fprintf(stderr, "%.17e\n", val);
  int okay = val != 0.0$4;
  return okay ? EXIT_SUCCESS : EXIT_FAILURE;]
)],[acx_cv_c_$1_denormal=yes],[acx_cv_c_$1_denormal=no])
      LIBS=$saved_LIBS
      AC_LANG_POP([C])])
   AS_IF([test $acx_cv_c_$1_denormal = yes],[$6],[$7])])

dnl
dnl ACX_C_STRTOD_DENORMAL([ACTION-IF-STRTOD-READS-DENORMAL],
dnl                       [ACTION-IF-STRTOD-FAILS-DENORMAL-READ])
AC_DEFUN([ACX_C_STRTOD_DENORMAL],
  [_ACX_C_STRTOLDF_DENORMAL([strtod],[double],[HUGE_VAL],[],
     [@%:@if DBL_MIN_EXP == -1021
  static const char str@<:@@:>@ = "0.20432305856979998E-307";
@%:@else
@%:@error "unexpected floating-point characteristics!"
@%:@endif],[$1],[$2])])
AC_DEFUN([ACX_C_STRTOF_DENORMAL],
  [_ACX_C_STRTOLDF_DENORMAL([strtof],[float],[HUGE_VALF],
     [@%:@if FLT_MIN_EXP == -125
  static const char str@<:@@:>@ = "5.87747175e-39";
@%:@else
@%:@error "unexpected floating-point characteristics!"
@%:@endif],[$1],[$2])])

dnl
dnl Local Variables:
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl mode: autoconf
dnl End:
