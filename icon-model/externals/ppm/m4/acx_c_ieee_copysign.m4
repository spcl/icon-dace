dnl acx_c_ieee_copysign.m4 --- test whether copysign reliably fetches the
dnl                            sign bit from -0.0
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
dnl
dnl ACX_C_IEEE_COPYSIGN([ACTION-IF-COPYSIGN-PRESERVES-NEGZERO],
dnl                     [ACTION-IF-NOT-COPYSIGN-PRESERVES-NEGZERO])
AC_DEFUN([ACX_C_IEEE_COPYSIGN],
  [AC_CACHE_CHECK([if copysign preserves negative zero],
     [acx_cv_c_ieee_copysign],
     [AC_LANG_PUSH([C])
      saved_LIBS=$LIBS
      LIBS="$LIBS -lm"
      AC_RUN_IFELSE([AC_LANG_PROGRAM([@%:@ifdef HAVE_INTTYPES_H
@%:@include <inttypes.h>
@%:@endif
@%:@ifdef HAVE_STDINT_H
@%:@include <stdint.h>
@%:@endif
@%:@include <stdio.h>
@%:@include <stdlib.h>
@%:@include <string.h>
@%:@include <math.h>
@%:@include <fenv.h>
],[  {
#pragma STDC FENV_ACCESS ON
    float sf, vf;
    double sd, vd;
    sf = -1.0f * 0.0f;
    if ((vf = copysignf(5.0f, sf)) != -5.0f)
      return EXIT_FAILURE;
    sd = -1.0 * 0.0;
    if ((vd = copysign(5.0, sd)) != -5.0)
      return EXIT_FAILURE;
    fprintf(stderr, "copysignf(5.0f, %f) = %f\n"
                    "copysign(5.0, %f) = %f\n", sf, vf, sd, vd);
  }]
)],[acx_cv_c_ieee_copysign=yes],[acx_cv_c_ieee_copysign=no])
      LIBS=$saved_LIBS
      AC_LANG_POP([C])
     ])
   AS_IF([test x$acx_cv_c_ieee_copysign = xyes],
     [AC_DEFINE([HAVE_IEEE_COPYSIGN],[1],
        [copysign(x, -0.0) is negative (for x != 0.0)])
      $1
     ],[$2])
  ])

dnl
dnl Local Variables:
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl mode: autoconf
dnl End:
