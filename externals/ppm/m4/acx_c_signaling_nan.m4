dnl acx_c_signaling_nan.m4 --- test whether we can produce an exception
dnl                            from 754 bit pattern
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
dnl ACX_C_SIGNALING_NAN([ACTION-IF-EXCEPTION-RAISED],[ACTION-IF-NOT-RAISED])
AC_DEFUN([ACX_C_SIGNALING_NAN],
  [AC_CACHE_CHECK([if platform supports signaling NaNs],
     [acx_cv_c_signaling_nan],
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
    float r;
    float f;
    uint32_t u;
    FILE *fp;
    fexcept_t fe_a, fe_b;
    if (fegetexceptflag(&fe_a, FE_ALL_EXCEPT))
      return EXIT_FAILURE;
    u = UINT32_C(0x7f800001);
    if (!(fp = fopen("conftest.data", "w+")))
      return EXIT_FAILURE;
    if (fwrite(&u, 1, sizeof(u), fp) != sizeof(u))
      return EXIT_FAILURE;
    if (fseek(fp, 0L, SEEK_SET))
      return EXIT_FAILURE;
    if (fread(&f, 1, sizeof(f), fp) != sizeof(u))
      return EXIT_FAILURE;
    r = sin(f);
    if (fegetexceptflag(&fe_b, FE_ALL_EXCEPT))
      return EXIT_FAILURE;
    if (fe_a == fe_b)
      return EXIT_FAILURE;
    fprintf(stderr, "%f, %d\n", r, memcmp(&fe_a, &fe_b, sizeof(fe_a)));
  }]
)],[acx_cv_c_signaling_nan=yes],[acx_cv_c_signaling_nan=no])
      LIBS=$saved_LIBS
      AC_LANG_POP([C])
     ])
   AS_IF([test x"$acx_cv_c_signaling_nan" = xyes],
     [AC_DEFINE([HAVE_IEEE_SIGNALING_NAN],,
        [Can we produce an NaN that raises an FPU exception?])
      $1
     ],[$2])
  ])

dnl
dnl Local Variables:
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl mode: autoconf
dnl End:
