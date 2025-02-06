dnl acx_c_limit.m4 --- check C compile-time limit
dnl
dnl Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Keywords:
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
dnl ACX_C_LIMIT(LIMIT_DEFINE,[ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
dnl
dnl
AC_DEFUN([ACX_C_LIMIT],
  [AS_VAR_PUSHDEF([acx_c_limit_name],[C_LIMIT_$1])
   AC_ARG_VAR(AS_TR_SH([acx_c_limit_name]),[C limit $1])
   AS_IF([test x"${acx_c_limit_name+set}" != xset],
     [AS_VAR_PUSHDEF([acx_cv_c_limit_name],[acx_cv_c_limit_]AS_TR_SH([$1]))
      AC_CACHE_CHECK([for invariant C limit $1],
        [acx_cv_c_limit_name],
        [AC_COMPUTE_INT([acx_cv_c_limit_name],[$1],[AC_INCLUDES_DEFAULT
@%:@include <limits.h>
])
        ],[$3])
      AS_VAR_COPY([acx_c_limit_name],[acx_cv_c_limit_name])
      AS_VAR_POPDEF([acx_cv_c_limit_name])
      $2
     ])
   AS_VAR_POPDEF([acx_c_limit_name])
  ])
dnl Local Variables:
dnl mode: autoconf
dnl license-default: "bsd"
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl End:
