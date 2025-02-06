dnl gcc_ac_c_char_bit.m4 --- probe number of bits in a byte
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
dnl Commentary:
dnl
dnl
dnl
dnl Code:
dnl
dnl Probe number of bits in a byte.
dnl Derived from gcc acinclude.m4 macro gcc_AC_C_CHAR_BIT.
dnl Note C89 requires CHAR_BIT >= 8.
dnl Defines CHAR_BIT if not found and sets shell variable acx_cv_c_char_bits
dnl to number of bits per character if determinable.
dnl
AC_DEFUN([ACX_C_CHAR_BITS],
[AC_CHECK_HEADERS([limits.h])
AC_CHECK_DECL([CHAR_BIT],
[ac_cv_have_decl_CHAR_BIT=yes],[ac_cv_have_decl_CHAR_BIT=no],
[[#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif]])
AC_CACHE_CHECK(number of bits in a byte, acx_cv_c_char_bits,
  [i=8
   acx_cv_c_char_bits=
   while test $i -lt 65; do
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,
       [switch(0) {
        case (unsigned char)((unsigned long)1 << $i)
             == ((unsigned long)1 << $i):
        case (unsigned char)((unsigned long)1<<($i-1))
             == ((unsigned long)1<<($i-1)):
        ; }])],
        [acx_cv_c_char_bits=$i; break])
     i=`expr $i + 1`
   done
   test -z "$acx_cv_c_char_bits" && acx_cv_c_char_bits=failed
  ])
AS_IF([test $acx_cv_c_char_bits = failed],
  [AC_MSG_ERROR(cannot determine number of bits in a byte)],
  [AS_IF([test $ac_cv_have_decl_CHAR_BIT = no],
      [AC_DEFINE_UNQUOTED(CHAR_BIT, $acx_cv_c_char_bits,
      [Define as the number of bits in a byte, if `limits.h' doesn't.])
  ])])])
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
