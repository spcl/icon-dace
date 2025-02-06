dnl acx_tls_xlc_retry.m4 --- check for TLS storage declarator and retry for
dnl                     IBM XL which might need an extra compiler option
dnl
dnl Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Keywords: configure configure.ac autoconf MPI mpirun mpiexec
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://www.dkrz.de/redmine/projects/show/scales-ppm
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
dnl _ACX_TLS_NEWFLAG run AX_TLS with compiler option
m4_define([_ACX_TLS_RETRY],
  [saved_CFLAGS=$CFLAGS
   $2
   AS_IF([test x"$CFLAGS" = x"$saved_CFLAGS"],
     [ac_cv_tls=none],
     [AC_MSG_NOTICE([retrying with $1 added to CFLAGS])
      AS_UNSET([ac_cv_tls])
      AX_TLS(,[CFLAGS=$saved_CFLAGS])])])
dnl
dnl ACX_TLS_XLC_RETRY([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl
dnl First runs AX_TLS and retries with xlc option if this doesn't work
dnl
dnl TODO: instead of assuming C language, perform test for active AC_LANG
AC_DEFUN([ACX_TLS_XLC_RETRY],
  [AX_TLS(,[ac_cv_tls=`$CC -qversion 2>&1 | sed -n '/^IBM XL C/{
n
s/^Version: \(@<:@0-9@:>@*\).*/\1/
t print
b
: print
p
}'`
     AS_IF([test x"$ac_cv_tls" = x],
       [ac_cv_tls=`$CC -V | sed -n '/^pgcc /{
s/^pgcc \([0-9][0-9.]*\).*/\1/
p
}'`
# pgcc 18.1 and newer support TLS if switched to C11 mode
        AS_VERSION_COMPARE([$ac_cv_tls],[18.9],
          [ac_cv_tls=none],[ac_cv_tls=none],
          [_ACX_TLS_RETRY([-c11],
             [AS_CASE([" $CFLAGS "],[ -c11 ],,[CFLAGS="$CFLAGS -c11"])])])],
       [test "$ac_cv_tls" -gt 7],
       [# unless the user already set the -qtls option, add it and retry test
        _ACX_TLS_RETRY([-qtls=initial-exec],[CFLAGS=`echo "$CFLAGS" | sed -n '/.*-qtls\(=@<:@^ @:>@*\)\{0,1\}/{
p
q
}
s/$/ -qtls=initial-exec/
p
q
'`])
],[ac_cv_tls=none])])
   m4_ifnblank([$1$2],
     [AS_IF([test "$ac_cv_tls" != "none"],
        [m4_ifnblank([$1],[$1])],
        [m4_ifnblank([$2],[$2])])])
])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/show/scales-ppm"
dnl license-default: "bsd"
dnl End:
