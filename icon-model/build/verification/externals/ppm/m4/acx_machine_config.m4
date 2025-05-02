dnl acx_machine_config.m4 --- load configuration file containing
dnl                           Makefile-like shell code
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
dnl Commentary:
dnl
dnl
dnl
dnl Code:
dnl
dnl
dnl AC_GET_MH
dnl
dnl Load machine specific compiler options
dnl
AC_DEFUN([AC_GET_MH],
[AC_MSG_CHECKING([for machine dependent configuration])
changequote(,)
cat > confsed <<EOF
:redo
/\\\\\$/ {
    s/\\\\\$//g
    h
    s/.*//
    n
    s/^[ 	]*//
    H
    x
    s/\\n//
    bredo
}
/^[ 	]*[A-Za-z_][0-9A-Za-z_]*[	 ]*=/ {
    s/[ 	]*=[ 	]*/="/
    /="/s/$/"/
    p
    s/^\\([ 	]*\\)\\([A-Za-z_][0-9A-Za-z_]*\\)[	 ]*=.*/\\1export \\2/
}
/^ *#/{
    d
}
p

:next
EOF
changequote([,])
sed -n -f confsed $1 > conftest
. ./conftest
/bin/rm -f confsed conftest
if test "$1" != 0 ; then
    AC_MSG_RESULT($1)
else
    AC_MSG_RESULT(unavailable)
fi
])dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
