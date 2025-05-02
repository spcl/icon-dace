dnl acx_prog_cc_posix.m4 --- test whether _POSIX_VERSION is defined
dnl
dnl Copyright  (C)  2017  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords:
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
dnl ACX_PROG_CC_POSIX([VERSION])
dnl Test if compiler is defining _POSIX_VERSION correctly to at least
dnl the numeral corresponding to VERSION. Supported values for
dnl version:
dnl     1988, 1990, 1996, 2001, 2008
dnl
AC_DEFUN([ACX_PROG_CC_POSIX],
  [m4_case([$1],
     [1988],[acx_prog_cc_posix_version_value=198808
             acx_prog_cc_posix_version_print="POSIX.1-1988"],
     [1990],[acx_prog_cc_posix_version_value=199009
             acx_prog_cc_posix_version_print="POSIX.1-1990"],
     [1996],[acx_prog_cc_posix_version_value=199506L
             acx_prog_cc_posix_version_print="POSIX.1-1996"],
     [2001],[acx_prog_cc_posix_version_value=200112
             acx_prog_cc_posix_version_print="POSIX.1-2001"],
     [2008],[acx_prog_cc_posix_version_value=200809
             acx_prog_cc_posix_version_print="POSIX.1-2008"],
     [m4_fatal([Unexpected POSIX version argument])])
   AC_LANG_PUSH([C])
   AS_VAR_PUSHDEF([acx_cv_cc_posix_support],[acx_cv_cc_posix_support]$1)
   AC_CACHE_CHECK([For conformance to ${acx_prog_cc_posix_version_print}.],
     [acx_cv_cc_posix_support],
     [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([@%:@include <unistd.h>],
        [  int n@<:@(_POSIX_VERSION >= ${acx_prog_cc_posix_version_value}L) * 2 - 1@:>@;])],
        [AS_VAR_SET([acx_cv_cc_posix_support],[yes])],
        [AS_VAR_SET([acx_cv_cc_posix_support],[no])])])
   AS_VAR_IF([acx_cv_cc_posix_support],[no],
     [AC_MSG_NOTICE(
[It seems your system does not define _POSIX_VERSION to a value
greater-or-equal ${acx_prog_cc_posix_version_value}. This is typically the case when the
compiler is instructed to make ISO C features available only,
e.g. when using gcc -std=c99])
      AC_MSG_FAILURE([${acx_prog_cc_posix_version_print} profile required])])])
