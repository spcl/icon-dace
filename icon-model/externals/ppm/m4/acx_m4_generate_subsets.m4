dnl acx_m4_generate_subsets.m4 --- generate all subsets of set of strings
dnl
dnl Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
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
dnl Helper function for recursion
m4_define([_ACX_M4_GENERATE_SUBSETS],
  [m4_ifval([$1],
     [_ACX_M4_GENERATE_SUBSETS(m4_cdr($1),m4_unquote(m4_cdr($@)),[$3])dnl
m4_ifval([$2],
[_ACX_M4_GENERATE_SUBSETS(m4_cdr($1),
m4_dquote(m4_unquote($2,m4_car($1))),[$3])],
[_ACX_M4_GENERATE_SUBSETS(m4_cdr($1),m4_dquote(m4_car($1)),[$3])])],
  [m4_ifval([$2],[,])m4_dquote(m4_join([$3],m4_unquote(m4_car(m4_shift($@)))))])])
dnl
dnl ACX_M4_GENERATE_SUBSETS(SET,[SEPARATOR])
dnl
dnl generates list of all subsets of SET, where SET is a
dnl comma-seperated list of elements like [[A],[B],[C]]
dnl
AC_DEFUN([ACX_M4_GENERATE_SUBSETS],dnl
[m4_dquote(_ACX_M4_GENERATE_SUBSETS([$1],[],[$2]))],dnl
)])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
