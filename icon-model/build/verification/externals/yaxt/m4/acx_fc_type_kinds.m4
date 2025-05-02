dnl acx_fc_type_kinds.m4 --- check availability of a fortran type (kind)
dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: configure configure.ac autotools fortran type check
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
# ACX_FORTRAN_CHECK_TYPE(TYPE,
#		     [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#		     [USEINCLUDE])
# ------------------------------------------------------------
# Check whether the type TYPE is supported by the system, maybe via
# the provided includes and module use statements.
#
# Analogous to C specific macro AC_CHECK_TYPE in the newer syntax version
AC_DEFUN([ACX_FORTRAN_CHECK_TYPE],
[
  AS_VAR_PUSHDEF([acx_fortran_Type], [acx_cv_fortran_type_$1])dnl
  AC_CACHE_CHECK([for Fortran type $1], [acx_fortran_Type],
    [AC_LANG_PUSH([Fortran])
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([$4],
[      $1 a(2, 2)
      a = reshape(a, (/ 2, 2 /))])],
       [AS_VAR_SET([acx_fortran_Type], [yes])],
       [AS_VAR_SET([acx_fortran_Type], [no])])
     AC_LANG_POP([Fortran])])
  m4_ifval([$2$3],dnl
    [AS_IF([test x]AS_VAR_GET([acx_fortran_Type])[ = xyes], [$2], [$3])])dnl
  AS_VAR_POPDEF([acx_fortran_Type])dnl
])
# ACX_FORTRAN_CHECK_TYPES(TYPES,
#		         [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#		         [INCLUDES = DEFAULT-INCLUDES])
# --------------------------------------------------------
# TYPES is an m4 list.  There are no ambiguities here, we mean the newer
# AC_CHECK_TYPE.
AC_DEFUN([ACX_FORTRAN_CHECK_TYPES],
[m4_foreach([AC_Type], [$1],
  [ACX_FORTRAN_CHECK_TYPE(AC_Type,
    [AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_FORTRAN_[]AC_Type), 1,
			[Define to 1 if the system has the type `]AC_Type['.])
$2],
    [$3],
    [$4])])])
dnl
dnl _ACX_CHECK_TYPE_DEFINE_HELPER
dnl
m4_define([_ACX_CHECK_TYPE_DEFINE_HELPER],
  [AS_VAR_PUSHDEF([shellVar], [$2])dnl
   ACX_FORTRAN_CHECK_TYPES([$1],
     [AS_VAR_SET(shellVar, 1)])
   AS_VAR_POPDEF([shellVar])])
dnl ACX_FORTRAN_USUAL_TYPE_KINDS tests to see if the following fortran datatypes are
dnl supported: INTEGER(1), INTEGER(2), INTEGER(4), INTEGER(8),
dnl            REAL(4), REAL(8), REAL(16),
dnl            DOUBLE_COMPLEX, COMPLEX(4), COMPLEX(8), COMPLEX(16)
dnl
AC_DEFUN([ACX_FORTRAN_USUAL_TYPE_KINDS],dnl
[
_ACX_CHECK_TYPE_DEFINE_HELPER([integer(kind=1)], [FORT_INT1])
_ACX_CHECK_TYPE_DEFINE_HELPER([integer(kind=2)], [FORT_INT2])
_ACX_CHECK_TYPE_DEFINE_HELPER([integer(kind=4)], [FORT_INT4])
_ACX_CHECK_TYPE_DEFINE_HELPER([integer(kind=8)], [FORT_INT8])
dnl
_ACX_CHECK_TYPE_DEFINE_HELPER([real(kind=4)], [FORT_REAL4])
_ACX_CHECK_TYPE_DEFINE_HELPER([real(kind=8)], [FORT_REAL8])
_ACX_CHECK_TYPE_DEFINE_HELPER([real(kind=16)], [FORT_REAL16])
dnl
_ACX_CHECK_TYPE_DEFINE_HELPER([double complex], [DOUBLE_COMPLEX])
_ACX_CHECK_TYPE_DEFINE_HELPER([complex(kind=4)], [FORT_COMPLEX8])
_ACX_CHECK_TYPE_DEFINE_HELPER([complex(kind=8)], [FORT_COMPLEX16])
_ACX_CHECK_TYPE_DEFINE_HELPER([complex(kind=16)], [FORT_COMPLEX32])
])dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
