dnl acx_fortran2003_kinds.m4 --- check availability of Fortran 2003/C
dnl                              interoperability type
dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: configure configure.ac autoconf Fortran ISO_C type kinds
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
dnl identical to AC_CHECK_SIZEOF but does not add a define
AC_DEFUN([ACX_CHECK_SIZEOF_NODEF],
  [# The cast to long int works around a bug in the HP C Compiler
# version HP92453-01 B.11.11.23709.GP, which incorrectly rejects
# declarations like `int a3[[(sizeof (unsigned char)) >= 0]];'.
# This bug is HP SR number 8606223364.
_AC_CACHE_CHECK_INT([size of $1], [AS_TR_SH([ac_cv_sizeof_$1])],
  [(long int) (sizeof ($1))],
  [AC_INCLUDES_DEFAULT([$3])],
  [if test "$AS_TR_SH([ac_cv_type_$1])" = yes; then
     AC_MSG_FAILURE([cannot compute sizeof ($1)], 77)
   else
     AS_TR_SH([ac_cv_sizeof_$1])=0
   fi])])

dnl
dnl
dnl ACX_FORTRAN_TYPE_KIND(FORTRAN_TYPE, FORTRAN_KIND, C_CORRESPONDENCE,
dnl                       [PROLOGUE], [C_PROLOGUE],[ACTION-IF-SUBST],
dnl                       [ACTION-IF-NOT-SUBST)
dnl
dnl First tries whether a fortran declaration
dnl FORTRAN_TYPE(kind=FORTRAN_KIND) compiles and if not tries to find
dnl a fortran kind constant that corresponds to C type
dnl C_CORRESPONDENCE.
dnl PROLOGUE could be something like 'use iso_c_binding'
dnl
AC_DEFUN([ACX_FORTRAN_TYPE_KIND],
 [m4_ifval([$2],[ACX_FORTRAN_CHECK_TYPE([$1(kind=$2)],
                   [acx_search_kind_subst=false],
                   [acx_search_kind_subst=true], [$4])],
           [acx_search_kind_subst=true])
  AS_IF([test $acx_search_kind_subst = true],
    [ACX_CHECK_SIZEOF_NODEF([$3],,[$5])
     AS_VAR_PUSHDEF([acx_c_size], [ac_cv_sizeof_$3])dnl
     AS_VAR_COPY([acx_temp3],[acx_c_size])
     m4_case([$1],
       [integer],
       [ac_kind_lo=1 ac_kind_hi=32
        while :; do
          ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE([$1(kind=$ac_kind_lo)])
          AS_VAR_PUSHDEF([acx_fortran_ikind_Type],
            [acx_cv_fortran_type_$1(kind=$ac_kind_lo)])dnl
          AS_VAR_PUSHDEF([acx_fortran_Sizeof],
            [acx_cv_fortran_sizeof_$1(kind=$ac_kind_lo)])dnl
          AS_VAR_COPY([acx_temp1],[acx_fortran_ikind_Type])
          AS_VAR_COPY([acx_temp2],[acx_fortran_Sizeof])
          AS_IF([test -n "$acx_temp1" -a "$acx_temp2" = "$acx_temp3"],
            [acx_fortran_kind_subst=$ac_kind_lo ac_kind_lo=
             ac_kind_hi= ; break],
            [ac_kind_lo=`expr $ac_kind_lo + 1`
             AS_IF([test $ac_kind_lo -gt $ac_kind_hi],
               [ac_kind_lo= ac_kind_hi= ; break])])
          ASX_VAR_UNSET([acx_temp2])
          ASX_VAR_UNSET([acx_temp1])
          AS_VAR_POPDEF([acx_fortran_Sizeof])
          AS_VAR_POPDEF([acx_fortran_ikind_Type])
        done],
       [real],
       [ac_kind_lo=1 ac_kind_hi=32
        while :; do
          ACX_FORTRAN_RUN_CHECK_SIZEOF([$1(kind=$ac_kind_lo)])
          AS_VAR_PUSHDEF([acx_fortran_rkind_Type],
            [acx_cv_fortran_type_$1(kind=$ac_kind_lo)])dnl
          AS_VAR_PUSHDEF([acx_fortran_Sizeof],
            [acx_cv_fortran_sizeof_$1(kind=$ac_kind_lo)])dnl
          AS_VAR_COPY([acx_temp1],[acx_fortran_rkind_Type])
          AS_VAR_COPY([acx_temp2],[acx_fortran_Sizeof])
          AS_IF([test -n "$acx_temp1" -a "$acx_temp2" = "$acx_temp3"],
            [acx_fortran_kind_subst="$ac_kind_lo" ac_kind_lo=
             ac_kind_hi= ; break],
            [ac_kind_lo=`expr $ac_kind_lo + 1`
             AS_IF([test $ac_kind_lo -gt $ac_kind_hi],
               [ac_kind_lo= ac_kind_hi= ; break])])
          ASX_VAR_UNSET([acx_temp2])
          ASX_VAR_UNSET([acx_temp1])
          AS_VAR_POPDEF([acx_fortran_Sizeof])
          AS_VAR_POPDEF([acx_fortran_rkind_Type])
        done],
       [AC_MSG_WARN([Cannot derive C type correspondence for Fortran type $1(kind=$2).])])
     ASX_VAR_UNSET([acx_temp3])
     AS_VAR_POPDEF([acx_c_size])
   ],
   [acx_fortran_kind_subst=$2])
  AS_IF([test "x$acx_fortran_kind_subst" != "x$2"],
    [AC_MSG_NOTICE([Substituting $1 kind $acx_fortran_kind_subst for $3])
     m4_ifval([$2],[AC_DEFINE_UNQUOTED(AS_TR_CPP($2), $acx_fortran_kind_subst,[type kind override])])
     AC_DEFINE_UNQUOTED(AS_TR_CPP([HAVE_FORTRAN_ISO_$2]),[$acx_fortran_kind_subst],[type kind override])m4_ifval([$6],[
     $6])],[$7])
  ASX_VAR_UNSET([acx_search_kind_subst])
  ASX_VAR_UNSET([acx_fortran_kind_subst])
])
dnl
AC_DEFUN([ACX_FORTRAN_C_INT],
  [ACX_FORTRAN_TYPE_KIND([integer], [c_int], [int],[      use iso_c_binding])])
AC_DEFUN([ACX_FORTRAN_C_INT64_T],
  [AC_REQUIRE([AC_TYPE_INT64_T])
   ACX_FORTRAN_TYPE_KIND([integer], [c_int64_t], [int64_t],
  [      use iso_c_binding])])
AC_DEFUN([ACX_FORTRAN_C_FLOAT],
  [ACX_FORTRAN_TYPE_KIND([real], [c_float], [float],[      use iso_c_binding])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
