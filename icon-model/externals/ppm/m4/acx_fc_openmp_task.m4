dnl acx_fc_openmp_task.m4 --- test whether the compiler supports omp task
dnl
dnl Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: OpenMP task Fortran
dnl Author: Thomas Jahns <Thomas.Jahns@gmx.net>
dnl Maintainer: Thomas Jahns <Thomas.Jahns@gmx.net>
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
dnl
m4_define([_ACX_FORTRAN_OPENMP_TASK],
  [AS_VAR_PUSHDEF([cv],[acx_cv_[]_AC_LANG_ABBREV(_AC_LANG)_omp_task_supported])dnl
AC_CACHE_CHECK([whether $[]_AC_CC(_AC_LANG) supports omp task],
     [cv],
     [AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[      SUBROUTINE f(p)
      REAL, INTENT(OUT) :: p(:)
]m4_case(_AC_LANG,[Fortran],[!],[Fortran 77],[C])[$omp task
      p = 0.0
]m4_case(_AC_LANG,[Fortran],[!],[Fortran 77],[C])[$omp end task
      END SUBROUTINE f])],
     [AS_VAR_SET([cv],[yes])],
     [AS_VAR_SET([cv],[no])])])
   AS_VAR_IF([cv],[yes],[$1],[$2])
   AS_VAR_POPDEF([cv])])

AC_DEFUN([ACX_LANG_OPENMP_TASK(Fortran)],
  [AC_REQUIRE([AC_PROG_FC])dnl
_ACX_FORTRAN_OPENMP_TASK($@)])

AC_DEFUN([ACX_LANG_OPENMP_TASK(Fortran 77)],
  [AC_REQUIRE([AC_PROG_F77])dnl
_ACX_FORTRAN_OPENMP_TASK($@)])

AC_DEFUN([ACX_FC_OPENMP_TASK],
  [AC_LANG_PUSH([Fortran])
   _AC_LANG_DISPATCH([ACX_LANG_OPENMP_TASK],_AC_LANG,$@)
   AC_LANG_POP([Fortran])])

dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
