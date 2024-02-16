dnl acx_fortran_check_include.m4 --- check module inclusion in compilation
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
dnl
# _ACX_FORTRAN_CHECK_MOD_IFELSE(MOD-FILE,
#   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#   [PREAMBLE])
# ------------------------------------------------------------------
# Check the compiler accepts MOD-FILE. The PREAMBLE might be defaulted.
m4_defun([_ACX_FORTRAN_CHECK_MOD_IFELSE],
  [AC_LANG_PUSH([Fortran])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,m4_ifval([$4],
[       $4
])[      use $1])],
     [AS_VAR_SET([acx_Mod], [yes])])
   AC_LANG_POP([Fortran])
   AS_VAR_SET_IF([acx_Mod],[$2],[$3])dnl
])# _ACX_FORTRAN_CHECK_MOD_IFELSE
# ACX_FORTRAN_CHECK_MOD_PATHS_IFELSE(MOD-FILE, PATH...
#   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
#   [MODS], [MOD-FLAGS], [TAG])
# ------------------------------------------------------------------
# Check the compiler loads MOD-FILE. The MODS might be defaulted.
# TAG defaults to extra
AC_DEFUN([ACX_FORTRAN_CHECK_MOD_PATHS_IFELSE],
  [AC_REQUIRE([AC_PROG_FC])
   AC_REQUIRE([ACX_SL_FC_CHECK_MOD_PATH_FLAG])
   AS_VAR_PUSHDEF([acx_Mod], [acx_cv_fortran_mod_$1])dnl
   AC_MSG_CHECKING([for $1 m4_default([$7],[extra]) module path])
   AC_CACHE_VAL([acx_Mod],
     [ac_mod_search_FCFLAGS_SAVE=$FCFLAGS
      AS_FOR([Ac_moddir],[ac_moddir],['' $2],
        [AS_IF([test -z "Ac_moddir"],
          [ac_res="none required"
           FCFLAGS="m4_ifval([$6],[$6 ])$ac_mod_search_FCFLAGS_SAVE"],
          [ac_res="$FC_MOD_FLAG]Ac_moddir["
           FCFLAGS="m4_ifval([$6],[$6 ])$ac_res $ac_mod_search_FCFLAGS_SAVE"])
        _ACX_FORTRAN_CHECK_MOD_IFELSE([$1],
          [AS_IF([test -z "Ac_moddir"],
             [AS_VAR_SET([acx_Mod],["$6"])],
             [AS_VAR_SET([acx_Mod],["]m4_ifval([$6],[$6 ])[$FC_MOD_FLAG]Ac_moddir["])])
           break],,[$5])])
      FCFLAGS=$ac_mod_search_FCFLAGS_SAVE])
   AS_VAR_SET_IF([acx_Mod],
     [AS_IF([test x"AS_VAR_GET([acx_Mod])" = x],
        [AC_MSG_RESULT([(none required)])],
        [AC_MSG_RESULT([AS_VAR_GET([acx_Mod])])])m4_ifval([$3],[
      $3])],
     [AC_MSG_RESULT([not found])m4_ifval([$4],[
      $4])])
   AS_VAR_POPDEF([acx_Mod])dnl
])# ACX_FORTRAN_CHECK_MOD_PATHS_IFELSE
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
