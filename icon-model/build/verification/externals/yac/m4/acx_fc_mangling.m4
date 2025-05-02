# Copyright (c) 2018-2024, MPI-M
#
# Author: Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# ACX_FC_MANGLING_GLOBAL()
# -----------------------------------------------------------------------------
# Detects the name mangling scheme for the global Fortran functions. The result
# is either "unknown" or a comma-separated pair of two strings each denoting
# the name mangling scheme of the global Fortran functions without an
# underscore ('_') in their original names and with underscores in their
# original names, respectively. Both strings start with either "name" or "NAME"
# depending on whether the Fortran compiler generates symbols that correspond
# to function names in the lowercase or in the uppercase. The suffix of each
# string is zero or more underscores that need to be added to the original
# names of the functions int the given letter case.
#
# The implementation implies that the Fortran compiler can link objects
# compiled with the C compiler.
#
# The result is cached in the acx_cv_fc_mangling_global variable.
#
AC_DEFUN([ACX_FC_MANGLING_GLOBAL],
  [AC_REQUIRE([AC_PROG_CC])dnl
   m4_pushdef([acx_cache_var], [acx_cv_fc_mangling_global])dnl
   AC_LANG_PUSH([Fortran])
   AC_CACHE_CHECK([for the name-mangling scheme for Fortran global functions],
     [acx_cache_var],
     [acx_cache_var=unknown
      acx_save_LIBS=$LIBS; LIBS="./conftest_c.$ac_objext $LIBS"
      AC_LANG_CONFTEST([AC_LANG_SOURCE(
[[      subroutine funcname
      end subroutine
      program main
      end program]])])
      AC_LANG_PUSH([C])
      acx_success=no
      for acx_name_case in 'funcname' 'FUNCNAME'; do
        for acx_underscore in '' '_'; do
          acx_func_name="$acx_name_case$acx_underscore"
          AC_COMPILE_IFELSE([AC_LANG_CALL(
            [[#define main conftest_c]], [$acx_func_name])],
            [mv ./conftest.$ac_objext ./conftest_c.$ac_objext
             AC_LANG_POP([C])
             AC_LINK_IFELSE([], [acx_success=yes])
             AC_LANG_PUSH([C])
             rm -f ./conftest_c.$ac_objext
             test "x$acx_success" = xyes && break])
        done
        test "x$acx_success" = xyes && break
      done
      AS_VAR_IF([acx_success], [yes],
        [AS_VAR_IF([acx_name_case], [funcname],
           [acx_name_case='func_name'
            acx_cache_var="name$acx_underscore,name$acx_underscore"],
           [acx_name_case='FUNC_NAME'
            acx_cache_var="NAME$acx_underscore,NAME$acx_underscore"])
         AC_LANG_POP([C])
         AC_LANG_CONFTEST([AC_LANG_SOURCE(
[[      subroutine func_name
      end subroutine
      program main
      end program]])])
         AC_LANG_PUSH([C])
         acx_success=no
         for acx_extra_underscore in '' '_'; do
           acx_func_name="$acx_name_case$acx_underscore$acx_extra_underscore"
           AC_COMPILE_IFELSE([AC_LANG_CALL(
            [[#define main conftest_c]], [$acx_func_name])],
            [mv ./conftest.$ac_objext ./conftest_c.$ac_objext
             AC_LANG_POP([C])
             AC_LINK_IFELSE([], [acx_success=yes])
             AC_LANG_PUSH([C])
             rm -f ./conftest_c.$ac_objext
             test "x$acx_success" = xyes && break])
         done
         AS_VAR_IF([acx_success], [yes],
           [AS_VAR_APPEND([acx_cache_var], ["$acx_extra_underscore"])],
           [acx_cache_var=unknown])])
      AC_LANG_POP([C])
      rm -f conftest.$ac_ext
      LIBS=$acx_save_LIBS])
   AC_LANG_POP([Fortran])
   m4_popdef([acx_cache_var])])

# ACX_FC_MANGLING_DEFINE([MACRO-PREFIX = FC]
#                        [ACTION-IF-SUCCESS],
#                        [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Defines C preprocessor macros <MACRO-PREFIX>_GLOBAL(name, NAME) and
# <MACRO-PREFIX>_GLOBAL_(name, NAME) to properly mangle the names of C/C++
# identifiers, and identifiers with underscores, respectively, so that they
# match the name-mangling scheme for the global functions used by the Fortran
# compiler. The default value for <MACRO-PREFIX> is "FC".
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
AC_DEFUN([ACX_FC_MANGLING_DEFINE],
  [AC_REQUIRE([ACX_FC_MANGLING_GLOBAL])dnl
   m4_pushdef([acx_scheme_var], [acx_cv_fc_mangling_global])dnl
   AS_VAR_IF([acx_scheme_var], [unknown], [m4_default([$3],
        [AC_MSG_FAILURE([unable to detect the name-mangling scheme for dnl
Fortran global functions])])],
     [AC_DEFINE_UNQUOTED(m4_default([$1], [FC])[_GLOBAL(name,NAME)],
        [`AS_ECHO(["$acx_scheme_var"]) | cut -d, -f1 | sed 's%__*% [##] &%'`],
        [Define to a macro mangling the given C identifier ]dnl
[(in lower and upper case), which must not contain underscores, for ]dnl
[linking Fortran global functions.])
      AC_DEFINE_UNQUOTED(m4_default([$1], [FC])[_GLOBAL_(name,NAME)],
        [`AS_ECHO(["$acx_scheme_var"]) | cut -d, -f2 | sed 's%__*% [##] &%'`],
        [As ]m4_default([$1], [FC])[_GLOBAL, but for identifiers ]dnl
[containing underscores.])
      $2])
   m4_popdef([acx_scheme_var])])

# ACX_FC_MANGLING_SHVAR(FUNC-NAME,
#                       [ACTION-IF-SUCCESS],
#                       [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Applies the name mangling scheme for the global Fortran functions to the
# function FUNC-NAME. The result is empty on error.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is stored in the acx_mangled_name variable and cached in the
# acx_cv_fc_mangling_shvar_[]AS_TR_CPP(FUNC-NAME) variable.
#
AC_DEFUN([ACX_FC_MANGLING_SHVAR],
  [AC_REQUIRE([ACX_FC_MANGLING_GLOBAL])dnl
   m4_pushdef([acx_cache_var], [acx_cv_fc_mangling_shvar_[]AS_TR_CPP([$1])])dnl
   AC_CACHE_CHECK([for the mangled name of the Fortran function $1],
     [acx_cache_var],
     [AS_VAR_SET([acx_cache_var])
      m4_pushdef([acx_scheme_var], [acx_cv_fc_mangling_global])dnl
      AS_VAR_IF([acx_scheme_var], [unknown], [],
        [acx_func_name=$1
         AS_CASE(["$acx_func_name"],
           [*_*],
           [acx_tmp=`AS_ECHO(["$acx_scheme_var"]) | cut -d, -f2`],
           [acx_tmp=`AS_ECHO(["$acx_scheme_var"]) | cut -d, -f1`])
         AS_CASE(["$acx_tmp"],
           [NAME*],
           [acx_tmp=`AS_ECHO(["$acx_tmp"]) | dnl
sed "s/^NAME/$acx_func_name/" | tr 'm4_cr_letters' 'm4_cr_LETTERS'`],
           [acx_tmp=`AS_ECHO(["$acx_tmp"]) | dnl
sed "s/^name/$acx_func_name/" | tr 'm4_cr_LETTERS' 'm4_cr_letters'`])
         AS_VAR_COPY([acx_cache_var], [acx_tmp])])
       m4_popdef([acx_scheme_var])])
   AS_VAR_COPY([acx_mangled_name], [acx_cache_var])
   m4_popdef([acx_cache_var])dnl
   AS_IF([test -n "$acx_mangled_name"], [$2], [m4_default([$3],
     [AC_MSG_FAILURE([unable to detect the mangled name of dnl
the Fortran function $1])])])])

# ACX_FC_MANGLING_MAIN([ACTION-IF-SUCCESS],
#                      [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Detects the symbol name of the Fortran entry point function (i.e. the real
# name of the main function expected by the Fortran compiler when linking an
# executable) and defines the following C preprocessor macros: FC_MAIN equaling
# the result (if the current language is Fortran), F77_MAIN equaling the result
# (if the current language Fortran 77), and FC_MAIN_EQ_F77 equaling to 1 (if
# the macro is called for both Fortran and Fortran 77 and the results equal to
# each other). The result is either "unknown" or the actual name of the Fortran
# entry point function. The mentioned C preprocessor macros are defined only if
# the result is successful.
#
# The implementation implies that the Fortran compiler can link objects
# compiled with the C compiler. For example, Intel C compiler injects the
# reference to the "__intel_new_feature_proc_init" symbol into an object file
# containing definition of the function "main", which is unknown Fortran
# compilers other than Intel Fortran compiler. Thus, if the name of the Fortran
# entry point function equals to "main" (e.g. GNU Fortran compiler) and the
# C compiler of choice is Intel C compiler, the result would be "unknown" if
# a workaround for this particular case had not been implemented.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_mangling_main variable.
#
AC_DEFUN([ACX_FC_MANGLING_MAIN],
  [_AC_FORTRAN_ASSERT()dnl
   AC_REQUIRE([AC_PROG_CC])dnl
   m4_pushdef([acx_cache_var], [acx_cv_fc_mangling_main])dnl
   AC_CACHE_CHECK([for the name of _AC_LANG entry point function],
     [acx_cache_var],
     [acx_cache_var=unknown
      acx_save_LIBS=$LIBS; LIBS="./conftest_c.$ac_objext $LIBS"
      AC_LANG_CONFTEST([AC_LANG_SOURCE(
[[      subroutine conftest
      end subroutine]])])
      AC_LANG_PUSH([C])
      for acx_func_name in main MAIN__ MAIN_ __main MAIN _MAIN __MAIN main_ \
main__ _main; do
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#ifdef __cplusplus
extern "C"
#endif
#ifdef __INTEL_COMPILER
__intel_new_feature_proc_init() {}
#endif
#define main $acx_func_name]])],
          [mv ./conftest.$ac_objext ./conftest_c.$ac_objext
           AC_LANG_POP([C])
           AC_LINK_IFELSE([], [acx_cache_var=$acx_func_name])
           AC_LANG_PUSH([C])
           rm -f conftest_c.$ac_objext])
        test "x$acx_cache_var" != xunknown && break
      done
      AC_LANG_POP([C])
      rm -f conftest.$ac_ext
      LIBS=$acx_save_LIBS])
   AS_VAR_IF([acx_cache_var], [unknown],
     [m4_default([$2], [AC_MSG_FAILURE([unable to detect the name of dnl
Fortran entry point function])])],
     [AC_DEFINE_UNQUOTED([FC_MAIN], [$acx_cache_var],
        [Define to the name of Fortran entry point function.])
      $1])
   m4_popdef([acx_cache_var])])
