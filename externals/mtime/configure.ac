# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

AC_PREREQ([2.69])
AC_INIT([libmtime], [1.1.3], [luis.kornblueh@mpimet.mpg.de])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([include/mtime_calendar.h])
AC_CONFIG_HEADERS([config/config.h])

AM_INIT_AUTOMAKE([1.16.1 foreign])
AM_MAINTAINER_MODE([disable])
AM_EXTRA_RECURSIVE_TARGETS([doc install-doc man install-man])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

LT_PREREQ([2.4.6])dnl
dnl Get rid of 'rm: conftest.dSYM: is a directory' when configuring on
dnl MacOS:
m4_pushdef([_AM_PROG_CC_C_O],m4_bpatsubst(m4_dquote(m4_defn([_AM_PROG_CC_C_O])),
  [rm -f\(.* conftest\*.*\)], [rm -rf\1]))dnl
LT_INIT([pic-only])
m4_popdef([_AM_PROG_CC_C_O])dnl
dnl Clean up after Libtool checks for C compiler (e.g. Intel compiler
dnl behind an MPI compiler wrapper):
rm -f a.out

dnl Fix a bug in Autoconf 2.70 (see https://savannah.gnu.org/support/?110396):
m4_if(m4_version_compare(m4_defn([AC_AUTOCONF_VERSION]), [2.70]), [0],
  [m4_pushdef([_AC_PROG_CC_C99],
        m4_bpatsubst(m4_dquote(m4_defn([_AC_PROG_CC_C99])),
          [\[ac_c_conftest_c89_program\]], [[ac_c_conftest_c99_program]]))dnl
   AC_PROG_CC
   m4_popdef([_AC_PROG_CC_C99])],
  [AC_PROG_CC])
dnl Expand the legacy macro for older version of Autoconf (2.69 and older):
m4_if(m4_version_compare(m4_defn([AC_AUTOCONF_VERSION]), [2.70]), [-1],
  [AC_PROG_CC_C99])
AS_VAR_IF([ac_cv_prog_cc_c99], [no],
  [AC_MSG_FAILURE([C compiler does not support ISO C99])])

AC_PROG_FC
dnl Clean up after Libtool checks for Fortran compiler (e.g. Intel
dnl compiler behind an MPI compiler wrapper):
rm -f a.out
AC_LANG_PUSH([Fortran])
dnl Automake does not trace AC_FC_PP_SRCEXT/ACX_FC_PP_SRCEXT and therefore
dnl does not add their results to the compilation rules. But AC_FC_SRCEXT is
dnl traced, therefore, we temporarily redefine AC_FC_SRCEXT so it would
dnl expand to ACX_FC_PP_SRCEXT:
m4_pushdef([AC_FC_SRCEXT], [ACX_FC_PP_SRCEXT($@)])
AC_FC_SRCEXT([f90])
m4_popdef([AC_FC_SRCEXT])
ACX_FC_MODULE_NAMING([AC_SUBST([FCMODEXT], ["$acx_cv_fc_module_naming_ext"])])
AM_CONDITIONAL([FCMODUC], [test x"$acx_cv_fc_module_naming_upper" = xyes])
ACX_FC_MODULE_IN_FLAG([AC_SUBST([FCMODINC], ["$acx_cv_fc_module_in_flag"])])
AC_LANG_POP([Fortran])

AM_PATH_PYTHON([2.6])

AC_SEARCH_LIBS([roundf], [m])

AC_ARG_ENABLE([openmp],
  [AS_HELP_STRING([--enable-openmp],
     [ensure compatibility with OpenMP applications @<:@default=no@:>@])], [],
  [enable_openmp=no])

AC_SUBST([OPENMP_FCFLAG], [''])
AS_VAR_IF([enable_openmp], [yes],
  [AC_LANG_PUSH([Fortran])
   ACX_LANG_OPENMP_FLAG([OPENMP_FCFLAG=$acx_cv_fc_openmp_flag])
   AC_LANG_POP([Fortran])])
AM_CONDITIONAL([ENABLE_OPENMP], [test x"$enable_openmp" = xyes])

AC_ARG_ENABLE([examples],
  [AS_HELP_STRING([--enable-examples],
     [build examples @<:@default=yes@:>@])], [],
  [enable_examples=yes])
AM_CONDITIONAL([ENABLE_EXAMPLES], [test x"$enable_examples" = xyes])

AC_ARG_ENABLE([check],
  [AS_HELP_STRING([--enable-check],
     [enable unit testing with check library @<:@default=auto@:>@])], [],
  [enable_check=auto])

ACX_LANG_PACKAGE_INIT([CHECK], [-I], [-L])

AS_IF([test x"$enable_check" != xno],
  [mtime_save_CFLAGS=$CFLAGS
   CHECK_CFLAGS=${CHECK_CFLAGS-$acx_c_CHECK_inc_search_args}
   CFLAGS="$CHECK_CFLAGS $CFLAGS"
   AC_CHECK_HEADER([check.h],
     [ACX_LANG_LIB_SEARCH([CHECK_CLIBS], [tcase_create],
        ["$acx_c_CHECK_lib_search_args -lcheck"],
        [CHECK_CLIBS=$acx_cv_c_lib_func_tcase_create
         enable_check=yes],
        [AS_VAR_IF([enable_check], [auto],
           [AC_MSG_WARN([cannot link to CHECK library, the tests are disabled])
            enable_check=no],
           [AC_MSG_FAILURE([cannot link to CHECK library])])])],
     [AS_VAR_IF([enable_check], [auto],
        [AC_MSG_WARN([check.h not found, the tests are disabled])
         enable_check=no],
        [AC_MSG_FAILURE([check.h not found])])], [//])
   CFLAGS=$mtime_save_CFLAGS])

AS_VAR_IF([enable_check], [yes],
  [AC_DEFINE([ENABLE_CHECK], [1],
     [Define to 1 if the unit testing is enabled])],
  [CHECK_CLIBS=;CHECK_CFLAGS=])
AM_CONDITIONAL([ENABLE_CHECK], [test x"$enable_check" = xyes])

AM_MISSING_PROG([RAGEL], [ragel])

DX_DOXYGEN_FEATURE([ON])
DX_MAN_FEATURE([ON])
DX_INIT_DOXYGEN([$PACKAGE_NAME], [$(builddir)/Doxyfile])
AM_CONDITIONAL([DX_FLAG_man], [test $DX_FLAG_man -eq 1])
AM_CONDITIONAL([DX_FLAG_pdf], [test $DX_FLAG_pdf -eq 1])
AM_CONDITIONAL([DX_FLAG_html], [test $DX_FLAG_html -eq 1])
AM_CONDITIONAL([DX_FLAG_ps], [test $DX_FLAG_ps -eq 1])

AC_CONFIG_FILES([
        Makefile
        doc/Makefile
        doc/Doxyfile
        examples/Makefile
        include/Makefile
        src/Makefile
        test/Makefile
])

AC_OUTPUT
