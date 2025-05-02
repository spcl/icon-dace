# Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#

AC_PREREQ([2.69])
LT_PREREQ([2.4.6])dnl

m4_pushdef([mtime_version], [1, 2, 2])dnl
AC_INIT([mtime],
  m4_join([],
    m4_argn(1, mtime_version).m4_argn(2, mtime_version).m4_argn(3, mtime_version),
    m4_argn([4], mtime_version)),
  [luis.kornblueh@mpimet.mpg.de])
AC_SUBST([MTIME_VERSION_MAJOR], ['m4_argn([1], mtime_version)'])dnl
AC_SUBST([MTIME_VERSION_MINOR], ['m4_argn([2], mtime_version)'])dnl
AC_SUBST([MTIME_VERSION_PATCH], ['m4_argn([3], mtime_version)'])dnl
m4_popdef([mtime_version])dnl

AC_CONFIG_MACRO_DIR([m4])
dnl Make sure all required m4 files are in place:
m4_pattern_forbid([^((A[SC]|D)X|_KPSE)_])dnl
m4_pattern_allow([^DX_(CLEAN(_|FILES)|(PDF|PS)_GOAL|SNIPPET_|V_)])dnl
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([include/mtime_calendar.h])
AC_CONFIG_HEADERS([config/config.h])

AM_INIT_AUTOMAKE([1.16.1 foreign])
AM_MAINTAINER_MODE([enable])
AM_EXTRA_RECURSIVE_TARGETS([doc install-doc man install-man])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

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

ACX_USE_LIBTOOL_CONFIGURATION([pic-only])
dnl Get extension of a dynamically loaded module file from Libtool:
{ module=yes; eval acx_modext=\"$shrext_cmds\"; unset module;}
AC_SUBST([acx_modext])dnl

AC_LANG_PUSH([Fortran])
dnl Automake does not trace AC_FC_PP_SRCEXT/ACX_FC_PP_SRCEXT and therefore
dnl does not add their results to the compilation rules. But AC_FC_SRCEXT is
dnl traced, therefore, we temporarily redefine AC_FC_SRCEXT so it would
dnl expand to ACX_FC_PP_SRCEXT:
m4_pushdef([AC_FC_SRCEXT], [ACX_FC_PP_SRCEXT($@)])
AC_FC_SRCEXT([f90])
m4_popdef([AC_FC_SRCEXT])
dnl Handle possibly long expansions of the preprocessor macro __FILE__:
ACX_FC_LINE_LENGTH([unlimited])
ACX_FC_MODULE_NAMING([AC_SUBST([FCMODEXT], ["$acx_cv_fc_module_naming_ext"])])
AM_CONDITIONAL([FCMODUC], [test x"$acx_cv_fc_module_naming_upper" = xyes])
ACX_FC_MODULE_IN_FLAG([AC_SUBST([FCMODINC], ["$acx_cv_fc_module_in_flag"])])
AC_LANG_POP([Fortran])

AC_SEARCH_LIBS([roundf], [m])

AC_ARG_ENABLE([examples],
  [AS_HELP_STRING([--enable-examples],
     [build examples @<:@default=yes@:>@])], [],
  [enable_examples=yes])
AM_CONDITIONAL([ENABLE_EXAMPLES], [test x"$enable_examples" = xyes])

AC_ARG_ENABLE([fortran-hl],
  [AS_HELP_STRING([--enable-fortran-hl],
     [enable the high-level Fortran interface @<:@default=yes@:>@])], [],
  [enable_fortran_hl=yes])
AM_CONDITIONAL([ENABLE_FORTRAN_HL], [test x"$enable_fortran_hl" = xyes])

AC_ARG_ENABLE([python],
  [AS_HELP_STRING([--enable-python],
     [enable the Python interface @<:@default=no@:>@])], [],
  [enable_python=no])
AS_VAR_IF([enable_python], [yes],
  [AS_IF([test "x$enable_shared" = xno || test "x$with_pic" = xno],
     [AC_MSG_ERROR([cannot provide the Python interface when the shared libraries or dnl
PIC objects are disabled])])])

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
           [enable_check=no],
           [AC_MSG_FAILURE([cannot link to CHECK library])])])],
     [AS_VAR_IF([enable_check], [auto],
        [enable_check=no],
        [AC_MSG_FAILURE([check.h not found])])],
     [//])
   CFLAGS=$mtime_save_CFLAGS])

AS_VAR_IF([enable_check], [yes],
  [AC_DEFINE([ENABLE_CHECK], [1],
     [Define to 1 if the unit testing is enabled])],
  [CHECK_CLIBS=;CHECK_CFLAGS=])
AM_CONDITIONAL([ENABLE_CHECK], [test x"$enable_check" = xyes])

AS_VAR_IF([enable_python], [yes], [AM_PATH_PYTHON([2.7])])
AM_CONDITIONAL([ENABLE_PYTHON], [test x"$enable_python" = xyes])

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
  doc/Doxyfile
  doc/Makefile
  examples/Makefile
  include/Makefile
  python/Makefile
  src/Makefile
  test/Makefile
])

AC_CONFIG_FILES([
  test/test_python.sh
], [chmod a+x "$ac_file"])

AS_VAR_IF([enable_examples], [yes],
  [AC_CONFIG_LINKS([
     examples/iconatm.nml:examples/iconatm.nml
     examples/iconoce.nml:examples/iconoce.nml
])])

AS_VAR_IF([enable_python], [yes],
  [AC_CONFIG_LINKS([
     python/mtime/__init__.py:python/mtime/__init__.py
     python/mtime/_mtime.py:python/mtime/_mtime.py
])])

AC_OUTPUT
