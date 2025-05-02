# ACX_FORTRAN_INCLUDE_FLAG([ACTION-IF-SUCCESS],
#                          [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Finds the compiler flag needed to specify search paths for the Fortran
# "INCLUDE" statement. The result is either "unknown" or the actual compiler
# flag, which may contain a significant trailing whitespace.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_[]_AC_LANG_ABBREV[]_ftn_include_flag variable.
#
AC_DEFUN([ACX_FORTRAN_INCLUDE_FLAG],
  [_AC_FORTRAN_ASSERT()dnl
   m4_pushdef([acx_cache_var], [acx_cv_[]_AC_LANG_ABBREV[]_ftn_include_flag])dnl
   AC_CACHE_CHECK([for _AC_LANG compiler flag needed to specify search paths dnl
for the "INCLUDE" statement], [acx_cache_var],
     [acx_cache_var=unknown
      AS_MKDIR_P([conftest.dir])
      AC_LANG_CONFTEST([AC_LANG_PROGRAM])
      mv conftest.$ac_ext conftest.dir/conftest.inc
      AC_LANG_CONFTEST([AC_LANG_SOURCE(
        [[      include "conftest.inc"]])])
      acx_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
      for acx_flag in -I '-I '; do
        _AC_LANG_PREFIX[]FLAGS="$acx_save_[]_AC_LANG_PREFIX[]FLAGS ${acx_flag}conftest.dir"
        AC_LINK_IFELSE([], [acx_cache_var=$acx_flag])
        test "x$acx_cache_var" != xunknown && break
      done
      []_AC_LANG_PREFIX[]FLAGS=$acx_save_[]_AC_LANG_PREFIX[]FLAGS
      rm -rf conftest.$ac_ext conftest.dir])
   AS_VAR_IF([acx_cache_var], [unknown], [m4_default([$2],
     [AC_MSG_FAILURE([unable to detect _AC_LANG compiler flag needed to dnl
specify search paths for the "INCLUDE" statement])])], [$1])
   m4_popdef([acx_cache_var])])
