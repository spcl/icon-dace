dnl configure.ac --- autoconf script for YAXT
dnl
dnl Copyright  (C)  2019  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: configure configure.ac autotools
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
dnl Commentary:
dnl
dnl
dnl
dnl Code:
dnl
dnl ######################################################################
dnl Basic autotools settings.
dnl ######################################################################
AC_PREREQ([2.64])
dnl Update the corresponding entry in doc/Doxyfile when changing this
AC_INIT([yaxt],[0.10.2],[hanke@dkrz.de],,
  [http://www.dkrz.de/redmine/projects/yaxt])
AC_CONFIG_SRCDIR([src/xt_idxlist.c])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
dnl forbid acx macro names from remaining unexpanded
m4_pattern_forbid([^ACX_])
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects parallel-tests tar-ustar])
AM_MAINTAINER_MODE
dnl ######################################################################
dnl Diagnose incompatible tool versions
dnl ######################################################################
m4_ifdef([LT_PACKAGE_VERSION],
  [m4_if(m4_version_compare(AC_AUTOCONF_VERSION,[2.67]),1,
    [m4_if(m4_version_compare(LT_PACKAGE_VERSION,[2.4.2]),[-1],
      [m4_fatal([autoconf versions 2.68 and newer require using libtool 2.4.2 or newer])])])],
  [m4_warn([syntax],[Please run 'libtoolize --install' before continuing.])])
dnl ######################################################################
dnl Checks for programs.
dnl ######################################################################
CHECK_GNU_MAKE
dnl we prefer GNU make because the System V make on many commercial
dnl platforms is slow
AC_CHECK_PROGS([MAKE], [$_cv_gnu_make_command gmake gnumake make])
dnl maintainers need GNU make for the regeneration rules, so enforce that
dnl
AS_IF([test "x${_cv_gnu_make_command}x$USE_MAINTAINER_MODE" = xxyes],
  [AC_MSG_FAILURE([maintainers need to use GNU make or compatible])])
dnl Only substitute MAKE in Makefiles if actually needed
AC_PROG_MAKE_SET
AM_SUBST_NOTMAKE([MAKE])
dnl Require an ISO C99 compliant C compiler
m4_pushdef(m4_case(m4_version_compare(AC_AUTOCONF_VERSION,[2.70]),
     [-1],[[_AC_PROG_CC_C89]],
     [0],[[_AC_PROG_CC_C11]],
     [1],[[_AC_PROG_CC_STDC_EDITION]]),
  [rm -rf conftest.dSYM # needed when configured on MacOS with CFLAGS='-g'
   m4_case(m4_cmp(m4_version_compare(AC_AUTOCONF_VERSION,[2.71]),[-1]),
     [0],[_AC_PROG_CC_C99(,[AC_MSG_FAILURE([CC=$CC is not an ISO C99 compiler])])],
     [1],[ac_prog_cc_stdc=no
_AC_PROG_CC_STDC_EDITION_TRY([99])
AS_IF([test "x$ac_prog_cc_stdc" != xc99],
  [AC_MSG_FAILURE([CC=$CC is not an ISO C99 compiler])])])])dnl
AC_PROG_CC([mpicc mpiicc mpxlc_r mpixlc_r gcc cc])
m4_popdef(m4_case(m4_version_compare(AC_AUTOCONF_VERSION,[2.70]),
  [-1],[[_AC_PROG_CC_C89]],
  [0],[[_AC_PROG_CC_C11]],
  [1],[[_AC_PROG_CC_STDC_EDITION]]))dnl

AX_PTHREAD([AC_DEFINE([HAVE_PTHREAD],1,[Define if you have POSIX threads libraries and header files.])
   LIBS="$PTHREAD_LIBS $LIBS"
   test -n "$PTHREAD_CFLAGS" || ASX_VAR_UNSET([PTHREAD_CFLAGS])
   CFLAGS="$CFLAGS${PTHREAD_CFLAGS+ $PTHREAD_CFLAGS}"
   CC="$PTHREAD_CC"])
ACX_TLS_XLC_RETRY
if test "X$FC" != Xno ; then
m4_pushdef([_AC_F95_FC],[mpifort mpiifort mpxlf_r mpxlf95_r mpxlf_r mpixlf95_r mpif90 mpif77 ]m4_dquote(m4_defn(
  [_AC_F95_FC])))dnl
AC_PROG_FC
m4_popdef([_AC_F95_FC])dnl
fi
m4_define([test_with_fortran],[test -n "$FC" && test "X$FC" != Xno])dnl
if test_with_fortran ; then
AC_FC_SRCEXT([f90])
AC_FC_FREEFORM
AC_LANG_PUSH([Fortran])
ACX_PROG_FC_CHECK_FPP
dnl Appending the preprocessor flag to FCFLAGS directly is wrong because
dnl it "must appear immediately before the source-code file name when
dnl compiling" (see the description of AC_FC_PP_SRCEXT in Autoconf manual).
dnl Moreover, the flag should not appear in the linking command because it
dnl might confuse the compiler, which is, for example, the case for Cray
dnl compiler when the preprocessor flag is equal to '-eZ'. Indeed, this looks
dnl more like a bug in the compiler but there might be more such cases.
dnl Meanwhile, the problem with Cray compiler is unexpectedly resolved by using
dnl '-eT' flag instead of '-eZ' (the original intention for switching the flag
dnl was to get rid of intermediate *.i files generated by the compiler).
FCFLAGS_f90="$FCFLAGS_f90 $FC_FPP_FLAG"
AC_PROG_FPP
AC_PROG_FC_C_O
AC_LANG_POP([Fortran])
fi
AM_CONDITIONAL([with_fortran], [test_with_fortran])

AC_PROG_MKDIR_P
AX_WITH_PERL([/bin/false], $PATH:$prefix/bin:/usr/bin:/usr/local/bin)
AX_PROG_PERL_VERSION([5.0005],[],
        AC_MSG_ERROR([Perl (version>=5.0005) is required]))
AC_CHECK_PROGS([M4],[m4 gm4 gnum4],[false],,)
AM_PROG_CC_C_O
# Set up libtool.
ACX_USE_LIBTOOL_CONFIGURATION
AC_SUBST([lt_cv_nm_interface])
acx_symprfx=`expr "$lt_cv_sys_global_symbol_pipe" : '.* prfx=\^\(@<:@_@:>@\{0,1\}\)@S|@'` \
  || acx_symprfx=`expr "$lt_cv_sys_global_symbol_pipe" : '.*\\\\1 \(@<:@_@:>@\{0,1\}\)\\\\2'`
AC_SUBST([acx_symprfx])
dnl ######################################################################
dnl Properties of the compiler
dnl ######################################################################
if test_with_fortran ; then
ACX_SL_FC_MOD_SUFFIX(,
  [AC_MSG_ERROR([Fortran module file extension must be known to install F90 modules.])])
AC_LANG_PUSH([Fortran])
ACX_SL_FC_CHECK_MOD_PATH_FLAG
AC_LANG_POP([Fortran])
ACX_FC_VOLATILE([XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define HAVE_VOLATILE"])
ACX_FC_HAVE_SIZE_KIND_ARGUMENT([XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define HAVE_SIZE_KIND_ARGUMENT"])
AC_FC_DUMMY_MAIN(,[:])
AC_FC_WRAPPERS
ACX_FC_LOGICAL_INTEROP([HAVE_FC_LOGICAL_INTEROP=yes ;dnl
 XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define HAVE_FC_LOGICAL_INTEROP"])
ACX_FC_LOGICAL_C_BOOL([XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define HAVE_FC_LOGICAL_C_BOOL"])
dnl
dnl check if the compiler supports pointer bounds remapping
ACX_FC_POINTER_BOUNDS_REMAP([HAVE_POINTER_REMAP=yes],[HAVE_POINTER_REMAP=no])
AS_VAR_IF([HAVE_POINTER_REMAP],[yes],
  [ACX_FC_POINTER_HUGE_BOUNDS_REMAP(,
  [AC_MSG_NOTICE([Disabling native low-bounds mapping, $FC has a known defective implementation.])
   HAVE_POINTER_REMAP=no])])
AS_IF([test $HAVE_POINTER_REMAP = yes],
  [XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define HAVE_FC_PTR_BOUND_REMAP"])
dnl Do we need helper functions for contiguous arrays?
ACX_FC_HAVE_IS_CONTIGUOUS([HAVE_FC_IS_CONTIGUOUS=yes
  XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define HAVE_FC_IS_CONTIGUOUS"
AC_DEFINE([HAVE_FC_IS_CONTIGUOUS],[1])])
ACX_FC_HAVE_C_LOC_OF_SLICE([HAVE_FC_C_LOC_OF_SLICE=yes
  XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define HAVE_FC_C_LOC_OF_SLICE"
  AC_DEFINE([HAVE_FC_C_LOC_OF_SLICE],[1])])
ACX_FC_INTEGER_EXP_IS_PRECISE(
  [XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define HAVE_FC_PRECISE_INTEGER_EXPONENTIATION"])
fi
AM_CONDITIONAL([FC_MOD_UPPERCASE],
  [test x$ac_fc_mod_uppercase = xyes])
AH_TEMPLATE([HAVE_FC_IS_CONTIGUOUS],
  [Defined to 1 iff Fortran compiler supports IS_CONTIGUOUS query])dnl
AH_TEMPLATE([HAVE_FC_C_LOC_OF_SLICE],
  [Defined to 1 iff Fortran compiler supports C_LOC of an array slice])dnl
dnl Check whether char is unsigned by default
AC_C_CHAR_UNSIGNED
AS_IF([test x$ac_cv_c_char_unsigned = xyes],
  [C_CHAR_IS_UNSIGNED=1],
  [C_CHAR_IS_UNSIGNED=0])
AC_SUBST([C_CHAR_IS_UNSIGNED])
dnl preprocessor implicit type conversions depends on the size of these, also
dnl the result is used in establishing Xt_long
AC_CHECK_SIZEOF([int])
AC_SUBST([C_SIZEOF_INT],[$ac_cv_sizeof_int])
AC_CHECK_SIZEOF([long])
AC_SUBST([C_SIZEOF_LONG],[$ac_cv_sizeof_long])
AC_CHECK_SIZEOF([long long])
AC_SUBST([C_SIZEOF_LONG_LONG],[$ac_cv_sizeof_long_long])
AC_CHECK_SIZEOF([unsigned long])
AC_CHECK_SIZEOF([__int128])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([__int128_t f(void) {
  __int128_t i = 0, j = i + 1;
  return j >= 0;
}])],
  [AC_CHECK_SIZEOF([__int128_t])],
  [ac_cv_sizeof___int128_t=0])
dnl test what type to use for int64_t (usually long or long long)
TJ_FIND_INTEGRAL_TYPE([int64_t],[DEFINED_INT64_T],[@%:@include <inttypes.h>])
AS_IF([test x"$DEFINED_INT64_T" = x"long long"],
  [AC_DEFINE([INT64_T_IS_LONG_LONG],[1],
     [defined if int64_t is long long in disguise])],
  [test x"$DEFINED_INT64_T" = x"long"],
  [AC_DEFINE([INT64_T_IS_LONG],[1],
     [defined if int64_t is long in disguise])])
dnl
dnl
dnl check for various builtins
AS_FOR([builtin_macro],[builtin],[__builtin_ctzl __builtin_clzl],
  [AS_VAR_PUSHDEF([builtin_cache],[acx_cv_have_decl_]builtin_macro)
   AC_CACHE_CHECK([whether ]builtin_macro[ is declared],
     [builtin_cache],
     [AC_LINK_IFELSE([AC_LANG_PROGRAM(,[  long lbz = builtin_macro][[(56UL)]])],
        [AS_VAR_SET([builtin_cache],[yes])],
        [AS_VAR_SET([builtin_cache],[no])])])
   AS_VAR_IF([builtin_cache],[yes],
     [AC_DEFINE_UNQUOTED([HAVE_DECL_]AS_TR_CPP([builtin_macro]),[1])],
     [AC_DEFINE_UNQUOTED([HAVE_DECL_]AS_TR_CPP([builtin_macro]),[0])])
   AS_VAR_POPDEF([builtin_cache])])
AH_TEMPLATE([HAVE_DECL___BUILTIN_CTZL],
  [Define to 1 if __builtin_ctzl is available, 0 if not])
AH_TEMPLATE([HAVE_DECL___BUILTIN_CLZL],
  [Define to 1 if __builtin_clzl is available, 0 if not])
AC_CACHE_CHECK([whether __builtin_clz is declared],
  [acx_cv_have_decl___builtin_clz],
  [AC_LINK_IFELSE([AC_LANG_PROGRAM(,[  int nlz = __builtin_clz(56U)])],
     [AS_VAR_SET([acx_cv_have_decl___builtin_clz],[yes])],
     [AS_VAR_SET([acx_cv_have_decl___builtin_clz],[no])])])
AC_CACHE_CHECK([whether __builtin_clzll is declared],
  [acx_cv_have_decl___builtin_clzll],
  [AC_LINK_IFELSE([AC_LANG_PROGRAM(,[  int nlz = __builtin_clzll(56ULL)])],
     [AS_VAR_SET([acx_cv_have_decl___builtin_clzll],[yes])],
     [AS_VAR_SET([acx_cv_have_decl___builtin_clzll],[no])])])
dnl
AC_CACHE_CHECK([whether __lzcnt16 is available],
  [acx_cv_have_decl___lzcnt16],
  [AC_LINK_IFELSE([AC_LANG_PROGRAM([@%:@include <intrin.h>],
       [  unsigned short nlz = __lzcnt16((unsigned short)56)])],
     [AS_VAR_SET([acx_cv_have_decl___lzcnt16],[yes])],
     [AS_VAR_SET([acx_cv_have_decl___lzcnt16],[no])])])
AC_DEFINE_UNQUOTED([HAVE_DECL___LZCNT16],
  [`test "$acx_cv_have_decl___lzcnt16" = yes ; expr 1 - $?`],
  [PGI/Nvidia style __lzcnt16 is available])
dnl
AC_CACHE_CHECK([whether __lzcnt is available],
  [acx_cv_have_decl___lzcnt],
  [AC_LINK_IFELSE([AC_LANG_PROGRAM([@%:@include <intrin.h>],
       [  unsigned nlz = __lzcnt(56U)])],
     [AS_VAR_SET([acx_cv_have_decl___lzcnt],[yes])],
     [AS_VAR_SET([acx_cv_have_decl___lzcnt],[no])])])
AC_DEFINE_UNQUOTED([HAVE_DECL___LZCNT],
  [`test "$acx_cv_have_decl___lzcnt" = yes ; expr 1 - $?`],
  [PGI/Nvidia style __lzcnt is available])
dnl
AC_CACHE_CHECK([whether __lzcnt64 is available],
  [acx_cv_have_decl___lzcnt64],
  [AC_LINK_IFELSE([AC_LANG_PROGRAM([@%:@include <intrin.h>],
       [  unsigned __int64 nlz = __lzcnt64(56ULL)])],
     [AS_VAR_SET([acx_cv_have_decl___lzcnt64],[yes])],
     [AS_VAR_SET([acx_cv_have_decl___lzcnt64],[no])])])
AC_DEFINE_UNQUOTED([HAVE_DECL___LZCNT64],
  [`test "$acx_cv_have_decl___lzcnt64" = yes ; expr 1 - $?`],
  [PGI/Nvidia style __lzcnt64 is available])
dnl
AC_CACHE_CHECK([whether __builtin_complex is available],
  [acx_cv_have_decl___builtin_complex],
  [AC_LINK_IFELSE([AC_LANG_PROGRAM([@%:@include <complex.h>],[   long double _Complex c =  __builtin_complex(1.0L, -2.5L)])],
     [AS_VAR_SET([acx_cv_have_decl___builtin_complex],[yes])],
     [AS_VAR_SET([acx_cv_have_decl___builtin_complex],[no])])])
dnl ######################################################################
dnl Test whether cfortran.h works correctly
dnl ######################################################################
if test_with_fortran ; then
ACX_FIND_CFORTRAN_DEF
ACX_XLF_QEXTNAME_ADD_APPENDUS
ACX_CHECK_CFORTRAN
dnl we want to copy strings to the stack for cf calls
AC_FUNC_ALLOCA
fi

dnl ######################################################################
dnl Checks for header files.
dnl ######################################################################
dnl
while : ; do[]m4_pushdef([cuda_includes],[AC_INCLUDES_DEFAULT
@%:@include <cuda.h>])m4_pushdef([no_cuda],[AC_MSG_NOTICE([$1, not building CUDA component])
acx_cuda_version=0 ; break])
AC_CHECK_HEADER([cuda.h],
  [AC_CHECK_DECL([CUDA_VERSION],
     [AC_COMPUTE_INT([acx_cuda_version],
        [CUDA_VERSION],[cuda_includes],[no_cuda([Cannot determine CUDA version])])
      AC_CHECK_HEADERS([dlfcn.h],,
        [no_cuda([cuda.h found but dlsym is unavailable])])],
     [no_cuda([CUDA_VERSION definition not found])],[cuda_includes])],
  [no_cuda([cuda header not found])])
  break
done[]m4_popdef([no_cuda])m4_popdef([cuda_includes])
AS_IF([test "$acx_cuda_version" -ge 6050],
  [AC_DEFINE([HAVE_CUDA],[1],
     [defined to 1 iff a cuda header of the 6.5 or newer version was found])])
AM_CONDITIONAL([HAVE_CUDA],[test "$acx_cuda_version" -ge 6050])
#  ----------------------------------------------------------------------
# Compile with MPI support
dnl
dnl ######################################################################
dnl Checks for MPI.
dnl ######################################################################
ACX_C_PACKAGE([MPI],[mpi.h],,,
  [AC_MSG_ERROR([Required header mpi.h not found or not compilable.])],
  [MPI_Waitall],[mpi mpich],,,
  [AC_MSG_ERROR([Cannot link C MPI programs.])])
if test_with_fortran ; then
ACX_F90_PACKAGE([MPI],[mpi],,,
  [AC_MSG_ERROR([Required Fortran 90 module not found or incompatible.])],
  [mpi_waitall],[mpi mpi_f90 mpi_f77 mpich],[[-lmpi_f77 -lmpi],[-lmpi]],,
  [AC_MSG_ERROR([Cannot link Fortran 90 MPI programs.])],,
  [      INTEGER :: req(1), stat(mpi_status_size, 1), ierror
],[(1, req, stat, ierror)])
fi

saved_CFLAGS=$CFLAGS
CFLAGS="$CFLAGS $MPI_C_INCLUDE"
TJ_FIND_INTEGRAL_TYPE([MPI_Fint],[MPI_FINT],[@%:@include <mpi.h>])
ACX_MPI_SEND_CONST_VOID_P_BUF_ARG([CAST_MPI_SEND_BUF="(buf)"
   XT_MPI_SEND_BUF_CONST=const],
  [CAST_MPI_SEND_BUF="((void *)(buf))"
   XT_MPI_SEND_BUF_CONST=])
AC_DEFINE_UNQUOTED([CAST_MPI_SEND_BUF(buf)],[$CAST_MPI_SEND_BUF],
  [Macro to cast arguments to MPI_Send-like functions (if necessary)])
AC_DEFINE_UNQUOTED([XT_MPI_SEND_BUF_CONST],[$XT_MPI_SEND_BUF_CONST],
    [Defined to const, iff buf argument to MPI_Send-like functions can be const])
TJ_FIND_INTEGRAL_TYPE([MPI_Aint],[xt_mpi_aint],[@%:@include <mpi.h>])
xt_mpi_aint_kind=c_[]AS_TR_SH([$xt_mpi_aint])
XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define XT_MPI_AINT_KIND $xt_mpi_aint_kind"
AC_CHECK_SIZEOF([MPI_Datatype],,[@%:@include <mpi.h>])
AS_IF([test $ac_cv_sizeof_MPI_Datatype -eq 0],
  [AC_MSG_FAILURE([Could not detect size of MPI_Datatype!])])
AC_COMPUTE_INT([MPI_VERSION], [MPI_VERSION], [@%:@include <mpi.h>],
  [AC_MSG_ERROR([cannot determine MPI version!])])
CFLAGS=$saved_CFLAGS
AS_CASE(["$MPI_FINT"],
  [int],[XT_MPI_FINT_MAX=INT_MAX XT_MPI_FINT_FC_KIND=c_int],
  [long],[XT_MPI_FINT_MAX=LONG_MAX XT_MPI_FINT_FC_KIND=c_long],
  ["long long"],[XT_MPI_FINT_MAX=LLONG_MAX XT_MPI_FINT_FC_KIND=c_long_long],
  [short],[XT_MPI_FINT_MAX=SHRT_MAX XT_MPI_FINT_FC_KIND=c_short])
AC_DEFINE_UNQUOTED([XT_MPI_FINT_MAX],[$XT_MPI_FINT_MAX],
  [Maximum value of MPI_Fint.])
# Match MPI_Fint in reverse to make sure Fortran/C-APIs match
XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define XT_MPI_FINT_FC_KIND $XT_MPI_FINT_FC_KIND"
# the following is needed for source parts from libraries with
# optional MPI dependencies
AC_DEFINE([USE_MPI],[1],[use MPI parallel process environment])
saved_CFLAGS=$CFLAGS
saved_FCFLAGS=$FCFLAGS
saved_LIBS=$LIBS
CFLAGS="$CFLAGS $MPI_C_INCLUDE"
FCFLAGS="$FCFLAGS $MPI_FC_MOD"
LIBS="$LIBS $MPI_C_LIB"
# check MPI defects
AC_ARG_WITH([regard-for-quality],
  [AS_HELP_STRING([--without-regard-for-quality],
     [ignore results of checks for known defects @<:@default: abort on error@:>@])],
  [],[with_regard_for_quality=yes])
ACX_MPIRUN(,[saved_MPI_LAUNCH=$MPI_LAUNCH],
  [AC_MSG_WARN([$MPI_LAUNCH_failMsg])
   AS_VAR_IF([with_regard_for_quality], [yes],
     [AC_MSG_FAILURE([unable to find a working MPI launch program, which is
required for checks for known MPI defects (see --without-regard-for-quality)])],
     [saved_MPI_LAUNCH=$MPI_LAUNCH ; MPI_LAUNCH=true])])
dnl do not remove conftest.err we need to inspect the output
m4_pushdef([_AC_RUN_LOG_LIMIT],
  m4_bpatsubst(m4_dquote(m4_defn([_AC_RUN_LOG_LIMIT])),
    [rm -f conftest.er1 conftest.err],[rm -f conftest.er1]))dnl
m4_pushdef([_AC_LINK_IFELSE],
  m4_bpatsubst(m4_dquote(m4_defn([_AC_LINK_IFELSE])),
    [rm -f core conftest.err],[rm -f core]))dnl
dnl common failure action for MPI defect tests
m4_define([mpi_defect_fail_action],
  [AS_VAR_IF([with_regard_for_quality],[yes],
     [AC_MSG_FAILURE([test for known defect $acx_subtestname failed,
re-configure with --without-regard-for-quality in case you must use the
given MPI implementation, but some tests might fail])],
     [AC_MSG_NOTICE([ignoring failed test for known defect $acx_subtestname,
because configure option --without-regard-for-quality is active])])])
AC_CONFIG_LIBOBJ_DIR([src/mpi-workarounds])
xt_disable_nb_a2a=no
m4_pushdef([AC_LANG(Fortran)],m4_bpatsubst(m4_dquote(m4_defn([AC_LANG(Fortran)])),
  [\(ac_link='.*\)\$LIBS],[\1$XT_LIBS4FC]))
XT_LIBS4FC="$LIBS $MPI_FC_LIB"
ACX_MPI_DEFECTS(,,
  [ACX_MPI_DEFECTS_DOCUMENT([${acx_subtestname}])
   AS_CASE(["$acx_subtestname"],
     [openmpi_datatype],
     [AC_MSG_WARN([known OpenMPI 2.0.3-4.0.1 problem detected, attempting interpositioning work-around!])
      ACX_OMPI_DT_WORKAROUND(,[openmpi_datatype.c],
        [xt_mpi_workaround_LIBS=$ac_cv_search_opal_output
         extrasub="${extrasub+$extrasub$as_nl}/^\# skip internal symbols for stubs created by xlf/i \\\\${as_nl}\# ignore opal symbols overriden by us\\\\${as_nl}/"'^\${exp_sym_prefix}'"${acx_symprfx}"'opal_datatype_(commit|add)/b'
         AS_CASE([$acx_opal_datatype_add_workaround],
           [*/xt_opal_datatype_optimize.c],
           [AC_LIBOBJ([xt_opal_datatype_optimize_wrap])
            XT_WORKAROUND_MPI_SYMBIND="${XT_WORKAROUND_MPI_SYMBIND+$XT_WORKAROUND_MPI_SYMBIND }XT_SYM(opal_datatype_commit)"],
           [*/xt_opal_datatype_add.c],
           [AC_LIBOBJ([xt_opal_datatype_add_wrap])
            XT_WORKAROUND_MPI_SYMBIND="${XT_WORKAROUND_MPI_SYMBIND+$XT_WORKAROUND_MPI_SYMBIND }XT_SYM(opal_datatype_commit)"],
           [AC_MSG_WARN([Unexpected Open MPI datatype bug work-around.])
            mpi_defect_fail_action])
         xt_temp=`expr "${acx_opal_datatype_add_workaround}" : '\(.*\)/@<:@^/@:>@*\.c$'`
         XT_WORKAROUND_MPI_CFLAGS=$OMPI_BUILD_CFLAGS
         XT_WORKAROUND_MPI_CPPFLAGS="-I$xt_temp $OMPI_MCA_C_INCLUDE ${XT_WORKAROUND_MPI_CPPFLAGS-$CPPFLAGS}"],
        [AC_MSG_WARN([Cannot apply Open MPI datatype bug work-around.])
         mpi_defect_fail_action])],
     [mpich_3.4_yaksa_hindexed],
     [AC_MSG_WARN([known MPICH 3.4-3.4.3 problem detected, attempting interpositioning work-around!])
      ACX_MPICH34_DT_WORKAROUND(,[mpich_3.4_yaksa_hindexed.c],
        [extrasub="${extrasub+$extrasub$as_nl}/^\# skip internal symbols for stubs created by xlf/i \\\\${as_nl}\# ignore yaksa symbols overriden by us\\\\${as_nl}/"'^\${exp_sym_prefix}'"${acx_symprfx}"'yaks[ai]_type_create_h?indexed/b'
         AC_LIBOBJ([xt_yaksa_indexed_wrap])dnl
         XT_WORKAROUND_MPI_SYMBIND="${XT_WORKAROUND_MPI_SYMBIND+$XT_WORKAROUND_MPI_SYMBIND }XT_SYM(yaksa_type_create_indexed)"
         XT_WORKAROUND_MPI_CFLAGS=$CFLAGS
         XT_WORKAROUND_MPI_CPPFLAGS="-I$ac_pwd/$acx_yaksa_add_workaround ${XT_WORKAROUND_MPI_CPPFLAGS-$CPPFLAGS}"],
        [AC_MSG_WARN([Cannot apply MPICH YAKSA datatype bug work-around.])
         mpi_defect_fail_action])],
     [mpi_fortran_startup],
     [workaround_successful=false
      AC_MSG_WARN([known problem detected, attempting work-around!])
      ACX_LIBXML2_WORKAROUND(,[mpi_fortran_startup.f90],
        [AC_LIBOBJ([xt_xmlXPathInit])
         workaround_successful=:],
        [AC_MSG_NOTICE(
           [Could not apply libxml2 SIGFPE bug work-around.])])
      AS_IF([$workaround_successful],,
        [ACX_XCODE_COMMONS_WORKAROUND(,[mpi_fortran_startup.f90],
          [LDFLAGS="${LDFLAGS+$LDFLAGS }${acx_workaround}"])])],
     [mpich_unweighted_topo_comm_dup],
     [AC_MSG_WARN([known MPICH 3.0-3.4.2 problem detected, activating work-around!])
      AC_DEFINE([NEED_MPICH_UNWEIGHTED_COMM_DUP_WORKAROUND],,
        [work around problem with MPI_Comm_dup])],
     [incorrect_mpi_abort_exitcode],
     [AC_MSG_WARN([Your installation drops the exit code of MPI_Abort!
Consider getting a better MPI/batch scheduler.])
      AS_IF([grep -v -i 'mpi_abort(@<:@^,()@:>@*, 2)' conftest.err >/dev/null],
        [AC_MSG_NOTICE([Implementing fragile work-around for tests!])
         AC_DEFINE([XT_NEED_MPI_ABORT_WORK_AROUND],[1],
           [To pass the value for MPI_Abort to the surrounding execution environment, a hack is needed.])],
        [grep -v -i 'mpi_abort(@<:@^,()@:>@*, 2)' conftest.error_exit.txt >/dev/null],
        [AC_MSG_NOTICE([Implementing very fragile work-around for tests!])
         AC_DEFINE([XT_NEED_MPI_ABORT_WORK_AROUND],[2])])],
     [mpich_neighbor_alltoallw],
     [AC_MSG_WARN([known MPICH 4.0-4.0.2 problem detected, disabling use of MPI_Neighbor_alltoall!])
      xt_disable_nb_a2a=yes
      mpi_defect_fail_action],
     [mpich_zero_size_contiguous],
     [AC_MSG_WARN([known MPICH 3.4.3, 4.0-4.0.2 problem detected, disabling tests of zero size datatypes!])
      AC_DEFINE([XT_CANNOT_SUPPORT_ZERO_SIZE_DT],[1],
        [defined to 1 if this MPI does not support zero-size datatypes])
      mpi_defect_fail_action],
     [mpi_defect_fail_action])])
m4_popdef([AC_LANG(Fortran)])dnl
m4_popdef([_AC_RUN_LOG_LIMIT])m4_popdef([_AC_LINK_IFELSE])dnl
rm -f conftest.err conftest.error_exit.txt
AC_SUBST([xt_mpi_workaround_LIBS])
AC_SUBST([XT_WORKAROUND_MPI_CFLAGS])
AC_SUBST([XT_WORKAROUND_MPI_CPPFLAGS])
AM_CONDITIONAL([USE_NB_A2A],[test $MPI_VERSION -ge 3 && test x"$xt_disable_nb_a2a" = xno])
AS_IF([test $MPI_VERSION -ge 3 && test x"$xt_disable_nb_a2a" = xno],
  [AC_DEFINE([XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL],[1],
     [defined if MPI_Neighbor_alltoall functions are actually working])])
AM_CONDITIONAL([XT_ADD_MPI_WORKAROUND],
  [test x${acx_opal_datatype_add_workaround+set} = xset \
     -o x${acx_yaksa_add_workaround+set} = xset])
AS_IF([test x${acx_opal_datatype_add_workaround+set} = xset \
         -o x${acx_yaksa_add_workaround+set} = xset],
  [AC_DEFINE_UNQUOTED([XT_WORKAROUND_MPI_SYMBIND],
     [$XT_WORKAROUND_MPI_SYMBIND],
     [symbols to add a dependency on for correct static builds])])
dnl
dnl test if MPI_INT16_T is available
AC_CHECK_DECL([MPI_INT16_T],,
  [TJ_FIND_INTEGRAL_TYPE([int16_t],[xt_mpi_dt_int16_t],[AC_INCLUDES_DEFAULT
@%:@include <mpi.h>],
    [XT_MPI_INT16_T=MPI_[]AS_TR_CPP([$xt_mpi_dt_int16_t])
     AC_MSG_NOTICE([adding macro definition of MPI_INT16_T as $XT_MPI_INT16_T])
     AC_DEFINE_UNQUOTED([MPI_INT16_T],[$XT_MPI_INT16_T],[work around missing define in older MPI implementations])],[:])],
  [AC_INCLUDES_DEFAULT
@%:@include <mpi.h>])
MPI_LAUNCH=$saved_MPI_LAUNCH
FCFLAGS=$saved_FCFLAGS
CFLAGS=$saved_CFLAGS
LIBS=$saved_LIBS

# enable user to select an alternative enumeration type in case
# 2^31 elements isn't enough or too much
AC_ARG_WITH([idxtype],
  AS_HELP_STRING([--with-idxtype],[select type to use for global index values @<:@default=int@:>@]),
  [XT_INT=$withval
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[[char a[((($XT_INT)-1)<0)];
char b[($XT_INT)((-($XT_INT)-0.5)+1.5)]]])],,
     [AC_MSG_ERROR([Selected index type not compilable or unsigned])])],
  [XT_INT=int])
dnl expand setting of variables to describe variable type
dnl vartype_setq(TYPE_NAME, TYPE_C_TYPE, FC_KIND, TYPE_MPI_DT,
dnl              TYPE_VAL_MAX, TYPE_VAL_MIN,
dnl              TYPE_C_PRINTF_FMT, TYPE_C_ABS, TYPE_C_CLZ)
m4_define([vartype_setq],[$1][_MPIDT=$4
     AC_CHECK_SIZEOF([$2])
     $1][_FC_KIND=$3
     if test_with_fortran ; then
     ACX_FORTRAN_TYPE_KIND([integer],[$3],[$2],
       [      use iso_c_binding],,$1[_FC_KIND=$acx_fortran_kind_subst])
     AS_IF([test x"$][$1][_FC_KIND" = x"$3"],
       [XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define $1][_FC_KIND_IN_ISO_C_BINDING"])
     fi
     $1][_MAX=$5
     $1][_MIN=$6
     $1][_FMT=$7
     $1][_abs=$8
     m4_foreach([clz_func],[$9],
       [AS_VAR_IF(clz_func,[yes],
          [$1][_clz="@{:@]AS_LITERAL_WORD_IF([clz_func],
            [m4_substr(clz_func,[17])],
            ["`expr "x]clz_func[" : 'x.\{17\}\(.*\)'`"])[@{:@@{:@m4_if([$2],
            [short], [@{:@unsigned short@:}@])v@:}@@:}@@:}@"])
])])dnl

dnl
dnl test if MPI_INTEGER[248] are available
if test_with_fortran ; then
saved_FCFLAGS=$FCFLAGS
AC_LANG_PUSH([Fortran])
FCFLAGS="$MPI_FC_MOD $FCFLAGS"
for mpi_int2test in MPI_INTEGER2 MPI_INTEGER4 MPI_INTEGER8; do
  AC_MSG_CHECKING([if Fortran MPI bindings define $mpi_int2test])
  AC_COMPILE_IFELSE(AC_LANG_PROGRAM([
        USE mpi
        IMPLICIT NONE
        INTEGER :: i
        i = $mpi_int2test
]),[AC_MSG_RESULT([yes])
  eval HAVE_FC_$mpi_int2test=1],
    [AC_MSG_RESULT([no])
     ASX_VAR_UNSET([HAVE_FC_$mpi_int2test])])
done
FCFLAGS=$saved_FCFLAGS
AC_LANG_POP([Fortran])
fi

m4_foreach([substtype],[[XT_INT]],
  [AS_CASE(["$substtype"],
     [int],
     [vartype_setq(substtype,[int],[c_int],[MPI_INT],
                   [INT_MAX],["(-INT_MAX)"],["\"d\""],['(abs((v)))'],
                   [[acx_cv_have_decl___builtin_clz],
                    [acx_cv_have_decl___lzcnt]])],
     [long],
     [AS_IF([test $ac_cv_sizeof_long -eq $ac_cv_sizeof_long_long],
        [acx_nvidia_intrinsic_presence=acx_cv_have_decl___lzcnt64],
        [test $ac_cv_sizeof_long -eq $ac_cv_sizeof_int],
        [acx_nvidia_intrinsic_presence=acx_cv_have_decl___lzcnt])
      vartype_setq(substtype,[long],[c_long],[MPI_LONG],
                   [LONG_MAX],["(-LONG_MAX)"],["\"ld\""],['(labs((v)))'],
                   [[acx_cv_have_decl___builtin_clzl],
                    [$acx_nvidia_intrinsic_presence]])],
     [short],
     [vartype_setq(substtype,[short],[c_short],[MPI_SHORT],
                   [SHRT_MAX],["(-SHRT_MAX)"],["\"hd\""],
                   ['((short)abs((v)))'],
                   [[acx_cv_have_decl___builtin_clz],
                    [acx_cv_have_decl___lzcnt16]])
      AS_IF([test x"$acx_cv_have_decl___builtin_clz" = xyes \
             && test $ac_cv_sizeof_int != $ac_cv_sizeof_short],
        [AS_IF([expr $ac_cv_sizeof_int = $ac_cv_sizeof_short \* 2 >/dev/null],
           substtype[_clz="($]substtype[_clz - (int)(CHAR_BIT * sizeof (short)))"],
           [AC_MSG_ERROR([Unexpected size relationship of int ($ac_cv_sizeof_int) and short ($ac_cv_sizeof_short)])])])],
     ["long long"],
     [vartype_setq(substtype,[long long],[c_long_long],[MPI_LONG_LONG],
                   [LLONG_MAX],["(-LLONG_MAX)"],["\"lld\""],['(llabs((v)))'],
                   [[acx_cv_have_decl___builtin_clzll],
                    [acx_cv_have_decl___lzcnt64]])],
     [AC_MSG_ERROR([Invalid type selected for substtype])])
   AS_IF([eval eval test x\\\"\\\${HAVE_FC_MPI_INTEGER[]\$ac_cv_sizeof_]AS_TR_SH([$substtype])[+set}\\\" = xset],
     [eval substtype[]_FC_MPIDT=\"MPI_INTEGER[]\$ac_cv_sizeof_]AS_TR_SH([$substtype])[\"],
     substtype[_FC_MPIDT=MPI_DATATYPE_NULL])
   XT_FC_FEATURE_DEFS="$XT_FC_FEATURE_DEFS
@%:@define substtype[]_FC_KIND $substtype[]_FC_KIND
@%:@define substtype[]_FC_MPIDT $substtype[]_FC_MPIDT"
   eval "XT_FC_FEATURE_DEFS=\"\$XT_FC_FEATURE_DEFS
@%:@define SIZEOF_[]substtype[] \$ac_cv_sizeof_]AS_TR_SH([$substtype])[\""
   AC_DEFINE_UNQUOTED([CONF_]substtype[_MAX],[$]substtype[_MAX],
     [Maximal value for ]substtype)
   AC_DEFINE_UNQUOTED([CONF_]substtype[_MIN],[$]substtype[_MIN],
     [Minimal value for ]substtype)
   AC_DEFINE_UNQUOTED([CONF_]substtype[_FMT],[$]substtype[_FMT],
     [printf(3) format string for ]subsstype)
   AC_DEFINE_UNQUOTED(substtype[_ABS(v)],[$]substtype[_abs],
     [return absolute value of v])
   AC_DEFINE_UNQUOTED([SIZEOF_]substtype,[SIZEOF_]AS_TR_CPP([$substtype]),
     [The size of `]substtype[', as computed by sizeof.])
   AS_IF([test -n "$]substtype[_clz"],
     [AC_DEFINE_UNQUOTED(substtype[_CLZ(v)],[$]substtype[_clz],
       [return number of leading zero-bits in v])])
  ])

AS_CASE([$XT_INT],
  [int],[XT_UINT_MPIDT=MPI_UNSIGNED],
  [short],[XT_UINT_MPIDT=MPI_UNSIGNED_SHORT],
  [long],[XT_UINT_MPIDT=MPI_UNSIGNED_LONG],
  ["long long"],[XT_UINT_MPIDT=MPI_UNSIGNED_LONG_LONG])
dnl find type to quickly do double width integer arithmetic on Xt_int,
dnl if available, because results of AC_COMPUTE_INT are limited to the range
dnl of long, sizeof must be used as a proxy for value range
AC_CACHE_CHECK([type to use for long multiplication of Xt_int],
  [acx_cv_XT_LONG],
  [eval acx_sizeof_xt_int=\"\$ac_cv_sizeof_[]AS_TR_SH($XT_INT)\"
   acx_cv_XT_LONG='(unavailable)'
   AS_FOR([itype_macro],[itype],[int long "long long" __int128 __int128_t],
     [eval acx_itype_sizeof=\"\$ac_cv_sizeof_[]AS_TR_SH([itype_macro])\"
      AS_IF([expr "$acx_sizeof_xt_int" \* 2 = "$acx_itype_sizeof" >/dev/null],
        [acx_cv_XT_LONG="itype_macro" ; break])])])
AS_IF([test x"$acx_cv_XT_LONG" != 'x(unavailable)'],
  [AC_DEFINE_UNQUOTED([XT_LONG],[$acx_cv_XT_LONG],[Type to use for fast long multiplication of Xt_int])
   AS_IF([test x"$acx_cv_XT_LONG" != x__int128_t],
     [acx_XT_ULONG="unsigned $acx_cv_XT_LONG"],
     [acx_XT_ULONG="__uint128_t"])
   AC_DEFINE_UNQUOTED([XT_ULONG],[$acx_XT_ULONG],[Type to use for fast long multiplication of Xt_int])],
  [AC_CACHE_CHECK([short type to use for manual long multiplication of Xt_int],
    [acx_cv_XT_SHORT],
    [acx_cv_XT_SHORT='(unavailable)'
     AS_FOR([itype_macro],[itype],[char int long],
       [AS_IF([test itype_macro = char],
          [acx_itype_sizeof=1],
          [eval acx_itype_sizeof=\"\$ac_cv_sizeof_[]AS_TR_SH([itype_macro])\"])
        AS_IF([expr "$acx_sizeof_xt_int" = "$acx_itype_sizeof" \* 2 >/dev/null],
          [acx_cv_XT_SHORT="itype_macro" ; break])])])])
AS_IF([test x"$acx_cv_XT_SHORT" != 'x(unavailable)'],
  [AC_DEFINE_UNQUOTED([XT_SHORT],[$acx_cv_XT_SHORT],[Type to use for manual long multiplication and division of Xt_int])])

AC_DEFINE_UNQUOTED([XT_INT],[$XT_INT],[Type to use for enumeration of distributed elements.])
AC_DEFINE_UNQUOTED([XT_INT_MPIDT],[$XT_INT_MPIDT],[MPI Datatype for type for enumeration of distributed elements.])
AC_DEFINE_UNQUOTED([XT_UINT_MPIDT],[$XT_UINT_MPIDT],[MPI Datatype for type for counts of distributed elements.])
dnl
dnl checks needed for tests only
AC_CHECK_TYPES([long double _Complex])
dnl
AC_ARG_WITH([on-demand-check-programs],
  [AS_HELP_STRING([--with-on-demand-check-programs],
     [only build test programs when calling make check (default: build during make all)])],
  [AS_CASE([$withval],
     [yes|no],,
     [n],[with_on_demand_check_programs=no],
     [y],[with_on_demand_check_programs=yes],
     [AC_MSG_FAILURE([unexpected value for --with-on-demand-check-programs: $withval])])],
  [with_on_demand_check_programs=no])
AM_CONDITIONAL([with_on_demand_check_programs],
  [test x$with_on_demand_check_programs = xyes])

AC_ARG_WITH([example-programs],
  [AS_HELP_STRING([--without-example-programs],
     [only build example programs when calling make examples in examples directory (default: build during make all)])],
  [AS_CASE([$withval],
     [yes|no],,
     [n],[with_example_programs=no],
     [y],[with_example_programs=yes],
     [AC_MSG_FAILURE([unexpected value for --with-example-programs: '$withval'])])],
  [with_example_programs=yes])
AM_CONDITIONAL([with_example_programs],
  [test x$with_example_programs = xyes])

AC_ARG_WITH([perf-programs],
  [AS_HELP_STRING([--without-perf-programs],
     [only build performance test programs when calling make perf in perf directory (default: build during make all)])],
  [AS_CASE([$withval],
     [yes|no],,
     [n],[with_perf_programs=no],
     [y],[with_perf_programs=yes],
     [AC_MSG_FAILURE([unexpected value for --with-perf-programs: '$withval'])])],
  [with_perf_programs=yes])
AM_CONDITIONAL([with_perf_programs],
  [test x$with_perf_programs = xyes])

m4_foreach([flags_var],[[CFLAGS],[FCFLAGS],[LDFLAGS]],
  [AC_ARG_VAR([BUILD_]flags_var,
    [append to ]flags_var[ during build but not in configure phase])
   AC_CONFIG_COMMANDS_PRE([AS_IF([test x"$BUILD_]]flags_var[[" != x],
     ]flags_var[[="$]]flags_var[[ ${BUILD_]]flags_var[[}"])])])
AC_ARG_VAR([BUILD_LIBS],
  [prepend to LIBS during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([AS_IF([test x"$BUILD_LIBS" != x],
  [LIBS="$BUILD_LIBS $LIBS"])])
AC_ARG_VAR([BUILD_CC],
  [replace CC with expansion of $BUILD_CC during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([AS_IF([test x"$BUILD_CC" != x],
  [CC=$BUILD_CC])])
AC_ARG_VAR([BUILD_FC],
  [replace FC with expansion of $BUILD_FC during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([AS_IF([test x"$BUILD_FC" != x],
  [FC=$BUILD_FC])])

AC_CONFIG_HEADERS([include/config.h])
AC_SUBST([XT_FC_FEATURE_DEFS])
AM_SUBST_NOTMAKE([XT_FC_FEATURE_DEFS])
AC_CONFIG_FILES([include/fc_feature_defs.inc])

AC_CONFIG_FILES([						\
	perf/perf_idxsection_get_positions_of_indices_run	\
	perf/perf_sparse_array_gather_run			\
	tests/test_exchanger_parallel_run			\
	tests/test_redist_p2p_parallel_run			\
	tests/test_redist_collection_parallel_run		\
	tests/test_redist_collection_static_parallel_run	\
	tests/test_redist_repeat_parallel_run			\
	tests/test_xmap_all2all_parallel_run			\
	tests/test_xmap_all2all_intercomm_parallel_run		\
	tests/test_xmap_dist_dir_parallel_run			\
	tests/test_xmap_dist_dir_intercomm_parallel_run		\
	tests/test_xmap_intersection_parallel_run		\
        tests/test_initialized_finalized_run			\
	tests/test_idxempty_run					\
	tests/test_idxvec_run					\
	tests/test_idxlist_collection_run			\
	tests/test_idxmod_run					\
	tests/test_idxsection_run				\
	tests/test_idxstripes_run				\
	tests/test_redist_single_array_base_run			\
	tests/test_redist_single_array_base_parallel_run	\
	tests/test_redist_collection_run			\
	tests/test_redist_collection_static_run			\
	tests/test_redist_repeat_run				\
	tests/test_redist_p2p_run				\
	tests/test_request_parallel_run				\
	tests/test_exported_symbols				\
	tests/test_xmap_dist_dir_run				\
	tests/test_xmap_all2all_run				\
	tests/test_xmap_all2all_fail_run			\
	tests/test_ut_run					\
        tests/test_perf_run					\
        tests/test_perf_stripes_run				\
        tests/test_sort_run					\
        tests/test_uid_run					\
        tests/test_misc_run					\
        tests/test_yaxt_run					\
	tests/test_mpi_generate_datatype_run			\
	tests/test_mpi_smartdedup_run				\
	tests/test_ddt_run					\
	examples/rrobin_run					\
	examples/row2col_run					\
	examples/tr_gp2fs_run					\
	examples/unstruct_halo_exchange_run			\
	util/serialrun],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([Makefile src/Makefile src/mpi-workarounds/Makefile tests/Makefile examples/Makefile \
  inst_headers/Makefile inst_headers/f90/Makefile \
  src/pkgconfig/yaxt.pc src/pkgconfig/yaxt_c.pc perf/Makefile])
# edit out the am--include-marker that's required for automake 1.16
# and later on older automake versions
AS_VERSION_COMPARE([$am__api_version],[1.16],[extrasub="${extrasub+$extrasub$as_nl}s/\# am--include-marker\$//"
])

dnl Create dummy Fortran dependency files for Automake version 1.16 or later,
dnl which don't do so automatically anymore, older versions
dnl created the files by the time this command is executed. Having
dnl created the dummy Fortran dependency files, we enforce their
dnl regeneration at the build time by updating the timestamps of the
dnl corresponding Makefiles, which the Fortran dependency file have as
dnl prerequisites.
AC_CONFIG_COMMANDS([force_dep_regeneration],[if test x"$with_fortran" = xyes; then
$MKDIR_P "tests/${DEPDIR}" "src/${DEPDIR}"
touch "tests/${DEPDIR}/FC.deps" "src/${DEPDIR}/FC.deps"
sleep 1
touch tests/Makefile src/Makefile
fi],
[with_fortran=`test -z "${with_fortran_TRUE}" && echo yes || echo no`
DEPDIR="${DEPDIR}"])

AC_OUTPUT
dnl
dnl Local Variables:
dnl mode: autoconf
dnl coding: utf-8
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
