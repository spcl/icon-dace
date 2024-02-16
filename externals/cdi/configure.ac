#  Process this file with autoconf to produce a configure script.
#
#  automake 1.16.1
#  autoconf 2.69
#  libtool  2.4.6

AC_PREREQ([2.69])
LT_PREREQ([2.4.6])

AC_INIT([cdi],[2.2.4],[https://mpimet.mpg.de/cdi])
AC_DEFINE_UNQUOTED(CDI, ["$PACKAGE_VERSION"], [CDI version])

AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
dnl forbid acx macro names from remaining unexpanded
m4_pattern_forbid([^ACX_])
AC_CANONICAL_HOST
AC_CANONICAL_BUILD

AM_INIT_AUTOMAKE([1.16.1 foreign])
AM_MAINTAINER_MODE([enable])
AM_EXTRA_RECURSIVE_TARGETS([examples])

m4_pushdef([_AC_PROG_CC_C89],
  [rm -rf conftest.dSYM # needed when configured on MacOS with CFLAGS='-g'
   _AC_PROG_CC_C99([], [AC_MSG_FAILURE([C compiler does not support ISO C99])])])dnl
AC_PROG_CC
m4_popdef([_AC_PROG_CC_C89])dnl
ACX_PROG_CC_POSIX([2001])
AC_C_RESTRICT

dnl Compile with MPI support (we declare the option here because its value
dnl affects values of the auto options below):
AC_ARG_ENABLE([mpi],
  [AS_HELP_STRING([--enable-mpi],
     [enable parallel I/O with MPI @<:@default=no@:>@])],
  [test "x$enableval" != xno && enable_mpi=yes],
  [enable_mpi=no])

AC_ARG_ENABLE([iso-c-interface],
  [AS_HELP_STRING([--enable-iso-c-interface],
     [create Fortran 90 interface using Fortran 2003 ISO_C_BINDING facility ]dnl
[@<:@default=no@:>@])],
  [test "x$enableval" != xno && enable_iso_c_interface=yes],
  [enable_iso_c_interface=no])
AM_CONDITIONAL([ENABLE_ISOC_INTERFACE],
  [test "x$enable_iso_c_interface" = xyes])

AC_ARG_ENABLE([cf-interface],
  [AS_HELP_STRING([--enable-cf-interface],
     [create Fortran 77 interface using cfortran.h @<:@default=auto@:>@])],
  [AS_IF([test "x$enableval" != xno && test "x$enableval" != xauto],
     [enable_cf_interface=yes])],
  [enable_cf_interface=auto])
AS_IF([test "x${enable_mpi}${enable_cf_interface}" = xnoauto],
  [enable_cf_interface=no])

dnl We need FC for the Fortran 90 and for the Fortran 77 interfaces. In the
dnl latter case, the compiler is used for testing. Note that we use neither
dnl AS_IF nor AS_VAR_IF here to prevent checking for FC when we do not need it:
if test "x$enable_iso_c_interface" = xyes || \
   test "x$enable_cf_interface" != xno; then
  if test "x$FC" != xno; then
    AC_PROG_FC
    test -z $FC 2>/dev/null && FC=no
    AS_IF([test "x$FC" != xno], [AC_FC_SRCEXT([f90], [], [FC=no])])
    AS_IF([test "x$FC" != xno],
      [AC_CACHE_CHECK([whether the Fortran compiler works],
         [acx_cv_fc_works],
         [acx_cv_fc_works=no
          AC_LANG_PUSH([Fortran])
          AC_RUN_IFELSE([AC_LANG_PROGRAM],
            [acx_cv_fc_works=yes], [], [acx_cv_fc_works=yes])
          AC_LANG_POP([Fortran])])
       AS_VAR_IF([acx_cv_fc_works], [no], [FC=no])])
dnl Additional tests with fallback values:
    AS_IF([test "x$FC" != xno],
      [AC_LANG_PUSH([Fortran])
       ACX_SL_FC_CHECK_MOD_PATH_FLAG([], [FC_MOD_FLAG=-I])
       ACX_FORTRAN_INCLUDE_FLAG(
        [AC_SUBST([FC_OPTINC], [$acx_cv_fc_ftn_include_flag])],
        [FC_OPTINC=-I])
       ACX_SL_FC_MOD_SUFFIX([], [FCMODEXT=mod])
dnl Find the Fortran compiler flag needed to specify a preprocessor macro
dnl definition for *.F90 source files:
       m4_pushdef([AC_FC_PP_SRCEXT])dnl
       saved_ac_fc_srcext=$ac_fc_srcext
       ac_fc_srcext='F90'
       AC_FC_PP_DEFINE([], [FC_DEFINE=-D])
       ac_fc_srcext=$saved_ac_fc_srcext
       m4_popdef([AC_FC_PP_DEFINE])
       AC_LANG_POP([Fortran])])
  fi
else
  FC=no
fi
AM_CONDITIONAL([FC_MOD_UPPERCASE], [test "x$FCMODCASE" = xuc])
dnl Additional checks required for the Fortran 2003 interface:
AS_VAR_IF([enable_iso_c_interface], [yes],
  [AS_IF([test "x$FC" != xno], [ACX_FC_CHECK_STRPTR_CONVERT([], [FC=no])])
   AS_VAR_IF([FC], [no],
     [AC_MSG_FAILURE([the Fortran 2003 interface is requested but the ]dnl
[Fortran compiler is disabled, missing or lacks the required features])])])
AM_CONDITIONAL([BUILD_FC_PROGRAMS], [test "x$FC" != xno])

dnl We need Fortran 77 if either FC is not available or the user has specified
dnl F77 explicitly to make sure that the generated Fortran 77 interface is
dnl compatible with it:
if test "x$enable_cf_interface" != xno; then
  AS_IF([test -n "$F77" && test "x$F77" != xno],
    [requested_F77=yes], [requested_F77=no])
  if test "x$FC" = xno || test "x$requested_F77" = xyes; then
    AC_PROG_F77
    test -z $F77 2>/dev/null && F77=no
    AS_IF([test "x$F77" != xno],
      [AC_CACHE_CHECK([whether the Fortran 77 compiler works],
         [acx_cv_f77_works],
         [acx_cv_f77_works=no
          AC_LANG_PUSH([Fortran 77])
          AC_RUN_IFELSE([AC_LANG_PROGRAM],
            [acx_cv_f77_works=yes], [], [acx_cv_f77_works=yes])
          AC_LANG_POP([Fortran 77])])
       AS_VAR_IF([acx_cv_f77_works], [no], [F77=no])])
  else
    F77=no
  fi
  AS_IF([test "x$F77" = xno && test "x$requested_F77" = xyes],
    [AC_MSG_FAILURE([the explicitly requested Fortran 77 compiler is ]dnl
[missing or lacks the required features])])
  AS_UNSET([requested_F77])
else
  F77=no
fi
dnl Additional checks required for the Fortran 77 interface:
AS_IF([test "x$enable_cf_interface" != xno],
  [dnl
dnl The following macro make sure that the user's request for the Fortran 77
dnl interface compatibility is fulfilled:
   ACX_FIND_CFORTRAN_DEF
   ACX_XLF_QEXTNAME_ADD_APPENDUS
   ACX_CHECK_CFORTRAN([$srcdir/src],
     [AC_DEFINE([HAVE_CF_INTERFACE], [1],
        [Defined to 1 if C / Fortran interface cfortran.h works])])
   AS_IF([test "x$F77" != xno],
     [AC_LANG_PUSH([Fortran 77])
      ACX_FORTRAN_INCLUDE_FLAG
      AC_LANG_POP([Fortran 77])
      AS_IF([test "x$FC" != xno && dnl
test "x$acx_cv_f77_ftn_include_flag" != "x$FC_OPTINC"],
        [AC_MSG_ERROR([Fortran and Fortran 77 compilers require different ]dnl
[flags needed to specify search paths for the "INCLUDE" statements])])])
   AS_IF([test "x$FC" = xno && test "x$F77" = xno],
     [AS_VAR_IF([enable_cf_interface], [auto], [enable_cf_interface=no],
        [AC_MSG_FAILURE([the Fortran 77 interface is requested but both ]dnl
[Fortran and Fortran 77 compilers are disabled, missing or lack the ]dnl
[required features])])])])
AS_VAR_IF([enable_cf_interface], [auto], [enable_cf_interface=yes])
AM_CONDITIONAL([ENABLE_CF_INTERFACE], [test "x$enable_cf_interface" = xyes])

AC_ARG_ENABLE([ruby-interface],
  [AS_HELP_STRING([--enable-ruby-interface],
     [create Ruby language interface (EXPERIMENTAL) @<:@default=no@:>@])],
  [test "x$enableval" != xno && enable_ruby_interface=yes],
  [enable_ruby_interface=no])

AC_ARG_ENABLE([python-interface],
  [AS_HELP_STRING([--enable-python-interface],
     [create Python language interface (EXPERIMENTAL) @<:@default=no@:>@])],
  [test "x$enableval" != xno && enable_python_interface=yes],
  [enable_python_interface=no])

if test "x$enable_ruby_interface" = xyes || \
   test "x$enable_python_interface" = xyes; then
  if test "x$CXX" != xno; then
    AC_PROG_CXX
    test -z $CXX 2>/dev/null && CXX=no
    AS_IF([test "x$CXX" != xno],
      [AC_CACHE_CHECK([whether the C++ compiler works],
         [acx_cv_cxx_works],
         [acx_cv_cxx_works=no
          AC_LANG_PUSH([C++])
dnl Do not use AC_RUN_IFELSE for now:
          AC_LINK_IFELSE([AC_LANG_PROGRAM], [acx_cv_cxx_works=yes])
          AC_LANG_POP([C++])])
       AS_VAR_IF([acx_cv_cxx_works], [no], [CXX=no])])
  fi
  AS_VAR_IF([CXX], [no],
    [AC_MSG_FAILURE([either the Ruby or the Python interface is requested ]dnl
[but the C++ compiler is disabled, missing or lacks the required features])])
else
  CXX=no
dnl The Automake conditional "am__fastdepCXX" is not set in this case because
dnl we skipped the check for CXX. Prevent the configure script from failing:
  AM_CONDITIONAL([am__fastdepCXX], [false])
fi

AC_PROG_INSTALL

dnl Override the option set by AC_OPENMP:
AC_ARG_ENABLE([openmp],
  [AS_HELP_STRING([--enable-openmp],
     [enable OpenMP support @<:@default=no@:>@])],
  [test "x$enableval" != xno && enable_openmp=yes],
dnl TODO: decide whether we want to compile with OpenMP support by default.
dnl Before this change, we ran AC_OPENMP but did not append OPENMP_CFLAGS to
dnl CFLAGS, which resulted in no OpenMP support. We preserve that behaviour:
  [enable_openmp=no])
m4_pushdef([AC_ARG_ENABLE])dnl
AC_OPENMP
m4_popdef([AC_ARG_ENABLE])dnl
AS_IF(
  [test "x$enable_openmp" = xno],
  [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#ifndef _OPENMP
 choke me
#endif]])],
     [AC_MSG_ERROR([OpenMP support is disabled but the compiler enables it dnl
by default: change CFLAGS to disable it])])],
  [test "x$ac_cv_prog_c_openmp" = xunsupported],
  [AC_MSG_FAILURE([cannot link C OpenMP programs])],
  [test -n "$OPENMP_CFLAGS"],
  [AS_VAR_APPEND([CFLAGS], [" $OPENMP_CFLAGS"])])

dnl Set up libtool:
AC_MSG_NOTICE([setting up libtool])
ACX_USE_LIBTOOL_CONFIGURATION([pic-only])

dnl Check endianess of system:
AC_C_BIGENDIAN
dnl Check large file support on 32 bit systems:
AC_SYS_LARGEFILE
dnl Checks for library functions:
AC_FUNC_MMAP
dnl Check for structures:
AC_CHECK_MEMBERS([struct stat.st_blksize])
dnl Checks for header files:
AC_CHECK_HEADERS([malloc.h unistd.h sys/time.h],,,[AC_INCLUDES_DEFAULT])
AX_EXECINFO
dnl Checks for the availability of functions:
AC_CHECK_FUNCS([mallinfo pwrite])
AC_REPLACE_FUNCS([getline])
dnl Check for UUID library:
ACX_UUID
AS_VAR_APPEND([CPPFLAGS], [" $UUID_C_INCLUDE"])
LIBS="$UUID_C_LIB $LIBS"
dnl Checks for the availability of ANSI-C99 functions:
AC_CHECK_DECLS([isnan], [], [],
  [AC_INCLUDES_DEFAULT
@%:@include <math.h>])
dnl Check for sysconf names:
AC_CHECK_DECLS(
  [_SC_LARGE_PAGESIZE, PAGESIZE, PAGE_SIZE, _SC_PAGE_SIZE, _SC_PAGESIZE, dnl
_PC_REC_XFER_ALIGN, POSIX_REC_XFER_ALIGN], [], [],
  [AC_INCLUDES_DEFAULT
@%:@include <limits.h>
@%:@include <unistd.h>])
dnl Check for math library:
AC_SEARCH_LIBS(floor, m)

dnl Check for non-standard builtin:
AS_FOR([builtin_macro], [builtin], [__builtin_ctz],
  [AS_VAR_PUSHDEF([builtin_cache], [acx_cv_have_decl_]builtin_macro)
   AC_CACHE_CHECK([whether ]builtin_macro[ is declared],
     [builtin_cache],
     [AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([], [  unsigned lbz = builtin_macro][[(56U)]])],
        [AS_VAR_SET([builtin_cache], [yes])],
        [AS_VAR_SET([builtin_cache], [no])])])
   AS_VAR_IF([builtin_cache], [yes],
     [AC_DEFINE_UNQUOTED([HAVE_DECL_]AS_TR_CPP([builtin_macro]), [1])
      break],
     [AC_DEFINE_UNQUOTED([HAVE_DECL_]AS_TR_CPP([builtin_macro]), [0])])])
AH_TEMPLATE([HAVE_DECL___BUILTIN_CTZ],
  [Define to 1 if __builtin_ctz is available, 0 if not])

dnl TODO: add comments on why we do not have all options in one place:
dnl Add basic configure options:
ACX_CDI_OPTIONS
AM_CONDITIONAL([ENABLE_NETCDF], [test "x$ENABLE_NETCDF" = xyes])
AC_CHECK_FUNCS([grib_get_length])

dnl TODO: check whether we really need all the following information. Is it just
dnl for 'app/cdi -d'? If so, we can skip this when we are not going to install
dnl the executable:
AC_DEFINE_UNQUOTED(COMPILER, ["$CC $CFLAGS"], [Compiler])
dnl TODO: replace the following with a macro that can do the detection better:
dnl Check compiler version:
COMP_VERSION=
AS_CASE([$CC],
  [pgcc*], [COMP_VERSION=`$CC -V | head -2 | tail -n 1`],
  [*gcc*], [COMP_VERSION=`$CC --version | head -n 1`],
  [g++*], [COMP_VERSION=`$CC --version | head -n 1`],
  [clang*], [COMP_VERSION=`$CC --version | head -n 1`],
  [sxc*], [COMP_VERSION=`$CC -V 2>&1   | tail -n 1`],
  [xlc*], [COMP_VERSION=`$CC -qversion 2>&1   | head -n 1`],
  [COMP_VERSION=`$CC -V 2>&1 | sed -n 1p`])
test -z "$COMP_VERSION" && COMP_VERSION=unknown
AC_DEFINE_UNQUOTED([COMP_VERSION], ["$COMP_VERSION"], [Compiler version])
dnl
dnl Checks for system type:
AC_DEFINE_UNQUOTED([SYSTEM_TYPE], ["$ac_cv_build"], [System type])
AC_SUBST([SYSTEM_TYPE], ["$ac_cv_build"])

dnl MPI-related options and checks:
AC_ARG_ENABLE([ppm-dist-array],
  [AS_HELP_STRING([--enable-ppm-dist-array],
     [enable usage of the PPM distributed array @<:@default=auto@:>@])],
  [AS_CASE([$enableval],
     [no|auto], [],
     [AS_VAR_IF([enable_mpi], [no],
        [AC_MSG_ERROR([usage of the PPM distributed array can be enabled dnl
(--enable-ppm-dist-array) only if the parallel I/O with MPI is enabled dnl
(--enable-mpi)])])
      enable_ppm_dist_array=yes])],
  [enable_ppm_dist_array=auto])

build_pio_fc_programs=no
have_ppm=no
AC_SUBST([HAVE_PPM_DIST_ARRAY], [no])
AM_SUBST_NOTMAKE([HAVE_PPM_DIST_ARRAY])
have_parallel_nc4=no
AS_VAR_IF([enable_mpi], [yes],
  [dnl
dnl Check for MPI C interface:
   ACX_C_PACKAGE([MPI], [mpi.h],,,
     [AC_MSG_FAILURE([required header mpi.h not found or not compilable])],
     [MPI_Waitall],[mpi mpich], [], [],
     [AC_MSG_FAILURE([cannot link C MPI programs])])

dnl Check for MPI_LAUNCH command:
   saved_CPPFLAGS=$CPPFLAGS
   saved_LIBS=$LIBS
   AS_VAR_APPEND([CPPFLAGS], [" $MPI_C_INCLUDE"])
   LIBS="$MPI_C_LIB $LIBS"
   ACX_MPIRUN(,,[AC_MSG_WARN([MPI launch command not found or not compilable])])
   CPPFLAGS=$saved_CPPFLAGS
   LIBS=$saved_LIBS

dnl Check for MPI Fortran interface (if needed):
   build_pio_fc_programs=$enable_cf_interface
   test "x$FC" = xno && build_pio_fc_programs=no
   AS_VAR_IF([build_pio_fc_programs], [yes],
     [ACX_F90_PACKAGE([MPI], [mpi],,,
     [AC_MSG_WARN([required Fortran module mpi not found or not compilable])
      build_pio_fc_programs=no],
     [mpi_waitall], [mpi mpi_f90 mpi_f77 mpich], [[-lmpi_f77 -lmpi],[-lmpi]],,
     [AC_MSG_WARN([cannot link Fortran MPI programs])
      build_pio_fc_programs=no],,
     [[INTEGER :: req(1), stat(mpi_status_size, 1), ierror]],
     [[(1, req, stat, ierror)]])])

dnl Check for YAXT C interface:
dnl
dnl Declare pkg-config precious variables:
   PKG_PROG_PKG_CONFIG
dnl Set precious variables following the naming convention of ACX_C_PACKAGE:
   AC_ARG_VAR([YAXT_C_INCLUDE],
     [specifically set flags to use when compiling sources using YAXT ]dnl
[includes.])dnl
   AC_ARG_VAR([YAXT_C_LIB],
     [specifically set flags to use when linking YAXT.])dnl
dnl Prevent PKG_CHECK_MODULES from declaring package-specific precious
dnl variables, which names do not follow our convention:
   m4_pushdef([AC_ARG_VAR])dnl
dnl PKG_CHECK_MODULES relies on yaxt_c_CFLAGS and yaxt_c_LIBS, which in our case
dnl have names YAXT_C_INCLUDE and YAXT_C_LIB, respectively:
   yaxt_c_CFLAGS=$YAXT_C_INCLUDE
   yaxt_c_LIBS=$YAXT_C_LIB
   PKG_CHECK_MODULES([yaxt_c],[yaxt_c],
     [YAXT_C_INCLUDE=$yaxt_c_CFLAGS
      YAXT_C_LIB=$yaxt_c_LIBS],
     [dnl
dnl There was no yaxt_c.pc before version 0.5.2, therefore we try with yaxt.pc.
dnl PKG_CHECK_MODULES relies on yaxt_CFLAGS and yaxt_LIBS, which in our case
dnl have names YAXT_C_INCLUDE and YAXT_C_LIB, respectively:
      yaxt_CFLAGS=$YAXT_C_INCLUDE
      yaxt_LIBS=$YAXT_C_LIB
      PKG_CHECK_MODULES([yaxt], [yaxt],
        [YAXT_VERSION=`$PKG_CONFIG --modversion yaxt`
         YAXT_MAJOR_VERSION=`expr X"$YAXT_VERSION" : "X\(@<:@0-9@:>@*\)"`
         YAXT_MINOR_VERSION=`echo "$YAXT_VERSION" | sed -e 's/^@<:@0-9@:>@*\.\(@<:@0-9@:>@*\).*/\1/'`
         AS_IF([test "$YAXT_MAJOR_VERSION" -gt 0 -o "$YAXT_MINOR_VERSION" -gt 4],
           [YAXT_C_INCLUDE=$yaxt_CFLAGS
            YAXT_C_LIB=`echo "$yaxt_LIBS" | sed -e 's/-lyaxt *$/-lyaxt_c/'`])],
        [:])])
   m4_popdef([AC_ARG_VAR])dnl
   saved_CPPFLAGS=$CPPFLAGS
   saved_LIBS=$LIBS
   AS_VAR_APPEND([CPPFLAGS], [" $MPI_C_INCLUDE"])
   LIBS="$MPI_C_LIB $LIBS"
dnl The following is extracted from ACX_GENERIC_PACKAGE, which is expanded by
dnl ACX_C_PACKAGE. We do not use the latter macro here to avoid introduction of
dnl the '--with-yaxt-XXX' configure options, which might be very confusing since
dnl whatever pkg-config finds will be used instead.
   ACX_LANG_CHECK_INCLUDE_PATHS_IFELSE([yaxt.h],,
     [YAXT_C_INCLUDE=$acx_cv_c_include_yaxt_h
      ACX_OPTION_SEARCH_LIBS_MULTI([xt_initialized], [yaxt_c yaxt],
        [YAXT_C_LIB=$acx_cv_option_search_xt_initialized_c],
        [AC_MSG_FAILURE([cannot link C YAXT programs])],,[$YAXT_C_LIB])],
     [AC_MSG_FAILURE([required header yaxt.h not found or not compilable])],,
     [$YAXT_C_INCLUDE],[[]])
dnl
   AH_TEMPLATE([YAXT_UID_DT],
     [Defined to MPI datatype to be used for Xt_uid])dnl
   defined_Xt_uid=no
   AC_CHECK_DECL([MPI_UINT64_T],
     [AC_DEFINE([YAXT_UID_DT], [MPI_UINT64_T])
      defined_Xt_uid=yes],
     [AS_VAR_APPEND([CPPFLAGS], [" $YAXT_C_INCLUDE"])
      LIBS="$YAXT_C_LIB $LIBS"
      AC_CHECK_SIZEOF([Xt_uid],,
        [AC_INCLUDES_DEFAULT
@%:@include <yaxt.h>])
      AC_CHECK_SIZEOF([unsigned long])
      AS_IF([test "$ac_cv_sizeof_unsigned_long" -eq "$ac_cv_sizeof_Xt_uid"],
        [AC_DEFINE([YAXT_UID_DT], [MPI_UNSIGNED_LONG])
         defined_Xt_uid=yes],
        [AC_CHECK_DECLS([MPI_UNSIGNED_LONG_LONG],
           [AC_CHECK_SIZEOF([unsigned long long])
            AS_IF([test "$ac_cv_sizeof_unsigned_long_long" -eq "$ac_cv_sizeof_Xt_uid"],
              [AC_DEFINE([YAXT_UID_DT],[MPI_UNSIGNED_LONG_LONG])
               defined_Xt_uid=yes])],,
           [AC_INCLUDES_DEFAULT
@%:@include <mpi.h>])])],
     [AC_INCLUDES_DEFAULT
@%:@include <mpi.h>])
   AS_VAR_IF([defined_Xt_uid], [no],
     [AC_MSG_FAILURE([no way to communicate Xt_uid found])])
   AS_UNSET([defined_Xt_uid])
   CPPFLAGS=$saved_CPPFLAGS
   LIBS=$saved_LIBS

dnl Check for YAXT Fortran interface (if needed):
   AS_VAR_IF([build_pio_fc_programs], [yes],
     [dnl
dnl Set precious variables following the naming convention of ACX_F90_PACKAGE:
      AC_ARG_VAR([YAXT_FC_MOD],
        [flags to enable 'USE YAXT' in Fortran program.])dnl
      AC_ARG_VAR([YAXT_FC_LIB],
        [specifically set flags to use when linking YAXT.])dnl
dnl Prevent PKG_CHECK_MODULES from declaring package-specific precious
dnl variables, which names do not follow our convention:
      m4_pushdef([AC_ARG_VAR])dnl
dnl PKG_CHECK_MODULES relies on yaxt_CFLAGS and yaxt_LIBS, which in our case
dnl have names YAXT_FC_MOD and YAXT_FC_LIB, respectively:
      yaxt_CFLAGS=$YAXT_FC_MOD
      yaxt_LIBS=$YAXT_FC_LIB
      PKG_CHECK_MODULES([yaxt],[yaxt],
        [YAXT_FC_MOD=$yaxt_CFLAGS
         YAXT_FC_LIB=$yaxt_LIBS], [:])
      m4_popdef([AC_ARG_VAR])dnl
      saved_FCFLAGS=$FCFLAGS
      saved_LIBS=$LIBS
      AS_VAR_APPEND([FCFLAGS], [" $MPI_FC_MOD"])
      LIBS="$MPI_FC_LIB $LIBS"
dnl The following is extracted from ACX_F90_PACKAGE. We do not use the
dnl aforementioned macro here to avoid introduction of the '--with-yaxt-XXX'
dnl configure options, which might be very confusing since whatever pkg-config
dnl finds will be used instead.
      AC_LANG_PUSH([Fortran])dnl
      ACX_FORTRAN_CHECK_MOD_PATHS_IFELSE([yaxt],,
        [YAXT_FC_MOD=$acx_cv_fortran_mod_yaxt
         AS_VAR_APPEND([FCFLAGS], [" $YAXT_FC_MOD"])
dnl Note that the following test must check for a function implemented in the
dnl Fortran library, as well as for a function implemented in the C library and
dnl available via C bindings:
         ACX_OPTION_SEARCH_LIBS_MULTI([xt_initialized],[yaxt],
           [YAXT_FC_LIB=$acx_cv_option_search_xt_initialized_fc],
           [AC_MSG_WARN([cannot link Fortran YAXT programs])
            build_pio_fc_programs=no],[[-lyaxt_c]],
           [$YAXT_FC_LIB],
           [[      use yaxt
      logical is_initialized]],
           [[      is_initialized = xt_initialized()
      call xt_finalize()]])],
        [AC_MSG_WARN([required Fortran module yaxt not found or not compilable])
         build_pio_fc_programs=no],,
        [$YAXT_FC_MOD])
      AC_LANG_POP([Fortran])dnl
      FCFLAGS=$saved_FCFLAGS
      LIBS=$saved_LIBS])

dnl Parallel netCDF support still requires ScalES-PPM and YAXT to re-arrange the
dnl data when running with more than one collector.
dnl
dnl TODO: check whether it makes sense to check for PPM when configuring with
dnl --disable-ppm-dist-array.
dnl
dnl Set precious variables following the naming convention of ACX_C_PACKAGE:
   AC_ARG_VAR([PPM_CORE_C_INCLUDE],
     [specifically set flags to use when compiling sources using PPM_CORE ]dnl
[includes.])dnl
   AC_ARG_VAR([PPM_CORE_C_LIB],
     [specifically set flags to use when linking PPM_CODE.])dnl
dnl Prevent PKG_CHECK_MODULES from declaring package-specific precious
dnl variables, which names do not follow our convention:
   m4_pushdef([AC_ARG_VAR])dnl
dnl PKG_CHECK_MODULES relies on ppm_core_CFLAGS and ppm_core_LIBS, which in our
dnl case have names PPM_CORE_C_INCLUDE and PPM_CORE_C_LIB, respectively:
   ppm_core_CFLAGS=$PPM_CORE_C_INCLUDE
   ppm_core_LIBS=$PPM_CORE_C_LIB
   PKG_CHECK_MODULES([ppm_core],[scales-ppm-core],
     [PPM_CORE_C_INCLUDE=$ppm_core_CFLAGS
      PPM_CORE_C_LIB=$ppm_core_LIBS], [:])
   m4_popdef([AC_ARG_VAR])dnl
   saved_CPPFLAGS=$CPPFLAGS
   saved_LIBS=$LIBS
   AS_VAR_APPEND([CPPFLAGS], [" $MPI_C_INCLUDE"])
   LIBS="$MPI_C_LIB $LIBS"
dnl The following is extracted from ACX_GENERIC_PACKAGE, which is expanded by
dnl ACX_C_PACKAGE. We do not use the latter macro here to avoid introduction of
dnl the '--with-ppm-core-XXX' configure options, which might be very confusing
dnl since whatever pkg-config finds will be used instead.
   ACX_LANG_CHECK_INCLUDE_PATHS_IFELSE([ppm/ppm.h],,
     [PPM_CORE_C_INCLUDE=$acx_cv_c_include_ppm_ppm_h
      ACX_OPTION_SEARCH_LIBS_MULTI([PPM_initialize], [scalesppmcore],
        [PPM_CORE_C_LIB=$acx_cv_option_search_PPM_initialize_c
         have_ppm=yes
         AC_DEFINE([HAVE_PPM_CORE],,[ScalES PPM C core library is available])
         AS_VAR_APPEND([CPPFLAGS], [" $PPM_CORE_C_INCLUDE"])
         AS_IF([test "x$enable_ppm_dist_array" != xno],
           [AC_CHECK_HEADERS([ppm/dist_array.h],
              [HAVE_PPM_DIST_ARRAY=yes], [], [AC_INCLUDES_DEFAULT])])],
        [AC_MSG_WARN([cannot link C SCALES-PPM-CORE programs])],,
        [$PPM_CORE_C_LIB])],
     [AC_MSG_WARN([required header ppm/ppm.h not found or not compilable])],,
     [$PPM_CORE_C_INCLUDE],[[]])
   CPPFLAGS=$saved_CPPFLAGS
   LIBS=$saved_LIBS

dnl NC_HAS_PNETCDF is useful when available to short-circuit some tests
   AS_VAR_IF([ENABLE_NETCDF], [yes],
     [AC_CHECK_HEADERS([netcdf_meta.h],,,
        [AC_INCLUDES_DEFAULT
@%:@include <netcdf.h>])
dnl The following test is relevant only for older versions of netCDF. Starting
dnl version 4.3.3, netCDF provides the NC_HAS_PNETCDF macro and its value
dnl overrides the result of the test.
      AC_CACHE_CHECK([parallel netcdf's pnetcdf support],
        [acx_cv_have_pnetcdf],
        [acx_cv_have_pnetcdf=no
         test "x$NC_CONFIG" != "x" && \
         test "x$($NC_CONFIG --has-pnetcdf)" = "xyes" && \
         acx_cv_have_pnetcdf=yes])
      AS_VAR_IF([acx_cv_have_pnetcdf], [yes],
        [AC_DEFINE([HAVE_NETCDF_PAR_PNETCDF], [1],
           [Defined to 1 if NetCDF parallel open supports NC_PNETCDF])])])

dnl If scales-ppm is not available, netcdf can only be used in serial mode:
   AS_IF([test "x$ENABLE_NC4" = xyes && test "x$have_ppm" = xyes],
     [saved_CPPFLAGS=$CPPFLAGS
      AS_VAR_APPEND([CPPFLAGS], [" $MPI_C_INCLUDE"])
      AC_CHECK_HEADERS([netcdf_par.h],,,
        [AC_INCLUDES_DEFAULT
@%:@include <netcdf.h>])
      AC_CHECK_DECL([MPI_Bcast],
        [have_parallel_nc4=yes
         AC_CHECK_DECLS([nc_inq_format_extended],,
           [AC_MSG_NOTICE([The supplied netCDF library does not support ]dnl
[nc_inq_format_extended.])],
           [AC_INCLUDES_DEFAULT
@%:@include <netcdf.h>])
         AC_DEFINE([HAVE_PARALLEL_NC4], [1],
           [netCDF library does support MPI parallel invocations])],
        [AC_MSG_WARN([the supplied netCDF library does not support MPI ]dnl
[parallel invocations])],
        [AC_INCLUDES_DEFAULT
@%:@include <netcdf.h>
@%:@ifdef HAVE_NETCDF_PAR_H
@%:@include <netcdf_par.h>
@%:@endif])
      CPPFLAGS=$saved_CPPFLAGS
dnl WARNING: the following macro overwrites saved_CFLAGS shell variable and
dnl modifies CFLAGS:
      ACX_TLS_XLC_RETRY])

dnl Check for librt:
   AC_SUBST([LIBRT], [''])
   AC_CHECK_DECL([_POSIX_ASYNCHRONOUS_IO],
     [saved_LIBS=$LIBS
      AC_SEARCH_LIBS([aio_suspend], [rt],
        [test "$ac_res" = "none required" || LIBRT=$ac_res])
      LIBS=$saved_LIBS])],
  [YAXT_C_INCLUDE=; YAXT_C_LIB=
   MPI_C_INCLUDE=; MPI_C_LIB=])
AS_CASE(["$enable_ppm_dist_array$HAVE_PPM_DIST_ARRAY"],
  [yesyes|no*], [],
  [auto*], [enable_ppm_dist_array=$HAVE_PPM_DIST_ARRAY],
  [AC_MSG_FAILURE([usage of the PPM distributed array is requested but the ]dnl
[required headers and libraries are not found])])
AS_VAR_IF([have_ppm], [no], [PPM_CORE_C_INCLUDE=; PPM_CORE_C_LIB=])
AS_VAR_IF([build_pio_fc_programs], [no],
  [YAXT_FC_MOD=; YAXT_FC_LIB=
   MPI_FC_MOD=; MPI_FC_LIB=])
AM_CONDITIONAL([BUILD_PIO_FC_PROGRAMS], [test "x$build_pio_fc_programs" = xyes])
AM_CONDITIONAL([ENABLE_MPI], [test "x$enable_mpi" = xyes])
AM_CONDITIONAL([HAVE_PARALLEL_NC4], [test "x$have_parallel_nc4" = xyes])

AM_MISSING_PROG([SWIG], [swig])

AS_VAR_IF([enable_ruby_interface], [yes],
  [AC_ARG_VAR([RUBY], [the Ruby interpreter])
   AC_PATH_PROGS([RUBY], [ruby])
   ruby_CPPFLAGS=`$RUBY $srcdir/config/interface.rb 2>/dev/null`
   AS_IF([test $? -ne 0],
     [AC_MSG_ERROR([Ruby interface is requested but the Ruby interpreter is ]dnl
[not found])])
   saved_CPPFLAGS=$CPPFLAGS
   AS_VAR_APPEND([CPPFLAGS], [" $ruby_CPPFLAGS"])
   AC_CHECK_HEADER([ruby.h], [],
     [AC_MSG_FAILURE([Ruby interface is requested but the required header ]dnl
[ruby.h not found or not compilable])],
     [AC_INCLUDES_DEFAULT])
   CPPFLAGS=$saved_CPPFLAGS],
  [AM_MISSING_PROG([RUBY], [ruby])])
AM_CONDITIONAL([ENABLE_RUBY_INTERFACE], [test "x$enable_ruby_interface" = xyes])

AS_VAR_IF([enable_python_interface], [yes],
  [dnl
dnl For consistency with the other commands (MPI_LAUNCH, SWIG, RUBY, etc.), we
dnl ignore $PYTHON if it is not an absolute path:
  AS_VAR_SET_IF([PYTHON],
    [AS_CASE([$PYTHON],
       [[[\\/]]* | ?:[[\\/]]*], [],
       [PYTHON=])])
  AM_PATH_PYTHON([], [], [:])
dnl AM_PATH_PYTHON does not check whether $PYTHON really works and
dnl AX_PYTHON_DEVEL fails with a misleading message if it doesn't:
   AS_IF([$PYTHON -c 'import sys' >/dev/null 2>&1], [],
     [AC_MSG_ERROR([Python interface is requested but the Python ]dnl
[interpreter is not found])])
dnl Do not check for extra flags by default:
   PYTHON_EXTRA_LDFLAGS=${PYTHON_EXTRA_LDFLAGS-' '}
   PYTHON_EXTRA_LIBS=${PYTHON_EXTRA_LIBS-' '}
   AX_PYTHON_DEVEL])
AM_CONDITIONAL([ENABLE_PYTHON_INTERFACE],
  [test "x$enable_python_interface" = xyes])

AC_ARG_WITH([on-demand-check-programs],
  [AS_HELP_STRING([--with-on-demand-check-programs],
     [only build test programs when calling 'make check' @<:@default: build ]dnl
[during 'make all'@:>@])],
  [test "x$withval" != xno && with_on_demand_check_programs=yes],
  [with_on_demand_check_programs=no])
AM_CONDITIONAL([with_on_demand_check_programs],
  [test "x$with_on_demand_check_programs" = xyes])

AC_ARG_WITH([concurrent-mpi-checks],
  [AS_HELP_STRING([--with-concurrent-mpi-checks],
     [allow for simultaneous runs of MPI tests @<:@default: run MPI test ]dnl
[programs one after another@:>@])],
  [test "x$withval" != xno && with_concurrent_mpi_checks=yes],
  [with_concurrent_mpi_checks=no])
AM_CONDITIONAL([with_concurrent_mpi_checks],
  [test "x$with_concurrent_mpi_checks" = xyes])

AC_ARG_WITH([example-programs],
  [AS_HELP_STRING([--without-example-programs],
     [only build example programs when calling 'make examples' ]dnl
[@<:@default: build during 'make all'@:>@])],
  [test "x$withval" != xno && with_example_programs=yes],
  [with_example_programs=yes])
AM_CONDITIONAL([with_example_programs], [test "x$with_example_programs" = xyes])

m4_foreach([build_flag_var],[[BUILD_CFLAGS],[BUILD_FCFLAGS],[BUILD_LDFLAGS],[BUILD_LIBS],[BUILD_MPI_C_LIB],[BUILD_MPI_FC_LIB]],
  [AC_ARG_VAR(build_flag_var,
     [append to ]m4_bpatsubst(build_flag_var, [BUILD_], [])[ during build but not in configure phase])dnl
AC_CONFIG_COMMANDS_PRE(m4_bpatsubst(build_flag_var, [BUILD_], [])[="$]m4_bpatsubst(build_flag_var, [BUILD_], [])[$]{build_flag_var:+ $build_flag_var[}"])dnl
AM_SUBST_NOTMAKE(build_flag_var)])dnl

m4_foreach([build_tool_var],[[BUILD_CC],[BUILD_CXX],[BUILD_FC],[BUILD_F77]],
  [AC_ARG_VAR(build_tool_var,
     [replace ]m4_bpatsubst(build_tool_var, [BUILD_], [])[ with expansion of $]build_tool_var[ during build but not in configure phase])dnl
AC_CONFIG_COMMANDS_PRE(m4_bpatsubst(build_tool_var, [BUILD_], [])[="$]{build_tool_var:-$[]m4_bpatsubst(build_tool_var, [BUILD_], [])}["])dnl
AM_SUBST_NOTMAKE(build_tool_var)])dnl

dnl Additional output variables for cdi.settings:
AS_VAR_IF([enable_mpi],[yes],
  [ENABLE_MPI=true],
  [ENABLE_MPI=false])
AC_SUBST([ENABLE_MPI])
AM_SUBST_NOTMAKE([ENABLE_MPI)])
AS_VAR_IF([have_parallel_nc4],[yes],
  [HAVE_PARALLEL_NC4=1],
  [HAVE_PARALLEL_NC4=0])
AC_SUBST([HAVE_PARALLEL_NC4])
AM_SUBST_NOTMAKE([HAVE_PARALLEL_NC4])
AS_VAR_IF([enable_ruby_interface],[yes],
  [ENABLE_RUBY_INTERFACE=true],
  [ENABLE_RUBY_INTERFACE=false])
AC_SUBST([ENABLE_RUBY_INTERFACE])
AM_SUBST_NOTMAKE([ENABLE_RUBY_INTERFACE])
AS_VAR_IF([enable_python_interface],[yes],
  [ENABLE_PYTHON_INTERFACE=true],
  [ENABLE_PYTHON_INTERFACE=false])
AC_SUBST([ENABLE_PYTHON_INTERFACE])
AM_SUBST_NOTMAKE([ENABLE_PYTHON_INTERFACE])

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([
  Makefile
  app/Makefile
  cdi.settings
  examples/Makefile
  examples/pio/Makefile
  interfaces/Makefile
  src/Makefile
  src/cmake/cdi/cdi-config-version.cmake
  src/cmake/cdi/cdi-config.cmake
  src/pkgconfig/cdi.pc
  src/pkgconfig/cdi_f2003.pc
  src/pkgconfig/cdipio.pc
  tests/Makefile
  tests/test_cksum_base
])
AC_CONFIG_FILES([
  interfaces/test_cpp_nc.run
  interfaces/test_python_grb.run
  interfaces/test_python_nc.run
  interfaces/test_ruby_grb.run
  interfaces/test_ruby_nc.run
  tables/gen_tableheaderfile
  tests/pio_cksum_asynch.run
  tests/pio_cksum_fpguard.run
  tests/pio_cksum_grb2.run
  tests/pio_cksum_mpi_fw_at_all.run
  tests/pio_cksum_mpi_fw_at_reblock.run
  tests/pio_cksum_mpi_fw_ordered.run
  tests/pio_cksum_mpinonb.run
  tests/pio_cksum_nc.run
  tests/pio_cksum_nc2.run
  tests/pio_cksum_nc4.run
  tests/pio_cksum_writer.run
  tests/pio_write.parallel.run
  tests/pio_write.run
  tests/pio_write_deco2d.parallel.run
  tests/pio_write_deco2d.run
  tests/pio_write_dist_array.run
  tests/test_byteswap.run
  tests/test_cdf_const.run
  tests/test_cdf_transformation.run
  tests/test_cksum_ext.run
  tests/test_cksum_grb.run
  tests/test_cksum_grb2.run
  tests/test_cksum_ieg.run
  tests/test_cksum_nc.run
  tests/test_cksum_nc2.run
  tests/test_cksum_nc4.run
  tests/test_cksum_nc_chunk.run
  tests/test_cksum_srv.run
  tests/test_f2003.run
  tests/test_grib.run
  tests/test_resource_copy.parallel.run
  tests/test_resource_copy.run
  tests/test_table.run
],
[chmod a+x "$ac_file"])

AC_OUTPUT

AC_MSG_NOTICE([CDI is configured with the following options:])
cat cdi.settings >&AS_MESSAGE_FD

AC_MSG_NOTICE([

  Configuration completed.

    You can now run 'make' to compile the CDI package and 'make install' to install it afterwards.
])
