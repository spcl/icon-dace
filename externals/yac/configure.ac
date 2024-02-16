AC_PREREQ([2.69])
AC_INIT([yac], [3.0.3], [hanke@dkrz.de], [],
  [https://dkrz-sw.gitlab-pages.dkrz.de/yac/])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([src/grid.c])

dnl Generate the private and public config headers and avoid C-style block
dnl comments in it to safely include it into a Fortran source file:
AC_CONFIG_HEADERS([src/config.h src/yac_config.h],
  [AS_IF([test x"$acx_cv_fc_pp_comments" != xyes],
     [sed 's%/\*.*\*/%%g;\%/\*%,\%\*/%d;/^$/d' "$ac_file" dnl
> "$tmp/tmp_config.h" || AC_MSG_ERROR([could not postprocess $ac_file])
      rm -f "$ac_file"
      mv "$tmp/tmp_config.h" "$ac_file" || dnl
AC_MSG_ERROR([could not postprocess $ac_file])])],
  [acx_cv_fc_pp_comments=$acx_cv_fc_pp_comments])
dnl Prevent yac_config from overriding macros defined in config.h when both
dnl headers are included:
AC_DEFINE([YAC_CONFIG_H], [], [do not include yac_config.h])
AC_DEFINE_UNQUOTED([YAC_VERSION], ["$PACKAGE_VERSION"], [YAC version])

AC_CANONICAL_HOST

dnl Enable silent rules by default with a workaround for the case when the
dnl version constraint above is relaxed:
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AM_INIT_AUTOMAKE([1.16.1 foreign subdir-objects])
AM_MAINTAINER_MODE

dnl Tests for the C compiler:
dnl   Make sure conftest.dSYM is removed when configured on MacOS with
dnl   CFLAGS='-g':
m4_pushdef([_AM_PROG_CC_C_O], m4_bpatsubst(m4_dquote(m4_defn([_AM_PROG_CC_C_O])),
  [rm -f\(.* conftest\*.*\)], [rm -rf\1]))dnl
dnl   Check for ISO C99 support
dnl   (Autoconf 2.69- and 2.70+ have different documented ways for this):
m4_case(m4_version_compare(AC_AUTOCONF_VERSION,[2.70]),
  [-1], [AC_PROG_CC_C99],
  [0],  [dnl
dnl   Fix a bug in Autoconf 2.70 (see https://savannah.gnu.org/support/?110396):
         m4_pushdef([_AC_PROG_CC_C99],
         m4_bpatsubst(m4_dquote(m4_defn([_AC_PROG_CC_C99])),
           [\[ac_c_conftest_c89_program\]], [[ac_c_conftest_c99_program]]))dnl
         AC_PROG_CC
         m4_popdef([_AC_PROG_CC_C99])],
  [1],  [AC_PROG_CC])
AS_VAR_IF([ac_cv_prog_cc_c99], [no],
  [AC_MSG_FAILURE([unable to detect C compiler flag needed to accept ISO C99])])
m4_popdef([AC_PROG_CC_C_O])dnl
ACX_COMPILER_CC_VENDOR
AC_C_RESTRICT

dnl Tests for the Fortran compiler:
AC_PROG_FC
AC_LANG_PUSH([Fortran])
AC_FC_SRCEXT([f90])
dnl Automake does not trace AC_FC_PP_SRCEXT/ACX_FC_PP_SRCEXT and therefore
dnl does not add their results to the compilation rules. But AC_FC_SRCEXT is
dnl traced, therefore, we temporarily redefine AC_FC_SRCEXT so it would
dnl expand to ACX_FC_PP_SRCEXT:
m4_pushdef([AC_FC_SRCEXT], [ACX_FC_PP_SRCEXT($@)])
AC_FC_SRCEXT([F90])
m4_popdef([AC_FC_SRCEXT])
dnl Handle possibly long expansions of the preprocessor macro __FILE__:
ACX_FC_LINE_LENGTH([unlimited])
ACX_FC_MODULE_NAMING([AC_SUBST([FCMODEXT], ["$acx_cv_fc_module_naming_ext"])])
AM_CONDITIONAL([FCMODUC], [test x"$acx_cv_fc_module_naming_upper" = xyes])
ACX_FC_MODULE_IN_FLAG([AC_SUBST([FCMODINC], ["$acx_cv_fc_module_in_flag"])])
ACX_FC_INCLUDE_FLAG_PP([AC_SUBST([FCPPINC], ["$acx_cv_fc_pp_include_flag"])])
ACX_LANG_MACRO_FLAG([AC_SUBST([FCPPDEF], ["$acx_cv_fc_macro_flag"])])
ACX_FC_MODULE_CHECK([ISO_C_BINDING])
ACX_FC_PP_COMMENTS([], [:])

AS_IF([test "`cd $srcdir && pwd`" != "`pwd`"],
  [AC_SUBST([yac_fcisrc], [" ${FCPPINC}\$(srcdir)"])dnl
   AM_SUBST_NOTMAKE([yac_fcisrc])])

AC_CACHE_CHECK([whether $FC can properly detect non-contiguous data],
  [yac_cv_fc_is_contiguous_works],
  [AC_RUN_IFELSE([AC_LANG_PROGRAM([], [[
      type dbl_ptr
      double precision, pointer :: ptr(:)
      end type
      double precision, contiguous, pointer :: mem(:,:)
      integer :: nproma = 8, nblks = 1000, psize
      type(dbl_ptr) :: ptrs(3)
      allocate(mem(nproma,nblks))
      psize=nproma*nblks
      ptrs(1)%ptr(1:psize) => mem
      ptrs(2)%ptr(1:nproma) => mem(:,2)
      ptrs(3)%ptr(1:nblks) => mem(2,:)
      call check_cont(ptrs)
      contains
      subroutine check_cont(ptrs)
      type(dbl_ptr), intent(in) :: ptrs(:)
      if (.not. is_contiguous(ptrs(1)%ptr) .or. &
      .not. is_contiguous(ptrs(2)%ptr) .or. &
      is_contiguous(ptrs(3)%ptr)) then
      stop 1
      end if
      end subroutine check_cont]])],
     [yac_cv_fc_is_contiguous_works=yes],
     [yac_cv_fc_is_contiguous_works=no],
     [yac_cv_fc_is_contiguous_works=unknown])])
AS_IF(
  [test "x$yac_cv_fc_is_contiguous_works" = xno],
  [AC_MSG_FAILURE([$FC cannot properly detect non-contiguous data])],
  [test "x$yac_cv_fc_is_contiguous_works" = xunknown],
  [AC_MSG_WARN([unable to check whether $FC can properly detect dnl
non-contiguous data: this might result into a non-working YAC build])])

AC_LANG_POP([Fortran])

dnl Tests for static library generation tools:
AC_ARG_VAR([AR], [archiver command])
AC_CHECK_TOOL([AR], [ar], [:])
AS_VAR_IF([AR], [:], [AC_MSG_ERROR([the archiver command not found])])
AC_ARG_VAR([ARFLAGS], [archiver flags])
ARFLAGS=${ARFLAGS-cr}
AC_ARG_VAR([RANLIB], [archive indexer command])
AC_PROG_RANLIB

AC_ARG_VAR([MPI_LAUNCH], [MPI launcher command])

AC_ARG_ENABLE([lib-only],
  [AS_HELP_STRING([--enable-lib-only],
     [omit building examples and utility programs @<:@default=no@:>@])], [],
  [enable_lib_only=no])
AM_CONDITIONAL([ENABLE_LIB_ONLY], [test x"$enable_lib_only" = xyes])

dnl Variables for pkg-config file generation:
AC_SUBST([YAC_PKGCONF_CLIBS], [''])
AM_SUBST_NOTMAKE([YAC_PKGCONF_CLIBS])

AC_ARG_WITH([external-lapack],
[  --with-external-lapack  specify the interface to be used for enabling LAPACK
                          functions. The value of the argument is either 'yes',
                          'no', 'auto' or a space-separated list of the
                          following:
                              'mkl' - Intel MKL LAPACKE interface;
                              'lapacke' - Netlib LAPACKE interface;
                              'atlas' - ATLAS CLAPACK interface;
                              'clapack' - Netlib CLAPACK interface;
                              'fortran' - Fortran LAPACK interface.
                          If set to 'no', the bundled (very limited) version of
                          CLAPACK library is used. @<:@default=auto@:>@], [],
[with_external_lapack=auto])

AC_ARG_WITH([external-mtime],
  [AS_HELP_STRING([--with-external-mtime],
     [use external mtime instead of the bundled one @<:@default=no@:>@])], [],
  [with_external_mtime=no])

dnl Check for required libraries. There are several considerations regarding
dnl how the detection of the required libraries should be implemented:
dnl   1. A library might directly or indirectly depend on many other libraries.
dnl      For example, NetCDF, depending on how it is configured, might require
dnl      a subset of the following libraries: -lhdf5_hl, -lhdf5, -ldf,
dnl      -lpthread, -lsz, -lm, -lcurl, -lssl, -lcrypto -lz, -ljpeg, -lmpi,
dnl      -lszip, -lpnetcdf, etc. This makes it very hard to guess how to link
dnl      to -lnetcdf correctly. Especially if some of the libraries are static.
dnl      Thus, we try only -lnetcdf, and libnetcdf.so is expected to exist and
dnl      to have all the NEEDED and RPATH entries required to link it
dnl      correctly. If it's not the case, the users are expected to set
dnl      LDFLAGS, LIBS or PACKAGE_LANG_LIB (e.g. NETCDF_C_LIB, MPI_FC_LIB,
dnl      etc.) accordingly.
dnl   2. The same is applicable for Fortran intrinsic and runtime libraries
dnl      (e.g. libgfortran.so) required for linking a C application to a
dnl      Fortran library (e.g. liblapack.so). There is a standard Autoconf
dnl      macro AC_FC_LIBRARY_LDFLAGS that could help with this but it does not
dnl      always work (e.g. for NAG compiler), it does not know about
dnl      vendor-specific flags (e.g. -pgf90libs), and it's result often contain
dnl      many redundant flags. Thus, the users are expected to provide the
dnl      required flags manually when they want to link to a Fortran library
dnl      that does not have the NEEDED and RPATH entries for Fortran intrinsic
dnl      and runtime Fortran libraries. The only exception implemented so far
dnl      is for PGI compiler.

dnl Get absolute paths to the source and build directories:
_AC_SRCDIRS([.])
yac_abs_top_srcdir=$ac_abs_top_srcdir
yac_abs_top_builddir=$ac_abs_top_builddir

dnl Check for libm.
AC_SEARCH_LIBS([cos], [m])

dnl Check for pthread.
dnl We build and run the multithreading tests if pthread is found and does not
dnl require switching the compiler executable (i.e. $PTHREAD_CC = $CC):
yac_have_pthread=no
AX_PTHREAD([AS_VAR_IF([PTHREAD_CC], ["$CC"], [yac_have_pthread=yes])])
AS_VAR_IF([yac_have_pthread], [no], [PTHREAD_LIBS=;PTHREAD_CFLAGS=])
AM_CONDITIONAL([TEST_PTHREAD], [test x"$yac_have_pthread" = xyes])

dnl Checks for MPI libraries.
AC_ARG_ENABLE([mpi-checks],
  [AS_HELP_STRING([--enable-mpi-checks],
     [enable configure-time checks of MPI library for known defects ]dnl
[@<:@default=yes@:>@])], [],
  [enable_mpi_checks=yes])

dnl Checks for C MPI headers and libraries.

AC_ARG_VAR([MPI_CFLAGS], [exact C compiler flags enabling MPI])
AC_ARG_VAR([MPI_CLIBS],
  [exact linker flags enabling MPI when linking with C compiler])

yac_save_CFLAGS=$CFLAGS; yac_save_LIBS=$LIBS

CFLAGS="$MPI_CFLAGS $CFLAGS"
AC_CHECK_HEADER([mpi.h],
  [ACX_LANG_LIB_SEARCH([MPI_CLIBS], [MPI_Init], [],
     [MPI_CLIBS=$acx_cv_c_lib_func_MPI_Init
      AS_IF([test -n "$MPI_CLIBS"],
        [YAC_PKGCONF_CLIBS="$MPI_CLIBS $YAC_PKGCONF_CLIBS"])
      LIBS="$MPI_CLIBS $LIBS"],
     [AC_MSG_FAILURE([cannot link to MPI C library])])],
  [AC_MSG_FAILURE([mpi.h not found])], [//])

dnl Check for MPI_LAUNCH command.
ACX_PROG_MPIRUN([MPI_LAUNCH],
  [MPI_LAUNCH=$acx_cv_prog_mpirun],
  [AS_VAR_IF([enable_mpi_checks], [yes],
     [AC_MSG_FAILURE([unable to find a valid MPI launch program, which dnl
is required for the MPI checks (see --enable-mpi-checks)])])
   MPI_LAUNCH=])
AM_CONDITIONAL([TEST_MPI], [AS_VAR_TEST_SET([acx_cv_prog_mpirun])])

TJ_FIND_TYPE([MPI_Fint], [], [[#include <mpi.h>]], [int long 'long long' short],
  [AS_CASE([$tj_cv_c_type_MPI_Fint],
     [int], [YAC_MPI_FINT_MAX=INT_MAX; YAC_MPI_FINT_FC_KIND=c_int],
     [long], [YAC_MPI_FINT_MAX=LONG_MAX; YAC_MPI_FINT_FC_KIND=c_long],
     ['long long'], [YAC_MPI_FINT_MAX=LLONG_MAX
                     YAC_MPI_FINT_FC_KIND=c_long_long],
     [short], [YAC_MPI_FINT_MAX=SHRT_MAX; YAC_MPI_FINT_FC_KIND=c_short],
     [AC_MSG_ERROR(
        [unexpected matching typedef '$tj_cv_c_type_MPI_Fint' for MPI_Fint])])
   AC_DEFINE_UNQUOTED([YAC_MPI_FINT_FC_KIND], [$YAC_MPI_FINT_FC_KIND],
     [Fortran kind type parameter for MPI_Fint])],
  [AC_MSG_FAILURE([unable to find a matching typedef for MPI_Fint])])

CFLAGS=$yac_save_CFLAGS; LIBS=$yac_save_LIBS

dnl Checks for Fortran MPI headers, modules and libraries.

AC_LANG_PUSH([Fortran])
AC_ARG_VAR([MPI_FCFLAGS], [exact Fortran compiler flags enabling MPI])
AC_ARG_VAR([MPI_FCLIBS],
  [exact linker flags enabling MPI when linking with Fortran compiler])

yac_save_FCFLAGS=$FCFLAGS; yac_save_LIBS=$LIBS

FCFLAGS="$MPI_FCFLAGS $FCFLAGS"
ACX_FC_MODULE_CHECK([MPI],
  [ACX_LANG_LIB_SEARCH([MPI_FCLIBS], [mpi_init], [],
     [MPI_FCLIBS=$acx_cv_fc_lib_func_MPI_INIT
      LIBS="$MPI_FCLIBS $LIBS"],
     [AC_MSG_FAILURE([cannot link to MPI Fortran library])],
     [AC_LANG_PROGRAM([],[[      use MPI, only : MPI_INIT
      implicit none
      integer :: ierror
      call MPI_INIT(ierror)]])])],
  [AC_MSG_FAILURE([MPI Fortran module is not available])])

dnl Check for compatibility of Fortran and C MPI libraries
yac_save_CFLAGS=$CFLAGS; yac_save_MPI_LAUNCH=$MPI_LAUNCH
CFLAGS="$MPI_CFLAGS $CFLAGS"
AS_VAR_SET_IF([acx_cv_prog_mpirun], [], [MPI_LAUNCH=true])
ACX_FC_C_COMPATIBLE_MPI([$MPI_LAUNCH], [],
  [AC_MSG_FAILURE([FC and C MPI libraries are not compatible])])
CFLAGS=$yac_save_CFLAGS; MPI_LAUNCH=$yac_save_MPI_LAUNCH

FCFLAGS=$yac_save_FCFLAGS; LIBS=$yac_save_LIBS
AC_LANG_POP([Fortran])

AC_ARG_ENABLE([concurrent-mpi-tests],
  [AS_HELP_STRING([--enable-concurrent-mpi-tests],
     [allow for simultaneous runs of MPI tests @<:@default=no@:>@])], [],
  [enable_concurrent_mpi_tests=no])
AM_CONDITIONAL([ENABLE_CONCURRENT_MPI_TESTS],
  [test x"$enable_concurrent_mpi_tests" = xyes])

dnl Checks for YAXT library.

ACX_LANG_PACKAGE_INIT([YAXT], [-I], [-L])

yac_save_CFLAGS=$CFLAGS
YAXT_CFLAGS=${YAXT_CFLAGS-$acx_c_YAXT_inc_search_args}
CFLAGS="$MPI_CFLAGS $YAXT_CFLAGS $CFLAGS"
AC_CHECK_HEADER([yaxt.h],
  [ACX_LANG_LIB_SEARCH([YAXT_CLIBS], [xt_initialized],
     [ASX_PREPEND_LDFLAGS([$acx_c_YAXT_lib_search_args], [-lyaxt_c])],
     [YAXT_CLIBS=$acx_cv_c_lib_func_xt_initialized
      AS_IF([test -n "$YAXT_CLIBS"],
        [YAC_PKGCONF_CLIBS="$YAXT_CLIBS $YAC_PKGCONF_CLIBS"])],
     [AC_MSG_FAILURE([cannot link to YAXT library])])],
  [AC_MSG_FAILURE([yaxt.h not found])], [//])
CFLAGS=$yac_save_CFLAGS

AC_LANG_PUSH([Fortran])
ACX_LANG_PACKAGE_INIT([YAXT], ["$FCMODINC"], [-L])
yac_have_yaxt_fc=no

yac_save_FCFLAGS=$FCFLAGS
YAXT_FCFLAGS=${YAXT_FCFLAGS-$acx_fc_YAXT_inc_search_args}
FCFLAGS="$YAXT_FCFLAGS $FCFLAGS"
ACX_FC_MODULE_CHECK([YAXT],
  [ACX_LANG_LIB_SEARCH([YAXT_FCLIBS], [xt_initialized],
     [ASX_PREPEND_LDFLAGS([$acx_c_YAXT_lib_search_args],
        [-lyaxt], [-lyaxt -lyaxt_c])],
     [YAXT_FCLIBS=$acx_cv_fc_lib_func_XT_INITIALIZED
      yac_have_yaxt_fc=yes],
     [AC_MSG_WARN([cannot link to YAXT Fortran library, dnl
the respective tests will be skipped])],
     [AC_LANG_PROGRAM([],
dnl Note that the following test must check for a function implemented in the
dnl Fortran library, as well as for a function implemented in the C library and
dnl available via C bindings:
[[      use YAXT, only : XT_INITIALIZED, XT_FINALIZE
      implicit none
      print *, XT_INITIALIZED()
      call XT_FINALIZE()]])])],
  [AC_MSG_WARN([YAXT Fortran module is not available, dnl
the respective tests will be skipped])])
FCFLAGS=$yac_save_FCFLAGS
AC_LANG_POP([Fortran])
AS_VAR_IF([yac_have_yaxt_fc], [no], [YAXT_FCLIBS=;YAXT_FCFLAGS=])
AM_CONDITIONAL([TEST_YAXT_FC], [test x"$yac_have_yaxt_fc" = xyes])

dnl Checks for MPI defects and features.

AS_VAR_IF([enable_mpi_checks], [yes],
  [yac_save_CFLAGS=$CFLAGS; yac_save_LIBS=$LIBS
   CFLAGS="$MPI_CFLAGS $YAXT_CFLAGS $CFLAGS"
   LIBS="$MPI_CLIBS $LIBS"
   ACX_LANG_MPI_CHECK([$MPI_LAUNCH],
     [config/checksrc], [],
     [ACX_LANG_MPI_CHECK_FAIL_MSG([yac_mpi_error],
        [$acx_prog_mpi_test_file], [config/checkdoc])
      AC_MSG_FAILURE([
Your MPI library failed test "$acx_prog_mpi_test_file".

$yac_mpi_error

It is strongly recommended to switch to a version of MPI that is not affected
with the issue.

Alternatively, you can disable this test with the configure argument
"--disable-mpi-checks" but this might result into a non-working YAC build.
])])
   CFLAGS=$yac_save_CFLAGS; LIBS=$yac_save_LIBS],
  [AC_MSG_NOTICE([skipping MPI checks])])

dnl Checks for NetCDF library.

AC_ARG_ENABLE([netcdf],
  [AS_HELP_STRING([--enable-netcdf],
     [enable NetCDF support @<:@default=auto@:>@])], [],
  [enable_netcdf=auto])

ACX_LANG_PACKAGE_INIT([NetCDF], [-I], [-L])

AS_IF([test x"$enable_netcdf" != xno],
  [yac_save_CFLAGS=$CFLAGS
   NETCDF_CFLAGS=${NETCDF_CFLAGS-$acx_c_NetCDF_inc_search_args}
   CFLAGS="$NETCDF_CFLAGS $CFLAGS"
   AC_CHECK_HEADER([netcdf.h],
     [ACX_LANG_LIB_SEARCH([NETCDF_CLIBS], [nc_open],
        [ASX_PREPEND_LDFLAGS([$acx_c_NetCDF_lib_search_args], [-lnetcdf])],
        [NETCDF_CLIBS=$acx_cv_c_lib_func_nc_open
         AS_IF([test -n "$NETCDF_CLIBS"],
           [YAC_PKGCONF_CLIBS="$NETCDF_CLIBS $YAC_PKGCONF_CLIBS"])
         AC_DEFINE([YAC_NETCDF_ENABLED], [], [Enable NetCDF support])
         AC_MSG_NOTICE([NetCDF support is enabled])
         enable_netcdf=yes],
        [AS_VAR_IF([enable_netcdf], [auto],
           [AC_MSG_WARN([cannot link to NetCDF library, dnl
NetCDF support is disabled])
            enable_netcdf=no],
           [AC_MSG_FAILURE([cannot link to NetCDF library])])])],
     [AS_VAR_IF([enable_netcdf], [auto],
        [AC_MSG_WARN([netcdf.h not found, NetCDF support is disabled])
         enable_netcdf=no],
        [AC_MSG_FAILURE([netcdf.h not found])])], [//])
   CFLAGS=$yac_save_CFLAGS])

AS_VAR_IF([enable_netcdf], [no], [NETCDF_CLIBS=;NETCDF_CFLAGS=])
AM_CONDITIONAL([HAVE_NETCDF], [test x"$enable_netcdf" = xyes])
AM_CONDITIONAL([TEST_NETCDF], [test x"$enable_netcdf" = xyes])

dnl Checks for fyaml.

ACX_LANG_PACKAGE_INIT([fyaml], [-I], [-L])

yac_save_CFLAGS=$CFLAGS
FYAML_CFLAGS=${FYAML_CFLAGS-$acx_c_fyaml_inc_search_args}
CFLAGS="$FYAML_CFLAGS $CFLAGS"
AC_CHECK_HEADER([libfyaml.h],
  [ACX_LANG_LIB_SEARCH([FYAML_CLIBS],
     [fy_node_mapping_lookup_scalar0_by_simple_key],
     [ASX_PREPEND_LDFLAGS([$acx_c_fyaml_lib_search_args], [-lfyaml])],
     [FYAML_CLIBS=$acx_cv_c_lib_func_fy_node_mapping_lookup_scalar0_by_simple_key
      AS_IF([test -n "$FYAML_CLIBS"],
        [YAC_PKGCONF_CLIBS="$FYAML_CLIBS $YAC_PKGCONF_CLIBS"])],
     [AC_MSG_FAILURE([cannot link to libfyaml library])])],
  [AC_MSG_FAILURE([libfyaml.h not found])], [//])
CFLAGS=$yac_save_CFLAGS

dnl Checks for mtime library.

ACX_LANG_PACKAGE_INIT([mtime], [-I], [-L])

AS_VAR_IF([with_external_mtime], [yes],
  [yac_save_CFLAGS=$CFLAGS
   MTIME_CFLAGS=${MTIME_CFLAGS-$acx_c_mtime_inc_search_args}
   CFLAGS="$MTIME_CFLAGS $CFLAGS"
   AC_CHECK_HEADER([mtime_datetime.h],
     [ACX_LANG_LIB_SEARCH([MTIME_CLIBS], [datetimeToString],
        [ASX_PREPEND_LDFLAGS([$acx_c_mtime_lib_search_args], [-lmtime])],
        [MTIME_CLIBS=$acx_cv_c_lib_func_datetimeToString
         AS_IF([test -n "$MTIME_CLIBS"],
           [YAC_PKGCONF_CLIBS="$MTIME_CLIBS $YAC_PKGCONF_CLIBS"])],
        [AC_MSG_FAILURE([cannot link to mtime library])])],
     [AC_MSG_FAILURE([mtime_datetime.h not found, enable usage of the dnl
bundled mtime library by re-configuring with the argument dnl
--without-external-mtime])], [//])
   CFLAGS=$yac_save_CFLAGS],
  [MTIME_CFLAGS="-I$yac_abs_top_srcdir/mtime/include"
   MTIME_CLIBS="$yac_abs_top_builddir/mtime/libyac_mtime.a"
   YAC_PKGCONF_CLIBS="-lyac_mtime $YAC_PKGCONF_CLIBS"])
AM_CONDITIONAL([EXTERNAL_MTIME], [test x"$with_external_mtime" = xyes])

dnl Checks for LAPACK library functions LAPACKE_dgels/dgels.
dnl There are several options for enabling LAPACK functions:
dnl   1. Intel MKL LAPACKE interface.
dnl   2. Netlib LAPACKE interface (Netlib LAPACK, OpenBLAS, etc.).
dnl   3. ATLAS CLAPACK interface (incompatible with Netlib CLAPACK).
dnl   4. Netlib CLAPACK interface.
dnl   5. Name-mangling interface to an external Fortran LAPACK library.
dnl   6. bundled Netlib CLAPACK library containing only those functions that
dnl      are required by YAC.
dnl Users can set the command line argument '--with-external-lapack' to specify
dnl which of the options the configure script should consider. The value of the
dnl argument is either one or a space-separated list of the following:
dnl   mkl - Intel MKL LAPACKE interface;
dnl   lapacke - Netlib LAPACKE interface;
dnl   atlas - ATLAS CLAPACK interface;
dnl   clapack - Netlib CLAPACK interface;
dnl   fortran - Fortran LAPACK interface.
dnl If '--with-external-lapack=yes', the configure script goes through the list
dnl of all possible values (in the mentioned order) and takes the first usable
dnl option. If none of the options is usable, the configuration fails. If
dnl '--with-external-lapack=no', the configure script enables usage of the
dnl bundled CLAPACK library. If '--with-external-lapack' is set to a single
dnl option, the configure script checks only that option and fails if the check
dnl fails. If '--with-external-lapack' is set to a list of options, the
dnl configure script goes through the list (respecting order) and enables the
dnl first usable option. If none of the options is usable, the configure script
dnl fails. By default, '--with-external-lapack=auto', which means the same as
dnl '--with-external-lapack=yes' but if none of the options is usable, the
dnl configure script enables usage of the bundled CLAPACK library.

AH_VERBATIM([YAC_LAPACK_INTERFACE_ID],
  [/* Define to the ID of LAPACK interface to be used:
   1 - for Intel MKL LAPACKE;
   2 - for Netlib LAPACKE;
   3 - for ATLAS CLAPACK;
   4 - for Netlib CLAPACK;
   5 - for Fortran LAPACK. */
#undef YAC_LAPACK_INTERFACE_ID])

ACX_LANG_PACKAGE_INIT([MKL], [-I], [-L])
AM_SUBST_NOTMAKE([MKL_CFLAGS])
AM_SUBST_NOTMAKE([MKL_CLIBS])
ACX_LANG_PACKAGE_INIT([LAPACKE], [-I], [-L])
AM_SUBST_NOTMAKE([LAPACKE_CFLAGS])
AM_SUBST_NOTMAKE([LAPACKE_CLIBS])
ACX_LANG_PACKAGE_INIT([ATLAS], [-I], [-L])
AM_SUBST_NOTMAKE([ATLAS_CFLAGS])
AM_SUBST_NOTMAKE([ATLAS_CLIBS])
ACX_LANG_PACKAGE_INIT([CLAPACK], [-I], [-L])
AM_SUBST_NOTMAKE([CLAPACK_CFLAGS])
AM_SUBST_NOTMAKE([CLAPACK_CLIBS])
ACX_LANG_PACKAGE_INIT([Fortran LAPACK], [], [-L])
AM_SUBST_NOTMAKE([FORTRAN_LAPACK_CLIBS])

AC_SUBST([LAPACK_CFLAGS])
AC_SUBST([LAPACK_CLIBS])
AC_SUBST([LAPACK_FCLIBS])

AS_CASE([$with_external_lapack],
  [yes|auto], [external_lapack_options='mkl lapacke atlas clapack fortran'],
  [no], [external_lapack_options=],
  [external_lapack_options=$with_external_lapack])

lapack_interface_found=no
yac_save_CFLAGS=$CFLAGS
for option in $external_lapack_options; do
  AS_CASE([$option],

    [mkl],
    [AS_UNSET([acx_cv_c_lib_func_LAPACKE_dgels_work])
     MKL_CFLAGS=${MKL_CFLAGS-$acx_c_MKL_inc_search_args}
     CFLAGS="$MKL_CFLAGS $CFLAGS"
     AC_CHECK_HEADER([mkl_lapacke.h],
       [ACX_LANG_LIB_SEARCH([MKL_CLIBS], [LAPACKE_dgels_work],
          [-mkl],
          [MKL_CLIBS=$acx_cv_c_lib_func_LAPACKE_dgels_work
           LAPACK_CFLAGS=$MKL_CFLAGS
           LAPACK_CLIBS=$MKL_CLIBS
           AC_DEFINE([YAC_LAPACK_INTERFACE_ID], [1])
           AC_MSG_NOTICE([using Intel MKL LAPACKE interface])
           lapack_interface_found=yes], [:])], [:], [//])],

    [lapacke],
    [AS_UNSET([acx_cv_c_lib_func_LAPACKE_dgels_work])
     LAPACKE_CFLAGS=${LAPACKE_CFLAGS-$acx_c_LAPACKE_inc_search_args}
     CFLAGS="$LAPACKE_CFLAGS $CFLAGS"
     AC_CHECK_HEADER([lapacke.h],
       [ACX_LANG_LIB_SEARCH([LAPACKE_CLIBS], [LAPACKE_dgels_work],
          [ASX_PREPEND_LDFLAGS([$acx_c_LAPACKE_lib_search_args],
             [-llapacke], [-lopenblas], [-llapack -pgf77libs],
             [-llapack -lblas -pgf77libs])],
          [LAPACKE_CLIBS=$acx_cv_c_lib_func_LAPACKE_dgels_work
           LAPACK_CFLAGS=$LAPACKE_CFLAGS
           LAPACK_CLIBS=$LAPACKE_CLIBS
           AC_DEFINE([YAC_LAPACK_INTERFACE_ID], [2])
           AC_MSG_NOTICE([using Netlib LAPACKE interface])
           lapack_interface_found=yes], [:])], [:], [//])],

    [atlas],
    [AS_UNSET([ac_cv_header_clapack_h])
     ATLAS_CFLAGS=${ATLAS_CFLAGS-$acx_c_ATLAS_inc_search_args}
     CFLAGS="$ATLAS_CFLAGS $CFLAGS"
     AC_CHECK_HEADER([clapack.h],
       [ACX_LANG_LIB_SEARCH([ATLAS_CLIBS], [clapack_dgels],
          [ASX_PREPEND_LDFLAGS([$acx_c_ATLAS_lib_search_args],
             [-lsatlas], [-ltatlas], [-latlas], [-llapack],
             [-llapack -lcblas -lf77blas -latlas])],
          [ATLAS_CLIBS=$acx_cv_c_lib_func_clapack_dgels
           LAPACK_CFLAGS=$ATLAS_CFLAGS
           LAPACK_CLIBS=$ATLAS_CLIBS
           AC_DEFINE([YAC_LAPACK_INTERFACE_ID], [3])
           AC_MSG_NOTICE([using ATLAS CLAPACK interface])
           lapack_interface_found=yes], [:])], [:], [//])],

    [clapack],
    [AS_UNSET([ac_cv_header_clapack_h])
     AS_UNSET([acx_cv_c_lib_func_dgels_])
     CLAPACK_CFLAGS=${CLAPACK_CFLAGS-$acx_c_CLAPACK_inc_search_args}
     CFLAGS="$CLAPACK_CFLAGS $CFLAGS"
     AC_CHECK_HEADER([clapack.h],
       [ACX_LANG_LIB_SEARCH([CLAPACK_CLIBS], [dgels_],
          [],
          [CLAPACK_CLIBS=$acx_cv_c_lib_func_dgels_
           LAPACK_CFLAGS=$CLAPACK_CFLAGS
           LAPACK_CLIBS=$CLAPACK_CLIBS
           AC_DEFINE([YAC_LAPACK_INTERFACE_ID], [4])
           AC_MSG_NOTICE([using Netlib CLAPACK interface])
           lapack_interface_found=yes], [:])], [:], [[#include <f2c.h>]])],

    [fortran],
    [AS_UNSET([acx_cv_c_lib_func_dgels_])
     ACX_FC_MANGLING_GLOBAL
     ACX_FC_MANGLING_SHVAR([dgels],
       [ACX_LANG_LIB_SEARCH([FORTRAN_LAPACK_CLIBS],
          [$acx_cv_fc_mangling_shvar_DGELS],
          [ASX_PREPEND_LDFLAGS([$acx_c_Fortran_LAPACK_lib_search_args],
             [-llapack], [-lopenblas], [-llapack -pgf77libs],
             [-llapack -lblas -pgf77libs])],
          [ACX_FC_MANGLING_DEFINE([YAC_FC],
             [LAPACK_CFLAGS=
              AS_VAR_COPY([LAPACK_CLIBS],
                [acx_cv_c_lib_func_$acx_cv_fc_mangling_shvar_DGELS])
              AC_DEFINE([YAC_LAPACK_INTERFACE_ID], [5])
              AC_MSG_NOTICE([using Fortran LAPACK interface])
              lapack_interface_found=yes], [:])], [:])], [:])],

    [yes|no|auto],
    [AC_MSG_ERROR(['yes', 'no' and 'auto' cannot be elements of the list dnl
provided as a value of the argument --with-external-lapack; valid elements dnl
are 'mkl', 'lapacke', 'atlas', 'fortran', and 'clapack'])],

    [AC_MSG_ERROR([unexpected value of the argument dnl
--with-external-lapack='$option'; valid values are either 'yes', 'no', dnl
'auto' or a space-separated list of the following: 'mkl', 'lapacke', dnl
'atlas', 'fortran', 'clapack'])])

  CFLAGS=$yac_save_CFLAGS
  AS_VAR_IF([lapack_interface_found], [yes], [break])
done

AS_VAR_IF([lapack_interface_found], [yes],
  [AS_IF([test -n "$LAPACK_CLIBS"],
     [YAC_PKGCONF_CLIBS="$LAPACK_CLIBS $YAC_PKGCONF_CLIBS"])],
  [AS_CASE([$with_external_lapack],

     [yes],
     [AC_MSG_FAILURE([none of the supported LAPACK interfaces is available])],

     [no|auto],
     [LAPACK_CFLAGS="-I$yac_abs_top_srcdir/clapack/INCLUDE"
      LAPACK_CLIBS="$yac_abs_top_builddir/clapack/libyac_clapack.a"
      YAC_PKGCONF_CLIBS="-lyac_clapack $YAC_PKGCONF_CLIBS"
      AC_DEFINE([YAC_LAPACK_INTERFACE_ID], [4])
      AC_MSG_NOTICE([using bundled Netlib CLAPACK])
      with_external_lapack=no],

     [AC_MSG_FAILURE(
        [none of the requested LAPACK interfaces is available])])])
AM_CONDITIONAL([EXTERNAL_LAPACK], [test x"$with_external_lapack" != xno])

dnl Checks for Python bindings.

AC_ARG_ENABLE([python-bindings],
  [AS_HELP_STRING([--enable-python-bindings],
     [enable Python bindings @<:@default=no@:>@])], [],
  [enable_python_bindings=no])
AM_CONDITIONAL([ENABLE_PYTHON_BINDINGS],
  [test x"$enable_python_bindings" = xyes])
yac_have_py_runtime_deps=$enable_python_bindings

AS_VAR_IF([enable_python_bindings], [yes],
  [ACX_SHLIB_CC_PIC_FLAG
   AC_SUBST([PIC_CFLAGS], [$acx_cv_c_pic_flag])
   ACX_SHLIB_PIC_FLAGS_CHECK([$PIC_CFLAGS], [], [PIC_CFLAGS=])
dnl
   AC_ARG_VAR([PYTHON], [Python interpreter command])dnl
   ACX_PROG_SEARCH([PYTHON], [python3 python],
     [yac_python_try=dnl
"$acx_candidate -c 'import sys; exit(sys.version_info@<:@0:2@:>@ < (3, 6))' dnl
>&AS_MESSAGE_LOG_FD"
      _AC_DO_VAR([yac_python_try])],
     [PYTHON=$acx_cv_prog_PYTHON],
     [AC_MSG_FAILURE([unable to find Python interpreter with version >= 3.6])])
dnl
   yac_have_py_runtime_deps=yes
   set 'setuptools' 'from setuptools import setup'       'build'     \
       'pip'        'import pip'                         'build'     \
       'cython'     'from Cython.Build import cythonize' 'build'     \
       'numpy'      'import numpy'                       'build'     \
       'mpi4py'     'import mpi4py'                      'run-extra' \
       'netCDF4'    'import netCDF4'                     'test'      \
       'matplotlib' 'import matplotlib.pyplot as plt'    'test-extra'
   while test $[]@%:@ -gt 1; do
     m4_pushdef([yac_cache_var], [yac_cv_has_pymodule_[]AS_TR_SH([$][1])])dnl
     AC_CACHE_CHECK([whether Python module '$[]1' is available],
       [yac_cache_var],
       [AS_VAR_SET([yac_cache_var], [no])
        yac_pymodule_try="$PYTHON -c '$[]2' >&AS_MESSAGE_LOG_FD"
        _AC_DO_VAR([yac_pymodule_try])
        AS_IF([test $? -eq 0], [AS_VAR_SET([yac_cache_var], [yes])])])
     AS_VAR_IF([yac_cache_var], [no],
       [AS_CASE([$[]3],
          [build],
          [AC_MSG_FAILURE(
             [Python module '$[]1' needed to build the Python bindings is dnl
not available])],
          [run],
          [AC_MSG_WARN([Python module '$[]1' needed to run the Python dnl
bindings is not available, all Python tests will be skipped])
           yac_have_py_runtime_deps=no
           break],
          [test],
          [AC_MSG_WARN([Python module '$[]1' needed to test certain features dnl
of the Python bindings is not available, the respective tests will be dnl
skipped])])])
     shift; shift; shift
     m4_popdef([yac_cache_var])dnl
   done
   set dummy; shift])
AM_CONDITIONAL([TEST_PYTHON_BINDINGS],
  [test x"$yac_have_py_runtime_deps" = xyes])

dnl Extract all paths provided with -L flag and prepend them to the
dnl LD_LIBRARY_PATH (or its analogue) environment variable when running tests:
dnl   variable for the shell scripts
AC_SUBST([SHLIB_ENVIRONMENT_SHELL], [''])
AM_SUBST_NOTMAKE([SHLIB_ENVIRONMENT_SHELL])
dnl   variable for the makefiles (see also Automake docs for the description of
dnl   the TESTS_ENVIRONMENT variable)
AC_SUBST([SHLIB_ENVIRONMENT_MAKE], [''])
AS_CASE([$host_os],
  [darwin*],
  [dnl
dnl We should not do this for MacOS for the following reasons:
dnl   1) the linker puts all the information required for locating the dynamic
dnl      libraries into the executables by default;
dnl   2) setting DYLD_LIBRARY_PATH (even if the System Integrity Protection is
dnl      enabled) might break the dynamic linking due to incompatibilities of
dnl      the system libraries and the libraries provided with MacPorts (e.g.
dnl      see https://trac.macports.org/ticket/57821).
  ],
  [yac_libdirs=
   ASX_EXTRACT_ARGS([yac_L_args],
     ["$LDFLAGS $YAC_PKGCONF_CLIBS $MPI_FCLIBS $LIBS"], ['-L@<:@ @:>@*'])
   for yac_L_arg in $yac_L_args; do
     if test -d "$yac_L_arg"; then
       ASX_VAR_APPEND_UNIQ([yac_libdirs], ["$yac_L_arg"], ["$PATH_SEPARATOR"])
     fi
   done
   AS_IF([test -n "$yac_libdirs"],
     [ACX_SHLIB_PATH_VAR
      AS_VAR_SET([SHLIB_ENVIRONMENT_SHELL],
        ["$acx_cv_shlib_path_var=dnl
\"$yac_libdirs:\$$acx_cv_shlib_path_var\"; export $acx_cv_shlib_path_var;"])
      AS_VAR_SET([SHLIB_ENVIRONMENT_MAKE],
        [`echo "$SHLIB_ENVIRONMENT_SHELL" | sed 's/\\$/$$/g'`])])])

AS_IF([test -n "$LDFLAGS"], [YAC_PKGCONF_CLIBS="$LDFLAGS $YAC_PKGCONF_CLIBS"])
AS_IF([test -n "$LIBS"], [AS_VAR_APPEND([YAC_PKGCONF_CLIBS], [" $LIBS"])])

AS_IF([test -n "$YAC_PKGCONF_CLIBS"],
  [YAC_PKGCONF_CLIBS="-lyac $YAC_PKGCONF_CLIBS"],
  [YAC_PKGCONF_CLIBS='-lyac'])

AC_CONFIG_FILES([
        examples/OASIS3_MCT_example/OASIS3_MCT_example.sh
        examples/toy_coupling/toy_coupling.sh
        examples/toy_icon_runoff/toy_icon_runoff.sh
        examples/toy_multi/toy_multi.sh
        examples/perf_toy/perf_toy.sh
        examples/toy_atm_ocn/toy_atm_ocn.sh
        examples/toy_dummy/dummy_c.sh
        examples/toy_dummy/dummy.sh
        examples/toy_callback/toy_callback.sh
        examples/dynamic_config/dynamic_config.sh
        examples/interpolation_benchmark/toy_scrip_GEN_MASK.sh
        examples/interpolation_benchmark/toy_scrip_ICOS_ICOH.sh
        examples/interpolation_benchmark/toy_scrip_NOGT_ICOH.sh
        examples/interpolation_benchmark/toy_scrip_NOICOH.sh
        python/tests/seq_plot_noise_test.sh
        python/tests/par_plot_noise_test.sh
        python/tests/seq_noise2nc_test.sh
        python/tests/par_noise2nc_test.sh
        python/tests/seq_nc2plot_test.sh
        python/tests/par_nc2plot_test.sh
        python/tests/multithreading.sh
        python/tests/test_interface.sh
        tests/test_abort.sh
        tests/test_component_config.sh
        tests/test_couple_config.sh
        tests/test_def_comps.sh
        tests/test_def_datetime.sh
        tests/test_def_grid.sh
        tests/test_def_mask.sh
        tests/test_def_points.sh
        tests/test_dist_grid_pair_parallel.sh
        tests/test_dummy_coupling.sh
        tests/test_dummy_coupling2_c.sh
        tests/test_dummy_coupling3.sh
        tests/test_dummy_coupling3_c.sh
        tests/test_dummy_coupling4_c.sh
        tests/test_dummy_coupling5.sh
        tests/test_dummy_coupling6.sh
        tests/test_dummy_coupling6_c.sh
        tests/test_dummy_coupling7.sh
        tests/test_dummy_coupling7_c.sh
        tests/test_dummy_coupling8_c.sh
        tests/test_dummy_coupling9.sh
        tests/test_group_comm.sh
        tests/test_init_comm_final.sh
        tests/test_init_final.sh
        tests/test_mpi_handshake.sh
        tests/test_mpi_handshake_c.sh
        tests/test_instance_parallel1.sh
        tests/test_instance_parallel2.sh
        tests/test_instance_parallel3.sh
        tests/test_instance_parallel4.sh
        tests/test_interp_grid_parallel.sh
        tests/test_interp_method_avg_parallel.sh
        tests/test_interp_method_callback_parallel.sh
        tests/test_interp_method_conserv_parallel.sh
        tests/test_interp_method_creep_parallel.sh
        tests/test_interp_method_file_parallel.sh
        tests/test_interp_method_hcsbb_parallel.sh
        tests/test_interp_method_nnn_parallel.sh
        tests/test_interp_method_nnn_parallel2.sh
        tests/test_interp_method_parallel.sh
        tests/test_interp_method_rbf_parallel.sh
        tests/test_interp_method_spmap_parallel.sh
        tests/test_interp_weights_parallel.sh
        tests/test_interpolation_exchange.sh
        tests/test_interpolation_parallel1.sh
        tests/test_interpolation_parallel2.sh
        tests/test_interpolation_parallel3.sh
        tests/test_interpolation_parallel4.sh
        tests/test_interpolation_parallel5.sh
        tests/test_interpolation_parallel6.sh
        tests/test_mpi_error.sh
        tests/test_proc_sphere_part_parallel.sh
        tests/test_read_fesom.sh
        tests/test_read_icon.sh
        tests/test_read_icon_parallel.sh
        tests/test_read_mpiom.sh
        tests/test_read_scrip.sh
        tests/test_read_woa_data.sh
        tests/test_redirstdout.sh
        tests/test_restart.sh
        tests/test_restart2.sh
        tests/test_version.sh
        tests/test_weights2vtk.sh
        tests/test_dynamic_config.sh
        tests/test_query_routines.sh
        tests/test_multithreading.sh
        tests/test_io_config.sh
],
[chmod a+x "$ac_file"])

AC_CONFIG_FILES([
        Makefile
        clapack/Makefile
        contrib/Makefile
        examples/Makefile
        mtime/Makefile
        python/Makefile
        python/setup.py
        src/Makefile
        src/pkgconfig/yac.pc
        tests/Makefile
])

AC_CONFIG_LINKS(m4_foreach_w([file],
        couple_config_test.yaml
        coupling_test.yaml
        coupling_test2.yaml
        coupling_test3_comp_a_comp_b.yaml
        coupling_test3_comp_a_comp_c.yaml
        coupling_test3_comp_a_intra.yaml
        coupling_test3_comp_b_comp_c.yaml
        coupling_test3_comp_b_intra.yaml
        coupling_test3_comp_c_intra.yaml
        coupling_test3_default.yaml
        coupling_test4.yaml
        coupling_test5.yaml
        coupling_test6.json
        coupling_test6_local.yaml
        coupling_test9.yaml
        cube_10x10.csv
        instance_test_1_1.yaml
        instance_test_1_10.yaml
        instance_test_1_2.yaml
        instance_test_1_3.yaml
        instance_test_1_4.yaml
        instance_test_1_5.yaml
        instance_test_1_6.yaml
        instance_test_1_7.yaml
        instance_test_1_8.yaml
        instance_test_1_9.yaml
        test_def_datetime.yaml
, [tests/file:input/file ]))

AC_CONFIG_LINKS([
        examples/GR30_lsm.nc:grids/GR30_lsm.nc
        examples/icon_grid_R02B03_G.nc:grids/icon_grid_R02B03_G.nc
        examples/icon_grid_R02B03_O.nc:grids/icon_grid_R02B03_O.nc
        python/yac.pyx:python/yac.pyx
        tests/GR30_lsm.nc:grids/GR30_lsm.nc
        tests/icon_grid_R02B00.nc:grids/icon_grid_R02B00.nc
        tests/icon_grid_R02B02.nc:grids/icon_grid_R02B02.nc
        tests/icon_grid_R02B03_G.nc:grids/icon_grid_R02B03_G.nc
])

AC_OUTPUT
