AC_PREREQ([2.69])
AC_INIT([ecrad], [1.1.0], [], [],
  [https://confluence.ecmwf.int/display/ECRAD])
AC_CONFIG_MACRO_DIR([m4])
m4_pattern_forbid([^A[SC]X_])dnl
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([radiation/radiation_flux.F90])

AC_ARG_ENABLE([silent-rules],
  [AC_HELP_STRING([--enable-silent-rules],
     [less verbose build output (undo: "make V=1") @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_silent_rules=yes])],
  [enable_silent_rules=yes])

AC_ARG_ENABLE([single-precision],
  [AC_HELP_STRING([--enable-single-precision],
     [switch to single precision @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_single_precision=yes])],
  [enable_single_precision=no])

AC_ARG_ENABLE([pgi-inlib],
  [AC_HELP_STRING([--enable-pgi-inlib],
     [enable PGI cross-file function inlining via an inline library (undo:
"make INLIB=0") @<:@default=auto@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_pgi_inlib=yes])],
  [enable_pgi_inlib=auto])

AC_PROG_FC
AC_LANG([Fortran])
ACX_FC_PP_SRCEXT([F90])
ACX_COMPILER_FC_VENDOR_SIMPLE
AC_SUBST([FC_VENDOR], [$acx_cv_fc_compiler_vendor])

ACX_FC_LINE_LENGTH([unlimited])

ACX_FC_MODULE_NAMING(
  [AC_SUBST([FCMODEXT], ["$acx_cv_fc_module_naming_ext"])
   AC_SUBST([FCMODUC], ["$acx_cv_fc_module_naming_upper"])])

ACX_FC_MODULE_IN_FLAG([AC_SUBST([FCMODINC], ["$acx_cv_fc_module_in_flag"])])
ACX_FC_MODULE_OUT_FLAG([AC_SUBST([FCMODOUT], ["$acx_cv_fc_module_out_flag"])])

ACX_FC_INCLUDE_FLAG([AC_SUBST([FCINCFLAG], ["$acx_cv_fc_ftn_include_flag"])])
ACX_FC_INCLUDE_ORDER(
  [AC_SUBST([FCINCORDER], ["$acx_cv_fc_ftn_include_order"])],
  [AC_SUBST([FCINCORDER], ['src,flg'])])

ACX_FC_INCLUDE_FLAG_PP(
  [AC_SUBST([FCINCFLAG_PP], ["$acx_cv_fc_pp_include_flag"])])
ACX_FC_INCLUDE_ORDER_PP(
  [AC_SUBST([FCINCORDER_PP], ["$acx_cv_fc_pp_include_order"])],
  [AC_SUBST([FCINCORDER_PP], ['inc,flg'])])

ACX_LANG_MACRO_FLAG([AC_SUBST([FCDEF_PP], ["$acx_cv_fc_macro_flag"])])

ACX_FC_ENDIANNESS_REAL(
dnl Kind of the double-precision real variables that are read from the binary
dnl files (JPRD from ifsaux/parkind1.F90):
  [SELECTED_REAL_KIND(13,300)],
  [AS_CASE([$acx_cv_fc_endianness_real],
dnl Big-endian:
     [replEhkm], [],
dnl Little-endian:
     [mkhElper],
     [AS_VAR_APPEND([FCFLAGS],
        [" ${acx_cv_fc_macro_flag}__ECRAD_LITTLE_ENDIAN"])],
dnl Unsupported mixed-endian:
     [AC_MSG_WARN([the endianness of the target system is not supported])])],
  [AC_MSG_WARN([unable to detect the endianness of the target system])])

AS_VAR_IF([enable_silent_rules], [yes],
  [AC_SUBST([DEFAULT_VERBOSITY], [0])],
  [AC_SUBST([DEFAULT_VERBOSITY], [1])])

AS_VAR_IF([enable_single_precision], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${acx_cv_fc_macro_flag}SINGLE_PRECISION"])])

AS_CASE([$enable_pgi_inlib],
  [auto],
  [AS_VAR_IF([acx_cv_fc_compiler_vendor], [portland],
     [AC_MSG_NOTICE([PGI cross-file function inlining via an inline library dnl
is enabled because $FC is recognized as the PGI compiler: disable the dnl
inlining if required (--disable-pgi-inlib)])
      enable_pgi_inlib=yes],
     [enable_pgi_inlib=no])],
  [yes],
  [AS_IF([test x"$acx_cv_fc_compiler_vendor" != xportland],
     [AC_MSG_WARN([PGI cross-file function inlining via an inline library is dnl
enabled but $FC is not recognized as the PGI compiler: it is strongly dnl
recommended to disable the compiler-specific feature dnl
(--disable-pgi-inlib)])])])

AS_VAR_IF([enable_pgi_inlib], [yes],
  [AC_SUBST([DEFAULT_INLIB], [1])],
  [AC_SUBST([DEFAULT_INLIB], [0])])

ACX_LANG_PACKAGE_INIT([NetCDF], ["$FCMODINC"], [-L])
AC_SUBST([DEPGEN_EXTERNAL_MODS], [netcdf])

dnl Check for NetCDF
acx_save_FCFLAGS=$FCFLAGS
NETCDF_FCFLAGS=${NETCDF_FCFLAGS-$acx_fc_NetCDF_inc_search_args}
FCFLAGS="$NETCDF_FCFLAGS $acx_save_FCFLAGS"
ACX_FC_MODULE_CHECK([netcdf],
  [ACX_LANG_LIB_SEARCH([NETCDF_FCLIBS], [nf90_open],
     [ASX_PREPEND_LDFLAGS([$acx_fc_NetCDF_lib_search_args],
        [-lnetcdff], [-lnetcdff -lnetcdf])],
     [NETCDF_FCLIBS=$acx_cv_fc_lib_func_NF90_OPEN
      acx_have_netcdf=yes], [],
     [AC_LANG_PROGRAM([],
[[      use netcdf, only: nf90_open, NF90_NOWRITE
      implicit none
      integer :: status, ncid
      status = nf90_open('conftest.nc', NF90_NOWRITE, ncid)]])])], [])
FCFLAGS=$acx_save_FCFLAGS

dnl Check for AR and ARFLAGS
AC_ARG_VAR([AR], [archiver command])
AC_ARG_VAR([ARFLAGS], [archiver flags])
ARFLAGS=${ARFLAGS-cr}
ACX_PROG_SEARCH([AR], [ar],
  [AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
     [acx_ar_try="$acx_candidate ${ARFLAGS} libconftest.a dnl
conftest.$ac_objext >&AS_MESSAGE_LOG_FD"
      _AC_DO_VAR([acx_ar_try])
      rm -f libconftest.a])
   test "$ac_status" -eq 0 >/dev/null 2>&1],
  [AR=$acx_cv_prog_AR])

dnl Check for PYTHON
AC_ARG_VAR([PYTHON], [Python interpreter command])dnl
ACX_PROG_SEARCH([PYTHON], [python python3],
  [acx_python_try="$acx_candidate $srcdir/mkhelper/depgen.py -h dnl
>&AS_MESSAGE_LOG_FD"
   _AC_DO_VAR([acx_python_try])],
  [PYTHON=$acx_cv_prog_PYTHON])

AC_PROG_INSTALL
AC_PROG_MKDIR_P

AC_CONFIG_FILES([Makefile])

dnl An additional step to generate dummy dependencies:
AC_CONFIG_COMMANDS([depfiles],
  [${MAKE-make} dummy-depend >/dev/null 2>&AS_MESSAGE_LOG_FD
   sleep 1
   touch Makefile])

AC_OUTPUT

