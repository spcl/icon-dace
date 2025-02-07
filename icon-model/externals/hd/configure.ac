AC_PREREQ([2.69])
AC_INIT([HD], [unknown])
AC_CONFIG_MACRO_DIR([build-aux/m4])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR([code/src/hd_driver.f90])
dnl Make sure all required m4 files are in place:
m4_pattern_forbid([^(A[SC]X|PKG)_])dnl

AC_PROG_FC
AC_LANG([Fortran])
ACX_FC_PP_SRCEXT([f90])
ACX_COMPILER_FC_VENDOR_SIMPLE
AC_SUBST([FC_VENDOR], [$acx_cv_fc_compiler_vendor])

ACX_FC_LINE_LENGTH([unlimited])

ACX_FC_MODULE_NAMING(
  [AC_SUBST([FC_MOD_FILE_UPPER], ["$acx_cv_fc_module_naming_upper"])
   AC_SUBST([FC_MOD_FILE_EXT], ["$acx_cv_fc_module_naming_ext"])])

ACX_FC_MODULE_IN_FLAG([AC_SUBST([FC_MOD_IN], ["$acx_cv_fc_module_in_flag"])])
ACX_FC_MODULE_OUT_FLAG([AC_SUBST([FC_MOD_OUT], ["$acx_cv_fc_module_out_flag"])])

ACX_FC_INCLUDE_FLAG([AC_SUBST([FC_INC_FLAG], ["$acx_cv_fc_ftn_include_flag"])])
ACX_FC_INCLUDE_ORDER(
  [AC_SUBST([FC_INC_ORDER], ["$acx_cv_fc_ftn_include_order"])],
  [AC_SUBST([FC_INC_ORDER], ['src,flg'])])

ACX_FC_INCLUDE_FLAG_PP(
  [AC_SUBST([FC_INC_FLAG_PP_f90], ["$acx_cv_fc_pp_include_flag"])])
ACX_FC_INCLUDE_ORDER_PP(
  [AC_SUBST([FC_INC_ORDER_PP_f90], ["$acx_cv_fc_pp_include_order"])],
  [AC_SUBST([FC_INC_ORDER_PP_f90], ['inc,flg'])])

ACX_LANG_MACRO_FLAG([AC_SUBST([FC_PP_DEF], ["$acx_cv_fc_macro_flag"])])

AC_ARG_VAR([PYTHON], [Python interpreter command])dnl
ACX_PROG_SEARCH([PYTHON], [python python3],
  [acx_python_try="$acx_candidate $srcdir/build-aux/mkhelper/depgen.py -h dnl
>&AS_MESSAGE_LOG_FD"
   _AC_DO_VAR([acx_python_try])],
  [PYTHON=$acx_cv_prog_PYTHON])

AC_PROG_INSTALL
AC_PROG_MKDIR_P

AC_ARG_ENABLE([rpaths],
  [AS_HELP_STRING([--enable-rpaths],
     [add directories specified with -L flags in LDFLAGS and LIBS to the
runtime library search paths (RPATH) @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_rpaths=yes])],
  [enable_rpaths=yes])

acx_lib_search_flags='-L'
AS_VAR_IF([enable_rpaths], [yes],
  [ACX_SHLIB_FC_RPATH_FLAG
   ASX_EXTRACT_ARGS([acx_L_args], ["$LDFLAGS $LIBS"], ['-L@<:@ @:>@*'])
   for acx_L_arg in $acx_L_args; do
     ASX_VAR_APPEND_UNIQ([LDFLAGS], ["$acx_cv_fc_rpath_flag$acx_L_arg"], [' '])
   done
   ACX_SHLIB_RPATH_FLAGS_CHECK
   AS_VAR_APPEND([acx_lib_search_flags], [" '$acx_cv_fc_rpath_flag'"])])

AC_ARG_ENABLE([silent-rules],
  [AS_HELP_STRING([--enable-silent-rules],
     [less verbose build output (undo: "make V=1") @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_silent_rules=yes])],
  [enable_silent_rules=yes])

AS_VAR_IF([enable_silent_rules], [yes],
  [AC_SUBST([DEFAULT_VERBOSITY], [0])],
  [AC_SUBST([DEFAULT_VERBOSITY], [1])])

AC_ARG_ENABLE([mpi],
  [AS_HELP_STRING([--enable-mpi],
     [enable MPI (parallelization) support @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_mpi=yes])],
  [enable_mpi=yes])

AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}STANDALONE ${FC_PP_DEF}__STANDALONE"])

AC_ARG_VAR([MPI_FCFLAGS], [exact Fortran compiler flags enabling MPI])
AC_ARG_VAR([MPI_FCLIBS], [exact linker flags enabling MPI when linking with
Fortran compiler])
AC_SUBST([DEPGEN_EXTERNAL_MODS], [mpi])

AS_VAR_IF([enable_mpi], [yes],
  [AS_VAR_SET_IF([MPI_FCFLAGS], [FCFLAGS="$MPI_FCFLAGS $FCFLAGS"])
   ACX_FC_MODULE_CHECK([MPI],
     [ACX_LANG_LIB_SEARCH([MPI_FCLIBS], [mpi_init], [],
        [MPI_FCLIBS=$acx_cv_fc_lib_func_MPI_INIT
         LIBS="$MPI_FCLIBS $LIBS"], [],
        [AC_LANG_PROGRAM([],[[      use MPI
      implicit none
      integer :: ierror
      call MPI_INIT(ierror)]])])])],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NOMPI"])])

AC_ARG_ENABLE([resolution],
[  --enable-resolution     build the model for the specific domain resolution.
                          The value of the argument must be one of the
                          following @<:@default=05deg@:>@:
                            0|05deg     0.5 degree global domain
                            1|5min      5min global domain
                            2|euro5min  5min European domain
                            3|aus5min   5min Australian domain],
[AS_CASE([$enableval],
   [0|05deg], [enable_resolution=05deg],
   [1|5min], [enable_resolution=5min],
   [2|euro5min], [enable_resolution=euro5min],
   [3|aus5min], [enable_resolution=aus5min],
   [AC_MSG_ERROR([unexpected value for the argument dnl
--enable-resolution='$enableval'; see the help message for the valid values dnl
of the argument])])],
[enable_resolution=05deg])

AC_SUBST([EXEC_FILE_NAME])
AS_CASE([$enable_resolution],
  [05deg],
  [EXEC_FILE_NAME='hd_05.exe'],
  [5min],
  [EXEC_FILE_NAME='hd_5min.exe'],
  [euro5min],
  [EXEC_FILE_NAME='hd_5min.exe'],
  [aus5min],
  [EXEC_FILE_NAME='hd_5min.exe'],
  [AC_MSG_ERROR([unexpected value for the domain resolution dnl
'$enable_resolution'; this is a bug in the configure script])])
AS_CASE([$enable_resolution],
  [*5min], [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HD_5MIN"])])

eval "set dummy $acx_lib_search_flags; shift"
ACX_LANG_PACKAGE_INIT([NetCDF], ["$FC_INC_FLAG" "$FC_MOD_IN"], ["$@"])
AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',netcdf'])
NETCDF_FCFLAGS=${NETCDF_FCFLAGS-$acx_fc_NetCDF_inc_search_args}
FCFLAGS="$NETCDF_FCFLAGS $FCFLAGS"
ACX_FC_INCLUDE_CHECK([netcdf.inc],
  [ACX_LANG_LIB_SEARCH([NETCDF_FCLIBS], [nf_open],
     [ASX_PREPEND_LDFLAGS([$acx_fc_NetCDF_lib_search_args],
        [-lnetcdff], [-lnetcdff -lnetcdf])],
     [ACX_FC_MODULE_CHECK([NETCDF],
        [NETCDF_FCLIBS=$acx_cv_fc_lib_func_NF_OPEN
         LIBS="$NETCDF_FCLIBS $LIBS"])])])

AC_ARG_ENABLE([coupling],
[  --enable-coupling       enable coupling using the specific backend
                          @<:@default=no@:>@:
                            yes|yac  use YAC coupler
                            oasis    use OASIS coupler (not yet supported)
                            no       disable coupling],
[AS_CASE([$enableval],
   [yes], [enable_coupling=yac],
   [yac|oasis|no], [],
   [AC_MSG_ERROR([unexpected value for the argument dnl
--enable-coupling='$enableval'; see the help message for the valid values dnl
of the argument])])],
[enable_coupling=no])
AS_CASE([$enable_coupling],
  [yac],
  [ACX_LANG_PACKAGE_INIT([YAC], ["$FC_MOD_IN"], [-L])
   AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',mo_yac_finterface'])
   YAC_FCFLAGS=${YAC_FCFLAGS-$acx_fc_YAC_inc_search_args}
   FCFLAGS="$YAC_FCFLAGS $FCFLAGS"
   ACX_FC_MODULE_CHECK([MO_YAC_FINTERFACE],
     [dnl
dnl If YAC_FCLIBS is not set, try getting linker flags from yac.pc. The problem
dnl is that the flags returned by pkg-config are meant for the C compiler and
dnl might not be compatible with the Fortran compiler. Therefore, we cannot
dnl rely on the result blindly and still need to check whether the flags work.
dnl The following variable will store the set of flags to undergo the checking:
      acx_pkg_YAC_FCLIBS=
      AS_VAR_SET_IF([YAC_FCLIBS], [],
        [PKG_PROG_PKG_CONFIG
dnl If the user has provided either --with-yac-root or --with-yac-lib, prepend
dnl the PKG_CONFIG_PATH environment variable:
         AS_VAR_SET_IF([with_yac_lib],
           [PKG_CONFIG_PATH="$with_yac_lib/pkgconfig:$PKG_CONFIG_PATH"
            export PKG_CONFIG_PATH])
dnl Prevent PKG_CHECK_MODULES from declaring package-specific precious
dnl variables, which names do not follow our convention:
         m4_pushdef([AC_ARG_VAR])dnl
dnl PKG_CHECK_MODULES relies on yac_CFLAGS and yac_LIBS, which in our case
dnl have names YAC_FCFAGS and YAC_FCLIBS, respectively. The value of the first
dnl one is already known and the second one we want to get:
         yac_CFLAGS=$YAC_FCFAGS
         yac_LIBS=
         PKG_CHECK_MODULES([yac], [yac],
           [acx_pkg_YAC_FCLIBS=$yac_LIBS
            AS_VAR_IF([enable_rpaths], [yes],
              [ASX_EXTRACT_ARGS([acx_L_args], ["$acx_pkg_YAC_FCLIBS"],
                 ['-L@<:@ @:>@*'])
               for acx_L_arg in $acx_L_args; do
                 ASX_VAR_APPEND_UNIQ([acx_pkg_YAC_FCLIBS],
                   ["$acx_cv_fc_rpath_flag$acx_L_arg"], [' '])
               done])], [:])
         m4_popdef([AC_ARG_VAR])])
     ACX_LANG_LIB_SEARCH([YAC_FCLIBS], [yac_fget_version],
        ["$acx_pkg_YAC_FCLIBS"],
        [YAC_FCFLAGS=$acx_cv_fc_lib_func_YAC_FGET_VERSION
         LIBS="$YAC_FCFLAGS $LIBS"], [],
        [AC_LANG_PROGRAM([],[[      use MO_YAC_FINTERFACE, ONLY: yac_fget_version
      implicit none
      print *, 'YAC Version: ' // TRIM(yac_fget_version())]])])])
   AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}COUP_YAC"])],
  [oasis],
  [eval "set dummy $acx_lib_search_flags; shift"
   ACX_LANG_PACKAGE_INIT([OASIS], ["$FC_MOD_IN"], ["$@"])
   AS_VAR_APPEND([DEPGEN_EXTERNAL_MODS], [',mod_oasis,mod_oasis_namcouple'])
   OASIS_FCFLAGS=${YAC_FCFLAGS-$acx_fc_OASIS_inc_search_args}
   FCFLAGS="$OASIS_FCFLAGS $FCFLAGS"
   ACX_FC_MODULE_CHECK([MOD_OASIS],
      [ACX_LANG_LIB_SEARCH([OASIS_FCLIBS], [oasis_init_comp],
         [ASX_PREPEND_LDFLAGS([$acx_fc_OASIS_lib_search_args],
            [-lpsmile.MPI1], [-lpsmile.MPI1 -lmct -lmpeu -lscrip])],
         [OASIS_FCFLAGS=$acx_cv_fc_lib_func_OASIS_INIT_COMP
          LIBS="$OASIS_FCFLAGS $LIBS"], [],
         [AC_LANG_PROGRAM([],[[      use MOD_OASIS, ONLY: oasis_init_comp
      implicit none
      integer :: ncomp_id, nerror
      call oasis_init_comp(ncomp_id, '', nerror)]])])])
   AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}COUP_OAS"])],
  [no], [],
  [AC_MSG_ERROR([unexpected value for the coupling backend
'$enable_coupling'; this is a bug in the configure script])])

AC_CONFIG_FILES([Makefile])

AC_CONFIG_COMMANDS([depfiles],
  [${MAKE-make} dummy-depend >/dev/null 2>&AS_MESSAGE_LOG_FD
   sleep 1
   touch Makefile])

AC_OUTPUT
