AC_DEFUN([ACX_CDI_OPTIONS],
[
#  ----------------------------------------------------------------------
#  Checks for multithreaded compiling + linking
AC_ARG_WITH([threads],
            [AS_HELP_STRING([--with-threads=<yes/no/directory>],
                            [Compile + link for multithreading [default=yes]])],
            [],
            [with_threads=yes])
THREADS_INCLUDE=''
THREADS_LIBS=''
AS_CASE([$with_threads],
        [no],[AC_MSG_CHECKING([multithreading])
              AC_MSG_RESULT([suppressed])],
        [yes],[AX_PTHREAD([AC_DEFINE([HAVE_LIBPTHREAD],[1],[Define 1 for multithread support])],[AC_MSG_ERROR([multithreaded settings NOT found])])
               LIBS="$PTHREAD_LIBS $LIBS"
               CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
               CC="$PTHREAD_CC"
               AS_ECHO(["CC:$CC CFLAGS:$CFLAGS LIBS:$LIBS"]) >&AS_MESSAGE_FD],
        [*],[THREADS_ROOT=$with_threads
             LDFLAGS="-L$THREADS_ROOT/lib $LDFLAGS"
             CPPFLAGS="-I$THREADS_ROOT/include $CPPFLAGS "
             AC_CHECK_HEADERS([pthread.h],,,[AC_INCLUDES_DEFAULT])
             AC_CHECK_LIB([pthread],[pthread_create])
             THREADS_LIBS=" -L$THREADS_ROOT/lib -lpthread"
             THREADS_INCLUDE=" -I$THREADS_ROOT/include"])
AC_SUBST([THREADS_INCLUDE])
AC_SUBST([THREADS_LIBS])
#  ----------------------------------------------------------------------
#  Compile application with FDB5 library
FDB5_INCLUDE=''
FDB5_LIBS=''
AC_ARG_WITH([fdb5],
            [AS_HELP_STRING([--with-fdb5=<yes|no|directory> (default=no)],[location of fdb5 library (optional)])],
            [AS_CASE(["$with_fdb5"],
                     [no],[AC_MSG_CHECKING([for fdb5 library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([fdb5/api/fdb_c.h],,
                                             [AC_MSG_ERROR([Could not find fdb5/api/fdb_c.h])],
                                             [AC_INCLUDES_DEFAULT])
                            AC_SEARCH_LIBS([fdb_initialise],
                                           [fdb5],
                                           [AC_DEFINE([HAVE_LIBFDB5],[1],[Define to 1 for FDB5 support])],
                                           [AC_MSG_ERROR([Could not link to fdb5])])
                            FDB5_LIBS=" -lfdb5"],
                     [*],[FDB5_ROOT=$with_fdb5
                          AS_IF([test -d "$FDB5_ROOT"],
                                [LDFLAGS="-L$FDB5_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$FDB5_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([fdb5/api/fdb_c.h],,
                                                  [AC_MSG_ERROR([Could not find fdb5/api/fdb_c.h])],
                                                  [AC_INCLUDES_DEFAULT])
                                 AC_SEARCH_LIBS([fdb_initialise],
                                                [fdb5],
                                                [AC_DEFINE([HAVE_LIBFDB5],[1],[Define to 1 for FDB5 support])],
                                                [AC_MSG_ERROR([Could not link to fdb5])])
                                 FDB5_LIBS=" -L$FDB5_ROOT/lib -lfdb5"
                                 FDB5_INCLUDE=" -I$FDB5_ROOT/include"],
                                [AC_MSG_NOTICE([$FDB5_ROOT is not a directory! FDB5 suppressed])])])],
            [AC_MSG_CHECKING([for fdb5 library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([FDB5_INCLUDE])
AC_SUBST([FDB5_LIBS])
#  ----------------------------------------------------------------------
#  Compile application with SZLIB library, needed for GRIB1
SZLIB_INCLUDE=''
SZLIB_LIBS=''
AC_ARG_WITH([szlib],
            [AS_HELP_STRING([--with-szlib=<yes|no|directory> (default=no)],[location of szlib library, optional for GRIB1 and NETCDF4 compression])],
            [AS_CASE(["$with_szlib"],
                     [no],[AC_MSG_CHECKING([for szlib library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([szlib.h],,
                                             [AC_MSG_ERROR([Could not find szlib.h])],
                                             [AC_INCLUDES_DEFAULT])
                            AC_SEARCH_LIBS([SZ_BufftoBuffCompress],
                                           [sz],
                                           [AC_DEFINE([HAVE_LIBSZ],[1],[Define to 1 for SZIP support])],
                                           [AC_MSG_ERROR([Could not link to szlib])])
                            SZLIB_LIBS=" -lsz"],
                     [*],[SZLIB_ROOT=$with_szlib
                          AS_IF([test -d "$SZLIB_ROOT"],
                                [LDFLAGS="-L$SZLIB_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$SZLIB_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([szlib.h],,
                                                  [AC_MSG_ERROR([Could not find szlib.h])],
                                                  [AC_INCLUDES_DEFAULT])
                                 AC_SEARCH_LIBS([SZ_BufftoBuffCompress],
                                                [sz],
                                                [AC_DEFINE([HAVE_LIBSZ],[1],[Define to 1 for SZIP support])],
                                                [AC_MSG_ERROR([Could not link to szlib])])
                                 SZLIB_LIBS=" -L$SZLIB_ROOT/lib -lsz"
                                 SZLIB_INCLUDE=" -I$SZLIB_ROOT/include"],
                                [AC_MSG_NOTICE([$SZLIB_ROOT is not a directory! SZLIB suppressed])])])],
            [AC_MSG_CHECKING([for szlib library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([SZLIB_INCLUDE])
AC_SUBST([SZLIB_LIBS])
#  ----------------------------------------------------------------------
#  Compile application with netcdf
NETCDF_ROOT=''
NETCDF_INCLUDE=''
NETCDF_LIBS=''
ENABLE_NETCDF=no
ENABLE_NC2=no
ENABLE_NC4=no
ENABLE_NC4HDF5=no
ENABLE_NC4SZLIB=no
AC_ARG_WITH([netcdf],
            [AS_HELP_STRING([--with-netcdf=<yes|no|directory> (default=no)],[location of NetCDF library (lib and include subdirs)])],
            [AS_CASE(["$with_netcdf"],
                     [no],[AC_MSG_CHECKING([for NetCDF library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([netcdf.h],,
                                             [AC_MSG_ERROR([Could not find netcdf.h])],
                                             [AC_INCLUDES_DEFAULT])
                            AC_SEARCH_LIBS([nc_open],
                                           [netcdf],
                                           [AC_DEFINE([HAVE_LIBNETCDF],[1],[Define to 1 for NetCDF support])
                                            ENABLE_NETCDF=yes],
                                           [AC_MSG_ERROR([Could not link to NetCDF library])])
                            NETCDF_LIBS=" -lnetcdf"
                            AC_CHECK_PROG(NC_CONFIG,nc-config,nc-config)],
                     [*],[AS_IF([test -d "$with_netcdf"],
                                [NETCDF_ROOT=$with_netcdf
                                 LDFLAGS="-L$NETCDF_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$NETCDF_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([netcdf.h],,
                                                  [AC_MSG_ERROR([Could not find netcdf.h])],
                                                  [AC_INCLUDES_DEFAULT])
                                 AC_SEARCH_LIBS([nc_open],
                                                [netcdf],
                                                [AC_DEFINE([HAVE_LIBNETCDF],[1],[Define to 1 for NetCDF support])
                                                 ENABLE_NETCDF=yes],
                                                [AC_MSG_ERROR([Could not link to NetCDF library])])
                                 NETCDF_LIBS=" -L$NETCDF_ROOT/lib -lnetcdf"
                                 NETCDF_INCLUDE=" -I$NETCDF_ROOT/include"
                                 AC_CHECK_PROG(NC_CONFIG,nc-config,[$NETCDF_ROOT/bin/nc-config],,["$NETCDF_ROOT/bin"])],
                                [AC_MSG_NOTICE([$with_netcdf is not a directory! NetCDF suppressed])])])],
            [AC_MSG_CHECKING([for NetCDF library])
             AC_MSG_RESULT([suppressed])])

AS_VAR_IF([ENABLE_NETCDF], [yes],
  [AS_VAR_IF([NC_CONFIG],,
             [AC_MSG_WARN([Could not find nc-config! go on with default configuration])])

   AC_CACHE_CHECK([netcdf's OpenDAP support],
                  [acx_cv_have_libnc_dap],
                  [acx_cv_have_libnc_dap=no
                   test "x$NC_CONFIG" != "x" && \
                   test "x$($NC_CONFIG --has-dap)" = "xyes" && \
                   acx_cv_have_libnc_dap=yes])
   AS_VAR_IF([acx_cv_have_libnc_dap], [yes],
             [AC_DEFINE([HAVE_LIBNC_DAP],[1],[Define to 1 for NetCDF OpenDAP])])

   AC_CACHE_CHECK([netcdf's Zarr support],
                  [acx_cv_have_nczarr],
                  [acx_cv_have_nczarr=no
                   test "x$NC_CONFIG" != "x" && \
                   test "x$($NC_CONFIG --has-nczarr)" = "xyes" && \
                   acx_cv_have_nczarr=yes])
   AS_VAR_IF([acx_cv_have_nczarr], [yes],
             [AC_DEFINE([HAVE_NCZARR],[1],[Define to 1 for NetCDF Zarr])])

   AC_CACHE_CHECK([netcdf's nc2 support],
                  [acx_cv_have_netcdf2],
                  [acx_cv_have_netcdf2=no
                   test "x$NC_CONFIG" != "x" && \
                   test "x$($NC_CONFIG --has-nc2)" = "xyes" && \
                   acx_cv_have_netcdf2=yes])
   AS_VAR_IF([acx_cv_have_netcdf2], [yes],
             [AC_DEFINE([HAVE_NETCDF2],[1],[Define to 1 for NetCDF2 support])
              ENABLE_NC2=yes])

   AC_CACHE_CHECK([netcdf's nc4 support],
                  [acx_cv_have_netcdf4],
                  [acx_cv_have_netcdf4=no
                   test "x$NC_CONFIG" != "x" && \
                   test "x$($NC_CONFIG --has-nc4)" = "xyes" && \
                   acx_cv_have_netcdf4=yes])
   AS_VAR_IF([acx_cv_have_netcdf4], [yes],
             [AC_DEFINE([HAVE_NETCDF4],[1],[Define to 1 for NetCDF4 support])
              ENABLE_NC4=yes])

   AC_CACHE_CHECK([netcdf's nc4/hdf5 support],
                  [acx_cv_have_nc4hdf5],
                  [acx_cv_have_nc4hdf5=no
                   test "x$NC_CONFIG" != "x" && \
                   test "x$($NC_CONFIG --has-hdf5)" = "xyes" && \
                   acx_cv_have_nc4hdf5=yes])
   AS_VAR_IF([acx_cv_have_nc4hdf5], [yes],
             [AC_DEFINE([HAVE_NC4HDF5],[1],[Define to 1 for NetCDF4/HDF5 support])
              ENABLE_NC4HDF5=yes])

   AC_CACHE_CHECK([netcdf's nc4/szlib support],
                  [acx_cv_have_nc4szlib],
                  [acx_cv_have_nc4szlib=no
                   test "x$NC_CONFIG" != "x" && \
                   test "x$($NC_CONFIG --has-szlib)" = "xyes" && \
                   acx_cv_have_nc4szlib=yes])
   AS_VAR_IF([acx_cv_have_nc4szlib], [yes],
             [AC_DEFINE([HAVE_NC4SZLIB],[1],[Define to 1 for NetCDF4/szlib support])
              ENABLE_NC4SZLIB=yes])])

AS_IF([test "x$ENABLE_NC4SZLIB" = "xyes"],
      [AC_SEARCH_LIBS([nc_def_var_szip], [netcdf],
               [AC_DEFINE([HAVE_NC_DEF_VAR_SZIP],[1],[Define to 1 for NetCDF4 nc_def_var_szip support])],,)])

AS_IF([test "x$ENABLE_NC4HDF5" = "xyes"],
      [AC_SEARCH_LIBS([H5get_libversion], [netcdf],
               [AC_DEFINE([HAVE_H5GET_LIBVERSION],[1],[Define to 1 for H5get_libversion support])],,[-lhdf5])])

AC_SUBST([ENABLE_NETCDF])
AC_SUBST([ENABLE_NC2])
AC_SUBST([ENABLE_NC4])
AC_SUBST([ENABLE_NC4HDF5])
AC_SUBST([ENABLE_NC4SZLIB])
AC_SUBST([NETCDF_ROOT])
AC_SUBST([NETCDF_INCLUDE])
AC_SUBST([NETCDF_LIBS])
#  ----------------------------------------------------------------------
#  Compile application with ECCODES library (for GRIB2 support)
ECCODES_INCLUDE=''
ECCODES_LIBS=''
AC_ARG_WITH([eccodes],
            [AS_HELP_STRING([--with-eccodes=<yes|no|directory>],
                            [location of ECCODES library for grib2 encoding/decoding (lib and include subdirs)])],
            [AS_CASE(["$with_eccodes"],
                     [no],[AC_MSG_CHECKING([for ECCODES library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([grib_api.h],,
                                             [AC_MSG_ERROR([Could not find grib_api.h])],
                                             [AC_INCLUDES_DEFAULT])
                            AC_SEARCH_LIBS([grib_get_message],
                                           [eccodes],
                                           [AC_DEFINE([HAVE_LIBGRIB_API],[1],[ECCODES library is present if defined to 1])],
                                           [AC_MSG_ERROR([Could not link to eccodes library])])],
                     [*],[ECCODES_ROOT=$with_eccodes
                          AS_IF([test -d "$ECCODES_ROOT"],
                                [LDFLAGS="-L$ECCODES_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$ECCODES_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([grib_api.h],,
                                                  [AC_MSG_ERROR([Could not find grib_api.h])],
                                                  [AC_INCLUDES_DEFAULT])
                                 AC_SEARCH_LIBS([grib_get_message],
                                                [eccodes],
                                                [AC_DEFINE([HAVE_LIBGRIB_API],[1],[ECCODES library is present if defined to 1])],
                                                [AC_MSG_ERROR([Could not link to eccodes library])])
                                 ECCODES_LIBS=" -L$ECCODES_ROOT/lib -leccodes"
                                 ECCODES_INCLUDE=" -I$ECCODES_ROOT/include"],
                                [AC_MSG_ERROR([$ECCODES_ROOT is not a directory! ECCODES suppressed])])])],
            [AC_MSG_CHECKING([for the ECCODES library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([ECCODES_INCLUDE])
AC_SUBST([ECCODES_LIBS])
# AM_CONDITIONAL([HAVE_LIBGRIB_API],[test "x$with_eccodes" != 'x' -a "x$with_eccodes" != 'xno' ])
#  ----------------------------------------------------------------------
#  Compile application with GRIB_API library (for GRIB2 support)
GRIB_API_INCLUDE=''
GRIB_API_LIBS=''
AC_ARG_WITH([grib_api],
            [AS_HELP_STRING([--with-grib_api=<yes|no|directory>],
                            [location of GRIB_API library for grib2 encoding/decoding (lib and include subdirs)])],
            [AS_CASE(["$with_grib_api"],
                     [no],[AC_MSG_CHECKING([for GRIB_API library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_HEADERS([grib_api.h],,
                                             [AC_MSG_ERROR([Could not find grib_api.h])],
                                             [AC_INCLUDES_DEFAULT])
                            AC_SEARCH_LIBS([grib_get_message],
                                           [grib_api],
                                           [AC_DEFINE([HAVE_LIBGRIB_API],[1],[GRIB_API library is present if defined to 1])],
                                           [AC_MSG_ERROR([Could not link to grib_api library])])],
                     [*],[GRIB_API_ROOT=$with_grib_api
                          AS_IF([test -d "$GRIB_API_ROOT"],
                                [LDFLAGS="-L$GRIB_API_ROOT/lib $LDFLAGS"
                                 CPPFLAGS="-I$GRIB_API_ROOT/include $CPPFLAGS"
                                 AC_CHECK_HEADERS([grib_api.h],,
                                                  [AC_MSG_ERROR([Could not find grib_api.h])],
                                                  [AC_INCLUDES_DEFAULT])
                                 AC_SEARCH_LIBS([grib_get_message],
                                                [grib_api],
                                                [AC_DEFINE([HAVE_LIBGRIB_API],[1],[GRIB_API library is present if defined to 1])],
                                                [AC_MSG_ERROR([Could not link to grib_api library])])
                                 GRIB_API_LIBS=" -L$GRIB_API_ROOT/lib -lgrib_api"
                                 GRIB_API_INCLUDE=" -I$GRIB_API_ROOT/include"],
                                [AC_MSG_ERROR([$GRIB_API_ROOT is not a directory! GRIB_API suppressed])])])],
            [AC_MSG_CHECKING([for the GRIB_API library])
             AC_MSG_RESULT([suppressed])])
AC_SUBST([GRIB_API_INCLUDE])
AC_SUBST([GRIB_API_LIBS])
AM_CONDITIONAL([HAVE_LIBGRIB_API],[test \( "x$with_grib_api" != 'x' -a "x$with_grib_api" != 'xno' \) -o \( "x$with_eccodes" != 'x' -a "x$with_eccodes" != 'xno' \) ])
#  ----------------------------------------------------------------------
#  Enable GRIB support
AC_MSG_CHECKING([for GRIB support])
AC_ARG_ENABLE([grib],
              [AS_HELP_STRING([--enable-grib],[GRIB support [default=yes]])],
              [AS_IF([test "x$enable_grib" != 'xno'],
                     [AC_DEFINE(HAVE_LIBGRIB, [1], [Define to 1 for GRIB support])
                      enable_grib=yes])],
              [AC_DEFINE(HAVE_LIBGRIB, [1], [Define to 1 for GRIB support])
               enable_grib=yes])
AC_MSG_RESULT([$enable_grib])
AC_SUBST([ENABLE_GRIB],[$enable_grib])
#  ----------------------------------------------------------------------
#  Enable ACROSS support
AC_MSG_CHECKING([for ACROSS support])
AC_ARG_ENABLE([across],
              [AS_HELP_STRING([--enable-across],[ACROSS support [default=yes]])],
              [AS_IF([test "x$enable_across" != 'xno'],
                     [AC_DEFINE(HAVE_ACROSS, [1], [Define to 1 for ACROSS support])
                      enable_across=yes])],
              [AC_DEFINE(HAVE_ACROSS, [1], [Define to 1 for ACROSS support])
               enable_across=yes])
AC_MSG_RESULT([$enable_across])
AC_SUBST([ENABLE_ACROSS],[$enable_across])
#  ----------------------------------------------------------------------
#  Compile interface with internal CGRIBEX library
AC_MSG_CHECKING([for CGRIBEX support])
AC_ARG_ENABLE([cgribex],
              [AS_HELP_STRING([--enable-cgribex],[Use the CGRIBEX library [default=yes]])],
              [AS_IF([test "x$enable_cgribex" != 'xno'],
                     [AC_DEFINE(HAVE_LIBCGRIBEX,[1],[Define to 1 for GRIB1 decoding/encoding with cgribex])
                      enable_cgribex=yes])],
              [AC_DEFINE(HAVE_LIBCGRIBEX,[1],[Define to 1 for GRIB1 decoding/encoding with cgribex])
               enable_cgribex=yes])
AC_MSG_RESULT([$enable_cgribex])
AC_SUBST([ENABLE_CGRIBEX],[$enable_cgribex])
#  ----------------------------------------------------------------------
#  Compile interface with internal SERVICE library
AC_MSG_CHECKING([for SERVICE support])
AC_ARG_ENABLE([service],
              [AS_HELP_STRING([--enable-service],[Use the service library [default=yes]])],
              [AS_IF([test "x$enable_service" != 'xno'],
                     [AC_DEFINE(HAVE_LIBSERVICE,[1],[Define to 1 for SERVICE interface])
                      enable_service=yes])],
              [AC_DEFINE(HAVE_LIBSERVICE,[1],[Define to 1 for SERVICE interface])
               enable_service=yes])
AC_MSG_RESULT([$enable_service])
AC_SUBST([ENABLE_SERVICE],[$enable_service])
#  ----------------------------------------------------------------------
#  Compile interface with internal EXTRA library
AC_MSG_CHECKING([for EXTRA support])
AC_ARG_ENABLE([extra],
              [AS_HELP_STRING([--enable-extra],[Use the extra library [default=yes]])],
              [AS_IF([test "x$enable_extra" != 'xno'],
                     [AC_DEFINE(HAVE_LIBEXTRA,[1],[Define to 1 for EXTRA interface])
                      enable_extra=yes])],
              [AC_DEFINE(HAVE_LIBEXTRA,[1],[Define to 1 for EXTRA interface])
               enable_extra=yes])
AC_MSG_RESULT([$enable_extra])
AC_SUBST([ENABLE_EXTRA],[$enable_extra])
#  ----------------------------------------------------------------------
#  Compile interface with internal IEG library
AC_MSG_CHECKING([for IEG support])
AC_ARG_ENABLE([ieg],
              [AS_HELP_STRING([--enable-ieg],[Use the ieg library [default=yes]])],
              [AS_IF([test "x$enable_ieg" != 'xno'],
                     [AC_DEFINE(HAVE_LIBIEG,[1],[Define to 1 for IEG interface])
                      enable_ieg=yes])],
              [AC_DEFINE(HAVE_LIBIEG,[1],[Define to 1 for IEG interface])
               enable_ieg=yes])
AC_MSG_RESULT([$enable_ieg])
AC_SUBST([ENABLE_IEG],[$enable_ieg])
#  ----------------------------------------------------------------------
# At the moment, there are two possible CDI bindings
# (default for CDO) linking directly to CDI convenience library with libtool
# (default for CDI) build and link to a shared CDI library
AS_IF([test "x$CDO_DISABLE_CDILIB" = "x1"],[enable_cdi_lib=no],[enable_cdi_lib=yes])
# save CDI binding mode for later automake use
AM_CONDITIONAL([ENABLE_CDI_LIB],[test x$enable_cdi_lib = 'xyes'])
# create shell variables for the representation of configure results
AS_IF([test x$enable_cdi_lib = 'xno'],[AC_SUBST([ENABLE_CDI_LIB],[false])],[AC_SUBST([ENABLE_CDI_LIB],[true])])
#  ----------------------------------------------------------------------
#  Build a static CDI
AC_MSG_CHECKING([for building an additional static CDI binary])
AC_ARG_ENABLE([all-static],
              [AS_HELP_STRING([--enable-all-static],[build a completely statically linked CDO binary [default=no]])],
              [AS_IF([test "x$enable_all_static" != "xno"],
                     [enable_all_static=yes],
                     [enable_all_static=no])],
              [enable_all_static=no])
AC_MSG_RESULT([$enable_all_static])
AM_CONDITIONAL([ENABLE_ALL_STATIC],[test x$enable_all_static = 'xyes'])
#  ----------------------------------------------------------------------
#  Build CDO with HIRLAM extensions
AC_MSG_CHECKING([for HIRLAM extensions])
AC_ARG_ENABLE([hirlam-extensions],
              [AS_HELP_STRING([--enable-hirlam-extensions],[HIRLAM extensions [default=no]])],
              [AS_IF([test "x$enable_hirlam_extensions" != "xno"],
                    [AC_DEFINE(HIRLAM_EXTENSIONS,[1],[Define to 1 for HIRLAM extensions])
                     enable_hirlam_extensions=yes],
                    [enable_hirlam_extensions=no])],
              [enable_hirlam_extensions=no])
AC_MSG_RESULT([$enable_hirlam_extensions])
AM_CONDITIONAL([ENABLE_HIRLAM_EXTENSIONS],[test x$enable_hirlam_extensions = 'xyes'])
# ----------------------------------------------------------------------
# Build CDI application
AC_ARG_ENABLE([cdi-app],
              [AS_HELP_STRING([--enable-cdi-app],[build and install CDI application [default=yes]])],
              [], [enable_cdi_app=yes])
AM_CONDITIONAL([ENABLE_CDI_APP], [test x$enable_cdi_app = 'xyes'])
#
])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
