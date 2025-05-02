dnl acx_uuid.m4 --- check for library to create UUIDs
dnl
dnl Copyright  (C)  2018  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Keywords: configure configure.ac autoconf MPI mpirun mpiexec
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
dnl
dnl ACX_UUID
dnl
dnl Find one out of the supported UUID packages, these are:
dnl util-linux libuuid, the most frequently found option on Linux
dnl OSSP UUID, used by CMOR
dnl DCE uuid_create is part of the standard C library on some commercial Unix
dnl                 systems
dnl
dnl Sets UUID_C_LIB to the flags needed for linking and UUID_C_INCLUDE to flags
dnl needed for corresponding header.
AC_DEFUN([ACX_UUID],
  [have_dce_uuid_c_bindings=no
   have_util_linux_uuid_c_bindings=no
   have_ossp_uuid_c_bindings=no
   AC_ARG_ENABLE([util-linux-uuid],
     [AS_HELP_STRING([--disable-util-linux-uuid],
        [Do not test for the util-linux UUID library, ]dnl
[use OSSP or DCE versions instead])],,
     [enable_util_linux_uuid=auto])
   AC_ARG_ENABLE([ossp-uuid],
     [AS_HELP_STRING([--disable-ossp-uuid],
        [Do not test for the OSSP UUID library, ]dnl
[use util-linux or DCE versions instead])],,
     [enable_ossp_uuid=auto])
   AC_ARG_ENABLE([dce-uuid],
     [AS_HELP_STRING([--disable-dce-uuid],
        [Do not test for the DCE UUID library, ]dnl
[use util-linux or OSSP versions instead])],,
     [enable_dce_uuid=auto])

   AC_ARG_WITH([util-linux-uuid],
     [AS_HELP_STRING([--with-util-linux-uuid@<:@=yes|no|directory>@:>@],
        [Use util-linux UUID library, do not test for OSSP UUID library.])],
     [AS_CASE([$with_util_linux_uuid],
        [yes],[enable_ossp_uuid=no enable_dce_uuid=no],
        [no],[enable_util_linux_uuid=no],
        [/*],
        [AS_IF([test x${with_util_linux_uuid_root+set} != xset],
           [with_util_linux_uuid_root=$with_util_linux_uuid],
           [test x"${with_util_linux_uuid_root}" != x"${with_util_linux_uuid}"],
           [AC_MSG_FAILURE([inconsistent directories specified for dnl
--with-util-linux-uuid and --with-util-linux-uuid-root])])
         enable_ossp_uuid=no enable_dce_uuid=no])])
   AC_ARG_WITH([ossp-uuid],
     [AS_HELP_STRING([--with-ossp-uuid@<:@=yes|no|directory>@:>@],
        [Use OSSP UUID library, do not test for the util-linux UUID library.])],
     [AS_CASE([$with_ossp_uuid],
        [yes],[enable_util_linux_uuid=no enable_dce_uuid=no],
        [no],[enable_ossp_uuid=no],
        [/*],
        [AS_IF([test x${with_ossp_uuid_root+set} != xset],
           [with_ossp_uuid_root=$with_ossp_uuid],
           [test x"${with_ossp_uuid_root}" != x"${with_ossp_uuid}"],
           [AC_MSG_FAILURE([inconsistent directories specified for dnl
--with-ossp-uuid and --with-ossp-uuid-root])])
         enable_util_linux_uuid=no enable_dce_uuid=no])])
   AC_ARG_WITH([dce-uuid],
     [AS_HELP_STRING([--with-dce-uuid@<:@=yes|no|directory>@:>@],
        [Use DCE UUID library, do not test for the util-linux UUID library.])],
     [AS_CASE([$with_dce_uuid],
        [yes],[enable_util_linux_uuid=no enable_ossp_uuid=no],
        [no],[enable_dce_uuid=no],
        [/*],
        [AS_IF([test x${with_dce_uuid_root+set} != xset],
           [with_dce_uuid_root=$with_dce_uuid],
           [test x"${with_dce_uuid_root}" != x"${with_dce_uuid}"],
           [AC_MSG_FAILURE([inconsistent directories specified for dnl
--with-dce-uuid and --with-dce-uuid-root])])
         enable_util_linux_uuid=no enable_ossp_uuid=no])])

   AS_IF([test x"$enable_util_linux_uuid" != xno],
     [ACX_C_PACKAGE([util-linux-uuid],[uuid/uuid.h],,,[],
        [uuid_generate],[[uuid]],,,[])
      AS_IF([test x"$have_util_linux_uuid_c_bindings" = xyes],
        [acx_save_CPPFLAGS=$CPPFLAGS
         CPPFLAGS="$CPPFLAGS $UTIL_LINUX_UUID_C_INCLUDE"
         AC_CHECK_HEADERS([uuid/uuid.h],
           [AC_CHECK_DECLS([uuid_generate],,
              [have_util_linux_uuid_c_bindings=no],[AC_INCLUDES_DEFAULT
@%:@include <uuid/uuid.h>])],[have_util_linux_uuid_c_bindings=no],[AC_INCLUDES_DEFAULT])
         CPPFLAGS=$acx_save_CPPFLAGS])],
     [have_util_linux_uuid_c_bindings=no])
   AS_IF([test x"$enable_ossp_uuid" != xno -a x"$have_util_linux_uuid_c_bindings" = xno],
     [ACX_C_PACKAGE([ossp-uuid],[uuid.h],[@%:@include <uuid.h>
AC_INCLUDES_DEFAULT],,[],
        [uuid_create],[[ossp-uuid] [uuid]],,,[],,[[/usr/include/ossp]])
      AS_IF([test x"$have_ossp_uuid_c_bindings" = xyes],
        [acx_save_CPPFLAGS=$CPPFLAGS
         CPPFLAGS="${CPPFLAGS+$CPPFLAGS }$OSSP_UUID_C_INCLUDE"
         AC_CHECK_HEADERS([uuid.h],
           [AC_CHECK_DECLS([uuid_create],
              [AC_CHECK_DECLS([UUID_MAKE_V5],,[have_ossp_uuid_c_bindings=no],
                 [@%:@include <uuid.h>
AC_INCLUDES_DEFAULT])],[have_ossp_uuid_c_bindings=no],
                 [@%:@include <uuid.h>
AC_INCLUDES_DEFAULT])],[have_ossp_uuid_c_bindings=no],[@%:@include <uuid.h>
AC_INCLUDES_DEFAULT])
         CPPFLAGS=$acx_save_CPPFLAGS])],
     [have_ossp_uuid_c_bindings=no])
   # check for DCE uuid_create if util-linux and OSSP variants cannot be found
   AS_IF([test x"$enable_dce_uuid" != xno -a x"$have_util_linux_uuid_c_bindings$have_ossp_uuid_c_bindings" = xnono],
     [ACX_C_PACKAGE([dce-uuid],[uuid.h],,,[],
        [uuid_create])
      AS_IF([test x"$have_dce_uuid_c_bindings" = xyes],
        [acx_save_CPPFLAGS=$CPPFLAGS
         CPPFLAGS="$CPPFLAGS $DCE_UUID_C_INCLUDE"
         AC_CHECK_HEADERS([uuid.h],
           [AC_CHECK_DECLS([uuid_create],
              [have_dce_uuid_c_bindings=yes],
              [have_dce_uuid_c_bindings=no],[AC_INCLUDES_DEFAULT
@%:@include <uuid.h>])],[have_dce_uuid_c_bindings=yes], [AC_INCLUDES_DEFAULT])
         CPPFLAGS=$acx_save_CPPFLAGS])])

AS_IF([test x"$have_util_linux_uuid_c_bindings" = xyes],
  [UUID_C_INCLUDE=$UTIL_LINUX_UUID_C_INCLUDE
   UUID_C_LIB=$UTIL_LINUX_UUID_C_LIB],
  [test x"$have_ossp_uuid_c_bindings" = xyes],
  [UUID_C_INCLUDE=$OSSP_UUID_C_INCLUDE
   UUID_C_LIB=$OSSP_UUID_C_LIB],
  [test x"$have_dce_uuid_c_bindings" = xyes],
  [UUID_C_INCLUDE=$DCE_UUID_C_INCLUDE
   UUID_C_LIB=$DCE_UUID_C_LIB])
AC_SUBST([UUID_C_INCLUDE])
AC_SUBST([UUID_C_LIB])])

dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
