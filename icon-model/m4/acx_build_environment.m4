# Copyright (c) 2018-2024, MPI-M
#
# Author: Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# ACX_BUILD_ENVIRONMENT()
# -----------------------------------------------------------------------------
# Declares a precious variable BUILD_ENV to be set to the initialization code
# of the building environment. If the variable is not empty (i.e. set by the
# user), tests whether the code is accepted by the $SHELL. If the test is
# successful, tries to run the code to initialize the shell of the configure
# script for the following checks. If either the test or the attempt to run the
# contents of BUILD_ENV inside the configure script fail, the expanded script
# fails the configuration.
#
# Additionally, the macro declares an output variable BUILD_ENV_MAKE, which
# holds a value of BUILD_ENV modified to conform with Makefile syntax (e.g.
# each dollar sign "$" is duplicated).
#
AC_DEFUN([ACX_BUILD_ENVIRONMENT],
  [AC_ARG_VAR([BUILD_ENV],
     [initialization code to set up the building environment (must be a ]dnl
[single line ending with a semicolon), e.g. ]dnl
['LICENSE=file.lic; export LICENSE;'])dnl
   AC_SUBST([BUILD_ENV_MAKE])
   AS_IF([test -n "$BUILD_ENV"],
     [AS_IF([AS_ECHO(["$BUILD_ENV"]) | grep ';@<:@ @:>@*$' >/dev/null 2>&1],
        [],
        [AC_MSG_ERROR(
           [\$BUILD_ENV does not end with a semicolon: '$BUILD_ENV'])])
      AS_CASE([$BUILD_ENV], [*${as_nl}*],
        [AC_MSG_ERROR([\$BUILD_ENV contains a newline: '$BUILD_ENV'])])
      BUILD_ENV_MAKE=`echo "$BUILD_ENV" | sed 's/\\$/$$/g'`
      AC_MSG_CHECKING([whether \$BUILD_ENV is accepted by '$SHELL -c'])
      acx_build_env_quoted=$BUILD_ENV
      ASX_ESCAPE_SINGLE_QUOTE([acx_build_env_quoted])
      _AS_ECHO_LOG([$SHELL -c '$acx_build_env_quoted'])
      eval \$SHELL -c "'$acx_build_env_quoted'" >&AS_MESSAGE_LOG_FD 2>&1
      AS_IF([test $? -eq 0],
        [AC_MSG_RESULT([yes])
dnl Check that $BUILD_ENV does not change variables that have been set on the
dnl command line:
         acx_build_env_vars_to_check=
         eval "set dummy $ac_configure_args"; shift
         for acx_arg; do
           AS_CASE([$acx_arg],
             [-*], [],
             [*=*], [dnl
dnl The configure script has already checked that all arguments matching
dnl pattern '*=*' have valid shell variable names on the left-hand side.
            acx_arg_name=`expr "x$acx_arg" : 'x\(@<:@^=@:>@*\)='`
            acx_arg_cmd_value=`expr "x$acx_arg" : '@<:@^=@:>@*=\(.*\)'`
            AS_VAR_COPY([acx_arg_${acx_arg_name}], [$acx_arg_name])
dnl Check only those variables that have not been modified since they were set
dnl on the command line (otherwise, it is responsibility of the configure
dnl script developers):
            AS_VAR_IF(
              [acx_arg_cmd_value], ["AS_VAR_GET([acx_arg_${acx_arg_name}])"],
              [AS_VAR_APPEND([acx_build_env_vars_to_check],
                 [" $acx_arg_name"])],
              [AS_UNSET([acx_arg_${acx_arg_name}])])])
         done
         AC_MSG_CHECKING(
           [whether \$BUILD_ENV is accepted by the current shell])
         _AS_ECHO_LOG([$BUILD_ENV])
         eval "$BUILD_ENV" >&AS_MESSAGE_LOG_FD 2>&1
         AS_IF([test $? -eq 0],
           [AC_MSG_RESULT([yes])
            AS_IF([test -n "$acx_build_env_vars_to_check"],
              [for acx_arg_name in $acx_build_env_vars_to_check; do
                 AS_IF([test x"AS_VAR_GET([$acx_arg_name])" != \
x"AS_VAR_GET([acx_arg_${acx_arg_name}])"],
                   [AC_MSG_WARN([\$BUILD_ENV has modified variable dnl
'$acx_arg_name', which was set on the command line to dnl
"AS_VAR_GET([acx_arg_${acx_arg_name}])": new value of the variable is dnl
"AS_VAR_GET([$acx_arg_name])"])])
                 AS_UNSET([acx_arg_${acx_arg_name}])
               done])],
           [AC_MSG_RESULT([no])])],
        [AC_MSG_RESULT([no])
         acx_failMsg="failed to initialize '$SHELL' with the provided dnl
BUILD_ENV='$BUILD_ENV'
A possible workaround for the problem is to re-run the configuration with dnl
the following command:
CONFIG_SHELL=\$SHELL \$SHELL $as_myself"
         for acx_config_arg in "$[@]"; do
           AS_VAR_APPEND([acx_failMsg], " '$acx_config_arg'")
         done
         AC_MSG_FAILURE([$acx_failMsg])])])])
