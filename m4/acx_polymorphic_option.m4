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

# ACX_POLYMORPHIC_OPTION(PREFIX,
#                        PATTERN,
#                        VARIABLE)
# -----------------------------------------------------------------------------
# Declares a polymorphic configure option --PREFIX-PATTERN. PREFIX must be
# either 'enable' or 'with'. PATTERN must be a string containing a single
# asterisk ('*'). A space-separated list of the polymorphic parts of the
# provided options is stored in the shell variable VARIABLE (must be provided
# as a literal string). Options with empty polymorphic parts are not
# recognized.
#
# For example, the following line will make the configure script accept options
# '-enable-group-A', '--enable-group-B=value', '-disable-group-B', etc., and
# the shell variable 'acx_groups' will be set to ' A B B':
#
#     ACX_POLYMORPHIC_OPTION([enable], [group-*], [acx_groups])
#
# The values of the options from the example above can be retrieved as follows
# (the asterisk in the PATTERN is replaced with a value from VARIABLE):
#     for acx_option in $acx_groups; do
#       acx_option_var=AS_TR_SH([enable-group-$acx_option])
#       acx_option_val=AS_VAR_GET([$acx_option_var])
#     done
#
# The macro must be expanded before AC_INIT.
#
AC_DEFUN([ACX_POLYMORPHIC_OPTION],
  [dnl
dnl Check that the macro is expanded before AC_INIT (_AC_INIT_SRCDIR is the only
dnl AC_DEFUNed macro axpanded with non-AC_DEFUNed macro AC_INIT):
   AC_PROVIDE_IFELSE([_AC_INIT_SRCDIR],
     [m4_fatal([$0 must be expanded before AC_INIT])])
dnl Check that that the line marker we need is present in _AC_INIT_PARSE_ARGS:
   m4_pushdef([acx_marker_string], [^if test -n "$ac_unrecognized_opts"; then$])
   m4_bmatch(
     m4_dquote(m4_defn([_AC_INIT_PARSE_ARGS])),
     acx_marker_string, [],
     [m4_fatal([$0 is not compatible with the version of Autoconf in use ]dnl
[(_AC_INIT_PARSE_ARGS does not have the expected marker string)])])
dnl Check the first argument is either 'enable' or 'with':
   m4_case([$1], [enable], [], [with], [],
     [m4_fatal(
        [unexpected option prefix: '$1' (must be either 'enable' or 'with')])])
dnl Check the second argument contains a single asterisk symbol:
   m4_bmatch([$2], [^[-0-9A-Za-z]*\*[-0-9A-Za-z]*$], [],
     [m4_fatal(
        [unexpected option pattern: '$2' (must consist of alphanumeric ]dnl
[characters, hyphens and contain a single asterisk)])])
dnl Check the third argument is a literal variable name:
   AS_LITERAL_WORD_IF([$3], [], [m4_fatal([unexpected variable name: '$3' ]dnl
[(must be a literal variable name)])])
dnl Monkey-patch _AC_INIT_PARSE_ARGS:
   m4_define([_AC_INIT_PARSE_ARGS],
     m4_bpatsubst(
       m4_dquote(m4_defn([_AC_INIT_PARSE_ARGS])),
       acx_marker_string,
       [AS_VAR_SET([$3])
acx_unrecognized_opts=
acx_save_IFS=$IFS
IFS=','
for acx_opt in $ac_unrecognized_opts; do
  IFS=$acx_save_IFS
  acx_opt=`AS_ECHO(["$acx_opt"]) | sed 's/^ //'`
  acx_opt_recognized=no
  AS_CASE([$acx_opt],
    [--]m4_if([$1], [enable], [disable], [without])[-$2 | --$1-$2],
    [m4_pushdef([acx_expr_pattern],
       [m4_bpatsubst([$2], [\*], [\\\\(.*\\\\)])])dnl
     acx_opt_name=`expr "$acx_opt" : dnl
'--]m4_if([$1], [enable], [@<:@^-@:>@*able], [with@<:@^-@:>@*])[]dnl
[-acx_expr_pattern'`
     m4_popdef([acx_expr_pattern])dnl
     AS_IF([test -n "$acx_opt_name"],
       [AS_VAR_APPEND([$3], [" $acx_opt_name"])
        acx_opt_recognized=yes])])
  AS_VAR_IF([acx_opt_recognized], [no],
    [AS_IF([test -n "$acx_unrecognized_opts"],
       [AS_VAR_APPEND([acx_unrecognized_opts], [", $acx_opt"])],
       [acx_unrecognized_opts=$acx_opt])])
done
IFS=$acx_save_IFS
ac_unrecognized_opts=$acx_unrecognized_opts
\&]))])
