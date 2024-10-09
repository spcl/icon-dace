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

# ACX_FC_PP_COMMENTS([ACTION-IF-SUCCESS],
#                    [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Checks whether the Fortran compiler handles C-style block comments as well as
# single line in the context of macro definitions. The result is either "yes"
# or "no".
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_pp_comments.
#
AC_DEFUN([ACX_FC_PP_COMMENTS],
  [AC_LANG_ASSERT([Fortran])dnl
   AC_PROVIDE_IFELSE([AC_FC_PP_SRCEXT], [],
     [m4_warn([syntax],
        [ACX_FC_PP_COMMENTS requires calling the Fortran compiler with a ]dnl
[preprocessor but no call to AC_FC_PP_SRCEXT is detected])])dnl
   AC_CACHE_CHECK([whether Fortran compiler supports C-style comments],
     [acx_cv_fc_pp_comments],
     [acx_cv_fc_pp_comments=no
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
[[/* block
comment */
#define CONFTEST_MACRO1 // single line comment
#define CONFTEST_MACRO2 /* block comment */
#ifndef CONFTEST_MACRO1
      choke me
#endif
#ifndef CONFTEST_MACRO2
      choke me
#endif]])],
        [acx_cv_fc_pp_comments=yes])])
   AS_VAR_IF([acx_cv_fc_pp_comments], [yes], [$1],
     [m4_default([$2],
        [AC_MSG_FAILURE(
           [Fortran compiler does not support C-style comments])])])])
