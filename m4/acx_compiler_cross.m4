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

# ACX_COMPILER_CROSS_CHECK_DELAY()
# -----------------------------------------------------------------------------
# Delays the standard cross-compilation check to let the user avoid a false
# positive result (e.g. extend LDFLAGS with RPATHs, so the executable that is
# run during the check would not cause an error while loading shared libraries
# provided in LIBS).
#
# Example:
#     ACX_COMPILER_CROSS_CHECK_DELAY
#     AC_PROG_FC
#     AC_LANG([Fortran])
#     AC_FC_PP_SRCEXT([f90])
#     ACX_SHLIB_FC_RPATH_FLAG
#     LDFLAGS=`AS_ECHO(["$LDFLAGS"]) | dnl
#     sed ['s%\(-L[ ]*\([^ ][^ ]*\)\)%\1 '"$acx_cv_fc_rpath_flag"'\2%g']`
#     ACX_COMPILER_CROSS_CHECK_NOW
#
AC_DEFUN([ACX_COMPILER_CROSS_CHECK_DELAY],
  [acx_save_cross_compiling=$cross_compiling
   acx_save_ac_tool_warned=$ac_tool_warned
   cross_compiling=yes
   ac_tool_warned=yes
   m4_pushdef([_AC_COMPILER_EXEEXT_CROSS])])

# ACX_COMPILER_CROSS_CHECK_NOW()
# -----------------------------------------------------------------------------
# Runs the cross-compilation check that has been delayed with
# ACX_COMPILER_CROSS_CHECK_DELAY.
#
AC_DEFUN([ACX_COMPILER_CROSS_CHECK_NOW],
  [AC_PROVIDE_IFELSE([ACX_COMPILER_CROSS_CHECK_DELAY], [],
     [m4_fatal([$0 must be called after ACX_COMPILER_CROSS_CHECK_DELAY])])dnl
   cross_compiling=$acx_save_cross_compiling
   ac_tool_warned=$acx_save_ac_tool_warned
   m4_popdef([_AC_COMPILER_EXEEXT_CROSS])dnl
   AC_LANG_CONFTEST([_AC_LANG_IO_PROGRAM])
   _AC_COMPILER_EXEEXT_CROSS
   rm -f conftest.$ac_ext conftest$ac_exeext conftest.out])
