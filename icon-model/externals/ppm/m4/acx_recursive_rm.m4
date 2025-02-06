dnl acx_recursive_rm.m4 --- fixes for autoconf macros that clean up
dnl                         incompletely, resulting in warnings
dnl
dnl Copyright  (C)  2020  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: configure configure.ac autotools
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
m4_pushdef([AC_EGREP_CPP],m4_bpatsubst(m4_dquote(m4_defn([AC_EGREP_CPP])),[
rm -f conftest\*
],[
rm -rf conftest*
]))dnl
m4_pushdef([_AC_PROG_FC_C_O],m4_bpatsubst(m4_dquote(m4_defn([_AC_PROG_FC_C_O])),[
rm -f conftest\*],[
rm -rf conftest*]))dnl
m4_pushdef([_AM_PROG_CC_C_O],m4_bpatsubst(m4_dquote(m4_defn([_AM_PROG_CC_C_O])),[
  rm -f core conftest\*],[
  rm -rf core conftest*]))dnl
m4_pushdef([LT_PATH_NM],m4_bpatsubst(m4_dquote(m4_defn([LT_PATH_NM])),[
  rm -f conftest\*],[
  rm -rf conftest*]))dnl
m4_pushdef([_LT_CMD_GLOBAL_SYMBOLS],m4_bpatsubst(m4_dquote(m4_defn([_LT_CMD_GLOBAL_SYMBOLS])),[pipe_works=no

  rm -f conftest\*],[pipe_works=no

  rm -rf conftest*]))dnl
m4_pushdef([_LT_PATH_MANIFEST_TOOL],m4_bpatsubst(m4_dquote(m4_defn([_LT_PATH_MANIFEST_TOOL])),[
  rm -f conftest\*],[
  rm -rf conftest*]))dnl
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
