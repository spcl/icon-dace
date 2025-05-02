dnl acx_fc_logical_interop.m4 --- test whether the compiler supports
dnl c_loc for arrays of LOGICALS
dnl
dnl Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords:
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
AC_DEFUN([ACX_FC_LOGICAL_INTEROP],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AC_CACHE_CHECK([whether $FC supports C_LOC for LOGICAL arrays],
     [acx_cv_fc_logical_interop],
     [AC_LANG_PUSH([Fortran])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [      USE iso_c_binding, ONLY: c_loc, c_ptr
      IMPLICIT NONE
      LOGICAL :: LL(10)
      CALL BAR(LL)
      CONTAINS
      SUBROUTINE BAR(L)
      LOGICAL, TARGET, INTENT(IN) :: L(*)
      TYPE(C_PTR) :: CPTR
      CPTR = C_LOC(L)
      END SUBROUTINE BAR
])],
        [acx_cv_fc_logical_interop=yes],
        [acx_cv_fc_logical_interop=no])])
   m4_ifval([$1$2],
     [AS_IF([test $acx_cv_fc_logical_interop = yes],
        [$1],[$2])])
   AC_LANG_POP([Fortran])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
