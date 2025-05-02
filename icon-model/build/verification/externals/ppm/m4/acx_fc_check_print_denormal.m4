dnl acx_fc_check_print_denormal.m4 --- Check if the Fortran compiler can
dnl                                 print denormal numbers
dnl
dnl Copyright  (C)  2013  Thomas Jahns <jahns@dkrz.de>
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
dnl ACX_FC_CHECK_PRINT_DENORMAL_IFELSE([ACTION-IF-SUPPORTED],[ACTION-IF-NOT])
dnl Test if programs produced by the Fortran compiler can print
dnl denormal values (This is a known bug in PGI compiler configurations)
dnl
AC_DEFUN([ACX_FC_CHECK_PRINT_DENORMAL_IFELSE],
  [AC_CACHE_CHECK([if Fortran programs can handle printing of denormals],
     [acx_cv_fc_print_denormal_works],
     [acx_cv_fc_print_denormal_works=no
      _AC_FORTRAN_ASSERT()dnl
      AC_LINK_IFELSE([AC_LANG_PROGRAM(,
[      INTEGER, PARAMETER :: pd = 12
      INTEGER, PARAMETER :: rd = 307
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd)
      REAL(dp) :: d_fraction, d
      INTEGER :: d_exponent
      d_fraction = -3.0512077779711766E-002_dp
      d_exponent = MINEXPONENT(d)
      d = SCALE(d_fraction, d_exponent)
      PRINT *, d])],
        [acx_cv_fc_print_denormal_works=`./conftest`
         _AC_RUN_LOG([echo "$acx_cv_fc_print_denormal_works" >&2],
            [_AS_ECHO_LOG([./conftest])])
         AS_IF([echo "$acx_cv_fc_print_denormal_works" | grep '^ \*\** *$' >/dev/null],
           [acx_cv_fc_print_denormal_works=no],
           [echo "$acx_cv_fc_print_denormal_works" | grep '^ *[[-0-9.]]*'  >/dev/null],
           [acx_cv_fc_print_denormal_works=yes],
           [acx_cv_fc_print_denormal_works=no])])],
     [acx_cv_fc_print_denormal_works=no])
   AS_IF([test x"$acx_cv_fc_print_denormal_works" = xyes],
     [$1],
     m4_dquote(m4_ifval([$2],[$2],
       [AC_MSG_FAILURE([Fortran programs cannot print denormals!
In case you are using the PGI Fortran compiler consider compiler flag -Mnodaz])])))])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
