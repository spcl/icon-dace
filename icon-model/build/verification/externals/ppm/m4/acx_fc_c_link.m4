dnl acx_fc_c_link.m4 --- transform library c flags into version
dnl                      that suits the fortran compiler
dnl
dnl Copyright  (C)  2011  Thomas Jahns <jahns@dkrz.de>
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
dnl ACX_FC_C_LINK([C_LIBS_VAR],[FC_LIBS_VAR],[OPTIONAL-TEST-SOURCE])
dnl Try to construct FC_LIBS_VAR from C_LIBS_VAR by checking known
dnl idiosyncrasies of Fortran compiler front-ends. Optionally compiles
dnl a user-provided program source given in OPTIONAL-TEST-SOURCE.
AC_DEFUN([ACX_FC_C_LINK],
  [AC_REQUIRE([AX_WITH_PERL])
   save_LIBS=$LIBS
   LIBS="$][$1 $LIBS"
   AC_LANG_PUSH([Fortran])
   AC_LINK_IFELSE(m4_ifval([$3],[$3],[AC_LANG_PROGRAM(,
[      REAL :: p
      p = 0.1])]),
     [$2=$][$1],
     [dnl try to massage linker flags to suit the Fortran compiler
      AC_MSG_CHECKING([how to pass flags $1 to Fortran compiler])
      $2=`echo "$][$1" | $PERL "$srcdir/util/naglinkhack.pl"`
      LIBS="$save_LIBS $][$2"
      AC_LINK_IFELSE(m4_ifval([$3],[$3],[AC_LANG_PROGRAM(,
[      REAL :: p
      p = 0.1])]),
        [AC_MSG_RESULT([$2=$][$2])],
        [AC_MSG_ERROR([Cannot successfully link Fortran programs with $1 value of $][$1])])
      ])
   AC_LANG_POP([Fortran])
   LIBS=$save_LIBS])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
