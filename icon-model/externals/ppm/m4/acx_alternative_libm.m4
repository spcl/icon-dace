dnl acx_alternative_libm.m4 --- construct alternative math library compiler flag
dnl
dnl Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: libm replacement
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
dnl ACX_ALTERNATIVE_LIBM([CANDIDATES])
dnl tries to find a suitable -lm substitute
dnl defines an output variable LIBM to use in tests potentially adding -lm
dnl otherwise
AC_DEFUN([ACX_ALTERNATIVE_LIBM],
  [AC_REQUIRE([AX_WITH_PERL])
   AC_LANG_PUSH([C])
   AC_MSG_CHECKING([for alternative math library])
   AS_VAR_PUSHDEF([ac_Search], [acx_cv_option_search_acosh_]_AC_LANG_ABBREV)dnl
   while :; do
     m4_foreach([libmSubst], [$1],
       [# try to remove -lm from LIBS first!
        LIBS_save=$LIBS
        LIBS=`echo "$LIBS" | $PERL -e 'while(<>) { s{\s\K-lm(?=\s|$)}{'"libmSubst"'}g ; print }'`
        _ACX_OPTION_SEARCH_LIBS([acosh],,[libmSubst])
        LIBS=$LIBS_save
        AS_VAR_SET_IF([ac_Search],[break])
       ])
     break
   done
   AS_VAR_SET_IF([ac_Search],
     [LIBM=AS_VAR_GET([ac_Search])
      LIBS=`echo "$LIBS" |  $PERL -e 'while(<>) { s{\s\K-lm(?=\s|$)}{'"$LIBM"'}g ; print }'`
      AC_MSG_RESULT([$LIBM])],
     [LIBM=-lm
      AC_MSG_RESULT([default -lm])])
   AS_VAR_POPDEF([ac_Search])dnl        _ACX_OPTION_SEARCH_LIBS([acosh],,
   AC_SUBST([LIBM])
   AC_LANG_POP([C])
  ])

dnl USE_ALTERNATIVE_LIBM([LIBVARS])
dnl substitute a libm alternative in LIBS FCLIBS output variables and
dnl in each element of [LIBVARS]
AC_DEFUN([ACX_USE_ALTERNATIVE_LIBM],
  [m4_foreach([ACX_LibSet],[[LIBS],[FCLIBS]],
     [ACX_LibSet=`echo "$ACX_LibSet" | $PERL -e 'while(<>) { s{\s\K-lm(?=\s|$)}{'"$LIBM"'}g ; print }'`
     ])
   m4_ifval($1,[m4_foreach([ACX_LibSet],[$1],
        [ACX_LibSet=`echo "$ACX_LibSet" | $PERL -e 'while(<>) { s{\s\K-lm(?=\s|$)}{'"$LIBM"'}g ; print }'`
     ])])
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
