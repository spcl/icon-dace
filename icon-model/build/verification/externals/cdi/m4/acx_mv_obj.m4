dnl acx_mv_obj.m4 --- shell function to rename a libtool .lo file
dnl
dnl Copyright  (C)  2019  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords:
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://swprojects.dkrz.de/redmine/projects/scales-ppm
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
dnl ACX_MV_OBJ([FROM],[TO],[ACTION-IF-SUCCEEDED],[ACTION-IF-FAILED])
dnl rename object file FROM.$ac_objext to TO.$ac_objext,
dnl in case libtool is used, adjusts .lo and moves
dnl .libs/FROM.$OBJEXT accordingly if necessary
dnl
AC_DEFUN([ACX_MV_OBJ],
  [AC_REQUIRE_SHELL_FN([acx_fn_mv_obj],
     [AS_FUNCTION_DESCRIBE([ac_fn_mv_obj], [LINENO FROM TO],
       [Rename FROM.$ac_objext to TO.$ac_objext, and return whether this succeeded.])],
     [AS_LINENO_PUSH([$[]1])
      acx_path_from=`echo "$[]2" | sed -e 's!/\{0,1\}@<:@^/@:>@*$[]!!'`
      test -z "$acx_path_from" && acx_path_from=.
      acx_fn_from=`echo "$[]2" | sed -e 's@^.*/@@'`
      acx_path_to=`echo "$[]3" | sed -e 's!/\{0,1\}@<:@^/@:>@*$[]!!'`
      test -z "$acx_path_to" && acx_path_to=.
      acx_fn_to=`echo "$[]3" | sed -e 's@^.*/@@'`
      AS_IF([expr "$ac_compile" : '.*/libtool --mode=compile' >/dev/null],
        [_AC_RUN_LOG([sed 's@\(pic_object='"'\)"'\(\(.libs/\)\{0,1\}\)'"$acx_fn_from"'\.o'"'"'@\1\2'"$acx_fn_to"'.o'"'"'@' "$][2.$ac_objext" >"$][3.$ac_objext" && rm "$][2.$ac_objext" \
        && if test -f "$][2.$OBJEXT" ; then mv "$][2.$OBJEXT" "$][3.$OBJEXT" ; fi \
        && if test -f "$acx_path_from/.libs/$acx_fn_from.$OBJEXT" ; then \
          as_dir="$acx_path_to/.libs" as_fn_mkdir_p ; \
          mv "$acx_path_from/.libs/$acx_fn_from.$OBJEXT" \
             "$acx_path_to/.libs/$acx_fn_to.$OBJEXT" ; fi],
            [_AS_ECHO_LOG([Renaming object file $][2.$ac_objext to $][3.$ac_objext.])])],
         [_AC_RUN_LOG([mv "$][2.$ac_objext" "$][3.$ac_objext"],
            [_AS_ECHO_LOG([Renaming object file $][2.$ac_objext to $][3.$ac_objext.])])])
       AS_LINENO_POP
       AS_SET_STATUS([$ac_status])])dnl
   m4_ifval([$3$4],[AS_IF([acx_fn_mv_obj "$LINENO" $1 $2], [$3], [$4])],
     [acx_fn_mv_obj "$LINENO" $1 $2])])
dnl
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
dnl
