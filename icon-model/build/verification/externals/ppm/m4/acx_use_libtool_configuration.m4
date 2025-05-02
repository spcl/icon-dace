dnl acx_use_libtool_configuration.m4 --- prevent problematic libtool build configurations
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
dnl _ACX_LT_FLAGS_MANGLE
m4_define([nag_filter_flag_var],
  [AS_IF([test x${$1+set} = xset],
     [$1=`echo "$][$1" | tr ' ' '\n' | sed -e '/^-W@<:@lc@:>@/{' \
-e 's/^\(-Wl\)/-XCClinker \1/;s/^\(-Wc,\)/-Xcompiler \1/' \
-e 's/^\(-Wc=.*\)/-Xcompiler \1 -XCClinker \1/' -e '}' \
   | tr '\n' ' ' | sed -e 's/ $//'`])])
AC_DEFUN([_ACX_LT_FORT_FLAGS_MANGLE],
  [_AC_FORTRAN_ASSERT
   AC_LANG_CASE([Fortran],
     [AS_VAR_PUSHDEF([acx_FC],[FC])dnl
      AS_VAR_PUSHDEF([acx_FCFLAGS],[FCFLAGS])dnl
      AS_VAR_PUSHDEF([acx_LDFLAGS],[FCLDFLAGS])],
     [Fortran 77],
     [AS_VAR_PUSHDEF([acx_FC],[F77])dnl
      AS_VAR_PUSHDEF([acx_FCFLAGS],[FFLAGS])dnl
      AS_VAR_PUSHDEF([acx_LDFLAGS],[F77LDFLAGS])])
   acx_temp=`$acx_FC -V 2>&1 | sed -n 1,5p`
dnl fix problems from NAG compiler
   AS_CASE(["$acx_temp"],
     [*NAG\ Fortran\ Compiler\ Release*],
     [nag_filter_flag_var([acx_FCFLAGS])
      nag_filter_flag_var([acx_LDFLAGS])])
dnl fix conflicting use of -module by libtool and ifort
   AC_LANG_CASE([Fortran],
     [AS_CASE(["x${FC_MODOUT}x"],
        [x'-module 'x|x'-mod 'x],
        [FC_MODOUT="-Xcompiler ${FC_MODOUT}-Xcompiler "])])dnl
   AS_VAR_POPDEF([acx_FC])dnl
   AS_VAR_POPDEF([acx_FCFLAGS])dnl
   AS_VAR_POPDEF([acx_LDFLAGS])])
dnl
dnl ACX_USE_LIBTOOL_CONFIGURATION([ARGS-TO-LT_INIT])
dnl Switch compiler to libtool wrapper and prevent occurrence of
dnl problematic setups
AC_DEFUN([ACX_USE_LIBTOOL_CONFIGURATION],
  [dnl before switching on libtool, identify compilers that prevent us from
   dnl certain build configurations
   ACX_LT_PROBLEMS
dnl add some monkey patching for older libtool versions that don't handle
dnl newer PGI or NAG configurations particularly well
   m4_if(m4_cmp(m4_version_compare(LT_PACKAGE_VERSION,[2.4.6]),1),-1,
     [m4_pushdef([_LT_COMPILER_PIC],m4_bpatsubst(m4_dquote(m4_bpatsubst(m4_dquote(m4_defn([_LT_COMPILER_PIC])),[;;
	\*Portland\\ ],[;; @%:@(
	*NAG\\ Fortran\\ Compiler*)
	  _LT_TAGVAR(lt_prog_compiler_wl, $][1)='-Wl,-Wl,,'
	  _LT_TAGVAR(lt_prog_compiler_pic, $][1)='-PIC'
	  _LT_TAGVAR(lt_prog_compiler_static, $][1)='-Bstatic'
	  ;;
	*PGI\\ Compilers\\ and\\ Tools*|*NVIDIA\\ Compilers\\ and\\ Tools*|*Port][land\\ ])),[sed 5q`],[sed -n 1,5p`]))dnl
      m4_pushdef([_LT_LANG_CXX_CONFIG],m4_bpatsubst(m4_dquote(m4_defn([_LT_LANG_CXX_CONFIG])),
        [sed 5q`],[sed -n 1,5p`]))dnl
      m4_pushdef([_LT_LINKER_SHLIBS],m4_bpatsubst(m4_dquote(
        m4_bpatsubst(m4_dquote(m4_bpatsubst(m4_dquote(
        m4_bpatsubst(m4_dquote(m4_bpatsubst(m4_dquote(m4_defn(
          [_LT_LINKER_SHLIBS])),[tmp_sharedflag='-shared'],
          [tmp_sharedflag='-shared'
	tmp_compiler_flags='$compiler_flags'])),
        [\$CC '"\$tmp_sharedflag""\$tmp_addflag"' \$libobjs \$deplibs \$compiler_flags \$wl-soname],
        [$CC '"$tmp_sharedflag""$tmp_addflag"' $libobjs $deplibs '"$tmp_compiler_flags"' $wl-soname])),
        [  tmp_sharedflag='-Wl,-shared'],
        [  tmp_sharedflag='-Wl,-shared'
	  tmp_compiler_flags='`echo \$compiler_flags | sed -e '"'"'s/ -W@<:@cl@:>@,-no-pie\\b//g'"'"'`'])),
        [\*Sun\\ F\*\(.\)[ 	]*# Sun Fortran 8\.3
[ 	]*tmp_sharedflag='-G' ;;
],[\&	*NAG\\ Fortran\\ Compiler*\1
	  tmp_sharedflag='-Wl,-shared'
	  tmp_compiler_flags='`echo \$compiler_flags | sed -e '"'"'s/ -W@<:@cl@:>@,-no-pie\\b//g'"'"'`' ;;
])),[sed 5q`],[sed -n 1,5p`]))dnl
      m4_pushdef([_LT_SYS_HIDDEN_LIBDEPS],[AS_UNSET([output_verbose_link_cmd])]
        m4_bpatsubst(m4_dquote(m4_bpatsubst(m4_dquote(
          m4_defn([_LT_SYS_HIDDEN_LIBDEPS])),[test x-\([LR]\) = "\$p"],
            [test x-\1 = x"$p"])),
          [test x-R = x"\$p"],[\& ||
	   test x-l = x"$p"]))])dnl
   m4_foreach([acx_ltconfig],[[[_LT_LANG_C_CONFIG]],[[_LT_LANG_F77_CONFIG]],[[_LT_LANG_FC_CONFIG]]],
     [m4_pushdef(acx_ltconfig,m4_bpatsubst(m4_dquote(
        m4_bpatsubst(m4_dquote(m4_defn(acx_ltconfig)),[_LT_TAG_COMPILER
],
       [\&  _LT_TAGDECL([with_nagfor], [acx_is_nagfor], [0], [Is the compiler the NAG Fortran compiler?])])),
       [_LT_CONFIG],
       [AC_MSG_CHECKING([whether this is the NAG Fortran compiler])
    $CC -V 2>&1 | grep '^NAG Fortran Compiler Release' >/dev/null 2>&1
    _lt_result=$?
    AS_IF([test $_lt_result -eq 0],[_lt_result=yes],[_lt_result=no])
    AC_MSG_RESULT([$_lt_result])
    _LT_TAGVAR([acx_is_nagfor], $][1)=$_lt_result
    \&]))])dnl
   LT_INIT([$1])
   m4_popdef([_LT_LANG_F77_CONFIG])dnl
   m4_popdef([_LT_LANG_FC_CONFIG])dnl
   m4_popdef([_LT_LANG_C_CONFIG])dnl
   m4_if(m4_cmp(m4_version_compare(LT_PACKAGE_VERSION,[2.4.6]),1),-1,
     [m4_popdef([_LT_COMPILER_PIC])m4_popdef([_LT_LINKER_SHLIBS])dnl
      m4_popdef([_LT_SYS_HIDDEN_LIBDEPS])])dnl
   dnl _KPSE_USE_LIBTOOL ensures libtool is also used for configure-time tests,
   dnl which deduces dependent libraries automatically
   _KPSE_USE_LIBTOOL
dnl substitute -shared-intel if present
   AS_FOR([acx_flag_var],[acx_flag_var_],[CFLAGS CXXFLAGS FCFLAGS F77FLAGS LDFLAGS FCLDFLAGS],
     [AS_IF([eval test x\$\{acx_flag_var+set\} = xset],
        [eval acx_temp="\" \$$acx_flag_var_ \""
         AS_CASE([$acx_temp],[*\ -shared-intel\ *|*\ -static-intel\ *],
           [acx_temp=`echo "$acx_temp" | sed -e 's/ \(-\(shared\|static\)-intel\)\b/ -Xcompiler \1 -XCClinker \1/g'`])
         AS_CASE([$acx_temp],[*\ -Qlocation,*\ *],
           [acx_temp=`echo "$acx_temp" | sed -e 's/ \(-Qlocation,@<:@^, @:>@*,@<:@^ @:>@*\)\b/ -Xcompiler \1 -XCClinker \1/g'`])
dnl take care of ifort/icc/icpc two-part options
         eval acx_flag_var=\"`echo "$acx_temp" | sed -e 's/ -\(align\|allow\|assume\|ccdefault\|check\|convert\|debug\|debug-parameters\|diag-type\|diag-enable\|diag-disable\|double-size\|dynamic-linker\|dyncom\|export-dir\|extend-source\|fp-model\|fpscomp\|gen-interfaces\|heap-arrays\|imacros\|integer-size\|iprefix\|iquote\|iwithprefixbefore\|module\|names\|opt-report\|opt-streaming-stores\|pch-dir\|pch-use\|prof-dir\|prof-file\|real-size\|reentrancy\|stand\|tcollect-filter\|tune\|warn\|watch\) \(@<:@^-@:>@@<:@^ @:>@*\)\b/ -Xcompiler -\1 -Xcompiler \2/g' -e 's/^ //;s/ $//'`\"])])
   AC_PROVIDE_IFELSE([AC_PROG_CC],
     [AS_IF([test -n "$CC" -a X"$CC" != Xno],
        [AC_LANG_PUSH([C])
         _KPSE_CHECK_LIBTOOL
         AC_LANG_POP([C])])])
   AC_PROVIDE_IFELSE([AC_PROG_FC],
     [AS_IF([test -n "$FC" -a X"$FC" != Xno],
        [AC_LANG_PUSH([Fortran])
         _ACX_LT_FORT_FLAGS_MANGLE
         _KPSE_CHECK_LIBTOOL
         AC_LANG_POP([Fortran])])])
   AC_PROVIDE_IFELSE([AC_PROG_F77],
     [AS_IF([test -n "$F77" -a X"$F77" != Xno],
        [AC_LANG_PUSH([Fortran 77])
         _ACX_LT_FORT_FLAGS_MANGLE
         _KPSE_CHECK_LIBTOOL
         AC_LANG_POP([Fortran 77])])])
   AC_PROVIDE_IFELSE([AC_PROG_CXX],
     [AS_IF([test -n "$CXX" -a X"$CXX" != Xno],
        [AC_LANG_PUSH([C++])
         _KPSE_CHECK_LIBTOOL
         AC_LANG_POP([C++])])])])dnl
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
