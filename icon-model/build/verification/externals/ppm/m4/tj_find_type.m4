dnl tj_find_type.m4 --- find built-in type for a typedef
dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
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
dnl Commentary:
dnl
dnl
dnl
dnl Code:
dnl
dnl This code is very similar to that found in the pthreads package
dnl from Christopher Provenzano, I copied it from the pmpthreads 1.8.9
dnl distribution. Thomas Jahns 2002/07/25
dnl
dnl The following has been omitted from the code:
dnl
dnl  define a C macro to expand to the type identified with
dnl  PTHREADS_FIND_TYPE
dnl
dnl because it's of so little value to have a Preprocessor macro
dnl duplicate a typedef
dnl
dnl TJ_FIND_TYPE:
dnl Determine proper typedef value for a typedef name and create a
dnl shell variable with that value .  If none of the specified types
dnl to try match, the macro is left undefined, and the shell variable
dnl empty.  If the typedef name cannot be found in the specified
dnl header files, this test errors out; perhaps it should be changed
dnl to simply leave the macro undefined...
dnl
dnl TJ_FIND_TYPE(typedefname,varname,includes,possible values...)
dnl
dnl
dnl  Usage: TJ_FIND_TYPE(system-typedef-name, new-macro-name,
dnl		         includes,
dnl			 list-of-types-to-try,
dnl                      action-if-found,
dnl                      action-if-failed)
dnl
dnl  TJ_FIND_INTEGRAL_TYPE automatically provides a set of integral
dnl  types, and does not permit specification of additional types.
dnl
dnl  The specified types must all be able to work as prefixes -- i.e., no
dnl  direct specification of array or function types.  If you need such
dnl  types, add typedefs for them to include/pthread/xtypes.h, and include
dnl  that in the set of header files.  For simple struct types, you can
dnl  try including the definition directly here, but it had better not
dnl  contain any commas or square brackets.
dnl
dnl  If necessary, you can include other preprocessing commands and such
dnl  in the `includes' portion.
dnl
dnl  Note:  For now, each of these needs a corresponding entry
dnl  in acconfig.h.
dnl
dnl
AC_DEFUN([TJ_FIND_TYPE],
  [AS_VAR_PUSHDEF([tj_type_equiv],[tj_cv_c_type_$1])dnl
   AC_CACHE_CHECK([type of $1],[tj_type_equiv],
     [AC_LANG_PUSH([C])dnl
      AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([$3], [extern $1 foo; ])],
        [AS_FOR([TRY_TYPE],[try_type],[$4],
           [AC_COMPILE_IFELSE(
             [AC_LANG_PROGRAM([$3],
               [extern $1 foo; extern $try_type foo;])],
             [AS_VAR_SET([tj_type_equiv],[$try_type]); break])])])
       AC_LANG_POP([C])])
   AS_VAR_SET_IF([tj_type_equiv],
     [$5],
     [AC_MSG_ERROR([Cannot find matching typedef for $1.])])
   m4_ifval([$2],[AS_VAR_COPY([$2],[tj_type_equiv])])
   AS_VAR_POPDEF([tj_type_equiv])])
dnl
dnl
dnl
dnl Like above, but the list of types to try is pre-specified.
dnl
AC_DEFUN([TJ_FIND_INTEGRAL_TYPE],
  [TJ_FIND_TYPE([$1], [$2], [$3],
	[int 'unsigned int' long 'unsigned long' short 'unsigned short' char 'unsigned char' 'long long' 'unsigned long long'])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
