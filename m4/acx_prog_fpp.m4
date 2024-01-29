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

# ACX_PROG_FPP([VARIABLE],
#              [ACTION-IF-SUCCESS],
#              [ACTION-IF-FAILURE = FAILURE],
#              [INCLUDE-FLAG = -I],
#              [MACRO-FLAG = -D],
#              [WRAPPER-SUBDIR = fpp-wrappers])
# -----------------------------------------------------------------------------
# Searches for the Fortran preprocessor command emitting its output to the
# standard output. If the argument VARIABLE is given and the shell variable
# VARIABLE is set, checks whether the value it stores is a valid Fortran
# preprocessor command. If VARIABLE is not set, iterates over the known
# commands (including the shell wrappers expected to reside in the
# WRAPPER-SUBDIR directory of ${srcdir}) and stops when the first valid command
# is found. The value of VARIABLE is never set or changed. The test requires
# the values of the compiler flag INCLUDE-FLAG (defaults to -I) needed to
# specify search paths for the quoted form of the preprocessor "#include"
# directive, and the compiler flag MACRO-FLAG (defaults to -D) needed to
# specify a preprocessor macro definition.
#
# Note that some compilers cannot output the result of the preprocessing to
# the standard output. For that reason, this macro is provided with a set of
# wrappers, which are supposed to circumvent this problem in some cases. You
# can find the wrappers in the 'fpp-wrappers' subdirectory residing in the
# route source directory of this project. If you decide to use this macro, copy
# the wrappers to any subdirectory of your project and specify the relative
# path of the directory containing them as the WRAPPER-SUBDIR argument of this
# macro (defaults to fpp-wrappers).
#
# It should be also noted that the aforementioned wrappers are written and
# maintained to work with the following limitations:
#   1) the list of input files for the preprocessing is the longest sequence of
#      arguments, each element of which is a path to an existing and readable
#      file, that is found in the end of the command line;
#   2) the input files must not have the same basename, e.g. './test.f90',
#      './somedir/test.f90'.
#
# An example of the full preprocessing command in a Makefile (assuming the
# variable FPP is set to the result of this macro):
#
# %.f90: %.F90
#   $(FPP) $(FCFLAGS) $< >$@
#
# Supported compilers:
#   gfortran: "$FC -E $ac_fcflags_srcext" (the variable $ac_fcflags_srcext
#             keeps the flag enabling the preprocessing at the compile time,
#             it either holds '-cpp' or is empty depending on whether the flag
#             '-cpp' was specified in the FCFLAGS by the user when calling the
#             configure script);
#
#   Intel: "$FC -E $ac_fcflags_srcext" (the same as for gfortran but the
#          variable $ac_fcflags_srcext is usually set to '-fpp');
#
#   NAGWare: "$SHELL ${acx_prog_fpp_wrapper_dir}/nag.sh $FC" or
#            "$FC -o - -Wp,-w,-P -F $ac_fcflags_srcext" (see comments in the
#             wrapper);
#
#   Cray: "$SHELL ${acx_prog_fpp_wrapper_dir}/cray.sh $FC" (uses the
#         wrapper);
#
#   PGI: "$FC -E" (the compiler does not need $ac_fcflags_srcext when called
#        with the flag '-E');
#
#   A chance for not yet supported compilers: "$FC -F" and
#                                             "$FC -F $ac_fcflags_srcext".
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# A positive result of this test is cached in the acx_cv_prog_fpp variable.
#
AC_DEFUN([ACX_PROG_FPP],
  [AC_LANG_ASSERT([Fortran])dnl
   AC_PROVIDE_IFELSE([AC_FC_PP_SRCEXT], [],
     [m4_warn([syntax],
        [ACX_PROG_FPP requires calling the Fortran compiler with a ]dnl
[preprocessor but no call to AC_FC_PP_SRCEXT is detected])])dnl
   AC_MSG_CHECKING([for Fortran preprocessor command])
   AC_CACHE_VAL([acx_cv_prog_fpp],
     [AS_MKDIR_P([conftest.dir/sub])
      cd conftest.dir
      cat > sub/conftest.fpp_inc <<_ACEOF
@%:@define CONFTEST_THREE
@%:@ifndef CONFTEST_ZERO
      integer conftest_zero
@%:@else
      choke me
@%:@endif
@%:@ifdef CONFTEST_ONE
      integer conftest_one
@%:@else
      choke me
@%:@endif
@%:@if CONFTEST_TWO != 42
      choke me
@%:@else
      integer conftest_two
@%:@endif
@%:@ifdef CONFTEST_THREE
      integer conftest_three
@%:@else
      choke me
@%:@endif

_ACEOF
      acx_save_FCFLAGS=$FCFLAGS
      AS_VAR_SET_IF([$1], [set dummy "AS_VAR_GET([$1])"],
        [_AC_SRCDIRS([.])
         acx_prog_fpp_wrapper_dir=dnl
"${ac_abs_top_srcdir}/m4_default([$6], [fpp-wrappers])"
         set dummy "$FC -E" "$FC -E $ac_fcflags_srcext" \
                   "$SHELL ${acx_prog_fpp_wrapper_dir}/nag.sh $FC" \
                   "$FC -o - -Wp,-w,-P -F $ac_fcflags_srcext" \
                   "$SHELL ${acx_prog_fpp_wrapper_dir}/cray.sh $FC" \
                   "$FC -F" "$FC -F $ac_fcflags_srcext"])
      shift
      acx_prog_fpp_include_flag=m4_default([$4], [-I])
      acx_prog_fpp_macro_flag=m4_default([$5], [-D])
      for acx_candidate in "$[@]"; do
        AC_LANG_CONFTEST([AC_LANG_PROGRAM([],[      implicit none
@%:@include "conftest.fpp_inc"
      conftest_zero = 0
      conftest_one = 1
      conftest_two = 2
      conftest_three = 3])])
        rm -f conftest.acx_prog_fpp
        acx_try_fpp="$acx_candidate \$FCFLAGS conftest.\$ac_ext dnl
>conftest.acx_prog_fpp"
        FCFLAGS="${acx_prog_fpp_include_flag}sub dnl
${acx_prog_fpp_macro_flag}CONFTEST_ONE dnl
${acx_prog_fpp_macro_flag}CONFTEST_TWO=42 $acx_save_FCFLAGS"
        _AC_DO_VAR([acx_try_fpp])
        AS_IF([test $? -eq 0 && test -f conftest.acx_prog_fpp],
          [acx_result=`grep -c 'integer conftest' dnl
conftest.acx_prog_fpp 2>/dev/null`
           AS_IF([test $? -eq 0 && test 4 -eq "$acx_result" 2>/dev/null],
             [mv conftest.acx_prog_fpp conftest.$ac_ext
              FCFLAGS=$acx_save_FCFLAGS
              AC_COMPILE_IFELSE([],
                [acx_cv_prog_fpp=$acx_candidate
                 break])])])
      done
      FCFLAGS=$acx_save_FCFLAGS
      cd ..
      rm -rf conftest.dir])
   AS_VAR_SET_IF([acx_cv_prog_fpp],
     [AC_MSG_RESULT([$acx_cv_prog_fpp])
      $2],
     [AC_MSG_RESULT([unknown])
      m4_default([$3],
        [AC_MSG_FAILURE(
           [unable to find a valid Fortran preprocessor command])])])])
