dnl acx_prog_fc_test_fpp.m4 --- test Fortran compiler is preprocessing
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
# ACX_PROG_FC_CHECK_FPP([ACTION-IF-FC-DOES-PREPROCESS],[ACTION-IF-FC-DOES-NOT-PREPROCESS])
# -------------------------
# check if Fortran source by name .f90 is preprocessed
AC_DEFUN([ACX_PROG_FC_CHECK_FPP],dnl
  [ACX_ASSERT_LANG_IS_FORTRAN_VARIANT
   m4_case(_AC_LANG,[Fortran],[AC_REQUIRE([AC_PROG_FC])],[Fortran 77],[AC_REQUIRE([AC_PROG_F77])])dnl
   AS_VAR_PUSHDEF([fpp_flag],[acx_cv_fc_fpp_flag_]_AC_LANG_ABBREV)dnl
   ASX_VAR_UNSET([fpp_flag])
   m4_pushdef([acx_flags_var],m4_case(_AC_LANG,[Fortran],[FCFLAGS],[Fortran 77],[FFLAGS]))dnl
   AC_CACHE_CHECK([for flag to enable preprocessing],[fpp_flag],
     []acx_flags_var[_save=$]acx_flags_var[
      for i in none -cpp -fpp -qsuffix=cpp=f90 -eT -eZ -Mpreprocess -x\ f95-cpp-input ; do
	AS_IF([test "x$i" != xnone],
	  []acx_flags_var[="$]acx_flags_var[ ${i}"])
	AC_COMPILE_IFELSE([_ACX_SL_LANG_PROGRAM_FPP_ONLY],
	  [AS_VAR_SET([fpp_flag],[$i]) ; break])
	]acx_flags_var[=$]acx_flags_var[_save
      done
      ]acx_flags_var[=$]acx_flags_var[_save
      AS_IF([test "x$i" = xnone],
	[AS_VAR_SET([fpp_flag],[])])
      AS_VAR_SET_IF([fpp_flag],
	[m4_default([$1],[:])],
	[m4_default([$2],[AC_MSG_ERROR([Cannot automatically find flag to preprocess Fortran sources.
Please consult your compiler documentation for a flag to make $FC
preprocess files with .$ac_ext suffix and add it to ]acx_flags_var[.
Also please consider reporting the flag/compiler combination to Thomas Jahns <jahns@dkrz.de>])])])])
   FC_FPP_FLAG=AS_VAR_GET([fpp_flag])
   m4_popdef([acx_flags_var])dnl
   AC_SUBST([FC_FPP_FLAG])
   AS_VAR_POPDEF([fpp_flag])dnl
])# ACX_PROG_FC_CHECK_FPP
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
