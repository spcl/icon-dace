dnl acx_fc_real_edit_descriptors.m4 --- check what edit descriptor
dnl                                     field widths to use for the kinds of
dnl                                     REAL used in Fortran programs
dnl
dnl Copyright  (C)  2018  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Keywords: configure configure.ac autoconf MPI mpirun mpiexec
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
AC_DEFUN([ACX_C_FLOAT_DIGITS],
  [dnl
dnl find number of decimal digits for string conversions
dnl 1. export C float.h limits to shell variables
AS_FOR([acx_c_limit_],[acx_c_limit],[FLT_RADIX FLT_MANT_DIG FLT_MAX_EXP FLT_MIN_EXP DBL_MANT_DIG DBL_MAX_EXP DBL_MIN_EXP],
  [AC_COMPUTE_INT([acx_c_limit_],[acx_c_limit_],[@%:@include <float.h>
@%:@include <math.h>],[AC_MSG_FAILURE([Cannot determine value of acx_c_limit_])])])
dnl 2. use perl to evaluate expression (C cannot compute log10 at compile time)
AC_MSG_CHECKING([number of decimal mantissa digits for identity conversions of C expressions of type float])
PPM_FLT_DECIMAL_DIG=`perl -e "use strict;
use warnings;
use POSIX qw(ceil log10);
print ceil(1.0 + ${FLT_MANT_DIG} * log10(${FLT_RADIX}));"`
AC_MSG_RESULT([$PPM_FLT_DECIMAL_DIG])
AC_DEFINE_UNQUOTED([PPM_CONF_FLT_DECIMAL_DIG],[$PPM_FLT_DECIMAL_DIG],[number of decimal digits for identity conversions of C expressions of type float])

AC_MSG_CHECKING([number of decimal exponent digits for identity conversions of C expressions of type float])
PPM_FLT_EXP_DECIMAL_DIG=`perl -e "use strict;
use warnings;
use POSIX qw(ceil log10);
sub max(\\\$\\\$) { my (\\\$x, \\\$y) = @_; return \\\$x > \\\$y ? \\\$x : \\\$y; }
print ceil(log10(max(abs(${FLT_MIN_EXP}),${FLT_MAX_EXP}) * log10(${FLT_RADIX})));"`
AC_MSG_RESULT([$PPM_FLT_EXP_DECIMAL_DIG])
AC_DEFINE_UNQUOTED([PPM_CONF_FLT_EXP_DECIMAL_DIG],[$PPM_FLT_EXP_DECIMAL_DIG],[number of decimal exponent digits for identity conversions of C expressions of type float])

AC_MSG_CHECKING([number of decimal mantissa digits for identity conversions of C expressions of type double])
PPM_DBL_DECIMAL_DIG=`perl -e "use strict;
use warnings;
use POSIX qw(ceil log10);
print ceil(1.0 + ${DBL_MANT_DIG} * log10(${FLT_RADIX}));"`
AC_MSG_RESULT([$PPM_DBL_DECIMAL_DIG])
AC_DEFINE_UNQUOTED([PPM_CONF_DBL_DECIMAL_DIG],[$PPM_DBL_DECIMAL_DIG],[number of decimal digits for identity conversions of C expressions of type double])

AC_MSG_CHECKING([number of decimal exponent digits for identity conversions of C expressions of type double])
PPM_DBL_EXP_DECIMAL_DIG=`perl -e "use strict;
use warnings;
use POSIX qw(ceil log10);
sub max(\\\$\\\$) { my (\\\$x, \\\$y) = @_; return \\\$x > \\\$y ? \\\$x : \\\$y; }
print ceil(log10(max(abs(${DBL_MIN_EXP}),${DBL_MAX_EXP}) * log10(${FLT_RADIX})));"`
AC_MSG_RESULT([$PPM_DBL_EXP_DECIMAL_DIG])
AC_DEFINE_UNQUOTED([PPM_CONF_DBL_EXP_DECIMAL_DIG],[$PPM_DBL_EXP_DECIMAL_DIG],[number of decimal exponent digits for identity conversions of C expressions of type float])
])

dnl ACX_FC_REAL_EDIT_DESCRIPTORS([PREFIX])
dnl Set PREFIX_de_g_sp, PREFIX_de_g_sp_width,
dnl PREFIX_de_g_dp,     PREFIX_de_g_dp_width
dnl to the edit descriptors and corresponding character length to
dnl print REAL variables of kinds sp and dp respectively in loss-less
dnl precision. PREFIX is optional.
AC_DEFUN([ACX_FC_REAL_EDIT_DESCRIPTORS],
  [AC_REQUIRE([ACX_C_FLOAT_DIGITS])
   AC_CACHE_CHECK([for Fortran sp and dp output field widths],
     [acx_cv_fc_real_edit_descriptors],
     [AS_IF([test $cross_compiling = yes],
        [acx_tmp1=`expr 2 + $PPM_FLT_DECIMAL_DIG + 2 + $PPM_FLT_EXP_DECIMAL_DIG + 1`
         acx_tmp2="${acx_tmp1}.$PPM_FLT_DECIMAL_DIG"
         AS_IF([test $PPM_FLT_EXP_DECIMAL_DIG -gt 2],
           [acx_tmp2="${acx_tmp2}e${PPM_FLT_EXP_DECIMAL_DIG}"])
         acx_tmp3=`expr 2 + $PPM_DBL_DECIMAL_DIG + 2 + $PPM_DBL_EXP_DECIMAL_DIG + 1`
         acx_tmp4="${acx_tmp3}.$PPM_DBL_DECIMAL_DIG"
         AS_IF([test $PPM_DBL_EXP_DECIMAL_DIG -gt 2],
           [acx_tmp4="${acx_tmp4}e${PPM_DBL_EXP_DECIMAL_DIG}"])
         acx_cv_fc_real_edit_descriptors="de_g_sp='g${acx_tmp2}' de_g_sp_width='${acx_tmp1}' de_g_dp='g${acx_tmp4}' de_g_dp_width=${acx_tmp3}"
         m4_ifval([$1],[acx_cv_fc_real_edit_descriptors=`echo "$acx_cv_fc_real_edit_descriptors" | sed -e 's/de_g_/$1_&/g'`])],
        [AC_LANG_PUSH([Fortran])
         acx_temp=`echo "$[]m4_default([$1],[PPM])_FC_FEATURE_DEFS" | sed -e '$b
s/$/\\\\/'`
         sed -e '/^ppm_real_sp_dp_edit_descriptor\.f90$/{
r '"$srcdir/src/core/ppm_real_sp_dp_edit_descriptor.f90"'
d
}
/^ppm_real_sp_dp.inc$/{
r '"$srcdir/src/core/ppm_real_sp_dp.inc"'
d
}' <<_ACEOF | sed -e '/^ *!@<:@ !>@:>@/d' -e '/^ *!$/d' \
   -e '/^@%:@include "fc_feature_defs.inc"$/{
i '"$acx_temp"'
d
}' \
   >conftest.${ac_fc_srcext-f}
      MODULE ppm_std_type_kinds
      IMPLICIT NONE
      PUBLIC
ppm_real_sp_dp.inc
      END MODULE ppm_std_type_kinds
ppm_real_sp_dp_edit_descriptor.f90
      PROGRAM conftest
        USE ppm_real_sp_dp_edit_descriptor, ONLY: &
        get_edit_descriptor_sp,  get_edit_descriptor_dp
        !> data edit descriptors for real kinds sp and dp
        CHARACTER(len=10) :: de_g_sp, de_g_dp
        !> character width needed for the corresponding data edit descriptor
        INTEGER :: de_g_sp_width, de_g_dp_width
        CALL get_edit_descriptor_sp(de_g_sp, de_g_sp_width)
        PRINT '(3a)', "de_g_sp='", TRIM(de_g_sp), "'"
        PRINT '(a,i0)', 'de_g_sp_width=', de_g_sp_width
        CALL get_edit_descriptor_dp(de_g_dp, de_g_dp_width)
        PRINT '(3a)', "de_g_dp='", TRIM(de_g_dp), "'"
        PRINT '(a,i0)', 'de_g_dp_width=', de_g_dp_width
      END PROGRAM conftest
_ACEOF
         AC_LINK_IFELSE(,
           [acx_run="./conftest$EXEEXT"
            AS_IF([expr "$ac_link" : '.*/libtool --mode=link' >/dev/null],
              [acx_run=`echo "$ac_link" | sed -e 's@\(.*/libtool --mode=\)link.*@\1@'`"execute $acx_run"])
	    acx_cv_fc_real_edit_descriptors=`$acx_run \
                | sed  -e ':a' -e 'N;$bt' -e 'ba' -e ':t' -e 's/\n/ /g' m4_ifval([$1],[-e 's/de_g_/$1_&/g'])`
            AS_IF([test $ac_fc_mod_uppercase = yes],
              [rm PPM_REAL_SP_DP_EDIT_DESCRIPTOR.$FCMODEXT PPM_STD_TYPE_KINDS.$FCMODEXT],
              [rm ppm_real_sp_dp_edit_descriptor.$FCMODEXT ppm_std_type_kinds.$FCMODEXT])])
         AC_LANG_POP([Fortran])])])
   AS_IF([test $cross_compiling = yes],
     [AC_MSG_NOTICE([derived Fortran REAL output field widths from C.])])
   AS_VAR_SET_IF([acx_cv_fc_real_edit_descriptors],
     [eval "${acx_cv_fc_real_edit_descriptors}"],
     [m4_default([$2],
        [AC_MSG_FAILURE([failed setup of real edit descriptors])])])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
