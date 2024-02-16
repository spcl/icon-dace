dnl Floating point precision checks
dnl
dnl This file contains configure-time checks for various floating
dnl point precision issues
dnl
dnl This file is released under public domain - or - in countries where this is
dnl not possible under the following license:
dnl
dnl    Permission is hereby granted, free of charge, to any person obtaining a
dnl    copy of this software, to deal in the software without restriction,
dnl    including without limitation the rights to use, copy, modify, merge,
dnl    publish, distribute, sublicense, and/or sell copies of the software,
dnl    and to permit persons to whom the software is furnished to do so, subject
dnl    to no condition whatsoever.
dnl
dnl    This software is provided AS IS, without warranty of any kind, express or
dnl    implied.

dnl FIXME:
dnl  This file was only built for x86 and x86_64 platforms but it does not check
dnl  for the platform since CMake does not provide a viable variable.

AC_DEFUN([CHECK_FLOAT_PRECISION],[
  acx_fpu_precision_control=unknown
  AC_MSG_CHECKING([for usable _FPU_SETCW])
  AC_LANG_CONFTEST([AC_LANG_SOURCE([[
@%:@include <stdio.h>
@%:@include <string.h>
@%:@include <fpu_control.h>

double div (double a, double b) {
  fpu_control_t fpu_oldcw, fpu_cw;
  volatile double result;

  _FPU_GETCW(fpu_oldcw);
  fpu_cw = (fpu_oldcw & ~_FPU_EXTENDED & ~_FPU_SINGLE) | _FPU_DOUBLE;
  _FPU_SETCW(fpu_cw);
  result = a / b;
  _FPU_SETCW(fpu_oldcw);
  return result;
}

int main (int argc, char **argv) {
  double d = div (2877.0, 1000000.0);
  char buf[255];
  sprintf(buf, "%.30f", d);
  fprintf(stderr, "%.30f\n", d);
  /* see if the result is actually in double precision */
  return strncmp(buf, "0.00287699", 10) == 0 ? 0 : 1;
}
]])])
  AC_RUN_IFELSE(,[ac_cfp_have__fpu_setcw=yes],[ac_cfp_have__fpu_setcw=no],
    [AC_LINK_IFELSE(,[ac_cfp_have__fpu_setcw=yes],[ac_cfp_have__fpu_setcw=no])])
  AS_IF([test "$ac_cfp_have__fpu_setcw" = "yes"],
    [acx_fpu_precision_control=yes
     AC_DEFINE([HAVE__FPU_SETCW],[1],[whether _FPU_SETCW is present and usable])
     AC_MSG_RESULT(yes)],
    [AC_MSG_RESULT(no)])

  AS_IF([test x"$acx_fpu_precision_control" = xunknown],
    [AC_MSG_CHECKING([for usable fpsetprec])
     AC_LANG_CONFTEST([AC_LANG_SOURCE([[
@%:@include <stdio.h>
@%:@include <string.h>
@%:@include <machine/ieeefp.h>

double div (double a, double b) {
  fp_prec_t fpu_oldprec;
  volatile double result;

  fpu_oldprec = fpgetprec();
  fpsetprec(FP_PD);
  result = a / b;
  fpsetprec(fpu_oldprec);
  return result;
}

int main (int argc, char **argv) {
  double d = div (2877.0, 1000000.0);
  char buf[255];
  sprintf(buf, "%.30f", d);
  /* see if the result is actually in double precision */
  return strncmp(buf, "0.00287699", 10) == 0 ? 0 : 1;
}
]])])
     AC_RUN_IFELSE(,[ac_cfp_have_fpsetprec=yes],[ac_cfp_have_fpsetprec=no],
       [AC_LINK_IFELSE(,[ac_cfp_have_fpsetprec=yes],
	  [ac_cfp_have_fpsetprec=no])])
     AS_IF([test x"$ac_cfp_have_fpsetprec" = xyes],
       [acx_fpu_precision_control=yes
	AC_DEFINE([HAVE_FPSETPREC],[1],
	  [whether fpsetprec is present and usable])
	AC_MSG_RESULT([yes])],
       [AC_MSG_RESULT([no])])
    ])

  AS_IF([test x"$acx_fpu_precision_control" = xunknown],
    [AC_MSG_CHECKING([for usable _controlfp])
     AC_LANG_CONFTEST([AC_LANG_SOURCE([[
@%:@include <stdio.h>
@%:@include <string.h>
@%:@include <float.h>

double div (double a, double b) {
  unsigned int fpu_oldcw;
  volatile double result;

  fpu_oldcw = _controlfp(0, 0);
  _controlfp(_PC_53, _MCW_PC);
  result = a / b;
  _controlfp(fpu_oldcw, _MCW_PC);
  return result;
}

int main (int argc, char **argv) {
  double d = div (2877.0, 1000000.0);
  char buf[255];
  sprintf(buf, "%.30f", d);
  /* see if the result is actually in double precision */
  return strncmp(buf, "0.00287699", 10) == 0 ? 0 : 1;
}
]])])
     AC_RUN_IFELSE(,[ac_cfp_have__controlfp=yes],[ac_cfp_have__controlfp=no],
       [AC_LINK_IFELSE(,[ac_cfp_have__controlfp=yes],
	  [ac_cfp_have__controlfp=no])])
     AS_IF([test x"$ac_cfp_have__controlfp" = xyes],
       [acx_fpu_precision_control=yes
	AC_DEFINE([HAVE__CONTROLFP],[1],[whether _controlfp is present usable])
	AC_MSG_RESULT([yes])],
       [AC_MSG_RESULT(no)])
   ])

  AS_IF([test x"$acx_fpu_precision_control" = xunknown],
    [AC_MSG_CHECKING([for usable _controlfp_s])
     AC_LANG_CONFTEST([AC_LANG_SOURCE([[
@%:@include <stdio.h>
@%:@include <string.h>
@%:@include <float.h>

double div (double a, double b) {
  unsigned int fpu_oldcw, fpu_cw;
  volatile double result;

  _controlfp_s(&fpu_cw, 0, 0);
  fpu_oldcw = fpu_cw;
  _controlfp_s(&fpu_cw, _PC_53, _MCW_PC);
  result = a / b;
  _controlfp_s(&fpu_cw, fpu_oldcw, _MCW_PC);
  return result;
}

int main (int argc, char **argv) {
  double d = div (2877.0, 1000000.0);
  char buf[255];
  sprintf(buf, "%.30f", d);
  /* see if the result is actually in double precision */
  return strncmp(buf, "0.00287699", 10) == 0 ? 0 : 1;
}
]])])
     AC_RUN_IFELSE(,[ac_cfp_have__controlfp_s=yes],
       [ac_cfp_have__controlfp_s=no],
       [AC_LINK_IFELSE(,[ac_cfp_have__controlfp_s=yes],
	  [ac_cfp_have__controlfp_s=no])])
     AS_IF([test $ac_cfp_have__controlfp_s = yes],
       [acx_fpu_precision_control=yes
	AC_DEFINE([HAVE__CONTROLFP_S],[1],
	  [whether _controlfp_s is present and usable])
	AC_MSG_RESULT([yes])],
       [AC_MSG_RESULT([no])])
   ])

  AS_IF([test x"$acx_fpu_precision_control" = xunknown],
    [AC_MSG_CHECKING([whether FPU control word can be manipulated by inline assembler])
     AC_LANG_CONFTEST([AC_LANG_SOURCE([[
@%:@include <stdio.h>
@%:@include <string.h>

double div (double a, double b) {
  unsigned int oldcw, cw;
  volatile double result;

  __asm__ __volatile__ ("fnstcw %0" : "=m" (*&oldcw));
  cw = (oldcw & ~0x0 & ~0x300) | 0x200;
  __asm__ __volatile__ ("fldcw %0" : : "m" (*&cw));

  result = a / b;

  __asm__ __volatile__ ("fldcw %0" : : "m" (*&oldcw));

  return result;
}

int main (int argc, char **argv) {
  double d = div (2877.0, 1000000.0);
  char buf[255];
  sprintf(buf, "%.30f", d);
  /* see if the result is actually in double precision */
  return strncmp(buf, "0.00287699", 10) == 0 ? 0 : 1;
}
]])])
     AC_RUN_IFELSE(,[ac_cfp_have_fpu_inline_asm_x86=yes],
       [ac_cfp_have_fpu_inline_asm_x86=no],
       [AC_LINK_IFELSE(,[ac_cfp_have_fpu_inline_asm_x86=yes],
	  [ac_cfp_have_fpu_inline_asm_x86=no])])
     AS_IF([test $ac_cfp_have_fpu_inline_asm_x86 = yes],
       [acx_fpu_precision_control=yes
	AC_DEFINE([HAVE_FPU_INLINE_ASM_X86],[1],
	  [whether FPU control word can be manipulated by inline assembler])
	AC_MSG_RESULT([yes])],
       [AC_MSG_RESULT([no])])
    ])
])

AC_DEFUN([ACX_CHECK_FLOAT_UNDERFLOW],
  [AC_MSG_CHECKING([for access to FTZ/DAZ control bits by inline assembler])
   save_LIBS=$LIBS
   LIBS="-lm $LIBS"
   AC_LANG_CONFTEST([AC_LANG_PROGRAM([
@%:@include <float.h>
@%:@include <inttypes.h>
@%:@include <math.h>
@%:@include <stdlib.h>
@%:@include <stdio.h>
@%:@include <string.h>

enum {
  PPM_FTZ_BIT = 15,
  PPM_DM_BIT = 8,
  PPM_DAZ_BIT = 6
};

@%:@pragma float_control (precise,on)
@%:@pragma fenv_access (on)

@%:@ifndef __GNUC__
@%:@  define  __attribute__(x)  /*NOTHING*/
@%:@endif
void add(double a, double b, double *result) {
  uint32_t oldmxcsr, mxcsr;

  __asm__ __volatile__ ("stmxcsr %0" : "=m" (*&oldmxcsr));
  /* get abrupt underflow result first */
  mxcsr = (oldmxcsr | 1 << PPM_FTZ_BIT | 1 << PPM_DAZ_BIT | 1 << PPM_DM_BIT);
  __asm__ __volatile__ ("ldmxcsr %0" : : "m" (*&mxcsr));
  /* result@<:@0@:>@ = a + b; with abrupt underflow */
  __asm__ __volatile__ ("ldmxcsr %1\n\t"
                        "vaddsd %2, %3, %%xmm2\n\t"
                        "vmovsd %%xmm2, %0"
                        : "=m" (result@<:@0@:>@) : "m" (*&mxcsr), "x"(a), "x"(b) : "xmm2");
  /* next compute denormal result */
  mxcsr = (oldmxcsr & ~(1 << PPM_FTZ_BIT | 1 << PPM_DAZ_BIT)) | (1 << PPM_DM_BIT);
  /* result@<:@1@:>@ = a + b; with denormal created */
  __asm__ __volatile__ ("ldmxcsr %1\n\t"
                        "vaddsd %2, %3, %%xmm2\n\t"
                        "vmovsd %%xmm2, %0"
                        : "=m" (result@<:@1@:>@): "m" (*&mxcsr), "x"(a), "x"(b) : "xmm2");
  __asm__ __volatile__ ("ldmxcsr %0" : : "m" (*&oldmxcsr));
}],
[  double d@<:@2@:>@, a = scalbln(-7.5, DBL_MIN_EXP),
    b = scalbln(7.75, DBL_MIN_EXP);
  add(a, b, d);
  char buf@<:@40@:>@;
  fprintf(stderr, "%.17e, %.17e\n", d@<:@0@:>@, d@<:@1@:>@);
  sprintf(buf, "%.17e", d@<:@0@:>@);
  int okay = !strcmp(buf, "0.00000000000000000e+00");
  sprintf(buf, "%.17e", d@<:@1@:>@);
  okay &= !strcmp(buf, "1.11253692925360069e-308");
  /* see if the result is actually in double precision */
  return okay ? EXIT_SUCCESS : EXIT_FAILURE;
])])
   AC_RUN_IFELSE(,[acx_cfp_have_mxcsr_inline_asm_x86=yes],
     [acx_cfp_have_mxcsr_inline_asm_x86=no],
     [AC_LINK_IFELSE(,[acx_cfp_have_mxcsr_inline_asm_x86=yes],
	[acx_cfp_have_mxcsr_inline_asm_x86=no])])
   LIBS=$save_LIBS
   AS_IF([test $acx_cfp_have_mxcsr_inline_asm_x86 = yes],
     [acx_fpu_underflow_control=yes
      AC_DEFINE([HAVE_MXCSR_INLINE_ASM_X86],[1],
	[whether FTZ/DAZ control bits can be manipulated by inline assembler])],
     [acx_fpu_underflow_control=no])
   AC_MSG_RESULT([$acx_cfp_have_mxcsr_inline_asm_x86])
  ])

dnl _ACX_TEST_KAHAN_SUM(ACTION-IF-WORKING,ACTION-IF-NON-WORKING)
AC_DEFUN([_ACX_TEST_KAHAN_SUM],
  [AC_REQUIRE([ACX_MV_OBJ])
   AS_IF([test $cross_compiling = no],
     [AC_LANG_PUSH([C])
      _AC_RUN_LOG([cat confdefs.h \
	  "$srcdir/src/core/ppm_math_extensions_ddp_c.c" \
	  >conftest.c],[_AS_ECHO_LOG([Creating source file.])])
      CPPFLAGS="$CPPFLAGS -I$srcdir/src/core -I$srcdir/include"
      AC_COMPILE_IFELSE(,
	[ACX_MV_OBJ([conftest],[conftest_ddp])
	 save_LIBS=$LIBS
	 LIBS="conftest_ddp.$ac_objext $LIBS"
	 AC_RUN_IFELSE([AC_LANG_PROGRAM(
	   [@%:@include <complex.h>
@%:@include <float.h>
@%:@include <math.h>
@%:@include <stdio.h>
@%:@include <stdlib.h>

double complex
PPM_ddp_sum_dp(size_t n, const double *a);
double complex
PPM_ddp_add_ddp_ddp(double complex a, double complex b);],
	   [static const double a@<:@@:>@ = { 0.65749795994724536E-305,
 0.48298396135057048E-308,
-0.65749795994724536E-305,
-0.48298396135057048E-308,
 0.51266255288712248E-305,
 0.0000000000000000,
 1.11253692925360069E-308,
 1.11253692925360069E-308,
-2.2250738585072014E-308,
-0.51266255288712248E-305 };

enum { asize = sizeof(a)/sizeof(a@<:@0@:>@) };

static const double b@<:@@:>@ = {
  0.73633538217922512,
  2.05676803375792914E-016,
 -0.73633538217922512,
 -2.05676803375792914E-016,
  0.67448228723532633,
  0.0000000000000000,
  4.33183514517286229E-016,
 -0.67448228723532633,
 -4.33183514517286229E-016 };

enum { bsize = sizeof(b)/sizeof(b@<:@0@:>@) };

static const double c@<:@@:>@ = {
  0.61825480456471904,
  2.75858288416734245E-016,
 -0.61825480456471904,
 -2.75858288416734245E-016,
  0.59473487038421935,
  0.0000000000000000,
  6.18759993964230386E-016,
 -0.59473487038421935,
 -6.18759993964230386E-016 };

enum { csize = sizeof(c)/sizeof(c@<:@0@:>@) };

static const double d@<:@@:>@ = {
-0.60013302811469493,
-8.3011177257981544e-17,
4.5323096146354345e-17,
-3.8125690077553134e-17,
0.73029831751529928,
8.9590614262320777e-17,
-0.73029831751529928,
-8.9590614262320777e-17,
3.8125690077553134e-17,
-7.7647018569815587e-17,
0.62334907629392666,
8.0275708681884551e-17,
-0.62334907629392666,
-8.0275708681884551e-17,
7.7647018569815587e-17,
-1.5532073730149895e-17,
0.68889771851855675,
5.6229556412310705e-17,
-0.68889771851855675,
-5.6229556412310705e-17,
1.5532073730149895e-17,
-1.5959170195851886e-17,
0.87746271427488143,
6.1191266179381476e-17,
-0.87746271427488143,
-6.1191266179381476e-17,
1.5959170195851886e-17,
1.5959170195851886e-17,
   };
enum { dsize = sizeof(d)/sizeof(d@<:@0@:>@) };

static const double e@<:@@:>@ = {
-4.5323096146354345e-17,
0.60013302811469493,
8.3011177257981544e-17,
};
enum { esize = sizeof(e)/sizeof(e@<:@0@:>@) };

double complex s = PPM_ddp_sum_dp(asize, a);

int passed = 1, passed_a = fabs(creal(s)) <= DBL_MIN;
passed &= passed_a;
if (!passed_a)
  fprintf(stderr, "Bad result from Kahan summation test a, "
	  "should have been zero: %24.17g\n", creal(s));

s = PPM_ddp_sum_dp(bsize, b);
int passed_b = fabs(creal(s)) <= DBL_MIN;
passed &= passed_b;
if (!passed_b)
  fprintf(stderr, "Bad result from Kahan summation test b, "
	  "should have been zero: %24.17g\n", creal(s));

s = PPM_ddp_sum_dp(csize, c);
int passed_c = fabs(creal(s)) <= DBL_MIN;
passed &= passed_c;
if (!passed_c)
  fprintf(stderr, "Bad result from Kahan summation test c, "
	  "should have been zero: %24.17g\n", creal(s));

s = PPM_ddp_add_ddp_ddp(PPM_ddp_sum_dp(dsize, d),
	                PPM_ddp_sum_dp(esize, e));
int passed_d = fabs(creal(s) - d@<:@dsize-1@:>@) < DBL_MIN;
passed &= passed_d;
if (!passed_d)
  fprintf(stderr, "Bad result from Kahan summation test d, "
	  "should have been close to %24.17g but is: %24.17g\n",
	  d@<:@dsize-1@:>@, creal(s));

if (!passed)
  return EXIT_FAILURE;
])],
	  [$1],[$2])
      LIBS=$save_LIBS])
      AC_LANG_POP([C])])])

dnl To yield the behaviour required for working Kahan-summation,
dnl
dnl 1. defines NEED_PRECISION_CONTROL if the FPU requires explicit
dnl    control of precision of evaluations, and
dnl 2. defines NEED_UNDERFLOW_CONTROL if the FPU must be configured
dnl    to prevent replacement of denormal quantities by zero
dnl
dnl PRECISION-CONTROL-EXTRA-ACTION extra code to execute in case
dnl NEED_PRECISION_CONTROL is defined
dnl
dnl UNDERFLOW-CONTROL-EXTRA-ACTION extra code to execute in case
dnl NEED_UNDERFLOW_CONTROL is defined
dnl
dnl
dnl ACX_CHECK_KAHAN_SUMMATION_SETTINGS([PRECISION-CONTROL-EXTRA-ACTION],
dnl                                    [UNDERFLOW-CONTROL-EXTRA-ACTION])
AC_DEFUN([ACX_CHECK_KAHAN_SUMMATION_SETTINGS],
  [CHECK_FLOAT_PRECISION
   ACX_CHECK_FLOAT_UNDERFLOW
   AC_CACHE_CHECK([for working Kahan summation],[acx_cv_kahan_summation_ok],
     [acx_cv_kahan_summation_ok=no
      AS_CASE([x"$acx_fpu_precision_control$acx_fpu_underflow_control"],
	[xyesyes],[acx_temp='oo ox xo xx'],
	[xyesno],[acx_temp='oo xo'],
	[xnoyes],[acx_temp='oo ox'],
	[acx_temp='oo'])
      save_CPPFLAGS=$CPPFLAGS
      AS_FOR([ACX_Temp],[acx_temp],[$acx_temp],
	[AS_CASE([ACX_Temp],[x?],
	   [CPPFLAGS="$save_CPPFLAGS -DNEED_PRECISION_CONTROL"])
	 AS_CASE([ACX_Temp],[?x],
	   [CPPFLAGS="$save_CPPFLAGS -DNEED_UNDERFLOW_CONTROL"])
	 _ACX_TEST_KAHAN_SUM([acx_cv_kahan_summation_ok=yes ; LIBS=$save_LIBS; break])])
      CPPFLAGS=$save_CPPFLAGS])
   AS_IF([test $acx_cv_kahan_summation_ok = yes],
     [AS_CASE([$acx_temp],[x?],
	[AC_DEFINE([NEED_PRECISION_CONTROL],[1],
	   [defined if precision control is needed])m4_ifval([$1],[$1])])
      AS_CASE([$acx_temp],[?x],
	[AC_DEFINE([NEED_UNDERFLOW_CONTROL],[1],
	   [defined if underflow control is needed])m4_ifval([$2],[$2])])])])

dnl Local Variables:
dnl mode: autoconf
dnl End:
