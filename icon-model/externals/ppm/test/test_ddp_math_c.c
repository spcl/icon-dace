/**
 * @file test_ddp_math_c.c
 * @brief test specific problem cases for Kahan summation from C
 *
 * @copyright Copyright  (C)  2010-2017  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/ppm_math_extensions.h"

int
main ()
{
static const double a[] = { 0.65749795994724536E-305,
 0.48298396135057048E-308,
-0.65749795994724536E-305,
-0.48298396135057048E-308,
 0.51266255288712248E-305,
 0.0000000000000000,
 1.11253692925360069E-308,
 1.11253692925360069E-308,
-2.2250738585072014E-308,
-0.51266255288712248E-305 };

enum { asize = sizeof(a)/sizeof(a[0]) };

static const double b[] = {
  0.73633538217922512,
  2.05676803375792914E-016,
 -0.73633538217922512,
 -2.05676803375792914E-016,
  0.67448228723532633,
  0.0000000000000000,
  4.33183514517286229E-016,
 -0.67448228723532633,
 -4.33183514517286229E-016 };

enum { bsize = sizeof(b)/sizeof(b[0]) };

static const double c[] = {
  0.61825480456471904,
  2.75858288416734245E-016,
 -0.61825480456471904,
 -2.75858288416734245E-016,
  0.59473487038421935,
  0.0000000000000000,
  6.18759993964230386E-016,
 -0.59473487038421935,
 -6.18759993964230386E-016 };

enum { csize = sizeof(c)/sizeof(c[0]) };

static const double d[] = {
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
enum { dsize = sizeof(d)/sizeof(d[0]) };

static const double e[] = {
-4.5323096146354345e-17,
0.60013302811469493,
8.3011177257981544e-17,
};
enum { esize = sizeof(e)/sizeof(e[0]) };

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
int passed_d = fabs(creal(s) - d[dsize-1]) < DBL_MIN;
passed &= passed_d;
if (!passed_d)
  fprintf(stderr, "Bad result from Kahan summation test d, "
          "should have been close to %24.17g but is: %24.17g\n",
          d[dsize-1], creal(s));

if (!passed)
  return EXIT_FAILURE;

  ;
  return EXIT_SUCCESS;
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
