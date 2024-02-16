/**
 * @file yarandom.h
 * @brief internal header for PRNG
 *
 * Changes for ScalES-PPM:
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Jamie Zawinski <jwz@jwz.org>
 */
/*
 * Keywords: pseudo random number generator
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
 */
/* xscreensaver, Copyright (c) 1997, 1998, 2003 by Jamie Zawinski <jwz@jwz.org>
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation.  No representations are made about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 */

#ifndef YARANDOM_H
#define YARANDOM_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <inttypes.h>
#include "core/core.h"

#undef random
#undef rand
#undef drand48
#undef srandom
#undef srand
#undef srand48
#undef frand
#undef RAND_MAX

uint32_t PPM_ya_random (void);
uint64_t PPM_ya_random64(void);
unsigned int PPM_ya_rand_init(MPI_Comm comm, int);
void PPM_ya_rand_finish(void);

#define RAND_MAX   0x7FFFFFFFL
#define random()   ((long)PPM_ya_random() & RAND_MAX)

/*#define srandom(i) ya_rand_init(0)*/

/* Define these away to keep people from using the wrong APIs in scales_ppm.
 */
#define rand          __ERROR_use_random_not_rand_in_scales_ppm__
#define drand48       __ERROR_use_frand_not_drand48_in_scales_ppm__
#define srandom       __ERROR_do_not_call_srandom_in_scales_ppm__
#define srand         __ERROR_do_not_call_srand_in_scales_ppm__
#define srand48       __ERROR_do_not_call_srand48_in_scales_ppm__

/**
 * PRNG function for double precision floating-point quantities
 * @return random number in range [0.0,1.0)
 **/
double
PPM_ya_frandom(void);

/**
 * PRNG function for double precision floating-point quantities
 * @return random number in range [0.0f,1.0f)
 **/
float
PPM_ya_frandomf(void);

/**
 * PRNG function for double precision floating-point quantities
 * @return random number in range (-1.0,1.0), including 0.0
 **/
double
PPM_ya_fsgrandom(void);

/**
 * PRNG function for double precision floating-point quantities
 * @return random number in range (-1.0f,1.0f), including 0.0f
 **/
float
PPM_ya_fsgrandomf(void);

#define frand(f) (PPM_ya_frandom() * (f))

#endif
/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
