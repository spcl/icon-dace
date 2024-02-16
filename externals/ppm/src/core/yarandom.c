/**
 * @file yarandom.c
 * @brief base random number generation routines from xscreensaver
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/* Some changes for:
 * - interoperability with OpenMP
 * - efficient generation of 64bit random numbers
 * Copyright (c) 2011 by Thomas Jahns <jahns@dkrz.de>
 */
/* yarandom.c -- Yet Another Random Number Generator.
 * Copyright (c) 1997, 1998, 2003 by Jamie Zawinski <jwz@jwz.org>
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation.  No representations are made about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 */

/* The unportable mess that is rand(), random(), drand48() and friends led me
   to ask Phil Karlton <karlton@netscape.com> what the Right Thing to Do was.
   He responded with this.  It is non-cryptographically secure, reasonably
   random (more so than anything that is in any C library), and very fast.

   I don't understand how it works at all, but he says "look at Knuth,
   Vol. 2 (original edition), page 26, Algorithm A.  In this case n=55,
   k=20 and m=2^32."

   So there you have it.

   ---------------------------
   Note: xlockmore 4.03a10 uses this very simple RNG:

        if ((seed = seed % 44488 * 48271 - seed / 44488 * 3399) < 0)
          seed += 2147483647;
        return seed-1;

   of which it says

        ``Dr. Park's algorithm published in the Oct. '88 ACM  "Random Number
          Generators: Good Ones Are Hard To Find" His version available at
          ftp://cs.wm.edu/pub/rngs.tar Present form by many authors.''

   Karlton says: ``the usual problem with that kind of RNG turns out to
   be unexepected short cycles for some word lengths.''

   Karlton's RNG is faster, since it does three adds and two stores, while the
   xlockmore RNG does two multiplies, two divides, three adds, and one store.

   Compiler optimizations make a big difference here:
       gcc -O:     difference is 1.2x.
       gcc -O2:    difference is 1.4x.
       gcc -O3:    difference is 1.5x.
       SGI cc -O:  difference is 2.4x.
       SGI cc -O2: difference is 2.4x.
       SGI cc -O3: difference is 5.1x.
   Irix 6.2; Indy r5k; SGI cc version 6; gcc version 2.7.2.1.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>  /* for getpid() */
#include <sys/time.h> /* for gettimeofday() */
#include <inttypes.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "core/core.h"
#include "core/yarandom.h"
# undef PPM_ya_rand_init

/* The following 'random' numbers are taken from CRC, 18th Edition, page 622.
   Each array element was taken from the corresponding line in the table,
   except that a[0] was from line 100. 8s and 9s in the table were simply
   skipped. The high order digit was taken mod 4.
 */
enum {
  VectorSize = 55
};

static const uint32_t a_source[VectorSize] = {
 035340171546, 010401501101, 022364657325, 024130436022, 002167303062, /*  5 */
 037570375137, 037210607110, 016272055420, 023011770546, 017143426366, /* 10 */
 014753657433, 021657231332, 023553406142, 004236526362, 010365611275, /* 14 */
 007117336710, 011051276551, 002362132524, 001011540233, 012162531646, /* 20 */
 007056762337, 006631245521, 014164542224, 032633236305, 023342700176, /* 25 */
 002433062234, 015257225043, 026762051606, 000742573230, 005366042132, /* 30 */
 012126416411, 000520471171, 000725646277, 020116577576, 025765742604, /* 35 */
 007633473735, 015674255275, 017555634041, 006503154145, 021576344247, /* 40 */
 014577627653, 002707523333, 034146376720, 030060227734, 013765414060, /* 45 */
 036072251540, 007255221037, 024364674123, 006200353166, 010126373326, /* 50 */
 015664104320, 016401041535, 016215305520, 033115351014, 017411670323  /* 55 */
};

#ifdef _OPENMP
static uint32_t (*a)[VectorSize];
static int (*vidx)[2];
#else
static uint32_t a[1][VectorSize];
static int vidx[1][2];
#endif


uint32_t
PPM_ya_random (void)
{
  int tid
#ifdef _OPENMP
    = omp_get_thread_num();
#else
    = 0;
#endif
  uint32_t ret;
  int i1 = vidx[tid][0], i2 = vidx[tid][1];
  ret = a[tid][i1] + a[tid][i2];
  a[tid][i1] = ret;
  if (++i1 >= VectorSize) i1 = 0;
  if (++i2 >= VectorSize) i2 = 0;
  vidx[tid][0] = i1;
  vidx[tid][1] = i2;
  return ret;
}

uint64_t
PPM_ya_random64(void)
{
  int tid
#ifdef _OPENMP
    = omp_get_thread_num();
#else
    = 0;
#endif
  uint32_t ret_low, ret_high;
  uint64_t ret;
  int i1 = vidx[tid][0], i2 = vidx[tid][1];
  ret_low = a[tid][i1] + a[tid][i2];
  a[tid][i1] = ret_low;
  if (++i1 >= VectorSize) i1 = 0;
  if (++i2 >= VectorSize) i2 = 0;
  ret_high = a[tid][i1] + a[tid][i2];
  a[tid][i1] = ret_high;
  if (++i1 >= VectorSize) i1 = 0;
  if (++i2 >= VectorSize) i2 = 0;
  vidx[tid][0] = i1;
  vidx[tid][1] = i2;
  ret = (uint64_t)ret_low | ((uint64_t)ret_high << 32);
  return ret;
}

#ifdef _OPENMP
#define PPM_NOOMP_UNUSED(x) x
#else
#define PPM_NOOMP_UNUSED(x) PPM_UNUSED(x)
#endif

unsigned
PPM_ya_rand_init(MPI_Comm PPM_NOOMP_UNUSED(comm), int seed_arg)
{
  unsigned seed = (unsigned)seed_arg;
  int i, tid;
#ifdef _OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif
  if (seed == 0)
  {
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    /* ignore overflow */
    seed
      = (unsigned)((999*tp.tv_sec) + (1001*tp.tv_usec) + (1003 * getpid()) + tid);
    seed_arg = (int)seed;
  }
#ifdef _OPENMP
#pragma omp master
  {
    int n = omp_get_num_threads();
    void *temp;
    if (!(temp = realloc(vidx, (unsigned)n * sizeof(*vidx))))
      free(vidx);
    vidx = temp;
    if (!(temp = realloc(a, (unsigned)n * sizeof(*a))))
      free(a);
    a = temp;
    if (!vidx || !a)
      PPM_abort(comm, "error in PPM_ya_rand_init", __FILE__, __LINE__);
  }
#pragma omp barrier
#endif
  /* allow to call PPM_ya_rand_init() multiple times */
  for (i = 0; i < VectorSize; i++)
    a[tid][i] = a_source[i];

  a[tid][0] += seed;
  for (i = 1; i < VectorSize; i++)
  {
    seed = a[tid][i-1]*1001 + seed*999;
    a[tid][i] += seed;
  }

  vidx[tid][0] = (int)(a[tid][0] % VectorSize);
  vidx[tid][1] = (int)((vidx[tid][0] + 024) % VectorSize);
  return (unsigned)seed_arg;
}

void
PPM_ya_rand_finish()
{
#ifdef _OPENMP
#pragma omp master
  {
    free(a);
    free(vidx);
  }
#endif
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
