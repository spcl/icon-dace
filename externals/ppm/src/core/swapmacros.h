/*
 * @file swapmacros.h
 * @brief efficient swapping macros used in sort routines
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * Keywords:
 * @author Thomas Jahns <jahns@dkrz.de>
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

#ifndef SWAPMACROS_H
#define SWAPMACROS_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stddef.h>

/*
 * swap macros from qsort
 */
#define SWAPINIT(a, es) do {                                     \
    swaptype                                                     \
      = (((unsigned char *)a - (unsigned char *)0)               \
         % (ptrdiff_t)sizeof (long)                              \
         || es % sizeof (long)) ?                                \
      2 : es == sizeof (long) ? 0 : 1;                           \
  } while (0)

#define swapcode(TYPE, parmi, parmj, n) do {            \
    size_t i = (n) / sizeof (TYPE);                     \
    TYPE *pi = (TYPE *) (parmi);                        \
    TYPE *pj = (TYPE *) (parmj);                        \
    do {                                                \
      TYPE t = *pi;                                     \
      *pi++ = *pj;                                      \
      *pj++ = t;                                        \
    } while (--i);                                      \
  } while (0)

static inline void
swapfunc(unsigned char *a, unsigned char *b, size_t n, int swaptype)
{
  if (swaptype <= 1)
    swapcode(long, a, b, n);
  else
    swapcode(unsigned char, a, b, n);
}

#undef swapcode

#define swap(a, b)                              \
  do {                                          \
    if (swaptype == 0) {                        \
      long *pa_ = (long *)(a);                  \
      long t = *pa_;                            \
      long *pb_ = (long *)(b);                  \
      *pa_ = *pb_;                              \
      *pb_ = t;                                 \
    }                                           \
    else                                        \
      swapfunc((a), (b), es, swaptype);         \
  } while (0)

#endif
/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
