/*
 * @file insertion_sort.c
 * @brief perform insertion sort on array
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>
#include "core/insertion_sort.h"
#include "core/swapmacros.h"
#include "core/fptr_api.h"

void
PPM_insertion_sort(void *a, size_t n, size_t es, void *data,
                   PPM_CompareWithData cmp)
{
  unsigned char *last = (unsigned char *)a + (n - 1) * es;
  /* initializing buf might not be needed on many
   * platforms, but is still required by the C standard which
   * prohibits accessing the contents of uninitilized memory */
  unsigned char buf[es];
  memset(buf, 0, es);
  int swaptype;
  size_t mv_size;
  SWAPINIT(a, es);
  for (unsigned char *q = (unsigned char *)a + es; q <= last; q += es)
  {
    swap(buf, q);
    unsigned char *p = q;
    do {
      p -= es;
    } while ((void *)p >= a && cmp(p, buf, data) >= 0);
    if ((mv_size = (size_t)(q - es - p)) > 0)
      memmove(p + 2 * es, p + es, mv_size);
    swap(p + es, buf);
  }
}


void
PPM_sorted_insertion(void *a, size_t n, size_t es, void *data,
                     PPM_CompareWithData cmp)
{
  unsigned char *p, *last = (unsigned char *)a + (n - 1) * es;
  unsigned char buf[es];
  int swaptype;
  size_t mv_size;
  SWAPINIT(a, es);
  memcpy(buf, last, es);
  p = last - es;
  while ((void *)p >= a && cmp(p, buf, data) > 0)
    p -= es;
  if ((mv_size = (size_t)(last - es - p)) > 0)
    memmove(p + 2 * es, p + es, mv_size);
  swap(p + es, buf);
}

void
PPM_insertion_sort_once(void *a, size_t n, size_t es, size_t inversion,
                        void *data, PPM_CompareWithData cmp)
{
  unsigned char *q = (unsigned char *)a + es * inversion,
    *last = (unsigned char *)a + (n - 1) * es;
  if (q < last && cmp(q, q + es, data) > 0)
  {
    /* search a[inversion+1..n-1] for insertion place */
    unsigned char *p;
    size_t mv_size;
    int swaptype;
    unsigned char buf[es];
    SWAPINIT(a, es);
    swap(q, buf);
    p = q + es;
    do {
      p += es;
    } while (p <= last && cmp(buf, p, data) > 0);
    mv_size = (size_t)(p - es - q);
    memmove(q, q + es, mv_size);
    swap(p - es, buf);
  }
  else if (q > (unsigned char *)a && cmp(q - es, q, data) > 0)
  {
    /* search a[0..inversion-1] for insertion place */
    unsigned char *p;
    size_t mv_size;
    int swaptype;
    unsigned char buf[es];
    SWAPINIT(a, es);
    swap(q, buf);
    p = q - es;
    do {
      p -= es;
    } while ((void *)p >= a && cmp(p, buf, data) > 0);
    mv_size = (size_t)(q - es - p);
    memmove(p + 2 * es, p + es, mv_size);
    swap(p + es, buf);
  }
  /* else no insert necessary, because a[inversion] is not really misplaced */
}

/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
