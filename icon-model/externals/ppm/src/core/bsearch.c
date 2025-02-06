/*
 * @file bsearch.c
 * @brief re-entrant binary search functions
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * Keywords: binary search
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
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "core/ppm_visibility.h"

#include "core/bsearch.h"

void PPM_DSO_API_EXPORT *
PPM_bsearch_r(const void *key, const void *a, size_t nmemb, size_t size,
              void *data, PPM_CompareWithData cmp)
{
  const char *p = a, *q;
  size_t lim;
  int cmp_res;
  for (lim = nmemb; lim != 0; lim >>= 1) {
    q = p + (lim >> 1) * size;
    cmp_res = cmp(key, q, data);
    if (cmp_res > 0)
    {
      /* key > q: move right */
      p = (const char *)q + size;
      lim--;
    }
    else if (cmp_res < 0) {/* left index doesn't move if key < q */}
    else                   /* in case key == q, return q */
      return (void *)q;
  }
  return NULL;
}

void PPM_DSO_API_EXPORT*
PPM_bsearch_el_r(const void *key, const void *a, size_t nmemb, size_t size,
                 void *data, PPM_CompareWithData cmp)
{
  const char *p = a, *q;
  size_t lim;
  int cmp_res;
  for (lim = nmemb; lim != 0; lim >>= 1) {
    q = p + (lim >> 1) * size;
    cmp_res = cmp(key, q, data);
    if (cmp_res > 0)
    {
      /* key > q: move right */
      p = (const char *)q + size;
      lim--;
    }
    else if (cmp_res < 0) {/* left index doesn't move if key < q */}
    else                   /* in case key == q, return q */
      return (void *)q;
  }
  cmp_res = cmp(key, p, data);
  if (cmp_res > 0)
    return (void *)p;
  else if (p != a)              /* cmp_res < 0 is implicit here
                                 * because cmp_res == 0 is caught above */
    return (void *)(p - size);
  else  /* (p == a) */
    return NULL;
}

/*
 * Local Variables:
 * license-markup: "doxygen"
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */
