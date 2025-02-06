/**
 * @file ppm_combinatorics_c.c
 * @brief combinatorial routines
 *
 * Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords: combinatorics
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

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "core/bittab.h"
#include "core/ppm_xfuncs.h"
#include "core/yarandom.h"
#include "core/ppm_visibility.h"
#include "core/ppm_combinatorics.h"

int PPM_DSO_API_EXPORT
PPM_prime_factorization_32(uint32_t n, uint32_t **factors)
{
  uint32_t lim, i, j, p;
  /* stores prime number candidates, where for each odd number i bit
   * i/2 is set if it's not a possible divisor, with the special case
   * that bit 0 encodes if 2 is not a divisor */
  PPM_Bittab *nprime_divs;
  int num_factors = 0;

  assert(n > 0);
  if (n == 1)
    return 0;
  lim = (uint32_t)ceil(sqrt(n));
  p = n;
  nprime_divs = PPM_bittab_new(lim/2 + 1);
  if (!(p % 2))
  {
    do {
      p /= 2;
      ++num_factors;
    } while (!(p % 2));
  }
  else
    PPM_bittab_set_bit(nprime_divs, 0);
  i = 3;
  while (i <= lim)
  {
    if (!(p % i))
    {
      do {
        p /= i;
        ++num_factors;
      } while (!(p % i));
    }
    else
      PPM_bittab_set_bit(nprime_divs, i/2);
    for (j = 3*i; j <= lim; j += 2*i)
      PPM_bittab_set_bit(nprime_divs, j/2);
    i = (uint32_t)PPM_bittab_get_next_clear_bitnum(nprime_divs, i/2);
    if (i == lim/2 + 1)
      break;
    i = 2 * i + 1;
  }
  /* assert(p == 1 || is_prime(p)); */
  if (p != 1)
    ++num_factors;
  if (factors)
  {
    uint32_t last_factor_pos = 0;
    i = 0;
    p = n;
    uint32_t *restrict factors_ = *factors;
    if (!factors_)
      *factors = factors_ = xmalloc(sizeof(**factors) * (size_t)num_factors);
    if (!PPM_bittab_bit_is_set(nprime_divs, 0))
      while (!(p % 2))
      {
        factors_[i] = 2;
        p /= 2;
        ++i;
      }
    while (i < (uint32_t)num_factors)
    {
      last_factor_pos
        = (uint32_t)PPM_bittab_get_next_clear_bitnum(nprime_divs,
                                                     last_factor_pos);
      if (last_factor_pos < lim/2 + 1)
      {
        uint32_t factor = last_factor_pos * 2 + 1;
        while (!(p % factor))
        {
          factors_[i] = factor;
          p /= factor;
          ++i;
        }
      }
      else
        break;
    }
    if (i < (uint32_t)num_factors)
    {
      assert(p /= 1);
      factors_[i] = p;
    }
  }
  PPM_bittab_delete(nprime_divs);
  return num_factors;
}

static inline size_t
PPM_rand_size_t(void)
{
#if SIZEOF_SIZE_T == 8
  uint64_t r = PPM_ya_random64();
#elif SIZEOF_SIZE_T == 4
  uint32_t r = PPM_ya_random();
#else
#error "Unexpected size of size_t"
#endif
  return r;
}

static inline void
PPM_permute_randomly_i8(uint64_t *a, size_t n)
{
  size_t left = n;
  for (size_t i = 0; i < n; ++i, --left)
  {
    uint64_t temp = a[i];
    size_t s = PPM_rand_size_t() % left + i;
    a[i] = a[s];
    a[s] = temp;
  }
}

static inline void
PPM_permute_randomly_i4(uint32_t *a, size_t n)
{
  size_t left = n;
  for (size_t i = 0; i < n; ++i, --left)
  {
    uint32_t temp = a[i];
    size_t s = PPM_rand_size_t() % left + i;
    a[i] = a[s];
    a[s] = temp;
  }
}


void PPM_DSO_API_EXPORT
PPM_permute_randomly(void *a, size_t esize, size_t n)
{
  if ((intptr_t)a % 4 == 0 && esize == 4)
    PPM_permute_randomly_i4(a, n);
  else if ((intptr_t)a % 8 == 0 && esize == 8)
    PPM_permute_randomly_i8(a, n);
  else
  {
    size_t left = n;
    unsigned char temp[esize];
    for (size_t i = 0; i < n; ++i, --left)
    {
      memcpy(temp, (unsigned char *)a + i * esize, esize);
      size_t s = PPM_rand_size_t() % left + i;
      memcpy((unsigned char *)a + i * esize,
             (unsigned char *)a + s * esize, esize);
      memcpy((unsigned char *)a + s * esize , temp, esize);
    }
  }
}


/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
