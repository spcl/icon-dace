/*
 * @file bittab.c
 * @brief genometools bit table class adapted for ScalES-PPM
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
 */
/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <assert.h>

#include "core/bittab.h"
#include "core/ppm_xfuncs.h"

struct PPM_Bittab {
  unsigned long *tabptr,
                tabsize,
                num_of_bits;
};

PPM_Bittab* PPM_bittab_new(unsigned long num_of_bits)
{
  PPM_Bittab *b;

  assert(num_of_bits);

  b = xmalloc(sizeof (PPM_Bittab));
  b->num_of_bits = num_of_bits;

  if (num_of_bits / (8UL * sizeof (unsigned long)))
    b->tabsize = 1 + ((num_of_bits - 1) / (8UL * sizeof (unsigned long)));
  else
    b->tabsize = 1UL;

  b->tabptr = xcalloc(b->tabsize, sizeof (unsigned long));

  return b;
}

void PPM_bittab_set_bit(PPM_Bittab *b, unsigned long bit)
{
  assert(b && bit < b->num_of_bits);
  b->tabptr[(bit >> 3) / sizeof (unsigned long)] |=
    1UL << (bit & (8UL * sizeof (unsigned long) - 1));
}

void PPM_bittab_unset_bit(PPM_Bittab *b, unsigned long bit)
{
  assert(b && bit < b->num_of_bits);
  b->tabptr[(bit >> 3) / sizeof (unsigned long)] &=
    ~(1UL << (bit & (8UL * sizeof (unsigned long) - 1)));
}

void PPM_bittab_complement(PPM_Bittab *dest, const PPM_Bittab *src)
{
  unsigned long i;

  assert(dest && src && dest->num_of_bits == src->num_of_bits);

  for (i = 0; i < src->tabsize - 1; i++)
    dest->tabptr[i] = ~src->tabptr[i];

  /* the ``last'' bittab gets special treatment to prevent that unused bits
     become set. this could disturb subsequent PPM_bittab_count_set_bits() calls.
   */
  dest->tabptr[src->tabsize - 1] = ~src->tabptr[src->tabsize - 1] &
                                   (~0UL >> (- src->num_of_bits +
                                             src->tabsize * 8UL *
                                             sizeof (unsigned long)));
}

void PPM_bittab_equal(PPM_Bittab *dest, const PPM_Bittab *src)
{
  unsigned long i;
  assert(dest && src && dest->num_of_bits == src->num_of_bits);
  for (i = 0; i < src->tabsize; i++)
    dest->tabptr[i] = src->tabptr[i];
}

void PPM_bittab_and(PPM_Bittab *dest, const PPM_Bittab *src1,
                   const PPM_Bittab *src2)
{
  unsigned long i;
  assert(dest && src1 && src2);
  assert(dest->num_of_bits == src1->num_of_bits);
  assert(dest->num_of_bits == src2->num_of_bits);
  for (i = 0; i < src1->tabsize; i++)
    dest->tabptr[i] = src1->tabptr[i] & src2->tabptr[i];
}

void PPM_bittab_or(PPM_Bittab *dest, const PPM_Bittab *src1, const PPM_Bittab *src2)
{
  unsigned long i;
  assert(dest && src1 && src2);
  assert(dest->num_of_bits == src1->num_of_bits);
  assert(dest->num_of_bits == src2->num_of_bits);
  for (i = 0; i < src1->tabsize; i++)
    dest->tabptr[i] = src1->tabptr[i] | src2->tabptr[i];
}

void PPM_bittab_nand(PPM_Bittab *dest,
                  const PPM_Bittab *minuend,
                  const PPM_Bittab *subtrahend)
{
  unsigned long i;
  assert(dest && minuend && subtrahend);
  assert(dest->num_of_bits == minuend->num_of_bits);
  assert(minuend->num_of_bits == subtrahend->num_of_bits);
  for (i = 0; i < dest->tabsize; i++)
    dest->tabptr[i] = minuend->tabptr[i] & ~subtrahend->tabptr[i];
}

void PPM_bittab_and_equal(PPM_Bittab *dest, const PPM_Bittab *src)
{
  unsigned long i;
  assert(dest && src);
  assert(dest->num_of_bits == src->num_of_bits);
  for (i = 0; i < dest->tabsize; i++)
    dest->tabptr[i] &= src->tabptr[i];
}

void PPM_bittab_or_equal(PPM_Bittab *dest, const PPM_Bittab *src)
{
  unsigned long i;
  assert(dest && src);
  assert(dest->num_of_bits == src->num_of_bits);
  for (i = 0; i < dest->tabsize; i++)
    dest->tabptr[i] |= src->tabptr[i];
}

void PPM_bittab_shift_left_equal(PPM_Bittab *b)
{
  unsigned long i, new_carry, old_carry = 0;
  assert(b);
  for (i = 0; i < b->tabsize; i++) {
    new_carry = b->tabptr[i] & (1UL << (8UL * sizeof (unsigned long) - 1));
    b->tabptr[i] = (b->tabptr[i] << 1) | old_carry;
    old_carry = new_carry >> (8UL * sizeof (unsigned long) - 1);
  }
}

void PPM_bittab_shift_right_equal(PPM_Bittab *b)
{
  unsigned long i, new_carry, old_carry = 0;
  assert(b);
  for (i = b->tabsize; i > 0; i--) {
    new_carry = b->tabptr[i-1] & 1UL;
    b->tabptr[i-1] = (b->tabptr[i-1] >> 1) | old_carry;
    old_carry = new_carry << (8UL * sizeof (unsigned long) - 1);
  }
}

void PPM_bittab_unset(PPM_Bittab *b)
{
  unsigned long i;
  assert(b);
  for (i = 0; i < b->tabsize; i++)
    b->tabptr[i] = 0;
}

void PPM_bittab_get_all_bitnums(const PPM_Bittab *b, PPM_Array *bitnums)
{
  unsigned long i;
  assert(b && bitnums);
  for (i = 0; i < b->num_of_bits; i++)
    if (PPM_bittab_bit_is_set(b, i)) PPM_array_add(bitnums, i);
}

bool PPM_bittab_bit_is_set(const PPM_Bittab *b, unsigned long bit)
{
  assert(b && bit < b->num_of_bits);
  if (b->tabptr[(bit >> 3) / sizeof (unsigned long)] &
      1UL << (bit & (8UL * sizeof (unsigned long) - 1))) {
    return true;
  }
  return false;
}

bool PPM_bittab_is_true(const PPM_Bittab *b)
{
  unsigned long i;
  assert(b);
  for (i = 0; i < b->tabsize; i++) {
    if (b->tabptr[i])
      return true;
  }
  return false;
}

bool PPM_bittab_cmp(const PPM_Bittab *b1, const PPM_Bittab *b2)
{
  unsigned long i;
  assert(b1 && b2 && b1->num_of_bits == b2->num_of_bits);
  for (i = 0; i < b1->tabsize; i++) {
    if (b1->tabptr[i] != b2->tabptr[i])
      return false;
  }
  return true;
}

unsigned long PPM_bittab_size(PPM_Bittab *b)
{
  assert(b);
  return b->num_of_bits;
}

unsigned long PPM_bittab_get_first_bitnum(const PPM_Bittab *b)
{
  unsigned long i, rval = (unsigned long)-1;
  assert(b);
  for (i = 0; i < b->num_of_bits; i++)
    if (PPM_bittab_bit_is_set(b, i)) {
      rval = i;
      break;
    }
  if (rval == (unsigned long)-1)
    return b->num_of_bits;
  return rval;
}

unsigned long PPM_bittab_get_last_bitnum(const PPM_Bittab *b)
{
  assert(b);
  return b->num_of_bits;
}

unsigned long PPM_bittab_get_next_bitnum(const PPM_Bittab *b,
                                        unsigned long curnum)
{
  unsigned long i, rval = (unsigned long)-1;

  assert(b);
  assert(curnum < b->num_of_bits);
  for (i = curnum + 1; i < b->num_of_bits; i++)
    if (PPM_bittab_bit_is_set(b, i)) {
      rval = i;
      break;
    }
  if (rval == (unsigned long)-1)
    return b->num_of_bits;
  return rval;
}

unsigned long
PPM_bittab_get_next_clear_bitnum(const PPM_Bittab *b,
                                 unsigned long curnum)
{
  unsigned long i, rval = (unsigned long)-1;

  assert(b);
  assert(curnum < b->num_of_bits);
  for (i = curnum + 1; i < b->num_of_bits; i++)
    if (!PPM_bittab_bit_is_set(b, i)) {
      rval = i;
      break;
    }
  if (rval == (unsigned long)-1)
    return b->num_of_bits;
  return rval;
}

unsigned long PPM_bittab_count_set_bits(const PPM_Bittab *b)
{
  static const unsigned char bits_in_char[256] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2,
                                                   2, 3, 2, 3, 3, 4, 1, 2, 2, 3,
                                                   2, 3, 3, 4, 2, 3, 3, 4, 3, 4,
                                                   4, 5, 1, 2, 2, 3, 2, 3, 3, 4,
                                                   2, 3, 3, 4, 3, 4, 4, 5, 2, 3,
                                                   3, 4, 3, 4, 4, 5, 3, 4, 4, 5,
                                                   4, 5, 5, 6, 1, 2, 2, 3, 2, 3,
                                                   3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                                                   2, 3, 3, 4, 3, 4, 4, 5, 3, 4,
                                                   4, 5, 4, 5, 5, 6, 2, 3, 3, 4,
                                                   3, 4, 4, 5, 3, 4, 4, 5, 4, 5,
                                                   5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
                                                   4, 5, 5, 6, 5, 6, 6, 7, 1, 2,
                                                   2, 3, 2, 3, 3, 4, 2, 3, 3, 4,
                                                   3, 4, 4, 5, 2, 3, 3, 4, 3, 4,
                                                   4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                                   2, 3, 3, 4, 3, 4, 4, 5, 3, 4,
                                                   4, 5, 4, 5, 5, 6, 3, 4, 4, 5,
                                                   4, 5, 5, 6, 4, 5, 5, 6, 5, 6,
                                                   6, 7, 2, 3, 3, 4, 3, 4, 4, 5,
                                                   3, 4, 4, 5, 4, 5, 5, 6, 3, 4,
                                                   4, 5, 4, 5, 5, 6, 4, 5, 5, 6,
                                                   5, 6, 6, 7, 3, 4, 4, 5, 4, 5,
                                                   5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                                                   4, 5, 5, 6, 5, 6, 6, 7, 5, 6,
                                                   6, 7, 6, 7, 7, 8 };
  unsigned long i, j, counter = 0;
  assert(b);
  for (i = 0; i < b->tabsize; i++)
    for (j = 0; j < sizeof (unsigned long); j++)
      counter += bits_in_char[((b->tabptr[i] >> (j * 8)) & 0xffu)];
  return counter;
}

void PPM_bittab_show(const PPM_Bittab *b, FILE *outfp)
{
  unsigned long i;
  assert(b && outfp);
  /* header line */
  for (i = 0; i < b->num_of_bits; i++)
    fprintf(outfp, "%lu", i % 10);
  xfputc('\n', outfp);
  /* actual bits */
  for (i = 0; i < b->num_of_bits; i++) {
    if (PPM_bittab_bit_is_set(b, i))
      xfputc('1', outfp);
    else
      xfputc('0', outfp);
  }
  xfputc('\n', outfp);
}

void PPM_bittab_delete(PPM_Bittab *b)
{
  if (!b) return;
  free(b->tabptr);
  free(b);
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
