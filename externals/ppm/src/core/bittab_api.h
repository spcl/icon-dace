/**
 * @file bittab_api.h
 * @brief genometools bit table class adapted for ScalES-PPM
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
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
 */
/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef BITTAB_API_H
#define BITTAB_API_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdbool.h>
#include <stdio.h>
#include "core/array_api.h"

/* Implements arbitrary-length bit arrays and various operations on them. */
typedef struct PPM_Bittab PPM_Bittab;

/* Return a new <PPM_Bittab> of length <num_of_bits>, initialised to 0. */
PPM_Bittab*     PPM_bittab_new(unsigned long num_of_bits);

/* Set bit <i> in <bittab> to 1. */
void          PPM_bittab_set_bit(PPM_Bittab *bittab, unsigned long i);

/* Set bit <i> in <bittab> to 0. */
void          PPM_bittab_unset_bit(PPM_Bittab *bittab, unsigned long i);

/* Set <bittab_a> to be the complement of <bittab_b>. */
void          PPM_bittab_complement(PPM_Bittab *bittab_a,
                                   const PPM_Bittab *bittab_b);

/* Set <bittab_a> to be equal to <bittab_b>. */
void          PPM_bittab_equal(PPM_Bittab *bittab_a, const PPM_Bittab *bittab_b);

/* Set <bittab_a> to be the bitwise AND of <bittab_b> and <bittab_c>. */
void          PPM_bittab_and(PPM_Bittab *bittab_a, const PPM_Bittab *bittab_b,
                            const PPM_Bittab *bittab_c);

/* Set <bittab_a> to be the bitwise OR of <bittab_b> and <bittab_c>. */
void          PPM_bittab_or(PPM_Bittab *bittab_a, const PPM_Bittab *bittab_b,
                           const PPM_Bittab *bittab_c);

/* Set <bittab_a> to be <bittab_b> NAND <bittab_c>. */
void          PPM_bittab_nand(PPM_Bittab *bittab_a, const PPM_Bittab *bittab_b,
                             const PPM_Bittab *bittab_c);

/* Set <bittab_a> to be the bitwise AND of <bittab_a> and <bittab_b>. */
void          PPM_bittab_and_equal(PPM_Bittab *bittab_a, const PPM_Bittab *bittab_b);

/* Set <bittab_a> to be the bitwise OR of <bittab_a> and <bittab_b>. */
void          PPM_bittab_or_equal(PPM_Bittab *bittab_a, const PPM_Bittab *bittab_b);

/* Shift <bittab> by one position to the left. */
void          PPM_bittab_shift_left_equal(PPM_Bittab *bittab);

/* Shift <bittab> by one position to the right. */
void          PPM_bittab_shift_right_equal(PPM_Bittab *bittab);

/* Set all bits in <bittab> to 0. */
void          PPM_bittab_unset(PPM_Bittab *bittab);

/* Output a representation of <bittab> to <fp>. */
void          PPM_bittab_show(const PPM_Bittab *bittab, FILE *fp);

/* Fill <array> with the indices of all set bits in <bittab>. */
void          PPM_bittab_get_all_bitnums(const PPM_Bittab *bittab, PPM_Array *array);

/* Return <true> if bit <i> is set in <bittab>. */
bool          PPM_bittab_bit_is_set(const PPM_Bittab *bittab, unsigned long i);

/* Return <true> if <bittab_a> and <bittab_b> are identical. */
bool          PPM_bittab_cmp(const PPM_Bittab *bittab_a, const PPM_Bittab *bittab_b);

/* Return the index of the first set bit in <bittab>. */
unsigned long PPM_bittab_get_first_bitnum(const PPM_Bittab *bittab);

/* Return the index of the last set bit in <bittab>. */
unsigned long PPM_bittab_get_last_bitnum(const PPM_Bittab *bittab);

/* Return the index of the next set bit in <bittab> with an index greater
   than <i>. */
unsigned long PPM_bittab_get_next_bitnum(const PPM_Bittab *bittab,
                                         unsigned long i);

/* Return the index of the next clear bit in <bittab> with an index greater
   than <i>. */
unsigned long
PPM_bittab_get_next_clear_bitnum(const PPM_Bittab *bittab,
                                 unsigned long i);

/* Return the number of set bits in <bittab>. */
unsigned long PPM_bittab_count_set_bits(const PPM_Bittab *bittab);

/* Return the total number of bits of <bittab>. */
unsigned long PPM_bittab_size(PPM_Bittab *bittab);

/* Delete <bittab>. */
void          PPM_bittab_delete(PPM_Bittab *bittab);

#endif
/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
