/**
 * @file xt_xmap_intersection_common.h
 *
 * @brief Utility functions shared by xt_xmap_intersection and
 *        xt_xmap_intersection_ext.
 *
 * @copyright Copyright  (C)  2019 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#ifndef XT_XMAP_INTERSECTION_COMMON_H
#define XT_XMAP_INTERSECTION_COMMON_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

struct pos_run {
  size_t len;
  int start, direction;
};

#include "core/core.h"
#include "xt_stripe_util.h"

/* how many pos values have monotonically either positively or
 * negatively consecutive values */
static inline struct pos_run get_pos_run_len(
  size_t num_pos, const int *restrict pos)
{
  size_t i = 0, j = 1;
  int direction = 0;
  int start = pos[0];
  if (j < num_pos) {
    direction = isign_mask(pos[1] - pos[0]);
    while (j < num_pos
           && pos[j] == start + (~direction & (int)(j - i)) +
                        (direction & -(int)(j - i)))
      ++j;
    direction = direction & ((j == 1) - 1);
  }
  return (struct pos_run){ .start = start, .len = j, .direction = direction };
}

/* compute number of position extents that would be required
   to represent positions array */
static size_t count_pos_ext(size_t num_pos, const int *restrict pos)
{
  size_t i = 0, num_pos_ext = 0;
  while (i < num_pos) {
    i += get_pos_run_len(num_pos - i, pos + i).len;
    ++num_pos_ext;
  }
  return num_pos_ext;
}

/*
 * generate position extents for given positions array
 */
static inline void generate_pos_ext(
  size_t num_pos, const int *restrict pos,
  size_t num_pos_ext, struct Xt_pos_ext *restrict pos_ext)
{
#ifdef NDEBUG
  (void)num_pos_ext;
#endif
  size_t i = 0, temp_num_pos_ext = 0;
  while (i < num_pos)
  {
    struct pos_run rn = get_pos_run_len(num_pos - i, pos + i);
    pos_ext[temp_num_pos_ext]
      = (struct Xt_pos_ext){ .start = rn.start,
                             .size = (rn.direction*2 + 1) * (int)rn.len };
    ++temp_num_pos_ext;
    i += rn.len;
  }
  assert(num_pos_ext == temp_num_pos_ext);
}

/* compute number of positions that would be required
   to represent position extents */
static inline size_t count_pos(
  size_t num_pos_ext, const struct Xt_pos_ext *restrict pos_ext)
{
  size_t num_pos = 0;
  for (size_t i = 0; i < num_pos_ext; ++i) num_pos += (size_t)(pos_ext[i].size);
  return num_pos;
}

static inline void generate_pos(
  size_t num_pos_ext, const struct Xt_pos_ext *restrict pos_ext,
  size_t num_pos, int *restrict pos)
{
#ifdef NDEBUG
  (void)num_pos;
#endif
  size_t ofs = 0;
  for (size_t i = 0; i < num_pos_ext; ++i) {
    int abssize = abs(pos_ext[i].size);
    int step = isign(pos_ext[i].size);
    for (int j = 0; j < abssize; ++j)
      pos[ofs + (size_t)j] = pos_ext[i].start + j * step;
    ofs += (size_t)abssize;
  }
  assert(ofs == num_pos);
}

static inline void
print_miss_msg(Xt_idxlist dst_idxlist, int missing_pos,
               MPI_Comm comm, const char *source, int line)
  __attribute__((noreturn));

static inline void
print_miss_msg(Xt_idxlist dst_idxlist, int missing_pos,
               MPI_Comm comm, const char *source, int line)
{
  Xt_int missing_index;
  xt_idxlist_get_index_at_position(dst_idxlist, missing_pos, &missing_index);
  static const char fmt[] = "ERROR: destination intersections do not match "
    "with destination index list (first missing index %lld "
    "at position %d)";
  char error_message[sizeof (fmt)
                     + (sizeof (long long) + sizeof (int)) * CHAR_BIT / 8 * 3];
  sprintf(error_message, fmt, (long long)missing_index, missing_pos);
  Xt_abort(comm, error_message, source, line);
}

#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
