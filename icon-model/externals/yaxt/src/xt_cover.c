/**
 * @file xt_cover.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>

#include "core/ppm_xfuncs.h"
#include "ensure_array_size.h"
#include "xt/xt_idxlist.h"
#include "xt_cover.h"

void
xt_cover_start(struct Xt_pos_ext_vec *restrict cover,
               size_t initial_size)
{
  cover->num_pos_ext = 0;
  cover->size_pos_ext = initial_size;
  cover->pos_ext = xmalloc(sizeof (*cover->pos_ext) * initial_size);
}

void
xt_cover_finish(struct Xt_pos_ext_vec *restrict cover)
{
  free(cover->pos_ext);
}

bool
xt_idxlist_pos_ext_is_full_cover(Xt_idxlist idxlist,
                              struct Xt_pos_ext_vec cover)
{
  int idxlist_size = xt_idxlist_get_num_indices(idxlist);
  if ((idxlist_size == 0) & (cover.num_pos_ext == 0))
    return true;
  if ((idxlist_size == 0) ^ (cover.num_pos_ext == 0))
    return false;
  /* at this point both idxlist and cover are non-empty */
  int after_last_end_pos = 0;
  int continuations_hold = 1;
  for (size_t i = 0; i < cover.num_pos_ext; ++i) {
    continuations_hold &= (cover.pos_ext[i].start == after_last_end_pos);
    after_last_end_pos += cover.pos_ext[i].size;
  }
  continuations_hold &= (after_last_end_pos == idxlist_size);
  return continuations_hold;
}

size_t
xt_cover_search(struct Xt_pos_ext_vec *restrict cover,
                struct Xt_pos_range query, bool forward,
                size_t search_start_pos)
{
  size_t i, num_pos_ext = cover->num_pos_ext;
  struct Xt_pos_ext *restrict covered_pos_ext = cover->pos_ext;
  if (num_pos_ext)
  {
    i = search_start_pos;
    if (forward) {
      while ((query.start >= (covered_pos_ext[i].start
                              + covered_pos_ext[i].size))
             & (i < num_pos_ext - 1))
        ++i;
    } else {
      while ((i > 0)
             & (query.end < covered_pos_ext[i - 1].start - 1))
        --i;
      if ((i == 0)
          & (query.start > covered_pos_ext[0].start + covered_pos_ext[0].size))
        i = 1;
    }
    if ((i >= num_pos_ext - 1)
        & (query.start >= (covered_pos_ext[num_pos_ext - 1].start
                           + covered_pos_ext[num_pos_ext - 1].size)))
      i = SIZE_MAX;
  }
  else
    i = SIZE_MAX;
  return i;
}

void
xt_cover_range_append(struct Xt_pos_ext_vec *restrict cover,
                      struct Xt_pos_ext range)
{
  size_t num_pos_ext = cover->num_pos_ext;
  struct Xt_pos_ext *restrict cover_pos_ext = cover->pos_ext;
  if (num_pos_ext > 0
      && (cover_pos_ext[num_pos_ext - 1].start
          + cover_pos_ext[num_pos_ext - 1].size
          == range.start)) {
    cover_pos_ext[num_pos_ext - 1].size += range.size;
  } else {
    ENSURE_ARRAY_SIZE(cover_pos_ext, cover->size_pos_ext, num_pos_ext + 1);
    cover->pos_ext = cover_pos_ext;
    cover_pos_ext[num_pos_ext] = range;
    ++cover->num_pos_ext;
  }
}


size_t
xt_cover_insert_or_overlap(struct Xt_pos_ext_vec *restrict cover,
                           struct Xt_pos_range range, bool forward,
                           size_t search_start_pos)
{
  size_t insert_pos = xt_cover_search(cover, range, forward, search_start_pos);
  struct Xt_pos_ext *restrict cover_pos_ext = cover->pos_ext;
  if (insert_pos != SIZE_MAX) {
    if ((range.start < (cover_pos_ext[insert_pos].start
                        + cover_pos_ext[insert_pos].size))
        & (range.end >= cover_pos_ext[insert_pos].start))
    {
      /* let caller handle overlap cases */
      return insert_pos;
    }
    else if (range.end + 1 < cover_pos_ext[insert_pos].start)
    {
      /* range precedes cover->pos_ext[insert_pos] with a hole
       * but might be a seemless extension of
       * cover->pos_ext[insert_pos-1] */
      if (insert_pos > 0
          && (cover_pos_ext[insert_pos - 1].start
              + cover_pos_ext[insert_pos - 1].size == range.start))
        cover_pos_ext[insert_pos - 1].size += range.end - range.start + 1;
      else
      {
        size_t num_pos_ext = cover->num_pos_ext;
        ENSURE_ARRAY_SIZE(cover_pos_ext,
                          cover->size_pos_ext,
                          num_pos_ext + 1);
        cover->pos_ext = cover_pos_ext;
        memmove(cover_pos_ext + insert_pos + 1, cover_pos_ext + insert_pos,
                sizeof (*cover_pos_ext)
                * (num_pos_ext - insert_pos));
        cover_pos_ext[insert_pos]
          = (struct Xt_pos_ext){ .start = range.start,
                                 .size = range.end - range.start + 1};
        cover->num_pos_ext = ++num_pos_ext;
      }
      insert_pos = SIZE_MAX;
    }
    else if (range.end + 1 == cover_pos_ext[insert_pos].start)
    {
      cover_pos_ext[insert_pos].start = range.start;
      cover_pos_ext[insert_pos].size += range.end - range.start + 1;
      if (insert_pos > 0
          && (cover_pos_ext[insert_pos].start
              == (cover_pos_ext[insert_pos - 1].start
                  + cover_pos_ext[insert_pos - 1].size)))
      {
        cover_pos_ext[insert_pos - 1].size
          += cover_pos_ext[insert_pos].size;
        memmove(cover_pos_ext + insert_pos, cover_pos_ext + insert_pos + 1,
                (--cover->num_pos_ext - insert_pos)
                * sizeof (*cover_pos_ext));
      }
      insert_pos = SIZE_MAX;
    }
  } else {
    /* range was not found -> append */
    xt_cover_range_append(cover, (struct Xt_pos_ext){ .start = range.start,
          .size = range.end - range.start + 1});
  }
  return insert_pos;
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
