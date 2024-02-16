/**
 * @file interval_tree.c
 *
 * @copyright Copyright  (C)  2014 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
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

#include <stdlib.h>

#include "ensure_array_size.h"
#include "interval_tree.h"
#include "utils.h"

static int
iv_compar(struct interval_node *a, struct interval_node *b)
{
  return (a->range.left > b->range.left) - (a->range.left < b->range.left);
}

static double
tree_part(struct interval_node intervals[], size_t num_nodes)
{
  size_t med = num_nodes / 2;
  double max = intervals[med].range.right;
  if (med)
  {
    double right;
    if ((right = tree_part(intervals, med)) > max)
      max = right;
  }
  if (med + 1 < num_nodes)
  {
    double right;
    if ((right = tree_part(intervals + med + 1, num_nodes - med - 1)) > max)
      max = right;
  }
  intervals[med].max = max;
  return max;
}

void
yac_generate_interval_tree(struct interval_node intervals[], size_t num_nodes)
{
  qsort(intervals, num_nodes, sizeof(intervals[0]),
        (int(*)(const void *, const void*))iv_compar);
  tree_part(intervals, num_nodes);
}

static inline void
overlap_push(struct overlaps *overlaps, size_t idx)
{
  size_t i = overlaps->num_overlaps;
  ENSURE_ARRAY_SIZE(overlaps->overlap_iv, overlaps->a_size, i + 1);
  overlaps->overlap_iv[i] = idx;
  overlaps->num_overlaps++;
}

static void
search_interval_tree_(struct interval_node tree[], size_t num_nodes,
                      struct interval query, struct overlaps *overlaps,
                      size_t ofs)
{
  /* while x is non-empty sub-trees */
  if (num_nodes)
  {
    size_t x = num_nodes/2;
    if (overlap_test(tree[x].range, query))
      overlap_push(overlaps, x + ofs);
    if (x && tree[x/2].max >= query.left)
      search_interval_tree_(tree, x, query, overlaps, ofs);
    if (x < num_nodes - 1
        && tree[x].range.left <= query.right
        && tree[x + 1 + (num_nodes - x - 1)/2].max >= query.left)
      search_interval_tree_(tree + x + 1, num_nodes - x - 1, query, overlaps,
                            ofs + x + 1);
  }
}

void
yac_search_interval_tree(struct interval_node tree[], size_t num_nodes,
                         struct interval query, struct overlaps *overlaps)
{
  search_interval_tree_(tree, num_nodes, query, overlaps, 0);
}


