// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>

#include "ensure_array_size.h"
#include "interval_tree.h"
#include "utils_core.h"

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


