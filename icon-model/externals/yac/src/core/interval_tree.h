// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

/** \example test_interval_tree.c
 * This contains a test of the interval_tree.
 */

#ifndef INTERVAL_TREE_H
#define INTERVAL_TREE_H

#include <stdlib.h>

struct interval
{
  double left, right;
};

static inline int
overlap_test(struct interval a, struct interval b)
{
  return (a.left <= b.left && a.right >= b.left) ||
    (a.left > b.left && a.left <= b.right);
}

struct interval_node
{
  struct interval range;
  double max;
  size_t value;
};

void
yac_generate_interval_tree(struct interval_node intervals[], size_t num_nodes);

struct overlaps
{
  size_t num_overlaps, a_size;
  size_t *overlap_iv;
};

void
yac_search_interval_tree(struct interval_node tree[], size_t num_nodes,
                         struct interval query, struct overlaps *overlaps);

#endif

