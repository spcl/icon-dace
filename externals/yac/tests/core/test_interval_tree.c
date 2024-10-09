// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "tests.h"
#include "interval_tree.h"

static void do_test(struct interval_node *tree, size_t num_nodes,
                    struct overlaps *overlaps);

enum {
  ntests = 100,
  ntrees = 100,
};


int main(void) {
  struct interval_node *tree;
  struct overlaps overlaps = (struct overlaps){ 0, 0, NULL };
  size_t nmax = 1550;
  tree = malloc(nmax * sizeof(tree[0]));
  if (tree)
  {
    size_t i;
    for (i = 0; i < ntrees; ++i)
    {
      size_t n = random()%(nmax + 1);
      do_test(tree, n, &overlaps);
    }
  }

  free(overlaps.overlap_iv);
  free(tree);

  return TEST_EXIT_CODE;
}

static inline struct interval rand_interval() {
  double x, y;
  struct interval iv;
  x = (double)(random() - RAND_MAX/2);
  y = (double)(x + random());
  iv.left = x;
  iv.right = y;
  return iv;
}


static void do_test(struct interval_node *tree, size_t num_nodes,
                    struct overlaps *overlaps)
{
  size_t j, noverlaps;
  unsigned i;
  for (i = 0; i < num_nodes; ++i)
    tree[i].range = rand_interval();
  yac_generate_interval_tree(tree, num_nodes);
  for (i = 0; i < ntests; ++i)
  {
    struct interval iv;
    iv = rand_interval();
    overlaps->num_overlaps = 0;
    yac_search_interval_tree(tree, num_nodes, iv, overlaps);
    for (j = 0; j < overlaps->num_overlaps; ++j)
      if (!overlap_test(iv, tree[overlaps->overlap_iv[j]].range))
        PUT_ERR("overlap result does not overlap\n");
    noverlaps = 0;
    for (j = 0; j < num_nodes; j++)
      noverlaps += overlap_test(iv, tree[j].range) != 0;
    if (noverlaps != overlaps->num_overlaps)
      PUT_ERR("overlap count mismatch\n");
  }
}

