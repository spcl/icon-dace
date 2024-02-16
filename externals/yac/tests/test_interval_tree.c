/**
 * @file test_interval_tree.c
 *
 * @copyright Copyright  (C)  2015 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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
#include "tests.h"
#include "interval_tree.h"
#include "utils.h"

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
  tree = xmalloc(nmax * sizeof(tree[0]));
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

