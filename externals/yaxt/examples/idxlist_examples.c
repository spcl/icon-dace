/**
 * @file idxlist_examples.c
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
#include <stdio.h>

#include <yaxt.h>

#include "print_index_list.h"

int main(void) {

  {
    //! [index list example]
    Xt_int indices[] = {0,3,32,26,44,14,48};
    int num_indices = sizeof(indices)/sizeof(indices[0]);

    Xt_idxlist idxvec = xt_idxvec_new(indices, num_indices);
    //! [index list example]

    puts("index vector example:");
    print_index_list(idxvec);

    xt_idxlist_delete(idxvec);
  }

  {
    //! [index stripes example]
    struct Xt_stripe stripes[] = {{.start = 0, .nstrides = 17, .stride = 3}};
    int num_stripes = sizeof(stripes)/sizeof(stripes[0]);

    Xt_idxlist idxstripes = xt_idxstripes_new(stripes, num_stripes);
    //! [index stripes example]

    puts("index stripes example:");
    print_index_list(idxstripes);

    xt_idxlist_delete(idxstripes);
  }

  {
    //! [index section example]
    Xt_int start = 0;
    int num_dimensions = 2;
    Xt_int global_size[2] = {5,10};
    int local_size[2] = {3,4};
    Xt_int local_start[2] = {2,3};

    Xt_idxlist idxsection
      = xt_idxsection_new(start, num_dimensions, global_size,
                          local_size, local_start);
    //! [index section example]

    puts("index section example:");
    print_index_list(idxsection);

    xt_idxlist_delete(idxsection);
  }

  {
    //! [index list collection example]
    Xt_int start = 0;
    int num_dimensions = 2;
    Xt_int global_size[2] = {5,10};
    int local_size[2][2] = {{2,2},{3,4}};
    Xt_int local_start[2][2] = {{0,0},{2,6}};

    Xt_idxlist idxlists[]
      = {xt_idxsection_new(start, num_dimensions, global_size,
                           local_size[0], local_start[0]),
         xt_idxsection_new(start, num_dimensions, global_size,
                           local_size[1], local_start[1])};
    int num_idxlists = sizeof(idxlists) / sizeof(idxlists[0]);

    Xt_idxlist idxcollection
      = xt_idxlist_collection_new(idxlists, num_idxlists);
    //! [index list collection example]

    puts("index list collection example:");
    print_index_list(idxcollection);

    for (int i = 0; i < num_idxlists; ++i)
      xt_idxlist_delete(idxlists[i]);
    xt_idxlist_delete(idxcollection);
  }

  {
    //![empty index list example]
    Xt_idxlist empty_list;

    empty_list = xt_idxempty_new();
    //![empty index list example]

    xt_idxlist_delete(empty_list);
  }

  {
    //![index list modifier example]
    Xt_int start = 0;
    int num_dimensions = 2;
    Xt_int global_size[2] = {5,10};
    int local_size[2] = {3,3};
    Xt_int local_start[2] = {2,0};

    Xt_idxlist source_list = xt_idxsection_new(start, num_dimensions,
                                               global_size, local_size,
                                               local_start);

    struct Xt_stripe extract_stripes[2] = {{.start = 0, .nstrides = 5, .stride = 10},
                                           {.start = 9, .nstrides = 5, .stride = 10}};
    struct Xt_stripe substitute_stripes[2] = {{.start = 8, .nstrides = 5, .stride = 10},
                                              {.start = 1, .nstrides = 5, .stride = 10}};
    struct Xt_modifier modifier[] = {{.extract = xt_idxstripes_new(extract_stripes, 2),
                                      .subst   = xt_idxstripes_new(substitute_stripes, 2)}};
    int num_modifier = sizeof(modifier) / sizeof(modifier[0]);

    Xt_idxlist target_list = xt_idxmod_new(source_list, modifier, num_modifier, NULL);
    //![index list modifier example]

    puts("index modifier example:");
    puts("    source_list:");
    print_index_list(source_list);
    puts("    target_list:");
    print_index_list(target_list);

    for (int i = 0; i < num_modifier; ++i) {
      xt_idxlist_delete(modifier[i].extract);
      xt_idxlist_delete(modifier[i].subst);
    }
    xt_idxlist_delete(source_list);
    xt_idxlist_delete(target_list);
  }

  return 0;
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
