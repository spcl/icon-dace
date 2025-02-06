/**
 * @file idxsection_examples.c
 *
 * @copyright Copyright  (C)  2012 Moritz Hanke <hanke@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

  { // test 2D section with negative global size

    for (int i = 0; i < 16; ++i) {


      Xt_int start = 0;
      enum { num_dimensions = 2 };
      static const Xt_int global_size[4][2] = {{5,10},{5,-10},{-5,10},{-5,-10}};
      static const int local_size [4][2] = {{3,4},{3,-4},{-3,4},{-3,-4}};
      static const Xt_int local_start[2] = {1,2};

      // create index section

      Xt_idxlist idxsection
        = xt_idxsection_new(start, num_dimensions, global_size[i >> 2],
                            local_size[i & 3], local_start);

      // testing

      printf("global size (x=%3d y=%2d) local size (x=%2d y=%2d): ",
             (int)global_size[i >> 2][1], (int)global_size[i >> 2][0],
             (int)local_size[i & 3][1], (int)local_size[i & 3][0]);
      print_index_list(idxsection);

      // clean up

      xt_idxlist_delete(idxsection);
    }
  }

  { // test 2D section with stride in x

    //! [custom stride in x]

    Xt_int start = 0;
    enum { num_dimensions = 3 };
    static const Xt_int global_size[3] = {5,5,2};
    static const int local_size [3] = {3,4,1};
    static const Xt_int local_start[3] = {2,0,1};

    // create index section

    Xt_idxlist idxsection
      = xt_idxsection_new(start, num_dimensions, global_size,
                          local_size, local_start);
    //! [custom stride in x]

    // testing

    printf("stride 2 in x: ");
    print_index_list(idxsection);

    // clean up

    xt_idxlist_delete(idxsection);
  }

  { // test 2D section with stride in x and y

    //! [custom stride in x and y]

    Xt_int start = 0;
    enum { num_dimensions = 4 };
    static const Xt_int global_size[4] = {3,2,5,2};
    static const int local_size [4] = {3,1,4,1};
    static const Xt_int local_start[4] = {0,1,1,0};

    // create index section

    Xt_idxlist idxsection
      = xt_idxsection_new(start, num_dimensions, global_size,
                          local_size, local_start);
    //! [custom stride in x and y]

    // testing

    printf("stride 2 in x and y: ");
    print_index_list(idxsection);

    // clean up

    xt_idxlist_delete(idxsection);
  }
}

/*
 * Local Variables:
 * coding: utf-8
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
 * license-default: "bsd"
 * End:
 */
