/**
 * @file perf_idxsection_get_positions_of_indices.c
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

#include <inttypes.h>
#include <limits.h>
#include <mpi.h>
#include <yaxt.h>
#include "core/ppm_xfuncs.h"

int main(void) {

#if INT_MAX >= INT32_C(512) * 256 * 32

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);
  xt_initialize(MPI_COMM_WORLD);

  int nstrides = 512 * 256 * 32;
  struct Xt_stripe stripe[1] = {{.start = 0, .stride = 1, .nstrides = nstrides}};
  Xt_idxlist selection_indices = xt_idxstripes_new(stripe, 1);

  int * positions = xmalloc((size_t)nstrides * sizeof(*positions));
  int i, j;
  int step = nstrides / 128;
  int single_match_only;
  int num_samples = 1;

  Xt_int global_size[3][5] = {{4, 64, 16, 256, 512},
                              {4, 16, 32, 128, 256},
                              {4, 8, 16, 64, 128}};
  int local_size [3][5] = {{4, 1, 16, 128, 512},
                           {4, 1, 32, 64, 256},
                           {4, 1, 16, 32, 128}};

  for (int test_id = 0; test_id < 3; ++test_id) {

    { // test with cached indices

      Xt_int start = 0;
      int num_dimension = 3;
      Xt_int local_start[] = {0,0,0};
      Xt_idxlist idxsection = xt_idxsection_new(start, num_dimension,
                                                global_size[test_id] + 2,
                                                local_size[test_id] + 2,
                                                local_start);

      // generate index cache in idxsection
      xt_idxlist_get_indices_const(idxsection);

      for (single_match_only = 0; single_match_only < 2;  ++single_match_only)
        for (i = step; i < nstrides; i += step) {
          if (i>1000000) break; //reduce workload
          for (j = 0; j < num_samples; ++j)
            xt_idxlist_get_positions_of_indices(idxsection,
                                                xt_idxlist_get_indices_const(selection_indices),
                                                i, positions, single_match_only);
        }
      xt_idxlist_delete(idxsection);
    }

    { // test without cached indices

      Xt_int start = 0;
      int num_dimension = 3;
      Xt_int local_start[] = {0,0,0};

      for (single_match_only = 0; single_match_only < 2;  ++single_match_only) {
        for (i = step; i < nstrides; i += step) {
          if (i>1000000) break; //reduce workload
          Xt_idxlist idxsection[num_samples];

          for (j = 0; j < num_samples; ++j)
            idxsection[j] = xt_idxsection_new(start, num_dimension,
                                              global_size[test_id] + 2,
                                              local_size[test_id] + 2,
                                              local_start);
          for (j = 0; j < num_samples; ++j)
            xt_idxlist_get_positions_of_indices(idxsection[j],
                                                xt_idxlist_get_indices_const(selection_indices),
                                                i, positions, single_match_only);
          for (j = 0; j < num_samples; ++j)
            xt_idxlist_delete(idxsection[j]);
        }
      }
    }
    for (int num_dimension = 3; num_dimension < 6; num_dimension += 2)
      { // test index enumeration without cached indices

      Xt_int local_start[] = {0,0,0,0,0};
      Xt_int *idxDump = NULL;
      size_t idxDumpSize = 0;
      num_samples = 52;
      Xt_idxlist idxsection[num_samples];

      for (single_match_only = 0; single_match_only < 2;  ++single_match_only) {

        for (j = 0; j < num_samples; ++j) {
          idxsection[j] = xt_idxsection_new(0, num_dimension,
                                            global_size[test_id],
                                            local_size[test_id],
                                            local_start);
        }
        for (j = 0; j < num_samples; ++j)
        {
          int numIdx = xt_idxlist_get_num_indices(idxsection[j]);
          if ((size_t)numIdx > idxDumpSize)
            idxDump = xrealloc(idxDump,
                               sizeof (idxDump[0])
                               * (idxDumpSize = (size_t)numIdx));
          xt_idxlist_get_indices(idxsection[j], idxDump);
        }
        for (j = 0; j < num_samples; ++j)
          xt_idxlist_delete(idxsection[j]);
      }
      num_samples = 1;
      free(idxDump);
    }
  }

  xt_idxlist_delete(selection_indices);
  free(positions);

  xt_finalize();
  xt_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);
#endif
  return 0;
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
