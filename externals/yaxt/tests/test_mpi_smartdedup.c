/**
 * @file test_mpi_smartdedup.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
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

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"

#include "../src/xt_mpi_internal.h"

static void comms_construct_destruct(MPI_Comm comm_world);


int main(void)
{
  // init mpi
  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_mpi_init();
  comms_construct_destruct(MPI_COMM_WORLD);
  xt_mpi_comm_mark_exclusive(MPI_COMM_WORLD);
  comms_construct_destruct(MPI_COMM_WORLD);
  xt_mpi_finalize();
  xt_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);

  return TEST_EXIT_CODE;
}

static void comms_construct_destruct(MPI_Comm comm_world)
{
  enum {
    num_comms = 128,
  };
  MPI_Comm comms[num_comms];
  int tag_offsets[num_comms];
  comms[0] = comm_world;
  tag_offsets[0] = 0;
  for (size_t i = 1; i < num_comms; ++i)
    comms[i] = xt_mpi_comm_smart_dup(comms[i-1], tag_offsets+i);
  for (size_t i = 1; i < num_comms; ++i)
    xt_mpi_comm_smart_dedup(comms + i, tag_offsets[i]);
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
