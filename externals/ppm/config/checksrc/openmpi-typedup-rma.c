/**
 * @file openmpi-typedup-rma.c
 * @brief detects problem some versions of OpenMPI have with
 * MPI_Type_dup and subsequent RMA get operations
 *
 * @copyright Copyright  (C)  2014 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
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
 *
 * acx_mpirun_num_tasks = 2
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stddef.h>
#include <stdlib.h>

#include <mpi.h>

#define xmpi(rc)                                \
  do {                                          \
    if (rc != MPI_SUCCESS) abort();             \
  } while (0)

static inline void *
xmalloc(size_t size)
{
  void *p = malloc(size);
  if (size && !p)
    abort();
  return p;
}

int main(int argc, char **argv)
{
  xmpi(MPI_Init(&argc, &argv));
  int world_size, world_rank;
  xmpi(MPI_Comm_size(MPI_COMM_WORLD, &world_size));
  xmpi(MPI_Comm_rank(MPI_COMM_WORLD, &world_rank));
  enum { ni = 16, nd = 8, nb = 2 };
  struct structExample { int i[ni];
    double d[nd];
  } *my_mem, *collected;
  collected = xmalloc(sizeof (*collected) * (size_t)world_size);
  xmpi(MPI_Alloc_mem(sizeof (*my_mem), MPI_INFO_NULL, &my_mem));
  MPI_Datatype sparse_dt;
  {
    static const int bl[nb] = { 1, 2 };
    static const MPI_Aint disp[nb] = { 0, offsetof(struct structExample, d) };
    MPI_Datatype t[2];
    xmpi(MPI_Type_dup(MPI_INT, t + 0));
    xmpi(MPI_Type_dup(MPI_DOUBLE, t + 1));
    xmpi(MPI_Type_create_struct(nb, (int*)bl, (MPI_Aint*)disp, t, &sparse_dt));
    xmpi(MPI_Type_commit(&sparse_dt));
    xmpi(MPI_Type_free(t + 1));
    xmpi(MPI_Type_free(t + 0));
  }
  MPI_Win win;
  xmpi(MPI_Win_create(my_mem, sizeof (*my_mem), 1, MPI_INFO_NULL,
                      MPI_COMM_WORLD, &win));

  for (int rank = 0; rank < world_size; ++rank)
  {
    xmpi(MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK, win));
    xmpi(MPI_Get(collected + rank, 1, sparse_dt, rank, 0, 1, sparse_dt, win));
    xmpi(MPI_Win_unlock(rank, win));
  }

  xmpi(MPI_Win_free(&win));
  xmpi(MPI_Type_free(&sparse_dt));
  xmpi(MPI_Free_mem(my_mem));
  free(collected);
  xmpi(MPI_Finalize());
  return EXIT_SUCCESS;
}
