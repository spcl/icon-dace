/**
 * @file test_mpi_rma.c
 * @brief assert that MPI_Get operations are indeed working and
 * complete in a reasonable amount of time
 *
 * @copyright Copyright  (C)  2017 Thomas Jahns <jahns@dkrz.de>
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

#define _XOPEN_SOURCE 600

#include <errno.h>
#include <setjmp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <signal.h>
#include <unistd.h>

#include <mpi.h>

#define xmpi(rc)                                \
  do {                                          \
    if ((rc) != MPI_SUCCESS) abort();           \
  } while (0)

static inline void *
xmalloc(size_t size)
{
  void *p = malloc(size);
  if (size && !p)
    abort();
  return p;
}

static sigjmp_buf errjmpbuf;

static void
handle_alarm(int sig)
{
  (void)sig;
  siglongjmp(errjmpbuf, 1);
}


int main(int argc, char **argv)
{
  xmpi(MPI_Init(&argc, &argv));
  int world_rank, world_size;
  xmpi(MPI_Comm_size(MPI_COMM_WORLD, &world_size));
  xmpi(MPI_Comm_rank(MPI_COMM_WORLD, &world_rank));

  /* establish signal handler for alarm */
  {
    if (sigsetjmp(errjmpbuf, 1))
    {
      fprintf(stderr, "%d: RMA transfers take too long, aborting attempts.\n",
              world_rank);
      _exit(EXIT_FAILURE);
    }
    struct sigaction act;
    act.sa_handler = handle_alarm;
    act.sa_flags = SA_RESETHAND;
#ifdef SA_INTERRUPT
    act.sa_flags |= SA_INTERRUPT;
#endif
    sigemptyset(&act.sa_mask);
    if (sigaction(SIGALRM, &act, NULL))
    {
      int errno_code = errno;
      fprintf(stderr, "%d: unexpected error setting signal handler: %d: %s\n",
              world_rank, errno_code, strerror(errno_code));
      xmpi(MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE));
      _exit(EXIT_FAILURE);
    }
  }
  /* set timer so this test cannot hang indefinitely */
  unsigned int t = alarm(30U);
  if (t != 0)
    fprintf(stderr, "warning: overriding previously set alarm for %u\n", t);

  enum { nbytes = 8192 };
  /* allocated exposed memory */
  void *my_mem;
  xmpi(MPI_Alloc_mem(nbytes, MPI_INFO_NULL, &my_mem));
  MPI_Win win;
  xmpi(MPI_Win_create(my_mem, nbytes, 1, MPI_INFO_NULL,
                      MPI_COMM_WORLD, &win));
  /* fill in data owned by rank */
  xmpi(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, world_rank, MPI_MODE_NOCHECK, win));
  {
    size_t rank_pattern = (size_t)(world_rank - world_size);
    for (size_t i = 0; i < nbytes; ++i)
      ((unsigned char *)my_mem)[i] = (unsigned char)(i ^ rank_pattern);
  }
  xmpi(MPI_Win_unlock(world_rank, win));
  /* wait till all tasks are ready */
  xmpi(MPI_Barrier(MPI_COMM_WORLD));
  /* fetch data from all ranks ... */
  unsigned char *collected = xmalloc(nbytes * (size_t)world_size);
  for (int rank = 0; rank < world_size; ++rank)
  {
    xmpi(MPI_Win_lock(MPI_LOCK_SHARED, rank, MPI_MODE_NOCHECK, win));
    xmpi(MPI_Get(collected + (size_t)rank * nbytes, nbytes, MPI_UNSIGNED_CHAR,
                 rank, 0, nbytes, MPI_UNSIGNED_CHAR, win));
    xmpi(MPI_Win_unlock(rank, win));
  }
  /* ... and verify the data is valid */
  unsigned long diff = 0;
  for (int rank = 0; rank < world_size; ++rank)
  {
    size_t rank_pattern = (size_t)(rank - world_size);
    for (size_t i = 0; i < nbytes; ++i)
      diff |= (collected[rank * nbytes + i] != (unsigned char)(i ^ rank_pattern));
  }
  /* next try put */
  /* disabled until it's actually used in PPM_dist_mult_array */
#if 0
  xmpi(MPI_Barrier(MPI_COMM_WORLD));
  {
    int target_rank = (world_rank + 1)%world_size,
      srank = (world_rank + world_size - 1)%world_size;
    xmpi(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, target_rank, MPI_MODE_NOCHECK, win));
    xmpi(MPI_Put(collected + (size_t)srank * nbytes, nbytes, MPI_UNSIGNED_CHAR,
                 target_rank, (MPI_Aint)0, nbytes, MPI_UNSIGNED_CHAR,
                 win));
    xmpi(MPI_Win_unlock(target_rank, win));
  }
  xmpi(MPI_Barrier(MPI_COMM_WORLD));
  {
    int srank = (world_rank + world_size - 2)%world_size;
    size_t rank_pattern = (size_t)(srank - world_size);
    for (size_t i = 0; i < nbytes; ++i)
      diff
        |= (((unsigned char *)my_mem)[i] != (unsigned char)(i ^ rank_pattern));
  }
#endif
  xmpi(MPI_Win_free(&win));
  xmpi(MPI_Free_mem(my_mem));
  free(collected);
  xmpi(MPI_Allreduce(MPI_IN_PLACE, &diff, 1, MPI_UNSIGNED_LONG, MPI_BOR,
                     MPI_COMM_WORLD));
  xmpi(MPI_Finalize());
  return diff ? EXIT_FAILURE : EXIT_SUCCESS;
}
