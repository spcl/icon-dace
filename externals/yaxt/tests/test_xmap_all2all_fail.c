/**
 * @file test_xmap_all2all_fail.c
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

#include <assert.h>
#include <string.h>
#ifdef XT_NEED_MPI_ABORT_WORK_AROUND
#  include <fcntl.h>
#  include <sys/stat.h>
#  include <sys/types.h>
#endif
#include <unistd.h>

#include <mpi.h>

#include <yaxt.h>

#define VERBOSE
#include "tests.h"
#include "ctest_common.h"
#include "test_xmap_common.h"

/* If we're not using GNU C, elide __attribute__ */
#ifndef __GNUC__
#  define  __attribute__(x)  /*NOTHING*/
#endif

static enum test_idxlist_size index_list_size = SMALL;

static void
parse_options(int *argc, char ***argv);

typedef void (*Xt_abort_func)(MPI_Comm comm, const char *msg,
                              const char *source, int line)
  __attribute__((noreturn));
extern Xt_abort_func Xt_abort;

static Xt_xmap
xmap_new_fail3(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist, MPI_Comm comm);



int main(int argc, char **argv) {

  MPI_Comm comm = MPI_COMM_WORLD;
  test_init_mpi(&argc, &argv, comm);

  xt_initialize(comm);

  int my_rank;
  MPI_Comm_rank(comm, &my_rank);

  parse_options(&argc, &argv);

  {
    // source index list
    struct Xt_stripe src_stripe;
    src_stripe.nstrides = (index_list_size == SMALL)?7:1023;
    src_stripe.start = (Xt_int)(1 + (Xt_int)my_rank * src_stripe.nstrides);
    src_stripe.stride = 1;

    // destination index list
    struct Xt_stripe dst_stripe;
    dst_stripe.nstrides = src_stripe.nstrides;
    dst_stripe.start = (Xt_int)(src_stripe.start + src_stripe.nstrides);
    dst_stripe.stride = -1;

    test_self_xmap_construct_idxstripes(&src_stripe, 1, &dst_stripe, 1,
                                        xmap_new_fail3, comm);
  }

  xt_finalize();
  xt_mpi_call(MPI_Finalize(), comm);

  return TEST_EXIT_CODE;
}

static void
xfail_abort(MPI_Comm comm, const char *msg, const char *source, int line)
  __attribute__((noreturn));

static Xt_xmap
xmap_new_fail3(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist, MPI_Comm comm)
{
  Xt_abort_func orig_xt_abort = Xt_abort;
  // test of exchange map constructor error handling
  Xt_abort = xfail_abort;
  // NOTE: this is the call which should fail
  Xt_xmap xmap
    = xt_xmap_all2all_new(src_idxlist, dst_idxlist, comm);
  /* this position should not be reached */
  Xt_abort = orig_xt_abort;
  MPI_Abort(MPI_COMM_WORLD, 1);
  return xmap;
}


static void
parse_options(int *argc, char ***argv)
{
  int opt;
  while ((opt = getopt(*argc, *argv, "s:")) != -1) {
    switch (opt) {
    case 's':
      if (!strcmp(optarg, "small"))
        index_list_size = SMALL;
      else if (!strcmp(optarg, "big"))
        index_list_size = BIG;
      else
      {
        fprintf(stderr, "Unknown data size \"%s\"\n", optarg);
        exit(EXIT_FAILURE);
      }
    }
  }
}

static void
xfail_abort(MPI_Comm comm, const char *msg, const char *source, int line)
{
  fprintf(stderr, "Fatal error in %s, line %d: %s\n", source, line, msg);
#ifdef XT_NEED_MPI_ABORT_WORK_AROUND
  static const char exit_msg[] = "MPI_Abort(0xdeadbeef, 3)\n";
  int fd;
#  if XT_NEED_MPI_ABORT_WORK_AROUND == 1
  fd = STDERR_FILENO;
#  elif XT_NEED_MPI_ABORT_WORK_AROUND == 2
  fd = open("test_xmap_all2all_fail.result.txt",
            O_WRONLY | O_TRUNC | O_CREAT | O_NOCTTY, 0777);
#  endif
  write(fd, exit_msg, sizeof (exit_msg)-1);
#  if XT_NEED_MPI_ABORT_WORK_AROUND == 2
  close(fd);
#  endif
#endif
  MPI_Abort(comm, 3);
  abort();
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
