/**
 * @file test_redirstdout_c.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "tests.h"
#include "test_common.h"

#include "yac_interface.h"

static int get_parallel(int argc, char* argv[]);
static void delete_output_files(int comm_size);
static void check_output_files(int parallel, int comm_size);

int main(int argc, char* argv[]) {

  MPI_Init(NULL, NULL);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  if (comm_rank == 0) delete_output_files(comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  int parallel = get_parallel(argc, argv);

  yac_redirstdout("test_redirstdout_c", parallel, comm_rank, comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  check_output_files(parallel, comm_size);

  MPI_Finalize();

  if (comm_rank == 0) delete_output_files(comm_size);

  return TEST_EXIT_CODE;
}

static int get_parallel(int argc, char* argv[]) {

  YAC_ASSERT(
    argc == 2, "ERROR(get_parallel): wrong number of arguments (has to be 1)")
  YAC_ASSERT(
    strlen(argv[1]) == 1,
    "ERROR(get_parallel): wrong argument length (has to be 1)")
  YAC_ASSERT(
    (argv[1][0] == 'T') || (argv[1][0] == 'F'),
    "ERROR(get_parallel): wrong argument (has to be \"T\" or \"F\")")

  return (argv[1][0] == 'T');
}

static void delete_output_files(int comm_size) {

  char filename[128];

  unlink("test_redirstdout_c.log");
  unlink("test_redirstdout_c.err");

  for (int i = 0; i < comm_size; ++i) {
    sprintf(filename, "test_redirstdout_c.%d", i);
    unlink(filename);
    sprintf(filename, "test_redirstdout_c.err.%d", i);
    unlink(filename);
  }
}

static void check_output_files(int parallel, int comm_size) {

  char filename[128];

  if (parallel == 0) {

    YAC_ASSERT(
      yac_file_exists("test_redirstdout_c.log"),
      "test_redirstdout_c.log does not exist")
    YAC_ASSERT(
      yac_file_exists("test_redirstdout_c.err"),
      "test_redirstdout_c.err does not exist")

    for (int i = 0; i < comm_size; ++i) {

      sprintf(filename, "test_redirstdout_c.%d", i);
      YAC_ASSERT(
        !yac_file_exists(filename), "test_redirstdout_c.* should not exist")
      sprintf(filename, "test_redirstdout_c.err.%d", i);
      YAC_ASSERT(
        !yac_file_exists(filename), "test_redirstdout_c.err.* should not exist")
    }

  } else {

    YAC_ASSERT(
      !yac_file_exists("test_redirstdout_c.log"),
      "test_redirstdout_c.log should not exist")
    YAC_ASSERT(
      !yac_file_exists("test_redirstdout_c.err"),
      "test_redirstdout_c.err should not exist")

    for (int i = 0; i < comm_size; ++i) {

      sprintf(filename, "test_redirstdout_c.%d", i);
      YAC_ASSERT(
        yac_file_exists(filename), "test_redirstdout_c.* does not exist")
      sprintf(filename, "test_redirstdout_c.err.%d", i);
      YAC_ASSERT(
        yac_file_exists(filename), "test_redirstdout_c.err.* does not exist")
    }
  }
}
