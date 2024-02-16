/**
 * @file test_abort_c.c
 *
 * @copyright Copyright  (C)  2021 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 *
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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include "tests.h"
#include "yac_interface.h"
#include "yac_mpi.h"

static void custom_error_handler(
  MPI_Comm comm, const char * msg, const char * source, int line);

char const * ref_msg = "reference message";
char const * ref_source = __FILE__;
int const ref_line = __LINE__;
MPI_Comm ref_comm_world;

int main (void) {

  MPI_Init(NULL, NULL);

  yac_mpi_call(MPI_Comm_dup(MPI_COMM_WORLD, &ref_comm_world), MPI_COMM_WORLD);

  if (yac_get_abort_handler() !=
      yac_get_default_abort_handler())
    PUT_ERR("error in yac_get_abort_handler/yac_get_default_abort_handler");

  yac_set_default_comm(ref_comm_world);
  yac_set_abort_handler((yac_abort_func)custom_error_handler);

  if (yac_get_abort_handler() != custom_error_handler)
    PUT_ERR("error in yac_get_abort_handler");

  yac_restore_default_abort_handler();

  if (yac_get_abort_handler() !=
      yac_get_default_abort_handler())
    PUT_ERR("error in yac_get_abort_handler/yac_get_default_abort_handler");

  yac_set_abort_handler((yac_abort_func)custom_error_handler);

  yac_abort_message(ref_msg, ref_source, ref_line);

  // test should never reach this point
  PUT_ERR("yac_abort_default did not abort program");

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void custom_error_handler(
  MPI_Comm comm, const char * msg, const char * source, int line) {

  int result;
  yac_mpi_call(
    MPI_Comm_compare(comm, ref_comm_world, &result),
    MPI_COMM_WORLD);
  if (result != MPI_IDENT) PUT_ERR("error in yac_abort_message (comm)");

  if (strcmp(msg, ref_msg)) PUT_ERR("error in yac_abort_message (msg)");

  if (strcmp(source, ref_source)) PUT_ERR("error in yac_abort_message (source)");

  if (line != ref_line) PUT_ERR("error in yac_abort_message (line)");

  // MPI_Abort may yield non-zero error codes of mpirun hence we
  // terminate the programm gracefully
  MPI_Finalize();
  exit(TEST_EXIT_CODE);
}
