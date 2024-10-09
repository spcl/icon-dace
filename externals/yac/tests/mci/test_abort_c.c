// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include "tests.h"
#include "yac.h"

static void custom_error_handler(
  MPI_Comm comm, const char * msg, const char * source, int line);

char const * ref_msg = "reference message";
char const * ref_source = __FILE__;
int const ref_line = __LINE__;
MPI_Comm ref_comm_world;

int main (void) {

  MPI_Init(NULL, NULL);

  MPI_Comm_dup(MPI_COMM_WORLD, &ref_comm_world);

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
  MPI_Comm_compare(comm, ref_comm_world, &result);
  if (result != MPI_IDENT) PUT_ERR("error in yac_abort_message (comm)");

  if (strcmp(msg, ref_msg)) PUT_ERR("error in yac_abort_message (msg)");

  if (strcmp(source, ref_source)) PUT_ERR("error in yac_abort_message (source)");

  if (line != ref_line) PUT_ERR("error in yac_abort_message (line)");

  // MPI_Abort may yield non-zero error codes of mpirun hence we
  // terminate the programm gracefully
  MPI_Finalize();
  exit(TEST_EXIT_CODE);
}
