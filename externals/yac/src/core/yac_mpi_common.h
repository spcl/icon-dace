// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef YAC_MPI_COMMON_H
#define YAC_MPI_COMMON_H

#include <mpi.h>

/**
 * check return code of MPI call and call abort function if needed
 */
#ifndef YAC_CODE_COVERAGE_TEST
#define yac_mpi_call(call, comm)                \
  do {                      \
    int error_code = (call);                    \
    if (error_code != MPI_SUCCESS)              \
      yac_mpi_error(error_code, comm);          \
  } while(0)
#else
#define yac_mpi_call(call, comm) do {(call);} while(0)
#endif

/**
 * report error return of MPI call
 *
 * @param[in] error_code return code of an MPI call
 * @param[in] comm       communicator which was used for the respective MPI call
 */

void yac_mpi_error(int error_code, MPI_Comm comm);

#endif // YAC_MPI_COMMON_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
