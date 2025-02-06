// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MPI_HANDSHAKE_H
#define MPI_HANDSHAKE_H

#include <mpi.h>

/* Collective algorithm that splits a given MPI communicator into n
 * communicators. Each new communicator contains all processes that
 * provided the same group name.
 * @param[in]  comm        MPI communicator that is to be split
 * @param[in]  n           number communicator that are to be generated
 * @param[in]  group_name  group names for the new communicators
 * @param[out] group_comms new communicators
 */
void yac_mpi_handshake(
  MPI_Comm comm, size_t n, char const** group_names,
  MPI_Comm * group_comms);

#endif // MPI_HANDSHAKE_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
