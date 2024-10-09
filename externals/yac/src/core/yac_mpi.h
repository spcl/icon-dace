// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef YAC_MPI_H
#define YAC_MPI_H

// YAC PUBLIC HEADER START

#include <mpi.h>

void yac_mpi_init();
void yac_yaxt_init(MPI_Comm comm);
int yac_mpi_is_initialised();
void yac_mpi_cleanup();
void yac_mpi_finalize();

// YAC PUBLIC HEADER STOP

#endif // YAC_MPI_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
