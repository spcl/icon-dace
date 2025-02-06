// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef IO_UTILS_H
#define IO_UTILS_H

#include "yac_config.h"

#ifdef YAC_NETCDF_ENABLED
#include <netcdf.h>
#endif

#include <mpi.h>

#include "yac_assert.h"

// YAC PUBLIC HEADER START

int yac_file_exists(const char * filename);

void yac_get_io_ranks(
  MPI_Comm comm, int * local_is_io, int ** io_ranks, int * num_io_ranks);

void yac_nc_open(const char * path, int omode, int * ncidp);
void yac_nc_create(const char * path, int cmode, int * ncidp);
void yac_nc_inq_dimid(int ncid, char const * name, int * dimidp);
void yac_nc_inq_varid(int ncid, char const * name, int * varidp);

#define YAC_HANDLE_ERROR(exp) \
  do { \
    int handle_error_status = (exp); \
    YAC_ASSERT_F( \
      (handle_error_status) == NC_NOERR, \
      "%s", nc_strerror(handle_error_status)) \
  } while(0)

// YAC PUBLIC HEADER STOP

#endif // IO_UTILS_H
