/**
 * @file io_utils.h
 *
 * @copyright Copyright  (C)  2020 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#ifndef IO_UTILS_H
#define IO_UTILS_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef YAC_NETCDF_ENABLED
#include <netcdf.h>
#endif

#include "yac_mpi.h"

void yac_get_io_ranks(
  MPI_Comm comm, int * local_is_io, int ** io_ranks, int * num_io_ranks);

void yac_nc_open(const char * path, int omode, int * ncidp);
void yac_nc_create(const char * path, int cmode, int * ncidp);
void yac_nc_inq_varid(int ncid, char const * name, int * varidp);

#ifdef YAC_NETCDF_ENABLED
#define HANDLE_ERROR(exp) \
  do { \
    int handle_error_status = (exp); \
    YAC_ASSERT_F( \
      (handle_error_status) == NC_NOERR, \
      "%s", nc_strerror(handle_error_status)) \
  } while(0)
#else
#define HANDLE_ERROR(exp) \
  YAC_ASSERT( \
    0, \
    "YAC was compiled without support for NETCDF"))
#endif

#endif // IO_UTILS_H
