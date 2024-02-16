/**
 * @file xt_mpi.h
 * @brief utility routines for MPI
 *
 * contains utility routines for handling MPI
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
 *
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

#ifndef XT_MPI_H
#define XT_MPI_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_core.h"

/** \example test_mpi_generate_datatype.c
 */

/**
 * check return code of MPI call and call abort function if needed
 */
#define xt_mpi_call(call, comm)                 \
  do {                                          \
    int error_code = (call);                    \
    if (error_code != MPI_SUCCESS)              \
      xt_mpi_error(error_code, comm);           \
  } while(0)

/**
 * report error return of MPI call
 *
 * @param[in] error_code return code of an MPI call
 * @param[in] comm       communicator which was used for the respective MPI call
 */
void xt_mpi_error(int error_code, MPI_Comm comm);

/**
 * generates an MPI datatype
 *
 * @param[in] displacements array of displacements
 * @param[in] count         number of elements
 * @param[in] old_type      base MPI datatype of all elements
 * @param[in] comm          MPI communicator
 * @return                  MPI datatype for the given data layout
 *
 * \remarks the returned datatype needs to be freed by the user (MPI_Type_free)
 */
MPI_Datatype xt_mpi_generate_datatype(int const * displacements, int count,
                                      MPI_Datatype old_type, MPI_Comm comm);

/**
 * generates an MPI datatype
 *
 * @param[in] displacements array of displacements
 * @param[in] blocklengths  array of block sizes
 * @param[in] count         number of blocks
 * @param[in] old_type      base MPI datatype of all elements
 * @param[in] comm          MPI communicator
 * @return                  MPI datatype for the given data layout
 *
 * \remarks the returned datatype needs to be freed by the user (MPI_Type_free)
 */
MPI_Datatype xt_mpi_generate_datatype_block(const int *displacements,
                                            const int *blocklengths,
                                            int count, MPI_Datatype old_type,
                                            MPI_Comm comm);

/**
 * generates an MPI datatype
 *
 * @param[in] v             array of displacement, stride and length
 *                          of sub-vectors
 * @param[in] count         number of sub-vectors
 * @param[in] old_type      base MPI datatype of all elements
 * @param[in] comm          MPI communicator
 * @return                  MPI datatype for the given data layout
 *
 * \remarks the returned datatype needs to be freed by the user (MPI_Type_free)
 */
MPI_Datatype xt_mpi_generate_datatype_stripe(const struct Xt_offset_ext *v,
                                             int count, MPI_Datatype old_type,
                                             MPI_Comm comm);


/**
 * Annotate communicator that is for exclusive use by YAXT.
 *
 * YAXT will trust communicators marked this way to not have active
 * communications other than those initiated by YAXT functions. This
 * is useful to prevent unnecessary MPI_Comm_dup() calls.
 *
 * @param[in,out] comm communicator that will not be used by
 *                    application code
 */
void
xt_mpi_comm_mark_exclusive(MPI_Comm comm);

#endif // XT_MPI_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
