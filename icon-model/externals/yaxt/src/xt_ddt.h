/**
 * @file xt_ddt.h
 * @brief utility routines for manual handling of MPI DDT's
 *
 * contains utility routines for handling manual MPI DDT's
 *
 * @copyright Copyright  (C)  2022 Jörg Behrens <behrens@dkrz.de>
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

#ifndef XT_DDT_H
#define XT_DDT_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "core/ppm_visibility.h"

typedef struct Xt_ddt_ *Xt_ddt;

/** \example test_ddt.c
 */

/**
 * gets the buffer size required to pack the data in mpi_ddt
 * @param[in] mpi_ddt MPI Datatype
 * @returns required packing buffer size
 */
PPM_DSO_INTERNAL size_t
xt_ddt_get_pack_size(MPI_Datatype mpi_ddt);

/**
 * packs the data from the source buffer into destination buffer, based on
 * the data layout described by MPI datatype
 * @param[in]  mpi_ddt MPI datatype
 * @param[in]  src     source buffer
 * @param[out] dst     destination buffer
 */
PPM_DSO_INTERNAL void
xt_ddt_pack(MPI_Datatype mpi_ddt, void const * src, void * dst);

/**
 * unpacks the data from the source buffer into destination buffer, based on
 * the data layout described by MPI datatype
 * @param[in]  mpi_ddt MPI datatype
 * @param[in]  src     source buffer
 * @param[out] dst     destination buffer
 */
PPM_DSO_INTERNAL void
xt_ddt_unpack(MPI_Datatype mpi_ddt, void const * src, void * dst);

#endif // XT_DDT_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
