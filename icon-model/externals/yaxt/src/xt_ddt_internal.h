/**
 * @file xt_ddt_internal.h
 * @brief internal utility routines for manual handling of MPI DDT's
 *
 * contains internal utility routines for handling manual MPI DDT's
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

#ifndef XT_DDT_INTERNAL_H
#define XT_DDT_INTERNAL_H

#include "xt_ddt.h"
#include "xt_gpu.h"

#include "core/ppm_visibility.h"

/**
 * Returns a xt_ddt object for the provided MPI datatype.
 * The routine will first check whether there is already a xt_ddt for this
 * MPI datatype and will retrieve it, if available. Otherwise a new xt_ddt
 * will be generated and cached in the MPI datatype.
 * @param[in] mpi_ddt MPI datatype
 * @return xt_ddt object
 * @remark the returned xt_ddt object is only valid as long as the mpi_ddt
 *         is not freed
 * @remark no \ref xt_ddt_delete needs to be called for the xt_ddt object
 *         returned by this routine
 */
PPM_DSO_INTERNAL Xt_ddt
xt_ddt_from_mpi_ddt(MPI_Datatype mpi_ddt);

/**
 * deletes a xt_ddt object
 * @param[in] ddt xt_ddt object
 */
PPM_DSO_INTERNAL void
xt_ddt_delete(Xt_ddt ddt);

/**
 * increases the internal reference counter for the given xt_ddt object
 * @param[in] ddt xt_ddt object
 * @remark for each reference of a xt_ddt object
 *         \ref xt_ddt_delete has to be called
 */
PPM_DSO_INTERNAL void
xt_ddt_inc_ref_count(Xt_ddt ddt);

/**
 * gets the buffer size required to pack the data in ddt
 * @param[in] ddt xt_ddt object
 * @returns required packing buffer size
 */
PPM_DSO_INTERNAL size_t
xt_ddt_get_pack_size_internal(Xt_ddt ddt);

/**
 * packs the data from the source buffer into destination buffer
 * @param[in]  ddt     xt_ddt object
 * @param[in]  src     source buffer
 * @param[out] dst     destination buffer
 * @param[in]  memtype type of source and destination memory
 * @remark both src and dst have to point to the same type of memory
 */
PPM_DSO_INTERNAL void xt_ddt_pack_internal(
  Xt_ddt ddt, void const * src, void * dst, enum xt_memtype memtype);

/**
 * unpacks the data from the source buffer into destination buffer
 * @param[in]  ddt     xt_ddt object
 * @param[in]  src     source buffer
 * @param[out] dst     destination buffer
 * @param[in]  memtype type of source and destination memory
 * @remark both src and dst have to point to the same type of memory
 */
PPM_DSO_INTERNAL void xt_ddt_unpack_internal(
  Xt_ddt ddt, void const * src, void * dst, enum xt_memtype memtype);

#endif // XT_DDT_INTERNAL_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
