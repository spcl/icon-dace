/**
 * @file xt_redist_p2p.h
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

#ifndef XT_REDIST_P2P_H
#define XT_REDIST_P2P_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_core.h"
#include "xt/xt_xmap.h"
#include "xt/xt_redist.h"
#include "xt/xt_config.h"

/** \example test_redist_p2p.c
 */
/** \example test_redist_p2p_parallel.c
 */

/**
 * constructor for a redistribution using point to point
 * communication for the exchange. Uses default settings.
 * @param[in] xmap     exchange map
 * @param[in] datatype MPI datatype of single element
 *                     in the data to be exchanged
 */
Xt_redist xt_redist_p2p_new(Xt_xmap xmap, MPI_Datatype datatype);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange. Uses custom settings.
 * @param[in] xmap     exchange map
 * @param[in] datatype MPI datatype of single element
 *                     in the data to be exchanged
 * @param[in] config   configuration object for custom settings
 */
Xt_redist
xt_redist_p2p_custom_new(Xt_xmap xmap, MPI_Datatype datatype, Xt_config config);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange. Uses default settings.
 * @param[in] xmap         exchange map
 * @param[in] src_offsets  array containing for all elements in
 *                         the source index list passed to the
 *                         exchange map the position of the respective
 *                         element in the input array passed to the
 *                         exchange routine
 * @param[in] dst_offsets  array containing for all elements in
 *                         the destination index list passed to the
 *                         exchange map the position of the respective
 *                         element in the output array passed to the
 *                         exchange routine
 * @param[in] datatype     MPI datatype of single element
 *                         in the data to be exchanged
 */
Xt_redist xt_redist_p2p_off_new(Xt_xmap xmap, const int *src_offsets,
                                const int *dst_offsets, MPI_Datatype datatype);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange. Uses custom settings.
 * @param[in] xmap         exchange map
 * @param[in] src_offsets  array containing for all elements in
 *                         the source index list passed to the
 *                         exchange map the position of the respective
 *                         element in the input array passed to the
 *                         exchange routine
 * @param[in] dst_offsets  array containing for all elements in
 *                         the destination index list passed to the
 *                         exchange map the position of the respective
 *                         element in the output array passed to the
 *                         exchange routine
 * @param[in] datatype     MPI datatype of single element
 *                         in the data to be exchanged
 * @param[in] config       configuration object for custom settings
 */
Xt_redist
xt_redist_p2p_off_custom_new(Xt_xmap xmap, const int *src_offsets,
                             const int *dst_offsets, MPI_Datatype datatype,
                             Xt_config config);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange. Uses default settings.
 * @param[in] xmap         exchange map
 * @param[in] num_src_ext  number of source extents
 * @param[in] src_extents  array of extents describing offsets for every
 *                         element of the index lists composing the xmap,
 *                         i.e. { 10, 5, 1 } denotes 5 offsets,
 *                         namely 10, 11, 12, 13, 14
 * @param[in] num_dst_ext  number of destination extents
 * @param[in] dst_extents  array of extents analogous to src_extents
 * @param[in] datatype     MPI datatype of single element
 *                         in the data to be exchanged
 */
Xt_redist xt_redist_p2p_ext_new(Xt_xmap xmap,
                                int num_src_ext,
                                const struct Xt_offset_ext src_extents[],
                                int num_dst_ext,
                                const struct Xt_offset_ext dst_extents[],
                                MPI_Datatype datatype);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange. Uses custom settings.
 * @param[in] xmap         exchange map
 * @param[in] num_src_ext  number of source extents
 * @param[in] src_extents  array of extents describing offsets for every
 *                         element of the index lists composing the xmap,
 *                         i.e. { 10, 5, 1 } denotes 5 offsets,
 *                         namely 10, 11, 12, 13, 14
 * @param[in] num_dst_ext  number of destination extents
 * @param[in] dst_extents  array of extents analogous to src_extents
 * @param[in] datatype     MPI datatype of single element
 *                         in the data to be exchanged
 * @param[in] config       configuration object for custom settings
 */
Xt_redist
xt_redist_p2p_ext_custom_new(Xt_xmap xmap,
                             int num_src_ext,
                             const struct Xt_offset_ext src_extents[],
                             int num_dst_ext,
                             const struct Xt_offset_ext dst_extents[],
                             MPI_Datatype datatype,
                             Xt_config config);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange. Uses default settings.
 * @param[in] xmap         exchange map
 * @param[in] num_src_ext  number of source extents
 * @param[in] src_extents  array of extents describing byte offsets for every
 *                         element of the index lists composing the xmap,
 *                         i.e. { 16, 5, 4 } denotes 5 offsets,
 *                         namely 16, 20, 24, 28, 32
 * @param[in] num_dst_ext  number of destination extents
 * @param[in] dst_extents  array of extents analogous to src_extents
 * @param[in] datatype     MPI datatype of single element
 *                         in the data to be exchanged
 */
Xt_redist xt_redist_p2p_aext_new(Xt_xmap xmap,
                                 int num_src_ext,
                                 const struct Xt_aoffset_ext src_extents[],
                                 int num_dst_ext,
                                 const struct Xt_aoffset_ext dst_extents[],
                                 MPI_Datatype datatype);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange. Uses custom settings.
 * @param[in] xmap         exchange map
 * @param[in] num_src_ext  number of source extents
 * @param[in] src_extents  array of extents describing byte offsets for every
 *                         element of the index lists composing the xmap,
 *                         i.e. { 16, 5, 4 } denotes 5 offsets,
 *                         namely 16, 20, 24, 28, 32
 * @param[in] num_dst_ext  number of destination extents
 * @param[in] dst_extents  array of extents analogous to src_extents
 * @param[in] datatype     MPI datatype of single element
 *                         in the data to be exchanged
 * @param[in] config       configuration object for custom settings
 */
Xt_redist
xt_redist_p2p_aext_custom_new(Xt_xmap xmap,
                              int num_src_ext,
                              const struct Xt_aoffset_ext src_extents[],
                              int num_dst_ext,
                              const struct Xt_aoffset_ext dst_extents[],
                              MPI_Datatype datatype,
                              Xt_config config);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange, special case: elements (which
 * correspond to each idxlist element) are blocks of variable length
 * with corresponding offsets, therefore src_block_num and
 * dst_block_num must match the lengths of the src/dst index lists
 * used for the construction of xmap
 *
 * @param[in] xmap              exchange map
 * @param[in] src_block_offsets array containing for all source index space positions
 *                              of xmap the offsets for blocks in data space
 * @param[in] src_block_sizes   source block lengths in unit of elements
 * @param[in] src_block_num     number of src blocks
 * @param[in] dst_block_offsets array containing for all destination index space positions
 *                              of xmap the offsets for blocks in data space
 * @param[in] dst_block_sizes   destination block lengths in unit of elements
 * @param[in] dst_block_num     number of dst blocks
 * @param[in] datatype          MPI datatype of a single element in data space,
 *                              all elements have the same
 * @remarks NULL offsets arguments mean: use zero based prefix sum of block sizes as
 *          effective offsets
 */
Xt_redist xt_redist_p2p_blocks_off_new(Xt_xmap xmap,
                                       const int *src_block_offsets,
                                       const int *src_block_sizes,
                                       int src_block_num,
                                       const int *dst_block_offsets,
                                       const int *dst_block_sizes,
                                       int dst_block_num,
                                       MPI_Datatype datatype);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange, special case: elements (which
 * correspond to each idxlist element) are blocks of variable length
 * with corresponding offsets, therefore src_block_num and
 * dst_block_num must match the lengths of the src/dst index lists
 * used for the construction of xmap. Uses custom settings.
 *
 * @param[in] xmap              exchange map
 * @param[in] src_block_offsets array containing for all source index space positions
 *                              of xmap the offsets for blocks in data space
 * @param[in] src_block_sizes   source block lengths in unit of elements
 * @param[in] src_block_num     number of src blocks
 * @param[in] dst_block_offsets array containing for all destination index space positions
 *                              of xmap the offsets for blocks in data space
 * @param[in] dst_block_sizes   destination block lengths in unit of elements
 * @param[in] dst_block_num     number of dst blocks
 * @param[in] datatype          MPI datatype of a single element in data space,
 *                              all elements have the same
 * @param[in] config   configuration object for custom settings
 * @remarks NULL offsets arguments mean: use zero based prefix sum of block sizes as
 *          effective offsets
 */
Xt_redist
xt_redist_p2p_blocks_off_custom_new(Xt_xmap xmap,
                                    const int *src_block_offsets,
                                    const int *src_block_sizes,
                                    int src_block_num,
                                    const int *dst_block_offsets,
                                    const int *dst_block_sizes,
                                    int dst_block_num,
                                    MPI_Datatype datatype,
                                    Xt_config config);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange, special case: blocks without
 * explicit offsets
 * @param[in] xmap         exchange map
 * @param[in] src_block_sizes  source block lengths in unit of elements
 * @param[in] src_block_num  number of src blocks
 * @param[in] dst_block_sizes  destination block lengths in unit of elements
 * @param[in] dst_block_num  number of dst blocks
 * @param[in] datatype     MPI datatype of a single element in data space,
 *                         all elements have the same
 * @remarks calls \ref xt_redist_p2p_blocks_off_new with NULL offsets
 */
Xt_redist xt_redist_p2p_blocks_new(Xt_xmap xmap,
                                   const int *src_block_sizes,
                                   int src_block_num,
                                   const int *dst_block_sizes,
                                   int dst_block_num,
                                   MPI_Datatype datatype);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange, special case: blocks without
 * explicit offsets
 * @param[in] xmap         exchange map
 * @param[in] src_block_sizes  source block lengths in unit of elements
 * @param[in] src_block_num  number of src blocks
 * @param[in] dst_block_sizes  destination block lengths in unit of elements
 * @param[in] dst_block_num  number of dst blocks
 * @param[in] datatype     MPI datatype of a single element in data space,
 *                         all elements have the same
 * @param[in] config   configuration object for custom settings
 * @remarks calls \ref xt_redist_p2p_blocks_off_new with NULL offsets
 */
Xt_redist
xt_redist_p2p_blocks_custom_new(Xt_xmap xmap,
                                const int *src_block_sizes,
                                int src_block_num,
                                const int *dst_block_sizes,
                                int dst_block_num,
                                MPI_Datatype datatype,
                                Xt_config config);

#endif // XT_REDIST_P2P_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
