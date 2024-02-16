/**
 * @file xt_idxlist.h
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

#ifndef XT_IDXLIST_H
#define XT_IDXLIST_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include "xt/xt_config.h"
#include "xt/xt_core.h"
#include "xt/xt_stripe.h"

/** \example test_idxlist_utils.c
 */
/** \example test_idxlist_utils_f.f90
 */
/** \example test_uid.c
 */

struct Xt_bounds {
   Xt_int start, size;
};

/**
 * \file xt_idxlist.h
 * \brief index list declaration
 *
 * contains declaration for the index list datatype
 */

/**
 * destructor for an index list
 *
 * @param[in] idxlist index list that is to be destroyed
 */
void xt_idxlist_delete(Xt_idxlist idxlist);

/**
 * computes the buffer size in byte that is required to pack the given
 * index list using MPI_Pack
 *
 * @param[in] idxlist index list for which the packed size needs to be computed
 * @param[in] comm    MPI communicator that will be used to send the buffer for
 *                    which the size is computed
 * @return            returns an upper bound on the size in bytes of the
 *                    packed index list
 *
 * @remark You need to provide the same buffer for all related calls of
 *         xt_idxlist_get_pack_size , \ref xt_idxlist_pack and
 *         \ref xt_idxlist_unpack
 *
 * @see xt_idxlist_pack
 * @see xt_idxlist_unpack
 */
size_t xt_idxlist_get_pack_size(Xt_idxlist idxlist, MPI_Comm comm);

/**
 * packs an index list into previously allocated buffer using MPI_Pack
 *
 * @param[in]     idxlist     index list, which needs be packed into the buffer
 * @param[in,out] buffer      previously allocated buffer into which the index
 *                            list will be packed
 * @param[in]     buffer_size size of buffer in bytes
 * @param[in,out] position    position in buffer at which the data of idxlist
 *                            will be copied
 *                            (position is automatically updated to the first
 *                            byte after the packed data)
 * @param[in]     comm        MPI communicator that will be used to send the
 *                            buffer
 *
 * @remark The provided buffer needs to have a size that is sufficient for the
 *         index list. You can determine the required size using
 *         \ref xt_idxlist_get_pack_size .
 * @remark You need to provide the same buffer for all related calls of
 *         \ref xt_idxlist_get_pack_size , xt_idxlist_pack and \ref xt_idxlist_unpack
 *
 * @see xt_idxlist_get_pack_size
 * @see xt_idxlist_unpack
 *
 */
void xt_idxlist_pack(Xt_idxlist idxlist, void* buffer, int buffer_size,
                     int* position, MPI_Comm comm);

/**
 * unpacks an index list from buffer
 *
 * @param[in]     buffer      buffer that contains the packed index list
 * @param[in]     buffer_size size of the buffer in bytes
 * @param[in,out] position    position in the buffer at which the unpacking
 *                            should start (will automatically be set to the
 *                            position after the unpacked data)
 * @param[in]     comm        MPI communicator that was used to receive the
 *                            buffer
 *
 * @remark You need to provide the same buffer for all related calls of
 *         \ref xt_idxlist_get_pack_size , \ref xt_idxlist_pack and
 *         \ref xt_idxlist_unpack
 *
 * @see xt_idxlist_pack
 */
Xt_idxlist xt_idxlist_unpack(void* buffer, int buffer_size,
                             int* position, MPI_Comm comm);

/**
 * computes the intersection between two index lists
 *
 * @param[in] idxlist_src index list for sender
 * @param[in] idxlist_dst index list for receiver
 * @return                return the intersection of idxlist_a and idxlist_b
 * @remark multiple occurrences of an element in idxlist_dst will result in
 *         multiple occurrences of the element in the intersection if
 *         idxlist_src contains it
 * @remark the elements in the resulting index list are sorted
 */
Xt_idxlist
xt_idxlist_get_intersection(Xt_idxlist idxlist_src, Xt_idxlist idxlist_dst);

/**
 * computes the intersection between two index lists with custom
 * configuration
 *
 * @param[in] idxlist_src index list for sender
 * @param[in] idxlist_dst index list for receiver
 * @param[in] config      custom configuration parameters
 * @return                return the intersection of idxlist_a and idxlist_b
 * @remark multiple occurrences of an element in idxlist_dst will result in
 *         multiple occurrences of the element in the intersection if
 *         idxlist_src contains it
 * @remark the elements in the resulting index list are sorted
 */
Xt_idxlist
xt_idxlist_get_intersection_custom(Xt_idxlist idxlist_src,
                                   Xt_idxlist idxlist_dst, Xt_config config);

/**
 * generates a copy of a given index list
 *
 * @param[in] idxlist index list that is to be copied
 * @return            copy of the given index list
 */
Xt_idxlist xt_idxlist_copy(Xt_idxlist idxlist);

/**
 * returns the number of indices stored in the given index list
 *
 * @param[in] idxlist index list for which the number of indices is required
 * @return            number of indices in the given index list
 */
int xt_idxlist_get_num_indices(Xt_idxlist idxlist);

/**
 * gets the indices stored in an index list
 *
 * @param[in]     idxlist index list for which the indices are to be returned
 * @param[in,out] indices array into which the indices are written
 */
void xt_idxlist_get_indices(Xt_idxlist idxlist, Xt_int *indices);

/**
 * gets a pointer of the constant indices stored in an index list
 *
 * @param[in]     idxlist index list for which the indices are to be returned
 * @return  pointer to const indices array, or NULL if
 *          xt_idxlist_get_num_indices(idxlist) == 0
 */
const Xt_int *xt_idxlist_get_indices_const(Xt_idxlist idxlist);

/**
 * returns the indices stored in an index list
 * the indices are returned in form of stripes (with stride 1)
 *
 * @param[in]  idxlist     index list for which the indices are to be returned
 * @param[out] stripes     array containing the stripes (user is responsible
 *                         for freeing the memory associated with the stripes)
 * @param[out] num_stripes number of stripes in array "stripes"
 */
void xt_idxlist_get_index_stripes(Xt_idxlist idxlist,
                                  struct Xt_stripe **stripes,
                                  int *num_stripes);

/**
 * gets the index stored at the specified position in an index list
 *
 * @param[in]  idxlist  index list
 * @param[in]  position position of the index that is to be retrieved
 * @param[out] index    index stored in the index list at the given position
 * @return              0 if index could be retrieved\n
 *                      1 if position is out of range
 */
int xt_idxlist_get_index_at_position(Xt_idxlist idxlist, int position,
                                     Xt_int * index);

/**
 * get indices stored at the specified positions in an index list
 *
 * @param[in]  idxlist    index list
 * @param[in]  positions  positions of the selected indices
 * @param[in]  num_pos    number of positions
 * @param[out] indices    selected indices
 * @param[in]  undef_idx  fallback value that is used if a position is invalid
 * @return                zero if no substitutions were necessary or\n
 *                        number of unresolved positions ( == number of undef
 *                        values assigned to the result)
 */
int xt_idxlist_get_indices_at_positions(Xt_idxlist idxlist,
                                        const int *positions,
                                        int num_pos, Xt_int *indices,
                                        Xt_int undef_idx);

/**
 * gets the position of the first occurrence of the given index
 * in the stored index list
 *
 * @param[in]  idxlist  index list in which the given index is to be searched
 * @param[in]  index    index that is to be searched
 * @param[out] position position of the first occurrence of index in idxlist,\n
 *                      or -1 if there is no match
 * @return              0 if index was found
 *                      1 if the index could not be found in idxlist
 */
int xt_idxlist_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                     int *position);

/**
 * gets the positions of the first occurrence of the given indices
 * in the stored index list
 *
 * @param[in]  idxlist     index list in which the given index is to be searched
 * @param[in]  indices     indices that are to be searched
 * @param[in]  num_indices number of indices in the array "indices"
 * @param[out] positions   positions of the first occurrence of the indices
 *                         in idxlist
 * @param[in]  single_match_only if true then do not consider a previous
 *                         matching position again
 *                         in the next index search - this is required on the
 *                         target side in order
 *                         to avoid multiple writes to the same memory position
 * @return                 0 if every index could be found
 *                         >0 number of indices which could not be found in
 *                         idxlist
 */
int xt_idxlist_get_positions_of_indices(Xt_idxlist idxlist,
                                        const Xt_int *indices,
                                        int num_indices, int *positions,
                                        int single_match_only);

/**
 * maps the positions of the first occurrence of the given index stripes
 * in the stored index list
 *
 * @param[in]  idxlist     index list in which the given index is to be searched
 * @param[in]  num_stripes number of stripes in the array \a stripes
 * @param[in]  stripes     stripes of indices that are to be searched
 * @param[out] num_ext     number of extents allocated in \a pos_exts
 * @param[out] pos_ext     ranges of positions for the first
 *                         occurrence of the indices from stripes in idxlist
 *                         the pointee is newly allocated and must be
 *                         freed by the caller
 * @param[in]  single_match_only if true then do not consider a previous
 *                         matching position again
 *                         in the next index search - this is required on the
 *                         target side in order
 *                         to avoid multiple writes to the same memory position
 * @return                 0 if every index could be found
 *                         >0 number of indices which could not be found in
 *                         idxlist
 */
int xt_idxlist_get_pos_exts_of_index_stripes(
  Xt_idxlist idxlist,
  int num_stripes,
  const struct Xt_stripe stripes[num_stripes],
  int *num_ext,
  struct Xt_pos_ext **pos_ext,
  int single_match_only);


/**
 * gets the position of the first occurrence of the given index
 * following the given offset in the stored index list idxlist
 *
 * @param[in]  idxlist  index list in which the given index is to be searched
 * @param[in]  index    index that is to be searched
 * @param[out] position position of the first occurrence of index in idxlist
 *                      with a position >= offset
 * @param[in]  offset   offset in the index list from which on the search
 *                      is to be conducted
 * @return              0 if index was found\n
 *                      1 if the index could not be found in idxlist
 */
int xt_idxlist_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                         int * position, int offset);

/**
 * gets the positions of the first occurrence of the given indices
 * following the given offsets in the stored index list
 *
 * @param[in]  idxlist     index list in which the given index is to be searched
 * @param[in]  indices     indices that are to be searched
 * @param[in]  num_indices number of indices in that array "indices"
 * @param[out] positions   positions of the first occurrence of the indices[i]
 *                         in idxlist with positions[i] >= offsets[i]
 * @param[in]  offsets     offsets in the index list from which on the search
 *                         is to be conducted
 * @return                 0 if index was found\n
 *                         1 if the index could not be found in idxlist
 *
 * \remarks a call to xt_idxlist_get_position_of_index_off with offset = NULL
 *          is the same as a call to xt_idxlist_get_positions_of_indices
 */
/* FIXME: this should probably return the number of failed lookups */
int xt_idxlist_get_positions_of_indices_off(Xt_idxlist idxlist,
                                            const Xt_int *indices,
                                            int num_indices, int * positions,
                                            int * offsets);

/**
 * gets the smallest index stored in the index list
 *
 * @param[in] idxlist index list
 * @return            smallest index in idxlist
 */
Xt_int xt_idxlist_get_min_index(Xt_idxlist idxlist);

/**
 * gets the largest index stored in the index list
 *
 * @param[in] idxlist index list
 * @return            biggest index in idxlist
 */
Xt_int xt_idxlist_get_max_index(Xt_idxlist idxlist);

/**
 * computes an n-dimensional bounding box around the given index list
 * @param[in]  idxlist            index list
 * @param[in]  ndim               number of dimension of the bounding box
 * @param[in]  global_size        global size of the n-dimensional index space
 *                                for which the bounding box is to be computed
 * @param[in]  global_start_index lowest index of the index space
 *                                (typically 0 or 1)
 * @param[out] bounds             bounds of the bounding box
 * @remarks the behaviour of the routine is undefined in case the given index
 *   space does not fit the indices
 * @remarks global_size has to be positive
 */
void xt_idxlist_get_bounding_box(Xt_idxlist idxlist, unsigned ndim,
                                 const Xt_int global_size[ndim],
                                 Xt_int global_start_index,
                                 struct Xt_bounds bounds[ndim]);

/**
 * return unique list id, where unique means no other index list will
 * return the same id within the same task, also all valid UIDs are non-zero
 * @param[in] idxlist index list of which to query unique id
 * @return unique id
 */
Xt_uid xt_idxlist_get_uid(Xt_idxlist idxlist);

#endif // XT_IDXLIST_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
