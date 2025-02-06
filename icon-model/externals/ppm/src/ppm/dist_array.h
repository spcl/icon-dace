/**
 * @file dist_array.h
 * @brief Distributed data structure of multiple global arrays
 *
 * @copyright Copyright  (C)  2014  Thomas Jahns <jahns@dkrz.de>
 *                                  Moritz Hanke <hanke@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
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
#ifndef DIST_ARRAY_H
#define DIST_ARRAY_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>

#include <mpi.h>

#include "core/ppm_extents.h"

enum {
  PPM_dma_max_rank = 7,
};

enum PPM_dma_sync_mode {
  /**
   *  in this mode calls to dist_mult_array_get will immediately
   * retrieve the requested value
   */
  PPM_dma_sync_mode_passive_target,
  /**
   * in this mode calls to dist_mult_array_get will result in
   * the passed variable to become defined only after the next call
   * to dist_mult_array_unexpose
   */
  PPM_dma_sync_mode_active_target,
  /**
   * in this mode, remote access is disabled at first but can later be
   * activated via @ref PPM_dist_mult_array_change_sync_mode
   */
  PPM_dma_sync_mode_local_only,
};

struct PPM_global_array_desc
{
  struct PPM_extent rect[PPM_dma_max_rank];
  MPI_Datatype element_dt;
  unsigned a_rank;
};

/**
 * @brief create distributed multi-array data structure
 *
 * The resulting data type represents a number of arrays distributed
 * over the ranks of the communicator passed to this function.
 *
 * @param num_sub_arrays number of distributed arrays
 *
 * @param sub_arrays contains array ranks, data types
 * and bounds of each distributed global array. The bounds are represented
 * by an @ref extent type that stores start and size.
 *
 * @param local_chunk
 * points to array of extents with dimensions
 * <tt>[num_sub_arrays][MAX(sub_arrays[].a_rank)]</tt> where
 * <tt>local_chunk[j][i]</tt> describes for dimension <tt>i</tt>
 * of sub_array <tt>j</tt> of the global arrays the local part available on
 * this MPI rank.\n Only contiguous local parts are possible.\n
 *
 * @param comm
 * MPI communicator for which this data structure is
 * collectively created
 *
 * @param cache_size number of ranks to cache remote local parts of or
 * automatically determined if 0
 *
 * @param sync_mode
 * switch synchronization modes, either PPM_sync_mode_passive_target
 * or sync_mode_active_target_mode.  For
 * sync_mode==PPM_dma_sync_mode_active_target_mode, execution of RMA is
 * deferred until the next synchronizing call
 * (dist_mult_array_unexpose or dist_mult_array_rma_sync).
 * (see @a PPM_dist_mult_array_get for example)
 *
 * @return initialized dist_mult_array structure
 * @remark This operation is collective for all MPI ranks in @a comm
 */
struct PPM_dist_mult_array *
PPM_dist_mult_array_new(size_t num_sub_arrays,
                        const struct PPM_global_array_desc *sub_arrays,
                        const struct PPM_extent *local_chunk,
                        MPI_Comm comm,
                        size_t cache_size,
                        enum PPM_dma_sync_mode sync_mode);

/**
 * @brief destroy dist_mult_array data type
 * @remark since an MPI window is freed, this routine is collective for
 * all processes which participated in the corresponding call to
 * @ref PPM_dist_mult_array_new
 * @param[in,out] dm_array distributed multi-array to delete
 */
void
PPM_dist_mult_array_delete(struct PPM_dist_mult_array *dm_array);

/**
 * @brief destroy dist_mult_array data type
 * @remark since an MPI window is freed, this routine is collective for
 * all processes which participated in the corresponding call to
 * @ref PPM_dist_mult_array_new
 * @param[in] dm_array distributed multi-array to delete
 */
struct PPM_dist_mult_array *
PPM_dist_mult_array_copy(struct PPM_dist_mult_array *dm_array);

/**
 * @brief make local data available to other ranks
 *
 * This call is collective for all ranks in the communicator used
 * to create @a dm_array
 * note: this enters an exposure epoch on the data shared via RMA
 * local data must not be changed while the array is in
 * exposed state
 */
void
PPM_dist_mult_array_expose(struct PPM_dist_mult_array *dm_array);

/**
 * @brief wait for all ranks to finish queries in current exposure epoch
 *
 * This call is collective for all ranks in the communicator used
 * to create @a dm_array
 *
 * note: local data can only be changed after the distributed array
 * was initially created or dist_mult_array_unexpose has been called
 */
void
PPM_dist_mult_array_unexpose(struct PPM_dist_mult_array *dm_array);

/**
 * @brief synchronize RMA updates only, ignore local updates
 *
 * This call is collective for all ranks in the communicator used
 * to create @a dm_array.
 * Also @a dm_array must be in exposed state.
 *
 * @param dm_array distributed multiple arrays to synchronize all RMA for
 */
void
PPM_dist_mult_array_rma_sync(struct PPM_dist_mult_array *dm_array);


/**
 * @brief query distributed array synchronization protocol
 *
 * @param dm_array distributed multiple arrays to query sync
 * mode for
 * @return sync mode set during initialization or the
 * last call to PPM_dist_mult_array_set_sync_mode for \a dm_array
 */
int
PPM_dist_mult_array_get_sync_mode(struct PPM_dist_mult_array *dm_array);

/**
 * @brief switch distributed array to different synchronization
 * protocol
 *
 * This call is collective for all ranks in the communicator used
 * to create @a dm_array
 *
 * note: this ends an exposure epoch on the data shared via RMA
 * note: in case the current sync mode is
 * PPM_dma_sync_mode_passive_target, the cache_size can still be
 * changed without actually switching the mode
 *
 * @param dm_array distributed multiple arrays to change access/sync
 * mode for
 * @param sync_mode new sync mode to set
 * @param cache_size in case @a sync_mode is
 * PPM_dma_sync_mode_passive_target, the number of ranks to cache data
 * for is set to this value (or determined automatically if 0)
  */
void
PPM_dist_mult_array_set_sync_mode(struct PPM_dist_mult_array *dm_array,
                                  enum PPM_dma_sync_mode sync_mode,
                                  size_t cache_size);

void *
PPM_dist_mult_array_local_ptr(struct PPM_dist_mult_array *dm_array,
                              size_t sub_array_idx);

/**
 * @brief Get value out of distributed multi-array, independent of
 * rank the data resides on.
 *
 * For example, to query the data assigned in the example for
 * dist_mult_array_local_ptr, use a query like this:
 * @code
 * struct PPMdist_mult_array *dma;
 * int i5;
 * PPM_dist_mult_array_expose(dma)
 * PPM_dist_mult_array_get(dma, 5, (/ i_l, j_l /), i5)
 * @endcode
 *
 * In case @a dma was initialized with sync_mode=sync_mode_active_target,
 * the value of i5 must not be accessed before calling a synchronizing
 * routine, either
 * @code
 * PPM_dist_mult_array_unexpose(dma);
 * @endcode
 * or
 * @code
 * PPM_dist_mult_array_rma_sync(dma);
 * @endcode
 */
void
PPM_dist_mult_array_get(struct PPM_dist_mult_array *dm_array,
                        size_t sub_array,
                        int32_t coord[],
                        void *v_out);

/**
 * Query array rank for sub-array of distributed multi-array.
 * @param[in] dm_array distributed multiple arrays to query
 * @param sub_array_idx index of sub-array to query
 * @return array rank in 1 to \a PPM_dma_max_rank
 */
unsigned
PPM_dist_mult_array_a_rank(struct PPM_dist_mult_array *dm_array,
                           size_t sub_array_idx);

/**
 * Query bounds of sub-array in distributed multi-array.
 * @param[in] dm_array distributed multiple arrays to query
 * @param sub_array_idx index of sub-array to query
 * @param rank MPI rank on communicator of \a dm_array to query chunk
 * for, own rank if less than 0
 * @param[out] rect this array of size
 * <tt>PPM_dist_mult_array_a_rank(dm_array, sub_array_idx)</tt> is set
 * to the local chunk size of \a rank
 */
void
PPM_dist_mult_array_rank_rect(struct PPM_dist_mult_array *dm_array,
                              size_t sub_array_idx,
                              int rank,
                              struct PPM_extent *rect);

MPI_Comm
PPM_dist_mult_array_comm(struct PPM_dist_mult_array *dm_array);

enum PPM_dma_transfer_mode {
  PPM_dma_transfer_mode_struct = 0,
  PPM_dma_transfer_mode_bytes = 1,
};

/**
 * set mode to transfer data for a distributed multi-array using
 * synchronization mode PPM_dma_sync_mode_passive_target.
 * @param[in,out] dm_array distributed multiple arrays to set mode for
 * @param[in] mode transfer data as struct or contiguous array of
 * bytes (the latter usually only works for MPI implementations on
 * homogeneous systems which don't feature different endianness or
 * alignments but is also often faster and this type of platform is
 * also the most common).
 */
void
PPM_dist_mult_array_set_transfer_mode(struct PPM_dist_mult_array *dm_array,
                                      int mode);

/**
 * get mode to transfer data for a distributed multi-array using
 * synchronization mode PPM_dma_sync_mode_passive_target.
 * @param[in] dm_array distributed multiple arrays to get mode for
 * @return currently active mode
 */
int
PPM_dist_mult_array_get_transfer_mode(struct PPM_dist_mult_array *dm_array);

#endif /* DIST_ARRAY_H */
/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
