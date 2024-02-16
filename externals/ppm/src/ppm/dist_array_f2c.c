/**
 * @file dist_array_f2c.c
 * @brief Distributed data structure of multiple global arrays,
 * Fortran to C interface
 *
 * @copyright Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stddef.h>

#include <mpi.h>

#define FCALLSC_QUALIFIER PPM_DSO_INTERNAL
#include "cfortran.h"

#include "core/ppm_visibility.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "ppm/dist_array.h"

struct PPM_global_array_desc_f
{
  MPI_Fint a_rank;
  struct PPM_extent rect[PPM_dma_max_rank];
  MPI_Fint element_dt;
};

static inline struct PPM_dist_mult_array *
handle2dma(MPI_Aint dma_handle)
{
  return (struct PPM_dist_mult_array *)(intptr_t)dma_handle;
}

static void
PPM_dist_mult_array_new_f2c(MPI_Aint *new_dma_handle,
                            MPI_Fint num_sub_arrays,
                            const struct PPM_global_array_desc_f *sub_arrays,
                            const struct PPM_extent *local_chunk,
                            MPI_Fint comm_f,
                            MPI_Fint cache_size,
                            MPI_Fint sync_mode)
{
  MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
  if (num_sub_arrays < 0)
    PPM_abort(comm_c, "invalid number of component arrays specified",
              __FILE__, __LINE__);
  size_t num_sub_arrays_ = (size_t)num_sub_arrays;
  struct PPM_global_array_desc *sub_arrays_
    = xmalloc(num_sub_arrays_ * sizeof (*sub_arrays_));
  size_t max_a_rank = 0;
  for (size_t i = 0; i < num_sub_arrays_; ++i)
  {
    size_t a_rank = (size_t)sub_arrays[i].a_rank;
    if (a_rank > max_a_rank)
      max_a_rank = a_rank;
    sub_arrays_[i].a_rank = (unsigned)a_rank;
    for (size_t j = 0; j < a_rank; ++j)
      sub_arrays_[i].rect[j] = sub_arrays[i].rect[a_rank - j - 1];
    sub_arrays_[i].element_dt = MPI_Type_f2c(sub_arrays[i].element_dt);
  }
  const struct PPM_extent (*local_chunk_)[max_a_rank]
    = (const struct PPM_extent (*)[max_a_rank])local_chunk;
  struct PPM_extent (*local_chunk_copy)[max_a_rank]
    = xmalloc(sizeof (*local_chunk_copy) * num_sub_arrays_);
  for (size_t i = 0; i < num_sub_arrays_; ++i)
  {
    size_t a_rank = (size_t)sub_arrays[i].a_rank;
    for (size_t j = 0; j < a_rank; ++j)
      local_chunk_copy[i][j] = local_chunk_[i][a_rank-j-1];
  }

  if (cache_size < 0)
    PPM_abort(comm_c, "invalid number of cached ranks specified",
              __FILE__, __LINE__);
  size_t cache_size_ = (size_t)cache_size;
  if (sync_mode > PPM_dma_sync_mode_local_only
      || sync_mode < PPM_dma_sync_mode_passive_target)
    PPM_abort(comm_c, "invalid synchronization mode specified",
              __FILE__, __LINE__);
  *new_dma_handle = (MPI_Aint)(intptr_t)PPM_dist_mult_array_new(
    num_sub_arrays_, sub_arrays_, (struct PPM_extent *)local_chunk_copy,
    comm_c, cache_size_, (enum PPM_dma_sync_mode)sync_mode);
  free(local_chunk_copy);
  free(sub_arrays_);
}

FCALLSCSUB7(PPM_dist_mult_array_new_f2c, PPM_DIST_MULT_ARRAY_NEW_F2C,
            ppm_dist_mult_array_new_f2c,
            PVOID, INT, PVOID, PVOID, INT, INT, INT)

static void
PPM_dist_mult_array_copy_f2c(MPI_Aint *dma_handle, MPI_Aint *copy_handle)
{
  *copy_handle = (MPI_Aint)(intptr_t)
    PPM_dist_mult_array_copy(handle2dma(*dma_handle));
}

FCALLSCSUB2(PPM_dist_mult_array_copy_f2c, PPM_DIST_MULT_ARRAY_COPY_F2C,
            ppm_dist_mult_array_copy_f2c, PVOID, PVOID)


static void
PPM_dist_mult_array_delete_f2c(MPI_Aint *dma_handle)
{
  PPM_dist_mult_array_delete(handle2dma(*dma_handle));
}

FCALLSCSUB1(PPM_dist_mult_array_delete_f2c, PPM_DIST_MULT_ARRAY_DELETE_F2C,
            ppm_dist_mult_array_delete_f2c, PVOID)


static void
PPM_dist_mult_array_local_ptr_f2c(MPI_Aint *dma_handle, int sub_array_idx,
                                  void **sub_array_ptr)
{
  struct PPM_dist_mult_array *dm_array = handle2dma(*dma_handle);
  if (sub_array_idx < 1)
    PPM_abort(PPM_dist_mult_array_comm(dm_array), "invalid sub-array specified",
              __FILE__, __LINE__);
  *sub_array_ptr
    = PPM_dist_mult_array_local_ptr(dm_array, (size_t)sub_array_idx - 1);
}

FCALLSCSUB3(PPM_dist_mult_array_local_ptr_f2c,
            PPM_DIST_MULT_ARRAY_LOCAL_PTR_F2C,
            ppm_dist_mult_array_local_ptr_f2c, PVOID, INT, PVOID)

static void
PPM_dist_mult_array_expose_f2c(MPI_Aint *dma_handle)
{
  PPM_dist_mult_array_expose(handle2dma(*dma_handle));
}

FCALLSCSUB1(PPM_dist_mult_array_expose_f2c,
            PPM_DIST_MULT_ARRAY_EXPOSE_F2C,
            ppm_dist_mult_array_expose_f2c, PVOID)

static void
PPM_dist_mult_array_unexpose_f2c(MPI_Aint *dma_handle)
{
  PPM_dist_mult_array_unexpose(handle2dma(*dma_handle));
}

FCALLSCSUB1(PPM_dist_mult_array_unexpose_f2c,
            PPM_DIST_MULT_ARRAY_UNEXPOSE_F2C,
            ppm_dist_mult_array_unexpose_f2c, PVOID)

static void
PPM_dist_mult_array_rma_sync_f2c(MPI_Aint *dma_handle)
{
  PPM_dist_mult_array_rma_sync(handle2dma(*dma_handle));
}

FCALLSCSUB1(PPM_dist_mult_array_rma_sync_f2c,
            PPM_DIST_MULT_ARRAY_RMA_SYNC_F2C,
            ppm_dist_mult_array_rma_sync_f2c, PVOID)

static int
PPM_dist_mult_array_get_sync_mode_f2c(MPI_Aint *dma_handle)
{
  return PPM_dist_mult_array_get_sync_mode(handle2dma(*dma_handle));
}

FCALLSCFUN1(INT, PPM_dist_mult_array_get_sync_mode_f2c,
            PPM_DIST_MULT_ARRAY_GET_SYNC_MODE_F2C,
            ppm_dist_mult_array_get_sync_mode_f2c, PVOID)

static void
PPM_dist_mult_array_set_sync_mode_f2c(MPI_Aint *dma_handle, int sync_mode,
                                      int cache_size)
{
  struct PPM_dist_mult_array *dm_array = handle2dma(*dma_handle);
  if (cache_size < 0)
    PPM_abort(PPM_dist_mult_array_comm(dm_array),
              "invalid number of cached ranks specified",
              __FILE__, __LINE__);
  PPM_dist_mult_array_set_sync_mode(handle2dma(*dma_handle),
                                    (enum PPM_dma_sync_mode)sync_mode,
                                    (size_t)cache_size);
}

FCALLSCSUB3(PPM_dist_mult_array_set_sync_mode_f2c,
            PPM_DIST_MULT_ARRAY_SET_SYNC_MODE_F2C,
            ppm_dist_mult_array_set_sync_mode_f2c, PVOID, INT, INT)

/* the Fortran version of PPM_dist_mult_array_get_f2c is directly
 * called by the user while all functions above are only called
 * through a wrapper in module ppm_distributed_array */
#undef FCALLSC_QUALIFIER
#define FCALLSC_QUALIFIER

static void
PPM_dist_mult_array_get_f2c(MPI_Aint *dma_handle,
                            int sub_array_idx,
                            const int coord[],
                            void *v_out)
{
  struct PPM_dist_mult_array *dm_array = handle2dma(*dma_handle);
  if (sub_array_idx < 1)
    PPM_abort(PPM_dist_mult_array_comm(dm_array),
              "invalid sub-array specified",
              __FILE__, __LINE__);
  size_t a_rank
    = PPM_dist_mult_array_a_rank(dm_array, (size_t)sub_array_idx - 1);
  int coord_[PPM_dma_max_rank];
  for (size_t i = 0; i < a_rank; ++i)
    coord_[i] = coord[a_rank - i - 1];
  PPM_dist_mult_array_get(dm_array, (size_t)sub_array_idx - 1,
                          coord_, v_out);
}

FCALLSCSUB4(PPM_dist_mult_array_get_f2c,
            PPM_DIST_MULT_ARRAY_GET_F2C,
            ppm_dist_mult_array_get_f2c, PVOID, INT, INTV, PVOID)

#undef FCALLSC_QUALIFIER
#define FCALLSC_QUALIFIER PPM_DSO_INTERNAL

static void
PPM_dist_mult_array_rank_rect_f2c(MPI_Aint *dma_handle,
                                  int sub_array_idx,
                                  int rank,
                                  struct PPM_extent *rect)
{
  struct PPM_dist_mult_array *dm_array = handle2dma(*dma_handle);
  if (sub_array_idx < 1)
    PPM_abort(PPM_dist_mult_array_comm(dm_array),
              "invalid sub-array specified",
              __FILE__, __LINE__);
  PPM_dist_mult_array_rank_rect(dm_array, (size_t)sub_array_idx - 1, rank,
                                rect);
  size_t a_rank
    = PPM_dist_mult_array_a_rank(dm_array, (size_t)sub_array_idx - 1);
  for (size_t i = 0; i < a_rank/2; ++i)
  {
    struct PPM_extent temp = rect[i];
    rect[i] = rect[a_rank - i - 1];
    rect[a_rank - i - 1] = temp;
  }
}

FCALLSCSUB4(PPM_dist_mult_array_rank_rect_f2c,
            PPM_DIST_MULT_ARRAY_RANK_RECT_F2C,
            ppm_dist_mult_array_rank_rect_f2c, PVOID, INT, INT, PVOID)


static MPI_Fint
PPM_dist_mult_array_comm_f2c(MPI_Aint *dma_handle)
{
  struct PPM_dist_mult_array *dm_array = handle2dma(*dma_handle);
  return MPI_Comm_c2f(PPM_dist_mult_array_comm(dm_array));
}

FCALLSCFUN1(INT, PPM_dist_mult_array_comm_f2c, PPM_DIST_MULT_ARRAY_COMM_F2C,
            ppm_dist_mult_array_comm_f2c, PVOID)

static void
PPM_dist_mult_array_set_transfer_mode_f2c(MPI_Aint *dma_handle, int mode)
{
  struct PPM_dist_mult_array *dm_array = handle2dma(*dma_handle);
  PPM_dist_mult_array_set_transfer_mode(dm_array, mode);
}

FCALLSCSUB2(PPM_dist_mult_array_set_transfer_mode_f2c,
            PPM_DIST_MULT_ARRAY_SET_TRANSFER_MODE_F2C,
            ppm_dist_mult_array_set_transfer_mode_f2c, PVOID, INT)

static int
PPM_dist_mult_array_get_transfer_mode_f2c(MPI_Aint *dma_handle)
{
  struct PPM_dist_mult_array *dm_array = handle2dma(*dma_handle);
  return PPM_dist_mult_array_get_transfer_mode(dm_array);
}

FCALLSCFUN1(INT, PPM_dist_mult_array_get_transfer_mode_f2c,
            PPM_DIST_MULT_ARRAY_GET_TRANSFER_MODE_F2C,
            ppm_dist_mult_array_get_transfer_mode_f2c, PVOID)

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
