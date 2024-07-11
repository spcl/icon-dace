/**
 * @file xt_ddt.c
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdbool.h>
#include <string.h>
#include <mpi.h>

#ifdef _OPENACC
#define STR(s) #s
#define xt_Pragma(args) _Pragma(args)
#define XtPragmaACC(args) xt_Pragma(STR(acc args))
#else
#define XtPragmaACC(args)
#endif

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_ddt.h"
#include "xt_ddt_internal.h"

/* attribute handle for attaching xt_ddt objects to MPI datatypes */
static int xt_ddt_internal_keyval = MPI_KEYVAL_INVALID;

/* attribute attached MPI_COMM_SELF which cleans up
 * xt_mpi_datatype_ddt_internal_keyval when finalize is called */
static int xt_ddt_cleanup_internal_keyval = MPI_KEYVAL_INVALID;

struct xt_ddt_data {
  size_t kernel_idx;   // index into the list of available pack/unpack kernels
  size_t displ_count;  // number of elements in displs
  ssize_t *displs[XT_MEMTYPE_COUNT];
                       // displacements in byte of elements (int does not work)
                       // (if required this array has to be available
                       //  in GPU memory)
};

struct Xt_ddt_ {

  int ref_count;

  size_t pack_size;

  int displs_available[XT_MEMTYPE_COUNT]; // determines in which memory type
                                          // the displacements are available

  // if exchange data is on GPU, one kernel per data entry is launched
  size_t count; // number of elements in data
  struct xt_ddt_data data[];
};

struct xt_ddt_tree_elem {
  MPI_Aint extent; // not necessarly the true extent
  enum {
    DTYPE,      // element is a basic MPI datatype
    SINGLE_SUB, // element contains a single sub element
    MULTI_SUB,  // element contains a dedicated sub element
                // for each displacement
  } type;       // type of this element
  size_t displ_count;    // number of displacement
  MPI_Aint displs[];
  // following these is one of
  // MPI_Datatype dtype; when type == DTYPE
  // struct xt_ddt_tree_elem *sub_elems[]; when type == SINGLE_SUB or MULTI_SUB
};

static void xt_ddt_pack_8(
  size_t count, ssize_t *restrict displs, const uint8_t *restrict src,
  uint8_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_16(
  size_t count, ssize_t *restrict displs, const uint16_t *restrict src,
  uint16_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_32(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_32_2(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_96(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_64(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_128(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_160(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_pack_256(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t (*restrict dst)[4], enum xt_memtype memtype);

static void xt_ddt_unpack_8(
  size_t count, ssize_t *restrict displs, const uint8_t *restrict src,
  uint8_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_16(
  size_t count, ssize_t *restrict displs, const uint16_t *restrict src,
  uint16_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_32(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_32_2(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[2],
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_96(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[3],
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_64(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_128(
  size_t count, ssize_t *restrict displs, const uint64_t (*restrict src)[2],
  uint64_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_160(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[5],
  uint32_t *restrict dst, enum xt_memtype memtype);
static void xt_ddt_unpack_256(
  size_t count, ssize_t *restrict displs, const uint64_t (*restrict src)[4],
  uint64_t *restrict dst, enum xt_memtype memtype);

typedef void (*kernel_func)
  (size_t, ssize_t*, const void*, void*, enum xt_memtype);

static struct xt_un_pack_kernels {
  size_t base_pack_size;
  size_t element_size;
  kernel_func pack;
  kernel_func unpack;
} valid_kernels[] = {
  {.base_pack_size = 1,
   .element_size = 1,
   .pack = (kernel_func)xt_ddt_pack_8,
   .unpack = (kernel_func)xt_ddt_unpack_8},
  {.base_pack_size = 2,
   .element_size = 2,
   .pack = (kernel_func)xt_ddt_pack_16,
   .unpack = (kernel_func)xt_ddt_unpack_16},
  {.base_pack_size = 4,
   .element_size = 4,
   .pack = (kernel_func)xt_ddt_pack_32,
   .unpack = (kernel_func)xt_ddt_unpack_32},
  {.base_pack_size = 4,
   .element_size = 8,
   .pack = (kernel_func)xt_ddt_pack_32_2,
   .unpack = (kernel_func)xt_ddt_unpack_32_2},
  {.base_pack_size = 8,
   .element_size = 8,
   .pack = (kernel_func)xt_ddt_pack_64,
   .unpack = (kernel_func)xt_ddt_unpack_64},
  {.base_pack_size = 4,
   .element_size = 12,
   .pack = (kernel_func)xt_ddt_pack_96,
   .unpack = (kernel_func)xt_ddt_unpack_96},
  {.base_pack_size = 8,
   .element_size = 16,
   .pack = (kernel_func)xt_ddt_pack_128,
   .unpack = (kernel_func)xt_ddt_unpack_128},
  {.base_pack_size = 4,
   .element_size = 20,
   .pack = (kernel_func)xt_ddt_pack_160,
   .unpack = (kernel_func)xt_ddt_unpack_160},
  {.base_pack_size = 8,
   .element_size = 32,
   .pack = (kernel_func)xt_ddt_pack_256,
   .unpack = (kernel_func)xt_ddt_unpack_256},
};
enum {NUM_VALID_KERNELS = sizeof(valid_kernels)/sizeof(valid_kernels[0])};

#define MAX(a,b) ((a) >= (b) ? (a) : (b))

static inline void *after_displs(struct xt_ddt_tree_elem *elem)
{
  return elem->displs + elem->displ_count;
}

static struct xt_ddt_tree_elem *
  mpi_ddt_2_xt_ddt_tree(MPI_Datatype mpi_ddt, int free_mpi_ddt);

static struct xt_ddt_tree_elem *
xt_ddt_tree_elem_named_new(MPI_Aint extent, MPI_Datatype mpi_ddt)
{
  struct xt_ddt_tree_elem *elem
    = xmalloc(sizeof(*elem) + sizeof (MPI_Aint) + sizeof (MPI_Datatype));
  elem->extent = extent;
  elem->type = DTYPE;
  elem->displ_count = 1;
  elem->displs[0] = 0;
  *(MPI_Datatype *)after_displs(elem) = mpi_ddt;
  return elem;
}

static void free_unnamed_mpi_ddt(MPI_Datatype mpi_ddt) {

  // get parameters for decoding of MPI datatype
  int num_ints, num_adds, num_dtypes, combiner;
  xt_mpi_call(
    MPI_Type_get_envelope(
      mpi_ddt, &num_ints, &num_adds, &num_dtypes, &combiner), Xt_default_comm);

  if (combiner != MPI_COMBINER_NAMED)
    xt_mpi_call(MPI_Type_free(&mpi_ddt), Xt_default_comm);
}

// generates a ddt tree element with count contiguous entries of mpi_ddt
static struct xt_ddt_tree_elem *
xt_ddt_tree_elem_cont_new(int count, MPI_Datatype mpi_ddt)
{

  // if the sub element count is zero
  if (count == 0) {
    free_unnamed_mpi_ddt(mpi_ddt);
    return NULL;
  }

  struct xt_ddt_tree_elem * sub_elem = mpi_ddt_2_xt_ddt_tree(mpi_ddt, 1);

  // if the sub element is empty
  if (sub_elem == NULL) return NULL;

  // if there is only a single sub element
  if (count == 1) return sub_elem;

  struct xt_ddt_tree_elem * elem;
  MPI_Aint sub_elem_extent = sub_elem->extent;

  // if the sub element is a basic MPI_Datatype
  if ((sub_elem->type == DTYPE) &&
      (sub_elem->displ_count == 1) &&
      (sub_elem->displs[0] == 0)) {

    elem = xrealloc(sub_elem,
                    sizeof (*elem) + (size_t)count * sizeof (MPI_Aint)
                    + sizeof (MPI_Datatype));
    MPI_Datatype dt = *(MPI_Datatype *)after_displs(elem);
    elem->displ_count = (size_t)count;
    *(MPI_Datatype *)after_displs(elem) = dt;
  } else {

    elem = xmalloc(sizeof(*elem) + (size_t)count * sizeof (MPI_Aint)
                   + sizeof (struct xt_ddt_tree_elem *));
    elem->type = SINGLE_SUB;
    elem->displ_count = (size_t)count;
    *(struct xt_ddt_tree_elem **)after_displs(elem) = sub_elem;
  }

  elem->extent = (MPI_Aint)count * sub_elem_extent;
  MPI_Aint *displs = elem->displs;
  for (int i = 0; i < count; ++i) displs[i] = (MPI_Aint)i * sub_elem_extent;

  return elem;
}

// generates a ddt tree element based on a hvector style memory layout
static struct xt_ddt_tree_elem *
xt_ddt_tree_elem_hvector_new(
  MPI_Aint extent, int count, int blocklength, MPI_Aint stride,
  MPI_Datatype mpi_ddt) {

  // if the sub element count is zero
  if (count == 0) {
    free_unnamed_mpi_ddt(mpi_ddt);
    return NULL;
  }

  struct xt_ddt_tree_elem *sub_elem =
    xt_ddt_tree_elem_cont_new(blocklength, mpi_ddt);

  // if the sub element is empty
  if (sub_elem == NULL) return NULL;

  // if there is only a single sub element
  if (count == 1) return sub_elem;

  struct xt_ddt_tree_elem *elem
    = xmalloc(sizeof(*elem) + (size_t)count * sizeof (MPI_Aint)
              + sizeof (struct xt_ddt_tree_elem *));

  elem->displ_count = (size_t)count;
  elem->extent = extent;
  elem->type = SINGLE_SUB;
  *(struct xt_ddt_tree_elem **)after_displs(elem) = sub_elem;

  MPI_Aint *displs = elem->displs;
  for (int i = 0; i < count; ++i) displs[i] = (MPI_Aint)i * stride;

  return elem;
}

// generates a ddt tree element based on a vector style memory layout
static struct xt_ddt_tree_elem *xt_ddt_tree_elem_vector_new(
  MPI_Aint extent, int count, int blocklength, int stride,
  MPI_Datatype mpi_ddt) {

  MPI_Aint lb, base_extent;
  xt_mpi_call(
    MPI_Type_get_extent(mpi_ddt, &lb, &base_extent), Xt_default_comm);

  MPI_Aint hstride = base_extent * (MPI_Aint)stride;

  return
    xt_ddt_tree_elem_hvector_new(
      extent, count, blocklength, hstride, mpi_ddt);
}

// generates a ddt tree element based on a hindexed style memory layout
static struct xt_ddt_tree_elem *xt_ddt_tree_elem_hindexed_new(
  MPI_Aint extent, int count, int *blocklengths,
  MPI_Aint *displacements, MPI_Datatype mpi_ddt) {

  size_t displ_count = 0;
  for (int i = 0; i < count; ++i)
    displ_count += (size_t)blocklengths[i];

  // if there are no entries
  if (displ_count == 0) {
    free_unnamed_mpi_ddt(mpi_ddt);
    return NULL;
  }

  MPI_Aint lb, base_extent;
  MPI_Type_get_extent(mpi_ddt, &lb, &base_extent);
  struct xt_ddt_tree_elem *sub_elem = mpi_ddt_2_xt_ddt_tree(mpi_ddt, 1);

  // if the sub element is empty
  if (sub_elem == NULL) return NULL;

  struct xt_ddt_tree_elem *elem
    = xmalloc(sizeof(*elem) + displ_count * sizeof (MPI_Aint)
              + sizeof (struct xt_ddt_tree_elem *));

  elem->displ_count = displ_count;
  elem->extent = extent;
  elem->type = SINGLE_SUB;
  MPI_Aint *displs = elem->displs;
  *(struct xt_ddt_tree_elem **)(displs + displ_count) = sub_elem;

  size_t idx = 0;
  for (int i = 0; i < count; ++i) {
    MPI_Aint block_displacement = displacements[i];
    for (int j = 0; j < blocklengths[i]; ++j, ++idx) {
      displs[idx] = block_displacement + (MPI_Aint)j * base_extent;
    }
  }

  return elem;
}

// generates a ddt tree element based on an indexed style memory layout
static struct xt_ddt_tree_elem *xt_ddt_tree_elem_indexed_new(
  MPI_Aint extent, int count, int *blocklengths,
  int *displacements, MPI_Datatype mpi_ddt) {

  MPI_Aint lb, base_extent;
  MPI_Type_get_extent(mpi_ddt, &lb, &base_extent);

  MPI_Aint *hdisplacements = xmalloc((size_t)count * sizeof(*hdisplacements));
  for (int i = 0; i < count; ++i)
    hdisplacements[i] = (MPI_Aint)displacements[i] * base_extent;

  struct xt_ddt_tree_elem *elem =
    xt_ddt_tree_elem_hindexed_new(
      extent, count, blocklengths, hdisplacements, mpi_ddt);

  free(hdisplacements);

  return elem;
}

// generates a ddt tree element based on an hindexed block style memory layout
static struct xt_ddt_tree_elem *xt_ddt_tree_elem_hindexed_block_new(
  MPI_Aint extent, int count, int blocklength, MPI_Aint *displacements,
  MPI_Datatype mpi_ddt) {

  // if there are no blocks or if the blocks are empty
  if ((count == 0) || (blocklength == 0)) {
    free_unnamed_mpi_ddt(mpi_ddt);
    return NULL;
  }

  struct xt_ddt_tree_elem *sub_elem =
    xt_ddt_tree_elem_cont_new(blocklength, mpi_ddt);

  // if the sub element is empty
  if (sub_elem == NULL) return NULL;

  struct xt_ddt_tree_elem *elem
    = xmalloc(sizeof(*elem) + (size_t)count * sizeof (MPI_Aint)
              + sizeof (struct xt_ddt_tree_elem *));

  elem->displ_count = (size_t)count;
  elem->extent = extent;
  elem->type = SINGLE_SUB;
  MPI_Aint *displs = elem->displs;
  *(struct xt_ddt_tree_elem **)(displs + count) = sub_elem;

  for (int i = 0; i < count; ++i) displs[i] = displacements[i];

  return elem;
}

// generates a ddt tree element based on an indexed block style memory layout
static struct xt_ddt_tree_elem *xt_ddt_tree_elem_indexed_block_new(
  MPI_Aint extent, int count, int blocklength, int *displacements,
  MPI_Datatype mpi_ddt) {

  MPI_Aint lb, base_extent;
  xt_mpi_call(
    MPI_Type_get_extent(mpi_ddt, &lb, &base_extent), Xt_default_comm);

  MPI_Aint *hdisplacements = xmalloc((size_t)count * sizeof(*hdisplacements));
  for (int i = 0; i < count; ++i)
    hdisplacements[i] = (MPI_Aint)displacements[i] * base_extent;

  struct xt_ddt_tree_elem *elem =
    xt_ddt_tree_elem_hindexed_block_new(
      extent, count, blocklength, hdisplacements, mpi_ddt);

  free(hdisplacements);

  return elem;
}

// generates a ddt tree element based on a struct style memory layout
static struct xt_ddt_tree_elem *xt_ddt_tree_elem_struct_new(
  MPI_Aint extent, int count, int *blocklengths, MPI_Aint *displacements,
  MPI_Datatype *mpi_ddts) {

  size_t count_ = count > 0 ? (size_t)count : (size_t)0;
  // remove empty blocks
  size_t new_count = 0;
  for (size_t i = 0; i < count_; ++i) {

    MPI_Aint true_lb, true_extent;
    xt_mpi_call(
      MPI_Type_get_true_extent(mpi_ddts[i], &true_lb, &true_extent),
      Xt_default_comm);

    if ((blocklengths[i] != 0) && (true_extent != 0)) {

      if (new_count != count_) {

        blocklengths[new_count] = blocklengths[i];
        displacements[new_count] = displacements[i];
        mpi_ddts[new_count] = mpi_ddts[i];
      }

      ++new_count;
    } else {
      free_unnamed_mpi_ddt(mpi_ddts[i]);
    }
  }
  if (new_count != count_) count_ = new_count;

  // if there are no blocks
  if (count_ == 0) return NULL;

  // if there is only a single block
  if (count_ == 1)
    return
      xt_ddt_tree_elem_hindexed_block_new(
        extent, 1, blocklengths[0], displacements, *mpi_ddts);

  struct xt_ddt_tree_elem *elem
    = xmalloc(sizeof(*elem) + count_ * sizeof (MPI_Aint)
              + count_ * sizeof(struct xt_ddt_tree_elem *));
  elem->displ_count = count_;
  elem->extent = extent;
  elem->type = MULTI_SUB;
  MPI_Aint *displs = elem->displs;
  struct xt_ddt_tree_elem **sub_elems
    = (struct xt_ddt_tree_elem **)(displs + count_);

  for (size_t i = 0; i < count_; ++i) {
    sub_elems[i] = xt_ddt_tree_elem_cont_new(blocklengths[i], mpi_ddts[i]);
    displs[i] = displacements[i];
  }

  return elem;
}

// generate the displacements for a subarray
// (this routine is recursive)
static void xt_ddt_tree_elem_subarray_get_displs(
  MPI_Aint *displs, MPI_Aint dim_displ,
  size_t ndim, int *sizes, int *sub_sizes, int *starts, int order,
  MPI_Aint base_extent) {

  // if the base element is reached
  if (ndim == 0) {
    *displs = dim_displ;
    return;
  }

  // number of entries in the current dimension and index of the first entry
  int dim_sub_size, dim_start;

  if (order == MPI_ORDER_FORTRAN) {
    dim_sub_size = sub_sizes[ndim - 1];
    dim_start = starts[ndim - 1];
  } else {
    dim_sub_size = sub_sizes[0];
    dim_start = starts[0];
    sizes += 1;
    sub_sizes += 1;
    starts += 1;
  }

  size_t sub_array_size = 1;
  MPI_Aint sub_array_extent = base_extent;
  for (size_t i = 0; i < ndim-1; ++i) {
    sub_array_size *= (size_t)(sub_sizes[i]);
    sub_array_extent *= sizes[i];
  }

  // set current offset to the first entry in the current dimension
  dim_displ += (MPI_Aint)dim_start * sub_array_extent;

  // for the number of entries in the current dimension
  for (int i = 0; i < dim_sub_size; ++i) {

    // recursive call for sub dimensions
    xt_ddt_tree_elem_subarray_get_displs(
      displs, dim_displ, ndim - 1,
      sizes, sub_sizes, starts, order, base_extent);

    // set to next entry in current dimension
    displs += sub_array_size;
    dim_displ += sub_array_extent;
  }
}

static struct xt_ddt_tree_elem *xt_ddt_tree_elem_subarray_new(
  MPI_Aint extent, int ndim, int *sizes, int *sub_sizes,
  int *starts, int order, MPI_Datatype mpi_ddt) {

  struct xt_ddt_tree_elem *sub_elem = mpi_ddt_2_xt_ddt_tree(mpi_ddt, 1);

  // if the sub element is empty
  if (sub_elem == NULL) return NULL;

  MPI_Aint base_extent = sub_elem->extent;

  size_t displ_count = 1;
  for (int i = 0; i < ndim; ++i)
    displ_count *= (size_t)sub_sizes[i];

  struct xt_ddt_tree_elem *elem
    = xmalloc(sizeof(*elem) + displ_count * sizeof (MPI_Aint)
              + sizeof (struct xt_ddt_tree_elem *));
  elem->displ_count = displ_count;
  elem->extent = extent;
  elem->type = SINGLE_SUB;
  MPI_Aint *displs = elem->displs;
  *(struct xt_ddt_tree_elem **)(displs + displ_count) = sub_elem;

  xt_ddt_tree_elem_subarray_get_displs(
    displs, 0, (size_t)ndim, sizes, sub_sizes, starts, order, base_extent);

  return elem;
}

static struct xt_ddt_tree_elem *xt_ddt_tree_elem_resized_new(
  MPI_Aint extent, MPI_Datatype mpi_ddt) {

  struct xt_ddt_tree_elem *elem = mpi_ddt_2_xt_ddt_tree(mpi_ddt, 1);
  elem->extent = extent;

  return elem;
}

// generates a ddt tree element from a MPI datatype
static struct xt_ddt_tree_elem *
  mpi_ddt_2_xt_ddt_tree(MPI_Datatype mpi_ddt, int free_mpi_ddt) {

  // get parameters for decoding of MPI datatype
  int num_ints, num_adds, num_dtypes, combiner;
  xt_mpi_call(
    MPI_Type_get_envelope(
      mpi_ddt, &num_ints, &num_adds, &num_dtypes, &combiner), Xt_default_comm);
  MPI_Datatype oldtype[1];

  // get the extent of the MPI datatype
  MPI_Aint lb, extent;
  xt_mpi_call(MPI_Type_get_extent(mpi_ddt, &lb, &extent), Xt_default_comm);

  struct xt_ddt_tree_elem *elem;

  // decode MPI datatype
  switch (combiner) {

    default: {
      Xt_abort(
        Xt_default_comm,
        "ERROR(mpi_ddt_2_xt_ddt_tree): unsupported combiner type",
        __FILE__, __LINE__);
      return NULL;
    }

    case MPI_COMBINER_NAMED:

      elem = xt_ddt_tree_elem_named_new(extent, mpi_ddt);
      break;

    case MPI_COMBINER_DUP: {

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          NULL, NULL, oldtype), Xt_default_comm);

      elem = mpi_ddt_2_xt_ddt_tree(oldtype[0], 1);
      break;
    }

    case MPI_COMBINER_CONTIGUOUS: {

      int array_of_ints[1];

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, NULL, oldtype), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_cont_new(
          array_of_ints[0], oldtype[0]);
      break;
    }

    case MPI_COMBINER_VECTOR: {

      int array_of_ints[3];

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, NULL, oldtype), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_vector_new(
          extent, array_of_ints[0], array_of_ints[1], array_of_ints[2],
          oldtype[0]);
      break;
    }

    case MPI_COMBINER_HVECTOR: {

      int array_of_ints[2];
      MPI_Aint array_of_adds[1];

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, oldtype), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_hvector_new(
          extent, array_of_ints[0], array_of_ints[1], array_of_adds[0],
          oldtype[0]);
      break;
    }
    case MPI_COMBINER_INDEXED:{

      int *array_of_ints = xmalloc((size_t)num_ints * sizeof(*array_of_ints));

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, NULL, oldtype), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_indexed_new(
          extent, array_of_ints[0],
          array_of_ints + 1, array_of_ints + 1 + 1 * array_of_ints[0],
          oldtype[0]);

      free(array_of_ints);

      break;
    }
    case MPI_COMBINER_HINDEXED:{

      MPI_Aint *array_of_adds =
        xmalloc((size_t)num_adds * sizeof(*array_of_adds)
                + (size_t)num_ints * sizeof(int));
      int *array_of_ints = (void *)(array_of_adds + num_adds);

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, oldtype), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_hindexed_new(
          extent, array_of_ints[0], array_of_ints + 1,
          array_of_adds, oldtype[0]);

      free(array_of_adds);

      break;
    }
    case MPI_COMBINER_INDEXED_BLOCK:{

      int *array_of_ints = xmalloc((size_t)num_ints * sizeof(*array_of_ints));

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, NULL, oldtype), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_indexed_block_new(
          extent, array_of_ints[0], array_of_ints[1], array_of_ints + 2,
          oldtype[0]);

      free(array_of_ints);

      break;
    }
#if MPI_VERSION >= 3
    case MPI_COMBINER_HINDEXED_BLOCK:{

      int array_of_ints[2];
      MPI_Aint *array_of_adds =
        xmalloc((size_t)num_adds * sizeof(*array_of_adds));

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, oldtype), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_hindexed_block_new(
          extent, array_of_ints[0], array_of_ints[1], array_of_adds,
          oldtype[0]);

      free(array_of_adds);

      break;
    }
#endif
    case MPI_COMBINER_STRUCT:{

      MPI_Aint *array_of_adds =
        xmalloc((size_t)num_adds * sizeof(*array_of_adds)
                + (size_t)num_dtypes * sizeof(MPI_Datatype)
                + (size_t)num_ints * sizeof(int));
      MPI_Datatype *array_of_dtypes = (void *)(array_of_adds + num_adds);
      int *array_of_ints = (void *)(array_of_dtypes + num_dtypes);

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_struct_new(
          extent, array_of_ints[0], array_of_ints + 1,
          array_of_adds, array_of_dtypes);

      free(array_of_adds);

      break;
    }
    case MPI_COMBINER_SUBARRAY:{

      int *array_of_ints =
        xmalloc((size_t)num_ints * sizeof(*array_of_ints));

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, NULL, oldtype), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_subarray_new(
          extent, array_of_ints[0],
          array_of_ints + 1 + 0 * array_of_ints[0],
          array_of_ints + 1 + 1 * array_of_ints[0],
          array_of_ints + 1 + 2 * array_of_ints[0],
          array_of_ints[1 + 3 * array_of_ints[0]],
          oldtype[0]);

      free(array_of_ints);

      break;
    }
    case MPI_COMBINER_RESIZED:{

      MPI_Aint array_of_adds[2];

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          NULL, array_of_adds, oldtype), Xt_default_comm);

      elem = xt_ddt_tree_elem_resized_new(extent, oldtype[0]);
      break;
    }
  };

  if ((combiner != MPI_COMBINER_NAMED) && (free_mpi_ddt))
    xt_mpi_call(MPI_Type_free(&mpi_ddt), Xt_default_comm);

  return elem;
}

static void xt_ddt_tree_elem_delete(struct xt_ddt_tree_elem *elem) {

  if (elem == NULL) return;

  size_t num_sub_elem = 0;
  switch (elem->type) {
  case DTYPE:
    num_sub_elem = 0;
    break;
  case SINGLE_SUB:
    num_sub_elem = 1;
    break;
  case MULTI_SUB:
    num_sub_elem = elem->displ_count;
    break;
  }
  struct xt_ddt_tree_elem **sub_elems
    = (struct xt_ddt_tree_elem **)(elem->displs + elem->displ_count);
  for (size_t i = 0; i < num_sub_elem; ++i)
    xt_ddt_tree_elem_delete(sub_elems[i]);
  free(elem);
}

static void xt_ddt_tree_delete(struct xt_ddt_tree_elem *tree) {

  if (tree != NULL) xt_ddt_tree_elem_delete(tree);
}

static inline size_t mpi_ddt_to_data_idx(MPI_Datatype mpi_ddt) {

  MPI_Aint true_lb, true_extent;
  xt_mpi_call(
    MPI_Type_get_true_extent(mpi_ddt, &true_lb, &true_extent),
    Xt_default_comm);

  size_t base_pack_size;
  if (mpi_ddt == MPI_FLOAT_INT)
    base_pack_size = MAX(sizeof(float), sizeof(int));
  else
    base_pack_size = (size_t)true_extent;
  for (size_t i = 0; i < NUM_VALID_KERNELS; ++i)
    if (valid_kernels[i].element_size == (size_t)true_extent
        && valid_kernels[i].base_pack_size == base_pack_size)
      return i;
  for (size_t i = 0; i < NUM_VALID_KERNELS; ++i)
    if (valid_kernels[i].element_size == (size_t)true_extent)
      return i;

  char err_msg[128];
  snprintf(
    err_msg, sizeof(err_msg),
    "ERROR(mpi_ddt_to_data_idx): unsupported datatype size (%d Byte)",
    (int)true_extent);
  Xt_abort(Xt_default_comm, err_msg, __FILE__, __LINE__);
  return SIZE_MAX;
}

static void xt_ddt_tree_elem_get_data_sizes(
  struct xt_ddt_tree_elem *elem, MPI_Datatype *prev_dtype,
  size_t *prev_data_idx, struct xt_ddt_data *data) {

  if (elem == NULL) return;

  void *dtype_data = after_displs(elem);
  size_t displ_count = elem->displ_count;
  switch (elem->type) {

    case(DTYPE): {

      MPI_Datatype dt = *(MPI_Datatype *)dtype_data;
      if (dt != *prev_dtype) {

        *prev_dtype = dt;
        *prev_data_idx = mpi_ddt_to_data_idx(dt);
      }

      data[*prev_data_idx].displ_count += displ_count;
      break;
    }
    case (SINGLE_SUB): {

      struct xt_ddt_tree_elem *sub_elem
        = *(struct xt_ddt_tree_elem **)dtype_data;

      for (size_t i = 0; i < displ_count; ++i)
        xt_ddt_tree_elem_get_data_sizes(
          sub_elem, prev_dtype, prev_data_idx, data);
      break;
    }
    case (MULTI_SUB): {

      struct xt_ddt_tree_elem **sub_elems
        = (struct xt_ddt_tree_elem **)dtype_data;

      for (size_t i = 0; i < displ_count; ++i)
        xt_ddt_tree_elem_get_data_sizes(
          sub_elems[i], prev_dtype, prev_data_idx, data);
      break;
    }
  };
}

static void xt_ddt_tree_get_data_sizes(
  struct xt_ddt_tree_elem *tree, struct xt_ddt_data *data) {

  MPI_Datatype prev_dtype = MPI_DATATYPE_NULL;
  size_t prev_data_idx = SIZE_MAX;

  if (tree != NULL)
    xt_ddt_tree_elem_get_data_sizes(
      tree, &prev_dtype, &prev_data_idx, data);
}

static void xt_ddt_tree_elem_to_data(
  struct xt_ddt_tree_elem *elem, struct xt_ddt_data *data, MPI_Aint displ,
  MPI_Datatype *prev_dtype, size_t *prev_data_idx) {

  if (elem == NULL) return;

  size_t elem_displ_count = elem->displ_count;
  MPI_Aint *elem_displs = elem->displs;
  void *dtype_data = after_displs(elem);

  switch (elem->type) {

    case(DTYPE): {

      MPI_Datatype dt = *(MPI_Datatype *)dtype_data;
      if (dt != *prev_dtype) {

        *prev_dtype = dt;
        *prev_data_idx = mpi_ddt_to_data_idx(dt);
      }

      size_t displ_count = data[*prev_data_idx].displ_count;
      ssize_t *displs = data[*prev_data_idx].displs[XT_MEMTYPE_HOST];
      for (size_t i = 0; i < elem_displ_count; ++i, ++displ_count)
        displs[displ_count] = (ssize_t)displ + (ssize_t)elem_displs[i];
      data[*prev_data_idx].displ_count = displ_count;

      break;
    }
    case (SINGLE_SUB): {

      struct xt_ddt_tree_elem *sub_elem
        = *(struct xt_ddt_tree_elem **)dtype_data;

      for (size_t i = 0; i < elem_displ_count; ++i)
        xt_ddt_tree_elem_to_data(
          sub_elem, data, displ + elem_displs[i],
          prev_dtype, prev_data_idx);
      break;
    }
    case (MULTI_SUB): {

      struct xt_ddt_tree_elem **sub_elems
        = (struct xt_ddt_tree_elem **)dtype_data;

      for (size_t i = 0; i < elem_displ_count; ++i)
        xt_ddt_tree_elem_to_data(
          sub_elems[i], data, displ + elem_displs[i],
          prev_dtype, prev_data_idx);
      break;
    }
  };
}

static void xt_ddt_tree_to_data(
  struct xt_ddt_tree_elem *tree, struct xt_ddt_data *data) {

  MPI_Aint displ = 0;
  MPI_Datatype prev_dtype = MPI_DATATYPE_NULL;
  size_t prev_data_idx = SIZE_MAX;

  if (tree != NULL)
    xt_ddt_tree_elem_to_data(
      tree, data, displ, &prev_dtype, &prev_data_idx);
}

static int compare_kernels(const void * a, const void * b) {
  struct xt_un_pack_kernels const * k_a = (struct xt_un_pack_kernels const *)a;
  struct xt_un_pack_kernels const * k_b = (struct xt_un_pack_kernels const *)b;

  int kcmp = (k_a->base_pack_size < k_b->base_pack_size)
    - (k_a->base_pack_size > k_b->base_pack_size);
  if (!kcmp)
    kcmp = (k_a->element_size < k_b->element_size)
      - (k_a->element_size > k_b->element_size);
  return kcmp;
}

static Xt_ddt xt_ddt_new(MPI_Datatype mpi_ddt) {

  static bool sort_kernels = true;
  if (sort_kernels) {
    // sort the supported kernels by their base_pack_size (this
    // avoids alignment issues when packing)
    qsort(valid_kernels, NUM_VALID_KERNELS, sizeof(*valid_kernels),
          compare_kernels);
    sort_kernels = false;
  }


  // parse the MPI datatype
  struct xt_ddt_tree_elem *tree = mpi_ddt_2_xt_ddt_tree(mpi_ddt, 0);

  struct xt_ddt_data data[NUM_VALID_KERNELS];
  for (size_t i = 0; i < NUM_VALID_KERNELS; ++i) {
    data[i].kernel_idx = i;
    data[i].displ_count = 0;
    for (int j = 0; j < XT_MEMTYPE_COUNT; ++j) data[i].displs[j] = NULL;
  }

  // determine the number of displacements for each elemental data size
  xt_ddt_tree_get_data_sizes(tree, data);

  // compute the total number of elemental data
  size_t total_displs_size = 0;
  size_t data_count = 0;
  size_t pack_size = 0;
  for (size_t i = 0; i < NUM_VALID_KERNELS; ++i) {
    total_displs_size += data[i].displ_count;
    if (data[i].displ_count > 0) {
      pack_size +=
        valid_kernels[data[i].kernel_idx].element_size * data[i].displ_count;
      ++data_count;
    }
  }

  Xt_ddt ddt = xmalloc(1 * sizeof(*ddt) + data_count * sizeof(ddt->data[0]));
  ddt->ref_count = 1;
  ddt->count = data_count;
  ddt->pack_size = pack_size;
  for (int i = 0; i < XT_MEMTYPE_COUNT; ++i) ddt->displs_available[i] = 0;
  ddt->displs_available[XT_MEMTYPE_HOST] = 1;

  if (total_displs_size > 0) {

    // allocate memory for all displacements
    ssize_t *displs =
      xt_gpu_malloc(total_displs_size * sizeof(*displs), XT_MEMTYPE_HOST);

    for (size_t i = 0, offset = 0; i < NUM_VALID_KERNELS; ++i) {
      data[i].displs[XT_MEMTYPE_HOST] = displs + offset;
      offset += data[i].displ_count;
      data[i].displ_count = 0;
    }

    // determine all displacements
    xt_ddt_tree_to_data(tree, data);

    for (int i = 0, j = 0; i < NUM_VALID_KERNELS; ++i)
      if (data[i].displ_count > 0)
        ddt->data[j++] = data[i];
  }

  // free internal representation of MPI datatype
  xt_ddt_tree_delete(tree);

  return ddt;
}

void xt_ddt_inc_ref_count(Xt_ddt ddt) {

  ddt->ref_count++;
}

void xt_ddt_delete(Xt_ddt ddt) {

  if (ddt == NULL) return;

  ddt->ref_count--;

  if (ddt->ref_count) return;

  // only the displacements of the first entry needs to be freed, all other
  // are part of the same allocation
  if (ddt->count > 0) {
    for (int i = 0; i < XT_MEMTYPE_COUNT; ++i) {
      if (ddt->displs_available[i]) {
        xt_gpu_free(ddt->data[0].displs[i], (enum xt_memtype)i);
      }
    }
  }
  free(ddt);
}

static int xt_ddt_internal_keyval_copy(
  MPI_Datatype XT_UNUSED(dtype), int XT_UNUSED(dtype_keyval),
  void *XT_UNUSED(extra_state), void *attribute_val_in,
  void *attribute_val_out, int *flag) {

  xt_ddt_inc_ref_count((Xt_ddt)attribute_val_in);
  *(Xt_ddt*)attribute_val_out = (Xt_ddt)attribute_val_in;
  *flag = 1;
  return MPI_SUCCESS;
}

static int xt_ddt_internal_keyval_delete(
  MPI_Datatype XT_UNUSED(dtype), int XT_UNUSED(dtype_keyval),
  void *attribute_val, void *XT_UNUSED(extra_state)) {

  xt_ddt_delete((Xt_ddt)attribute_val);
  return MPI_SUCCESS;
}

static int xt_ddt_cleanup_internal_keyval_delete(
  MPI_Comm comm, int XT_UNUSED(dtype_keyval),
  void *XT_UNUSED(attribute_val), void *XT_UNUSED(extra_state))
{

  if (xt_ddt_internal_keyval != MPI_KEYVAL_INVALID) {
    xt_mpi_call(
      MPI_Type_free_keyval(&xt_ddt_internal_keyval), comm);
    xt_ddt_internal_keyval = MPI_KEYVAL_INVALID;
  }

  xt_mpi_call(
    MPI_Comm_free_keyval(&xt_ddt_cleanup_internal_keyval), comm);
  xt_ddt_cleanup_internal_keyval = MPI_KEYVAL_INVALID;

  return MPI_SUCCESS;
}

Xt_ddt xt_ddt_from_mpi_ddt(MPI_Datatype mpi_ddt) {

  XT_GPU_INSTR_PUSH(xt_ddt_from_mpi_ddt);

  // if xt_ddt_internal_keyval has not yet been created
  if (xt_ddt_internal_keyval == MPI_KEYVAL_INVALID) {

    // create keyval that will store xt_ddt's in MPI_Datatype's
    xt_mpi_call(
      MPI_Type_create_keyval(
        // MPI_TYPE_NULL_COPY_FN,
        xt_ddt_internal_keyval_copy,
        xt_ddt_internal_keyval_delete,
        &xt_ddt_internal_keyval, NULL),
      Xt_default_comm);

    // register a callback in MPI_Finalize (via an attribute to MPI_COMM_SELF)
    // to clean up xt_ddt_internal_keyval
    xt_mpi_call(
      MPI_Comm_create_keyval(
        MPI_COMM_NULL_COPY_FN,
        xt_ddt_cleanup_internal_keyval_delete,
        &xt_ddt_cleanup_internal_keyval, NULL),
      Xt_default_comm);
    xt_mpi_call(
      MPI_Comm_set_attr(
        MPI_COMM_SELF, xt_ddt_cleanup_internal_keyval, NULL),
      Xt_default_comm);
  }

  // get xt_ddt from MPI datatype, if available
  void *attr;
  int flag;
  xt_mpi_call(
    MPI_Type_get_attr(
      mpi_ddt, xt_ddt_internal_keyval, &attr, &flag), Xt_default_comm);

  // if the MPI datatype had already a xt_ddt attached to itself
  if (flag) {XT_GPU_INSTR_POP; return (Xt_ddt)attr;}

  // generate a xt_ddt from the MPI datatype
  Xt_ddt ddt = xt_ddt_new(mpi_ddt);

  // attach this xt_ddt to the MPI datatype for later use
  xt_mpi_call(
    MPI_Type_set_attr(
      mpi_ddt, xt_ddt_internal_keyval, ddt), Xt_default_comm);

  XT_GPU_INSTR_POP; return ddt;
}

size_t xt_ddt_get_pack_size_internal(Xt_ddt ddt) {

  return (ddt == NULL)?0:(ddt->pack_size);
}

size_t xt_ddt_get_pack_size(MPI_Datatype mpi_ddt) {

  return xt_ddt_get_pack_size_internal(xt_ddt_from_mpi_ddt(mpi_ddt));
}

static void xt_ddt_copy_displs(Xt_ddt ddt, enum xt_memtype memtype) {

  // count total number of displacements
  size_t total_displs_size = 0, count = ddt->count;
  for (size_t i = 0; i < count; ++i)
    total_displs_size += ddt->data[i].displ_count;

  // allocate displacements in specified memory type
  ssize_t *displs;
  size_t buffer_size = total_displs_size * sizeof(*displs);
  displs = xt_gpu_malloc(buffer_size, memtype);

  // copy displacements from host to specified memory type
  xt_gpu_memcpy(
    displs, ddt->data[0].displs[XT_MEMTYPE_HOST],
    buffer_size, memtype, XT_MEMTYPE_HOST);

  // set displacements for all data entries
  for (size_t i = 0, offset = 0; i < count; ++i) {
    ddt->data[i].displs[memtype] = displs + offset;
    offset += ddt->data[i].displ_count;
  }

  ddt->displs_available[memtype] = 1;
}

#define add_rhs_byte_displ(rtype,ptr,disp) \
  ((const rtype *)(const void *)((const unsigned char *)(ptr) + (disp)))

static void xt_ddt_pack_8(
  size_t count, ssize_t *restrict displs, const uint8_t *restrict src,
  uint8_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i)
    dst[i] = *add_rhs_byte_displ(uint8_t, src, displs[i]);
}

static void xt_ddt_pack_16(
  size_t count, ssize_t *restrict displs, const uint16_t *restrict src,
  uint16_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i)
    dst[i] = *add_rhs_byte_displ(uint16_t, src, + displs[i]);
}

static void xt_ddt_pack_32(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i)
    dst[i] = *add_rhs_byte_displ(uint32_t, src, + displs[i]);
}

static void xt_ddt_pack_32_2(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype) {
  uint32_t (*restrict dst_)[2] = (uint32_t(*)[2])dst;
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst_, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    const uint32_t *src_32 = add_rhs_byte_displ(uint32_t, src, displs[i]);
XtPragmaACC(loop independent)
    for (int j = 0; j < 2; ++j) dst_[i][j] = src_32[j];
  }
}

static void xt_ddt_pack_96(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype) {
  uint32_t (*restrict dst_)[3] = (uint32_t(*)[3])dst;
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst_, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    const uint32_t *src_32 = (const void *)((const unsigned char *)src + displs[i]);
XtPragmaACC(loop independent)
    for (int j = 0; j < 3; ++j) dst_[i][j] = src_32[j];
  }
}

static void xt_ddt_pack_64(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i)
    dst[i] = *add_rhs_byte_displ(uint64_t, src, displs[i]);
}

static void xt_ddt_pack_128(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype) {
  uint64_t (*restrict dst_)[2] = (uint64_t(*)[2])dst;
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst_, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    const uint64_t *src_64 = (const void *)((const unsigned char *)src + displs[i]);
XtPragmaACC(loop independent)
    for (int j = 0; j < 2; ++j) dst_[i][j] = src_64[j];
  }
}

static void xt_ddt_pack_160(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype) {
  uint32_t (*restrict dst_)[5] = (uint32_t(*)[5])dst;
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst_, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    const uint32_t *src_32 = (const void *)((const unsigned char *)src + displs[i]);
XtPragmaACC(loop independent)
    for (int j = 0; j < 5; ++j) dst_[i][j] = src_32[j];
  }
}

static void xt_ddt_pack_256(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t (*restrict dst)[4], enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    const uint64_t *src_64 = (const void *)((const unsigned char *)src + displs[i]);
XtPragmaACC(loop independent)
    for (int j = 0; j < 4; ++j) dst[i][j] = src_64[j];
  }
}

void xt_ddt_pack_internal(
  Xt_ddt ddt, const void *src, void *dst, enum xt_memtype memtype) {

  XT_GPU_INSTR_PUSH(xt_ddt_pack_internal);

  size_t dst_offset = 0;

  // if the displacements are not avaible in the required memory type
  if (!ddt->displs_available[memtype]) xt_ddt_copy_displs(ddt, memtype);

  size_t count = ddt->count;

  // for all sections with the same elemental datatype extent
  for (size_t i = 0; i < count; ++i) {

    struct xt_un_pack_kernels * kernel =
      &valid_kernels[ddt->data[i].kernel_idx];
    size_t displ_count = ddt->data[i].displ_count;
    kernel->pack(
      displ_count, ddt->data[i].displs[memtype], src,
      (unsigned char *)dst + dst_offset, memtype);

    dst_offset += displ_count * kernel->element_size;
  }
  XT_GPU_INSTR_POP;
}

void xt_ddt_pack(MPI_Datatype mpi_ddt, const void *src, void *dst) {

  XT_GPU_INSTR_PUSH(xt_ddt_pack);
  XT_GPU_INSTR_PUSH(xt_ddt_pack:initialise);

  enum xt_memtype src_memtype = xt_gpu_get_memtype(src);
  enum xt_memtype dst_memtype = xt_gpu_get_memtype(dst);

  size_t pack_size;
  void *orig_dst;
  /* pacify buggy -Wmaybe-uninitialized */
#if defined __GNUC__ && __GNUC__ <= 11
  pack_size = 0;
  orig_dst = NULL;
#endif

  // if the source and destination are in different memory types
  if (src_memtype != dst_memtype) {
    pack_size = xt_ddt_get_pack_size(mpi_ddt);
    orig_dst = dst;
    dst = xt_gpu_malloc(pack_size, src_memtype);
  }

  XT_GPU_INSTR_POP; //xt_ddt_pack:initialise

  xt_ddt_pack_internal(
    xt_ddt_from_mpi_ddt(mpi_ddt), src, dst, src_memtype);

  XT_GPU_INSTR_PUSH(xt_ddt_pack:finalise);

  // if the source and destination are in different memory types
  if (src_memtype != dst_memtype) {
    xt_gpu_memcpy(orig_dst, dst, pack_size, dst_memtype, src_memtype);
    xt_gpu_free(dst, src_memtype);
  }

  XT_GPU_INSTR_POP; // xt_ddt_pack:finalise
  XT_GPU_INSTR_POP; // xt_ddt_pack
}

static void xt_ddt_unpack_8(
  size_t count, ssize_t *restrict displs, const uint8_t *restrict src,
  uint8_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i)
    dst[displs[i]] = src[i];
}


static void xt_ddt_unpack_16(
  size_t count, ssize_t *restrict displs, const uint16_t *restrict src,
  uint16_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint16_t *dst_ = (void *)((unsigned char *)dst + displs[i]);
    dst_[0] = src[i];
  }
}

static void xt_ddt_unpack_32(
  size_t count, ssize_t *restrict displs, const uint32_t *restrict src,
  uint32_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint32_t *dst_ = (void *)((unsigned char *)dst + displs[i]);
    dst_[0] = src[i];
  }
}

static void xt_ddt_unpack_32_2(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[2],
  uint32_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint32_t *dst_32 = (void *)((unsigned char *)dst + displs[i]);
    dst_32[0] = src[i][0];
    dst_32[1] = src[i][1];
  }
}

static void xt_ddt_unpack_96(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[3],
  uint32_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint32_t *dst_32 = (void *)((unsigned char *)dst + displs[i]);
    dst_32[0] = src[i][0];
    dst_32[1] = src[i][1];
    dst_32[2] = src[i][2];
  }
}

static void xt_ddt_unpack_64(
  size_t count, ssize_t *restrict displs, const uint64_t *restrict src,
  uint64_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint64_t *dst_ = (void *)((unsigned char *)dst + displs[i]);
    dst_[0] = src[i];
  }
}

static void xt_ddt_unpack_128(
  size_t count, ssize_t *restrict displs, const uint64_t (*restrict src)[2],
  uint64_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint64_t *dst_64 = (void *)((unsigned char *)dst + displs[i]);
    dst_64[0] = src[i][0];
    dst_64[1] = src[i][1];
  }
}

static void xt_ddt_unpack_160(
  size_t count, ssize_t *restrict displs, const uint32_t (*restrict src)[5],
  uint32_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint32_t *dst_32 = (void *)((unsigned char *)dst + displs[i]);
    dst_32[0] = src[i][0];
    dst_32[1] = src[i][1];
    dst_32[2] = src[i][2];
    dst_32[3] = src[i][3];
    dst_32[4] = src[i][4];
  }
}

static void xt_ddt_unpack_256(
  size_t count, ssize_t *restrict displs, const uint64_t (*restrict src)[4],
  uint64_t *restrict dst, enum xt_memtype memtype) {
#ifndef _OPENACC
  (void)memtype;
#endif
XtPragmaACC(
  parallel loop independent deviceptr(src, dst, displs)
  if (memtype != XT_MEMTYPE_HOST))
  for (size_t i = 0; i < count; ++i) {
    uint64_t *dst_64 = (void *)((unsigned char *)dst + displs[i]);
    dst_64[0] = src[i][0];
    dst_64[1] = src[i][1];
    dst_64[2] = src[i][2];
    dst_64[3] = src[i][3];
  }
}

void xt_ddt_unpack_internal(
  Xt_ddt ddt, const void *src, void *dst, enum xt_memtype memtype) {

  XT_GPU_INSTR_PUSH(xt_ddt_unpack_internal);

  size_t src_offset = 0;

  // if the displacements are not avaible in the required memory type
  if (!ddt->displs_available[memtype]) xt_ddt_copy_displs(ddt, memtype);

  size_t count = ddt->count;

  // for all sections with the same elemental datatype extent
  for (size_t i = 0; i < count; ++i) {

    struct xt_un_pack_kernels * kernel =
      &valid_kernels[ddt->data[i].kernel_idx];
    size_t displ_count = ddt->data[i].displ_count;
    kernel->unpack(
      displ_count, ddt->data[i].displs[memtype],
      (unsigned char *)src + src_offset, dst, memtype);

    src_offset += displ_count * kernel->element_size;
  }
  XT_GPU_INSTR_POP;
}

void xt_ddt_unpack(MPI_Datatype mpi_ddt, const void *src, void *dst) {

  XT_GPU_INSTR_PUSH(xt_ddt_unpack);
  XT_GPU_INSTR_PUSH(xt_ddt_unpack:initialise);

  enum xt_memtype src_memtype = xt_gpu_get_memtype(src);
  enum xt_memtype dst_memtype = xt_gpu_get_memtype(dst);

  void *src_;
  // if the source and destination are in different memory types
  if (src_memtype != dst_memtype) {
    size_t pack_size = xt_ddt_get_pack_size(mpi_ddt);
    src_ = xt_gpu_malloc(pack_size, dst_memtype);
    xt_gpu_memcpy(src_, src, pack_size, dst_memtype, src_memtype);
  } else
    src_ = (void *)src;

  XT_GPU_INSTR_POP; // xt_ddt_unpack:initialise

  xt_ddt_unpack_internal(
    xt_ddt_from_mpi_ddt(mpi_ddt), src_, dst, dst_memtype);

  XT_GPU_INSTR_PUSH(xt_ddt_unpack:finalise);

  // if the source and destination are in different memory types
  if (src_memtype != dst_memtype) xt_gpu_free(src_, dst_memtype);

  XT_GPU_INSTR_POP; // xt_ddt_unpack:finalise
  XT_GPU_INSTR_POP; // xt_ddt_unpack
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
