/**
 * @file xt_mpi.c
 *
 * @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
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

#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_core.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"

#if ! (HAVE_DECL___BUILTIN_CTZL || HAVE_DECL___BUILTIN_CLZL)       \
  && (HAVE_DECL___LZCNT && SIZEOF_LONG == SIZEOF_INT               \
      || HAVE_DECL___LZCNT64 && SIZEOF_LONG == 8 && CHAR_BIT == 8)
#include <intrin.h>
#endif

//! COMPACT_DT enables the anlysis of displacements in order to give a
//! more compact description to the datatype generators of MPI. For strong
//! enough MPI implementations this not required. Then you can undefine
//! COMPACT_DT and save some prcessing time within yaxt without losing communication
//! performance.
#define COMPACT_DT

static MPI_Datatype
xt_mpi_generate_compact_datatype_block(const int *disp, const int *blocklengths,
                                       int count, MPI_Datatype old_type,
                                       MPI_Comm comm);
static MPI_Datatype
xt_mpi_generate_compact_datatype(int const *disp, int disp_len,
                                 MPI_Datatype old_type, MPI_Comm comm);


//taken from http://beige.ucs.indiana.edu/I590/node85.html
void xt_mpi_error(int error_code, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  char error_string[MPI_MAX_ERROR_STRING];
  int length_of_error_string, error_class;

  MPI_Error_class(error_code, &error_class);
  MPI_Error_string(error_class, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Error_string(error_code, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Abort(comm, error_code);
}

#ifndef COMPACT_DT
static MPI_Datatype copy_mpi_datatype(MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_dup(old_type, &datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_simple(int displacement, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_create_indexed_block(1, 1, &displacement, old_type,
                                                &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_contiguous(int displacement, int blocklength,
                            MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  if (displacement == 0)
    xt_mpi_call(MPI_Type_contiguous(blocklength, old_type, &datatype),
                    comm);
  else
    xt_mpi_call(MPI_Type_create_indexed_block(1, blocklength,
                                                  &displacement, old_type,
                                                  &datatype), comm);

  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;

}

static MPI_Datatype
gen_mpi_datatype_vector(int count, int blocklength, int stride,
                        int offset, MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_vector(count, blocklength, stride, old_type,
                              &datatype), comm);
  if (offset != 0) {

    MPI_Datatype datatype_;
    int hindexed_blocklength = 1;
    MPI_Aint old_type_size, old_type_lb;

    xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                    &old_type_size), comm);

    MPI_Aint displacement = offset * old_type_size;

    xt_mpi_call(MPI_Type_create_hindexed(1, &hindexed_blocklength,
                                         &displacement, datatype, &datatype_),
                comm);
    xt_mpi_call(MPI_Type_free(&datatype), comm);
    datatype = datatype_;
  }
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_indexed_block(int const * displacements, int blocklength,
                               int count, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_create_indexed_block(count, blocklength,
                                                (void *)displacements,
                                                old_type, &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_indexed(const int *displacements, const int *blocklengths,
                         int count, MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_indexed(count, (int*)blocklengths, (void*)displacements,
                                   old_type, &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static inline int
check_for_vector_type(const int *displacements, const int *blocklengths,
                      int count) {

  int blocklength = blocklengths[0];

  for (int i = 1; i < count; ++i)
    if (blocklengths[i] != blocklength)
      return 0;

  int stride = displacements[1] - displacements[0];

  for (int i = 1; i + 1 < count; ++i)
    if (displacements[i+1] - displacements[i] != stride)
      return 0;

  return 1;
}

static inline int check_for_indexed_block_type(const int *blocklengths,
                                               int count) {

  int blocklength = blocklengths[0];

  for (int i = 1; i < count; ++i)
    if (blocklengths[i] != blocklength)
      return 0;

  return 1;
}
#endif

MPI_Datatype
xt_mpi_generate_datatype_block(const int *displacements,
                               const int *blocklengths,
                               int count, MPI_Datatype old_type,
                               MPI_Comm comm) {

#ifdef COMPACT_DT
  return xt_mpi_generate_compact_datatype_block(displacements, blocklengths,
                                                count, old_type, comm);
#else
  MPI_Datatype datatype;

  if (count == 0)
    datatype = MPI_DATATYPE_NULL;
  else if (count == 1 && blocklengths[0] == 1 && displacements[0] == 0)
    datatype = copy_mpi_datatype(old_type, comm);
  else if (count == 1 && blocklengths[0] == 1)
    datatype = gen_mpi_datatype_simple(displacements[0], old_type, comm);
  else if (count == 1)
    datatype = gen_mpi_datatype_contiguous(displacements[0], blocklengths[0],
                                           old_type, comm);
  else if (check_for_vector_type(displacements, blocklengths, count))
    datatype = gen_mpi_datatype_vector(count, blocklengths[0],
                                       displacements[1] - displacements[0],
                                       displacements[0], old_type, comm);
  else if (check_for_indexed_block_type(blocklengths, count))
    datatype = gen_mpi_datatype_indexed_block(displacements, blocklengths[0],
                                              count, old_type, comm);
  else
    datatype = gen_mpi_datatype_indexed(displacements, blocklengths, count,
                                        old_type, comm);

  return datatype;
#endif
}

MPI_Datatype
xt_mpi_generate_datatype(int const * displacements, int count,
                         MPI_Datatype old_type, MPI_Comm comm)
{
  if (count <= 0)
    return MPI_DATATYPE_NULL;

#ifdef COMPACT_DT
  return xt_mpi_generate_compact_datatype(displacements, count, old_type, comm);
#else
  int * blocklengths = xmalloc((size_t)count * sizeof(*blocklengths));
  int new_count = 0;
  {
    int i = 0;
    do {
      int j = 1;
      while (i + j < count && displacements[i] + j == displacements[i + j])
        ++j;
      blocklengths[new_count++] = j;
      i += j;
    } while (i < count);
  }

  int * tmp_displ = NULL;
  const int *displ;

  if (new_count != count) {

    tmp_displ = xmalloc((size_t)new_count * sizeof(*tmp_displ));

    int offset = 0;

    for (int i = 0; i < new_count; ++i) {

      tmp_displ[i] = displacements[offset];
      offset += blocklengths[i];
    }

    displ = tmp_displ;
  } else
    displ = displacements;

  MPI_Datatype datatype;

  datatype = xt_mpi_generate_datatype_block(displ, blocklengths, new_count,
                                            old_type, comm);

  free(blocklengths);

  free(tmp_displ);

  return datatype;
#endif
}

size_t
xt_disp2ext_count(size_t disp_len, const int *disp)
{
  if (!disp_len) return 0;
  size_t i = 0;
  int cur_stride = 1, cur_size = 1;
  int last_disp = disp[0];
  for (size_t p = 1; p < disp_len; ++p) {
    int new_disp = disp[p];
    int new_stride = new_disp - last_disp;
    if (cur_size == 1) {
      cur_stride = new_stride;
      cur_size = 2;
    } else if (new_stride == cur_stride) {
      // cur_size >= 2:
      cur_size++;
    } else if (cur_size > 2 || (cur_size == 2 && cur_stride == 1) ) {
      // we accept small contiguous vectors (nstrides==2, stride==1)
      i++;
      cur_stride = 1;
      cur_size = 1;
    } else { // cur_size == 2, next offset doesn't match current stride
      // break up trivial vec:
      i++;
      cur_size = 2;
      cur_stride = new_stride;
    }
    last_disp = new_disp;
  }
  // tail cases:
  if (cur_size > 2 || (cur_size == 2 && cur_stride == 1)) {
    i++;
  } else if (cur_size == 2) {
    i+=2;
  } else { // cur_size == 1
    i++;
  }

  return i;
}

size_t
xt_disp2ext(size_t disp_len, const int *disp,
            struct Xt_offset_ext *restrict v)
{
  if (disp_len<1) return 0;

  int cur_start = disp[0], cur_stride = 1, cur_size = 1;
  int last_disp = cur_start;
  size_t i = 0;
  for (size_t p = 1; p < disp_len; ++p) {
    int new_disp = disp[p];
    int new_stride = new_disp - last_disp;
    if (cur_size == 1) {
      cur_stride = new_stride;
      cur_size = 2;
    } else if (new_stride == cur_stride) {
      // cur_size >= 2:
      cur_size++;
    } else if (cur_size > 2 || (cur_size == 2 && cur_stride == 1) ) {
      // we accept small contiguous vectors (nstrides==2, stride==1)
      v[i] = (struct Xt_offset_ext){ .start = cur_start, .stride = cur_stride,
                                     .size = cur_size };
      i++;
      cur_start = new_disp;
      cur_stride = 1;
      cur_size = 1;
    } else { // cur_size == 2, next offset doesn't match current stride
      // break up trivial vec:
      v[i].start = cur_start;
      v[i].size = 1;
      v[i].stride = 1;
      i++;
      cur_start += cur_stride;
      cur_size = 2;
      cur_stride = new_stride;
    }
    last_disp = new_disp;
  }
  // tail cases:
  if (cur_size > 2 || (cur_size == 2 && cur_stride == 1)) {
    v[i] = (struct Xt_offset_ext){ .start = cur_start, .stride = cur_stride,
                                   .size = cur_size };
    i++;
  } else if (cur_size == 2) {
    v[i].start = cur_start;
    v[i].size = 1;
    v[i].stride = 1;
    i++;
    v[i].start = cur_start + cur_stride;
    v[i].size = 1;
    v[i].stride = 1;
    i++;
  } else { // cur_size == 1
    v[i].start = cur_start;
    v[i].size = 1;
    v[i].stride = 1;
    i++;
  }

  return i;
}

#define XT_MPI_STRP_PRS_PREFIX
#define XT_MPI_STRP_PRS_UNITSTRIDE 1
#define XT_MPI_STRP_PRS_AOFS_TYPE int
#define XT_MPI_STRP_PRS_DISP_ADJUST(val) ((val) * old_type_extent)
#define XT_MPI_STRP_PRS_BLOCK_VEC_CREATE MPI_Type_vector
#define XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE MPI_Type_create_indexed_block
#define XT_MPI_STRP_PRS_INDEXED_CREATE MPI_Type_indexed
#define XT_MPI_STRP_PRS_FALLBACK_NEEDS_OLD_TYPE_EXTENT
#include "xt_mpi_stripe_parse.c"
#undef XT_MPI_STRP_PRS_PREFIX
#undef XT_MPI_STRP_PRS_UNITSTRIDE
#undef XT_MPI_STRP_PRS_AOFS_TYPE
#undef XT_MPI_STRP_PRS_DISP_ADJUST
#undef XT_MPI_STRP_PRS_BLOCK_VEC_CREATE
#undef XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE
#undef XT_MPI_STRP_PRS_INDEXED_CREATE
#undef XT_MPI_STRP_PRS_FALLBACK_NEEDS_OLD_TYPE_EXTENT

#if MPI_VERSION < 3
static inline int
XtMPI_Type_create_hindexed_block(int count, int blocklength,
                                 const MPI_Aint array_of_displacements[],
                                 MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  size_t count_ = count > 0 ? (size_t)count : 0;
  MPI_Datatype *restrict oldtypes = xmalloc(count_ * sizeof (*oldtypes)
                                   + count_ * sizeof (int));
  int *restrict blocklengths = (int *)(oldtypes + count_);
  for (size_t i = 0; i < count_; ++i) {
    blocklengths[i] = blocklength;
    oldtypes[i] = oldtype;
  }
  int rc = MPI_Type_create_struct(count, blocklengths,
                                  CAST_MPI_SEND_BUF(array_of_displacements),
                                  oldtypes, newtype);
  free(oldtypes);
  return rc;
}

#define MPI_Type_create_hindexed_block XtMPI_Type_create_hindexed_block
#endif

#define XT_MPI_STRP_PRS_PREFIX a
#define XT_MPI_STRP_PRS_UNITSTRIDE old_type_extent
#define XT_MPI_STRP_PRS_AOFS_TYPE MPI_Aint
#define XT_MPI_STRP_PRS_DISP_ADJUST(val) (val)
#define XT_MPI_STRP_PRS_BLOCK_VEC_CREATE MPI_Type_create_hvector
#define XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE MPI_Type_create_hindexed_block
#define XT_MPI_STRP_PRS_INDEXED_CREATE MPI_Type_create_hindexed
#include "xt_mpi_stripe_parse.c"


static MPI_Datatype
xt_mpi_generate_compact_datatype_block(const int *disp, const int *blocklengths,
                                       int count, MPI_Datatype old_type,
                                       MPI_Comm comm)
{
  size_t count_ = (size_t)0;
  for (int i=0; i<count; ++i)
    count_ += (size_t)(blocklengths[i] > 0);
  if (count_ < 1) return MPI_DATATYPE_NULL;
  struct Xt_offset_ext *restrict v = xmalloc(sizeof(*v) * count_);
  size_t j=0;
  for (size_t i=0; i<(size_t)count; ++i) {
    v[j].start = disp[i];
    v[j].stride = 1;
    int bl = blocklengths[i];
    v[j].size = bl;
    j += (size_t)(bl > 0);
  }
  MPI_Datatype dt = parse_stripe(v, count_, old_type, comm);
  free(v);
  return dt;
}

static MPI_Datatype
xt_mpi_generate_compact_datatype(const int *disp, int disp_len,
                                 MPI_Datatype old_type, MPI_Comm comm)
{
  if (disp_len < 1) return MPI_DATATYPE_NULL;

  size_t vlen = xt_disp2ext_count((size_t)disp_len, disp);
  struct Xt_offset_ext *v = xmalloc(sizeof(*v) * vlen);
  xt_disp2ext((size_t)disp_len, disp, v);
  MPI_Datatype dt = parse_stripe(v, vlen, old_type, comm);
  free(v);
  return dt;
}

/* functions to handle optimizations on communicators */
static int xt_mpi_comm_internal_keyval = MPI_KEYVAL_INVALID;

typedef unsigned long used_map_elem;

enum {
  used_map_elem_bits = sizeof (used_map_elem) * CHAR_BIT,
};

struct xt_mpi_comm_internal_attr {
  int refcount;
  unsigned used_map_size;
  used_map_elem used_map[];
};

static int
xt_mpi_comm_internal_keyval_copy(
  MPI_Comm XT_UNUSED(oldcomm), int XT_UNUSED(keyval),
  void *XT_UNUSED(extra_state), void *XT_UNUSED(attribute_val_in),
  void *attribute_val_out, int *flag)
{
  struct xt_mpi_comm_internal_attr *new_comm_attr
    = malloc(sizeof (struct xt_mpi_comm_internal_attr)
             + sizeof (used_map_elem));
  int retval;
  if (new_comm_attr)
  {
    new_comm_attr->refcount = 1;
    new_comm_attr->used_map_size = 1;
    new_comm_attr->used_map[0] = 1U;
    *(void **)attribute_val_out = new_comm_attr;
    *flag = 1;
    retval = MPI_SUCCESS;
  } else {
    *flag = 0;
    retval = MPI_ERR_NO_MEM;
  }
  return retval;
}

static int
xt_mpi_comm_internal_keyval_delete(
  MPI_Comm XT_UNUSED(comm), int XT_UNUSED(comm_keyval),
  void *attribute_val, void *XT_UNUSED(extra_state))
{
  free(attribute_val);
  return MPI_SUCCESS;
}

static int xt_mpi_tag_ub_val;

void
xt_mpi_init(void) {
  assert(xt_mpi_comm_internal_keyval == MPI_KEYVAL_INVALID);
  xt_mpi_call(MPI_Comm_create_keyval(xt_mpi_comm_internal_keyval_copy,
                                     xt_mpi_comm_internal_keyval_delete,
                                     &xt_mpi_comm_internal_keyval, NULL),
              Xt_default_comm);
  void *attr;
  int flag;
  xt_mpi_call(MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &attr, &flag),
              MPI_COMM_WORLD);
  assert(flag);
  xt_mpi_tag_ub_val = *(int *)attr;
}

void
xt_mpi_finalize(void) {
  assert(xt_mpi_comm_internal_keyval != MPI_KEYVAL_INVALID);
  xt_mpi_call(MPI_Comm_free_keyval(&xt_mpi_comm_internal_keyval),
              Xt_default_comm);
}

static struct xt_mpi_comm_internal_attr *
xt_mpi_comm_get_internal_attr(MPI_Comm comm)
{
  int attr_found;
  void *attr_val;
  assert(xt_mpi_comm_internal_keyval != MPI_KEYVAL_INVALID);
  xt_mpi_call(MPI_Comm_get_attr(comm, xt_mpi_comm_internal_keyval,
                                &attr_val, &attr_found),
              comm);
  return attr_found ? attr_val : NULL;
}

#if HAVE_DECL___BUILTIN_CTZL
#define ctzl(v) (__builtin_ctzl(v))
#elif HAVE_DECL___BUILTIN_CLZL                                     \
  || HAVE_DECL___LZCNT && SIZEOF_LONG == SIZEOF_INT                \
  || HAVE_DECL___LZCNT64 && SIZEOF_LONG == 8 && CHAR_BIT == 8
static inline int
ctzl(unsigned long v) {
  enum {
    ulong_bits = sizeof (unsigned long) * CHAR_BIT,
  };
  /* clear all but lowest 1 bit */
  v = v & ~(v - 1);
  int c = ulong_bits - 1 - (int)
#if HAVE_DECL___BUILTIN_CTZL
    __builtin_clzl(v)
#elif HAVE_DECL___LZCNT && SIZEOF_LONG == SIZEOF_INT
    __lzcnt(v)
#else
    __lzcnt64(v)
#endif
    ;
  return c;
}
#else
static inline int
ctzl(unsigned long v) {
  enum {
    ulong_bits = sizeof (unsigned long) * CHAR_BIT,
  };
  // c will be the number of zero bits on the right
  unsigned int c = ulong_bits;
  v &= (unsigned long)-(long)v;
  if (v) c--;
#if SIZEOF_UNSIGNED_LONG * CHAR_BIT == 64
  if (v & UINT64_C(0x00000000ffffffff)) c -= 32;
  if (v & UINT64_C(0x0000ffff0000ffff)) c -= 16;
  if (v & UINT64_C(0x00ff00ff00ff00ff)) c -= 8;
  if (v & UINT64_C(0x0f0f0f0f0f0f0f0f)) c -= 4;
  if (v & UINT64_C(0x3333333333333333)) c -= 2;
  if (v & UINT64_C(0x5555555555555555)) c -= 1;
#elif SIZEOF_UNSIGNED_LONG * CHAR_BIT == 32
  if (v & 0x0000FFFFUL) c -= 16;
  if (v & 0x00FF00FFUL) c -= 8;
  if (v & 0x0F0F0F0FUL) c -= 4;
  if (v & 0x33333333UL) c -= 2;
  if (v & 0x55555555UL) c -= 1;
#else
  error "Unexpected size of long.\n"
#endif
  return (int)c;
}
#endif

MPI_Comm
xt_mpi_comm_smart_dup(MPI_Comm comm, int *tag_offset)
{
  MPI_Comm comm_dest;
  struct xt_mpi_comm_internal_attr *comm_xt_attr_val
    = xt_mpi_comm_get_internal_attr(comm);
  size_t position = 0;
  int refcount = comm_xt_attr_val ? comm_xt_attr_val->refcount : 0;
  if (comm_xt_attr_val
      && (refcount + 1) < xt_mpi_tag_ub_val / xt_mpi_num_tags) {
    comm_dest = comm;
    comm_xt_attr_val->refcount = ++refcount;
    size_t used_map_size = comm_xt_attr_val->used_map_size;
    while (position < used_map_size
           && comm_xt_attr_val->used_map[position] == ~(used_map_elem)0)
      ++position;
    if (position >= used_map_size) {
      /* sadly, we need to recreate the value to enlarge it */
      struct xt_mpi_comm_internal_attr *new_comm_xt_attr_val
        = xmalloc(sizeof (*new_comm_xt_attr_val)
                  + (used_map_size + 1) * sizeof (used_map_elem));
      new_comm_xt_attr_val->refcount = refcount;
      new_comm_xt_attr_val->used_map_size = (unsigned)(used_map_size + 1);
      for (size_t i = 0; i < used_map_size; ++i)
        new_comm_xt_attr_val->used_map[i] = comm_xt_attr_val->used_map[i];
      new_comm_xt_attr_val->used_map[used_map_size] = 1U;
      position *= used_map_elem_bits;
      assert(xt_mpi_comm_internal_keyval != MPI_KEYVAL_INVALID);
      xt_mpi_call(MPI_Comm_set_attr(comm_dest, xt_mpi_comm_internal_keyval,
                                    new_comm_xt_attr_val), comm_dest);
    } else {
      /* not all bits are set, find first unset position and insert */
      used_map_elem used_map_entry = comm_xt_attr_val->used_map[position],
        unset_lsb = ~used_map_entry & (used_map_entry + 1),
        bit_pos = (used_map_elem)ctzl(unset_lsb);
      comm_xt_attr_val->used_map[position] = used_map_entry | unset_lsb;
      position = position * used_map_elem_bits + (size_t)bit_pos;
    }
  } else {
    struct xt_mpi_comm_internal_attr *comm_attr
      = xmalloc(sizeof (*comm_attr) + sizeof (used_map_elem));
    comm_attr->refcount = 1;
    comm_attr->used_map_size = 1;
    comm_attr->used_map[0] = 1U;
    xt_mpi_call(MPI_Comm_dup(comm, &comm_dest), comm);
    assert(xt_mpi_comm_internal_keyval != MPI_KEYVAL_INVALID);
    xt_mpi_call(MPI_Comm_set_attr(comm_dest, xt_mpi_comm_internal_keyval,
                                  comm_attr), comm_dest);
  }
  *tag_offset = (int)(position * xt_mpi_num_tags);
  return comm_dest;
}

void
xt_mpi_comm_smart_dedup(MPI_Comm *comm, int tag_offset)
{
  struct xt_mpi_comm_internal_attr *comm_xt_attr_val
    = xt_mpi_comm_get_internal_attr(*comm);
  int refcount = comm_xt_attr_val ? --(comm_xt_attr_val->refcount) : 0;
  if (refcount < 1) {
    xt_mpi_call(MPI_Comm_free(comm), MPI_COMM_WORLD);
    *comm = MPI_COMM_NULL;
  } else {
    size_t position = (size_t)tag_offset / xt_mpi_num_tags,
      map_elem = position / used_map_elem_bits,
      in_elem_bit = position % used_map_elem_bits;
    comm_xt_attr_val->used_map[map_elem] &= ~((used_map_elem)1 << in_elem_bit);
  }
}

void
xt_mpi_comm_mark_exclusive(MPI_Comm comm) {
  struct xt_mpi_comm_internal_attr *comm_attr
    = xmalloc(sizeof (*comm_attr) + sizeof (used_map_elem));
  comm_attr->refcount = 1;
  comm_attr->used_map_size = 1;
  comm_attr->used_map[0] = 1U;
  assert(xt_mpi_comm_internal_keyval != MPI_KEYVAL_INVALID);
  xt_mpi_call(MPI_Comm_set_attr(comm, xt_mpi_comm_internal_keyval,
                                comm_attr), comm);
}

bool
xt_mpi_test_some(int *restrict num_req,
                 MPI_Request *restrict req,
                 int *restrict ops_completed, MPI_Comm comm)
{
  int done_count;
  size_t num_req_ = (size_t)*num_req;

#if __GNUC__ >= 11 && __GNUC__ <= 13
  /* GCC 11 has no means to specify that the special value pointer
   * MPI_STATUSES_IGNORE does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
  xt_mpi_call(MPI_Testsome(*num_req, req, &done_count, ops_completed,
                           MPI_STATUSES_IGNORE), comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13
#pragma GCC diagnostic pop
#endif

  if (done_count != MPI_UNDEFINED) {
    if (num_req_ > (size_t)done_count) {
      for (size_t i = 0, j = num_req_;
           i < (size_t)done_count && j >= num_req_ - (size_t)done_count;
           ++i)
        if (ops_completed[i] < (int)num_req_ - done_count) {
          while (req[--j] == MPI_REQUEST_NULL);
          req[ops_completed[i]] = req[j];
        }
      num_req_ -= (size_t)done_count;
    }
    else
      num_req_ = 0;
  }
  *num_req = (int)num_req_;
  return num_req_ == 0;
}

#ifdef _OPENMP
bool
xt_mpi_test_some_mt(int *restrict num_req,
                    MPI_Request *restrict req,
                    int *restrict ops_completed, MPI_Comm comm)
{
  int done_count;
  size_t num_req_ = (size_t)*num_req;

  size_t num_threads = (size_t)omp_get_num_threads(),
    tid = (size_t)omp_get_thread_num();
  size_t start_req = (num_req_ * tid) / num_threads,
    nreq_ = (num_req_ * (tid+1)) / num_threads - start_req;

  for (size_t i = start_req; i < start_req + nreq_; ++i)
    ops_completed[i] = -1;
#if __GNUC__ >= 11 && __GNUC__ <= 13
  /* GCC 11 has no means to specify that the special value pointer
   * MPI_STATUSES_IGNORE does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
  xt_mpi_call(MPI_Testsome((int)nreq_, req+start_req, &done_count,
                           ops_completed+start_req, MPI_STATUSES_IGNORE), comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13
#pragma GCC diagnostic pop
#endif
  if (done_count == MPI_UNDEFINED)
    done_count = 0;
#pragma omp barrier
#pragma omp atomic
  *num_req -= done_count;
#pragma omp barrier
  done_count = (int)num_req_ - *num_req;
#pragma omp single
  {
    if (num_req_ > (size_t)done_count) {
      for (size_t i = 0, j = 0; i < num_req_; ++i)
        if (req[i] != MPI_REQUEST_NULL)
          req[j++] = req[i];
    }
    *num_req = (int)num_req_ - done_count;
  }
  num_req_ -= (size_t)done_count;
  return num_req_ == 0;
}
#endif


/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
