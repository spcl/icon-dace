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


static size_t
scan_stripe(const int *disp, size_t disp_len, struct Xt_offset_ext *restrict v)
{
  if (disp_len<1) return 0;

  struct Xt_offset_ext x = (struct Xt_offset_ext){ disp[0], 1, 1 };
  size_t i = 0;
  for (size_t p = 1; p < disp_len; ++p) {
    int new_stride = disp[p] - disp[p-1];
    if (x.size == 1) {
      x.stride = new_stride;
      x.size = 2;
    } else if (new_stride == x.stride) {
      // x.size >= 2:
      x.size++;
    } else if (x.size > 2 || (x.size == 2 && x.stride == 1) ) {
      // we accept small contiguous vectors (nstrides==2, stride==1)
      v[i]= x;
      i++;
      x = (struct Xt_offset_ext){ disp[p], 1, 1 };
    } else { // x.size == 2, next offset doesn't match current stride
      // break up trivial vec:
      v[i].start = x.start;
      v[i].size = 1;
      v[i].stride = 1;
      i++;
      x.start += x.stride;
      x.size = 2;
      x.stride = new_stride;
    }
  }
  // tail cases:
  if (x.size > 2 || (x.size == 2 && x.stride == 1)) {
    v[i] = x;
    i++;
  } else if (x.size == 2) {
    v[i].start = x.start;
    v[i].size = 1;
    v[i].stride = 1;
    i++;
    v[i].start = x.start + x.stride;
    v[i].size = 1;
    v[i].stride = 1;
    i++;
  } else { // x.size == 1
    v[i].start = x.start;
    v[i].size = 1;
    v[i].stride = 1;
    i++;
  }

  return i;
}

static bool
match_simple_vec(size_t *pstart_, const struct Xt_offset_ext *v, size_t vlen,
                 MPI_Datatype old_type, MPI_Aint old_type_extent,
                 MPI_Aint *disp, MPI_Datatype *dt,
                 MPI_Comm comm) {
  // we only accept non-trivial matches (nsteps>2) with stride /= 1
  // using only one vector from v
  size_t p = *pstart_;
  if (p >= vlen) return false;
  int nstrides = v[p].size;
  int stride = v[p].stride;
  if (nstrides < 2 || stride == 1 ) return false;

  *pstart_ = p + 1;

  int disp_ = vlen > 1 ? v[p].start : 0;
  *disp = disp_ * old_type_extent;

  xt_mpi_call(MPI_Type_vector(nstrides, 1, stride, old_type, dt), comm);

  int start = v[p].start - disp_;
  if (start) {
    MPI_Datatype dt1 = *dt;

    // (start != 0) => add offset:
    MPI_Aint displacement = start * old_type_extent;
    xt_mpi_call(MPI_Type_create_hindexed(1, &(int){1}, &displacement, dt1, dt),
                comm);
    xt_mpi_call(MPI_Type_free(&dt1), comm);
  }
  return nstrides != 0;
}

/**
 * @return true if matched, false if not matched
 */
static bool
match_block_vec(size_t *pstart_, const struct Xt_offset_ext *v, size_t vlen,
                MPI_Datatype old_type, MPI_Aint old_type_extent,
                MPI_Aint *disp, MPI_Datatype *dt,
                MPI_Comm comm) {
  // using at least 3 vectors
  size_t p = *pstart_, pstart = p;
  if (p+2 >= vlen || v[p].stride != 1 || v[p+1].stride != 1 ) return false;
  int bl = v[p].size;
  assert(bl > 0);
  if (v[p+1].size != bl) return false;

  int vstride = v[p+1].start - v[p].start;

  p += 2;
  while( p < vlen && v[p].stride == 1 && v[p].size == bl &&
         v[p].start - v[p-1].start == vstride ) {
    p++;
  }
  size_t n = p - pstart;
  if (n<3) return false;
  *pstart_ = p;

  int disp_ = n == vlen ? 0 : v[pstart].start;
  *disp = disp_ * old_type_extent;

  xt_mpi_call(MPI_Type_vector((int)n, bl, vstride, old_type, dt), comm);

  int start = v[pstart].start - disp_;

  if (start) {
    MPI_Datatype dt1 = *dt;
    // (start != 0) => add offset:
    MPI_Aint displacement = start * old_type_extent;
    xt_mpi_call(MPI_Type_create_hindexed(1, &(int){1}, &displacement, dt1, dt),
                comm);
    xt_mpi_call(MPI_Type_free(&dt1), comm);
  }
  return n != 0;
}

static bool
match_contiguous(size_t *pstart_, const struct Xt_offset_ext *v, size_t vlen,
                 MPI_Datatype old_type, MPI_Aint old_type_extent,
                 MPI_Aint *restrict disp, MPI_Datatype *dt,
                 MPI_Comm comm) {
  size_t p = *pstart_;
  if (p >= vlen || v[p].stride != 1 || v[p].size < 2) return false;

  int disp_ = vlen > 1 ? v[p].start : 0;
  *disp = disp_ * old_type_extent;
  int d = v[p].start - disp_;

  if (!d)
    xt_mpi_call(MPI_Type_contiguous(v[p].size, old_type, dt), comm) ;
  else
    xt_mpi_call(MPI_Type_create_indexed_block(1, v[p].size, &d, old_type, dt),
                comm);

  *pstart_ = p+1;
  return v[p].size != 0;
}

static bool
match_indexed(size_t *pstart_, const struct Xt_offset_ext *v, size_t vlen,
              MPI_Datatype old_type, MPI_Aint old_type_extent,
              MPI_Aint *disp, MPI_Datatype *dt,
              MPI_Comm comm) {
  // we only accept non-trivial matches
  size_t p = *pstart_, pstart = p;
  if (p >= vlen || v[p].stride != 1 || v[p].size < 2) return false;

  do
    ++p;
  while (p < vlen && v[p].stride == 1);

  size_t n = p - pstart;

  if (n < 2) return false;
  *pstart_ = p;

  int start = n == vlen ? 0 : v[pstart].start;
  *disp = start * old_type_extent;
  int *restrict bl = xmalloc(2 * n * sizeof (*bl)),
    *restrict d = bl + n;
  bool hom_bl = true;
  d[0] = v[pstart].start - start;
  int bl0 = bl[0] = v[pstart].size;
  for (size_t i = 1; i < n; i++) {
    size_t iv = pstart + i;
    d[i] = v[iv].start - start;
    bl[i] = v[iv].size;
    hom_bl &= (bl[i] == bl0);
  }

  if (hom_bl) {
    xt_mpi_call(MPI_Type_create_indexed_block((int)n, bl0, d, old_type, dt),
                comm);
  } else {
    xt_mpi_call(MPI_Type_indexed((int)n, bl, d, old_type, dt), comm);
  }

  free(bl);
  return n != 0;
}

static void
gen_fallback_type(size_t set_start, size_t set_end,
                  const struct Xt_offset_ext *v, size_t vlen,
                  MPI_Datatype old_type, MPI_Aint old_type_extent,
                  MPI_Aint *disp,
                  MPI_Datatype *dt, MPI_Comm comm) {
  size_t ia = set_start;
  size_t ib = set_end;
  if (ib <= ia || ib > vlen) return;

  int n = 0;
  for (size_t i=ia; i < ib; i++)
    n += v[i].size;

  /* todo: given the guarantees for v that fceb584 introduced,
   * this check should never fire */
  assert(n>0);

  // generate absolute datatype if ia == 0 && ib == vlen,
  // else generate relative datatype that gets embedded by the caller
  int start = (ia == 0 && ib == vlen) ? 0 : v[ia].start;

  *disp = start * old_type_extent;

  int *restrict d = xmalloc(sizeof (*d) * (size_t)n);
  size_t p=0;
#ifndef NDEBUG
  /* did any element of v have non-positive size? */
  bool found_np = false;
#endif

  for (size_t i=ia; i < ib; i++) {
#ifndef NDEBUG
    found_np |= v[i].size <= 0;
#endif
    size_t v_i_size = (size_t)(v[i].size > 0 ? v[i].size : 0);
    for (size_t k=0; k < v_i_size; k++) {
      d[p] = v[i].start + (int)k * v[i].stride - start;
      p++;
    }
  }
  assert(!found_np);

  if (n==1 && d[0] == 0) {
    *dt = old_type;
  } else {
    xt_mpi_call(MPI_Type_create_indexed_block(n, 1, d, old_type, dt),
                comm);
  }
  free(d);
}

static MPI_Datatype
parse_stripe(const struct Xt_offset_ext *v, size_t vlen, MPI_Datatype old_type,
             MPI_Comm comm)
{
  /* [set_start,set_end) describes the prefix of non-matching
   * elements in v that then need to be handled with gen_fallback_type */
  size_t set_start = 0, set_end = 0;
  MPI_Aint old_type_lb, old_type_extent;
  xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                  &old_type_extent), comm);
  MPI_Aint *restrict wdisp
    = xmalloc(sizeof(MPI_Datatype) * (size_t)vlen
              + sizeof (MPI_Aint) * (size_t)vlen);
  MPI_Datatype *restrict wdt = (MPI_Datatype *)(wdisp + vlen);
  /* [p,vlen) is the part of v that still needs matching performed */
  /* m is the index of the next datatype and displacements to write
   * to wdt and wdisp respectively */
  size_t p = 0, m = 0;
  while (p<vlen) {
    /* depending on whether there is a non-empty prefix, the datatype
     * and displacement corresponding to a match need to be written
     * to wdt[m+1] and wdisp[m+1] or wdt[m] and wdisp[m] respectively */
    size_t mm = m + (set_start < set_end);
    if (match_block_vec(&p, v, vlen, old_type, old_type_extent,
                        wdisp+mm, wdt+mm, comm)
        || match_indexed(&p, v, vlen, old_type, old_type_extent,
                         wdisp+mm, wdt+mm, comm)
        || match_simple_vec(&p, v, vlen, old_type, old_type_extent,
                            wdisp+mm, wdt+mm, comm)
        || match_contiguous(&p, v, vlen, old_type, old_type_extent,
                            wdisp+mm, wdt+mm, comm) ) {
      /* in case a match is found, generate fallback datatype for
       * non-matching, preceding extents */
      if (set_start < set_end) {
        gen_fallback_type(set_start, set_end, v, vlen, old_type,
                          old_type_extent, wdisp+m, wdt+m, comm);
        m++;
      }
      m++;
      set_start = p;
    } else {
      /* assign ext investigated last to prefix */
      set_end = ++p;
    }
  }
  if (set_start <  set_end) {
    gen_fallback_type(set_start, set_end, v, vlen, old_type, old_type_extent,
                      wdisp+m, wdt+m, comm);
    m++;
  }
  size_t wlen = m;
  MPI_Datatype result_dt;
  if (wlen == 1 ) {
    assert(wdisp[0] == 0);
    if (wdt[0] == old_type)
      xt_mpi_call(MPI_Type_dup(old_type, wdt), comm);
    result_dt = wdt[0];
  } else {
    int *restrict wblocklength
      = wlen * sizeof (int) <= (vlen - wlen) * sizeof (*wdt)
      ? (void *)(wdt + wlen) : xmalloc(wlen * sizeof (*wblocklength));
    for(size_t i=0; i<wlen; i++)
      wblocklength[i] = 1;
    xt_mpi_call(MPI_Type_create_struct((int)wlen, wblocklength, wdisp,
                                       wdt, &result_dt), comm);
    if (wlen * sizeof (int) > (vlen - wlen) * sizeof (*wdt))
      free(wblocklength);
    for (size_t i = 0; i < wlen; i++)
      if (wdt[i] != old_type)
        xt_mpi_call(MPI_Type_free(wdt+i), comm);
  }
  xt_mpi_call(MPI_Type_commit(&result_dt), comm);
  free(wdisp);
  return result_dt;
}

MPI_Datatype
xt_mpi_generate_datatype_stripe(const struct Xt_offset_ext *v,
                                int count, MPI_Datatype old_type,
                                MPI_Comm comm)
{
  size_t count_ = (size_t)0;
  for (int i=0; i<count; ++i)
    count_ += (size_t)(v[i].size > 0);
  if (count_ < 1) return MPI_DATATYPE_NULL;
  struct Xt_offset_ext *v_comp;
  if ((size_t)count != count_) {
    v_comp = xmalloc(count_ * sizeof (*v_comp));
    for (size_t i=0, j=0; i<(size_t)count; ++i) {
      v_comp[j] = v[i];
      j+= v[i].size > 0;
    }
  } else
    v_comp = (struct Xt_offset_ext *)v;
  MPI_Datatype dt = parse_stripe(v_comp, count_, old_type, comm);
  if ((size_t)count != count_)
    free(v_comp);
  return dt;
}


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

  struct Xt_offset_ext *v = xmalloc(sizeof(*v) * (size_t)disp_len);
  size_t vlen = scan_stripe(disp, (size_t)disp_len, v);
  MPI_Datatype dt = parse_stripe(v, vlen, old_type, comm);
  free(v);
  return dt;
}

/* functions to handle optimizations on communicators */
static int xt_mpi_comm_internal_keyval;

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
  xt_mpi_call(MPI_Comm_free_keyval(&xt_mpi_comm_internal_keyval),
              Xt_default_comm);
}

static struct xt_mpi_comm_internal_attr *
xt_mpi_comm_get_internal_attr(MPI_Comm comm)
{
  int attr_found;
  void *attr_val;
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

#if __GNUC__ == 11
  /* GCC 11 has no means to specify that the special value pointer
   * MPI_STATUSES_IGNORE does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
  xt_mpi_call(MPI_Testsome(*num_req, req, &done_count, ops_completed,
                           MPI_STATUSES_IGNORE), comm);
#if __GNUC__ == 11
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


/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
