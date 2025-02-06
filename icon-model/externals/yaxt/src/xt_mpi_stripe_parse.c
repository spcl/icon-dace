/**
 * @file xt_mpi_stripe_parse.c
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


#define XT_TOKEN_PASTE2_(a,b) a##b
#define XT_TOKEN_PASTE2(a,b) XT_TOKEN_PASTE2_(a,b)
#define XT_TOKEN_PASTE3_(a,b,c) a##b##c
#define XT_TOKEN_PASTE3(a,b,c) XT_TOKEN_PASTE3_(a,b,c)

#define XT_MPI_OFFSET_EXT \
  struct XT_TOKEN_PASTE3(Xt_,XT_MPI_STRP_PRS_PREFIX,offset_ext)
#define XT_MPI_STRP_PRS_MATCH_BLOCK_VEC \
  XT_TOKEN_PASTE2(XT_MPI_STRP_PRS_PREFIX,match_block_vec)
#define XT_MPI_STRP_PRS_MATCH_INDEXED \
  XT_TOKEN_PASTE2(XT_MPI_STRP_PRS_PREFIX,match_indexed)
#define XT_MPI_STRP_PRS_MATCH_SIMPLE_VEC \
  XT_TOKEN_PASTE2(XT_MPI_STRP_PRS_PREFIX,match_simple_vec)
#define XT_MPI_STRP_PRS_MATCH_CONTIGUOUS \
  XT_TOKEN_PASTE2(XT_MPI_STRP_PRS_PREFIX,match_contiguous)
#define XT_MPI_STRP_PRS_GEN_FALLBACK_TYPE                \
  XT_TOKEN_PASTE2(XT_MPI_STRP_PRS_PREFIX,gen_fallback_type)
#define XT_MPI_STRP_PRS_ENTRY \
  XT_TOKEN_PASTE3(parse_,XT_MPI_STRP_PRS_PREFIX,stripe)
#define XT_MPI_STRP_PRS_DRIVER \
  XT_TOKEN_PASTE3(xt_mpi_generate_datatype_,XT_MPI_STRP_PRS_PREFIX,stripe)

/**
 * @return true if matched, false if not matched
 */
static bool
XT_MPI_STRP_PRS_MATCH_BLOCK_VEC
  (size_t *pstart_,
   const XT_MPI_OFFSET_EXT *v,
   size_t vlen,
   MPI_Datatype old_type, MPI_Aint old_type_extent,
   MPI_Aint *disp, MPI_Datatype *dt,
   MPI_Comm comm) {
  // using at least 3 vectors
  size_t p = *pstart_, pstart = p;
  if (p+2 >= vlen || v[p].stride != XT_MPI_STRP_PRS_UNITSTRIDE
      || v[p+1].stride != XT_MPI_STRP_PRS_UNITSTRIDE ) return false;
  int bl = v[p].size;
  assert(bl > 0);
  if (v[p+1].size != bl) return false;

  XT_MPI_STRP_PRS_AOFS_TYPE vstride = v[p+1].start - v[p].start;

  p += 2;
  while( p < vlen && v[p].stride == XT_MPI_STRP_PRS_UNITSTRIDE
         && v[p].size == bl && v[p].start - v[p-1].start == vstride ) {
    p++;
  }
  size_t n = p - pstart;
  if (n<3) return false;
  *pstart_ = p;

  XT_MPI_STRP_PRS_AOFS_TYPE disp_ = n == vlen ? 0 : v[pstart].start;
  *disp = XT_MPI_STRP_PRS_DISP_ADJUST(disp_);

  xt_mpi_call(
    XT_MPI_STRP_PRS_BLOCK_VEC_CREATE((int)n, bl, vstride, old_type, dt), comm);

  XT_MPI_STRP_PRS_AOFS_TYPE start = v[pstart].start - disp_;

  if (start) {
    MPI_Datatype dt1 = *dt;
    // (start != 0) => add offset:
    xt_mpi_call(MPI_Type_create_hindexed(
                  1, &(int){1},
                  &(MPI_Aint){XT_MPI_STRP_PRS_DISP_ADJUST(start)}, dt1, dt),
                comm);
    xt_mpi_call(MPI_Type_free(&dt1), comm);
  }
  return n != 0;
}

static bool
XT_MPI_STRP_PRS_MATCH_INDEXED(
  size_t *pstart_,
  const XT_MPI_OFFSET_EXT *v,
  size_t vlen,
  MPI_Datatype old_type, MPI_Aint old_type_extent,
  MPI_Aint *disp, MPI_Datatype *dt,
  MPI_Comm comm) {
  // we only accept non-trivial matches
  size_t p = *pstart_, pstart = p;
  if (p >= vlen || v[p].stride != XT_MPI_STRP_PRS_UNITSTRIDE || v[p].size < 2)
    return false;

  do
    ++p;
  while (p < vlen && v[p].stride == XT_MPI_STRP_PRS_UNITSTRIDE);

  size_t n = p - pstart;

  if (n < 2) return false;
  *pstart_ = p;

  XT_MPI_STRP_PRS_AOFS_TYPE start = n == vlen ? 0 : v[pstart].start;
  *disp = XT_MPI_STRP_PRS_DISP_ADJUST(start);
  XT_MPI_STRP_PRS_AOFS_TYPE *restrict d
    = xmalloc(n * sizeof (int) + n * sizeof (*d));
  int *restrict bl = (int * restrict)(d + n);
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
    xt_mpi_call(XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE(
                  (int)n, bl0, d, old_type, dt), comm);
  } else {
    xt_mpi_call(XT_MPI_STRP_PRS_INDEXED_CREATE(
                  (int)n, bl, d, old_type, dt), comm);
  }

  free(d);
  return n != 0;
}

static bool
XT_MPI_STRP_PRS_MATCH_SIMPLE_VEC(
  size_t *pstart_,
  const XT_MPI_OFFSET_EXT *v,
  size_t vlen,
  MPI_Datatype old_type, MPI_Aint old_type_extent,
  MPI_Aint *disp, MPI_Datatype *dt,
  MPI_Comm comm) {
  // we only accept non-trivial matches (nsteps>2) with stride /= 1
  // using only one vector from v
  size_t p = *pstart_;
  if (p >= vlen) return false;
  int nstrides = v[p].size;
  XT_MPI_STRP_PRS_AOFS_TYPE stride = v[p].stride;
  if (nstrides < 2 || stride == XT_MPI_STRP_PRS_UNITSTRIDE ) return false;

  *pstart_ = p + 1;

  XT_MPI_STRP_PRS_AOFS_TYPE disp_ = vlen > 1 ? v[p].start : 0;
  *disp = XT_MPI_STRP_PRS_DISP_ADJUST(disp_);

  xt_mpi_call(XT_MPI_STRP_PRS_BLOCK_VEC_CREATE(
                nstrides, 1, stride, old_type, dt), comm);

  XT_MPI_STRP_PRS_AOFS_TYPE start = v[p].start - disp_;
  if (start) {
    MPI_Datatype dt1 = *dt;

    // (start != 0) => add offset:
    xt_mpi_call(MPI_Type_create_hindexed(
                  1, &(int){1},
                  &(MPI_Aint){XT_MPI_STRP_PRS_DISP_ADJUST(start)}, dt1, dt),
                comm);
    xt_mpi_call(MPI_Type_free(&dt1), comm);
  }
  return nstrides != 0;
}

static bool
XT_MPI_STRP_PRS_MATCH_CONTIGUOUS(
  size_t *pstart_,
  const XT_MPI_OFFSET_EXT *v,
  size_t vlen,
  MPI_Datatype old_type, MPI_Aint old_type_extent,
  MPI_Aint *restrict disp, MPI_Datatype *dt,
  MPI_Comm comm) {
  size_t p = *pstart_;
  if (p >= vlen || v[p].stride != XT_MPI_STRP_PRS_UNITSTRIDE || v[p].size < 2)
    return false;

  XT_MPI_STRP_PRS_AOFS_TYPE disp_ = vlen > 1 ? v[p].start : 0;
  *disp = XT_MPI_STRP_PRS_DISP_ADJUST(disp_);
  XT_MPI_STRP_PRS_AOFS_TYPE d = v[p].start - disp_;

  if (!d)
    xt_mpi_call(MPI_Type_contiguous(v[p].size, old_type, dt), comm);
  else
    xt_mpi_call(XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE(
                  1, v[p].size, &d, old_type, dt), comm);

  *pstart_ = p+1;
  return true;
}

static void
XT_MPI_STRP_PRS_GEN_FALLBACK_TYPE(
  size_t set_start, size_t set_end,
  const XT_MPI_OFFSET_EXT *v,
  size_t vlen,
  MPI_Datatype old_type,
#ifdef XT_MPI_STRP_PRS_FALLBACK_NEEDS_OLD_TYPE_EXTENT
  MPI_Aint old_type_extent,
#endif
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
  XT_MPI_STRP_PRS_AOFS_TYPE start = (ia == 0 && ib == vlen) ? 0 : v[ia].start;

  *disp = XT_MPI_STRP_PRS_DISP_ADJUST(start);

  XT_MPI_STRP_PRS_AOFS_TYPE *restrict d = xmalloc(sizeof (*d) * (size_t)n);
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
      d[p] = v[i].start + (XT_MPI_STRP_PRS_AOFS_TYPE)k * v[i].stride - start;
      p++;
    }
  }
  assert(!found_np);

  if (n==1 && d[0] == 0) {
    *dt = old_type;
  } else {
    xt_mpi_call(XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE(n, 1, d, old_type, dt),
                comm);
  }
  free(d);
}

static MPI_Datatype
XT_MPI_STRP_PRS_ENTRY(
  const XT_MPI_OFFSET_EXT *v,
  size_t vlen,
  MPI_Datatype old_type,
  MPI_Comm comm)
{
  /* [set_start,set_end) describes the prefix of non-matching
   * elements in v that then need to be handled with gen_fallback_type */
  size_t set_start = 0, set_end = 0;
  MPI_Aint old_type_lb, old_type_extent;
  xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                  &old_type_extent), comm);
  MPI_Aint *restrict wdisp
    = xmalloc(sizeof(MPI_Datatype) * vlen + sizeof (MPI_Aint) * vlen);
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
    if ((XT_MPI_STRP_PRS_MATCH_BLOCK_VEC(
           &p, v, vlen, old_type, old_type_extent,
           wdisp+mm, wdt+mm, comm))
        || (XT_MPI_STRP_PRS_MATCH_INDEXED(
              &p, v, vlen, old_type, old_type_extent,
              wdisp+mm, wdt+mm, comm))
        || (XT_MPI_STRP_PRS_MATCH_SIMPLE_VEC(
              &p, v, vlen, old_type, old_type_extent,
              wdisp+mm, wdt+mm, comm))
        || (XT_MPI_STRP_PRS_MATCH_CONTIGUOUS(
              &p, v, vlen, old_type, old_type_extent,
              wdisp+mm, wdt+mm, comm)) ) {
      /* in case a match is found, generate fallback datatype for
       * non-matching, preceding extents */
      if (set_start < set_end) {
        XT_MPI_STRP_PRS_GEN_FALLBACK_TYPE(
          set_start, set_end, v, vlen, old_type,
#ifdef XT_MPI_STRP_PRS_FALLBACK_NEEDS_OLD_TYPE_EXTENT
          old_type_extent,
#endif
          wdisp+m, wdt+m, comm);
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
    XT_MPI_STRP_PRS_GEN_FALLBACK_TYPE(
      set_start, set_end, v, vlen, old_type,
#ifdef XT_MPI_STRP_PRS_FALLBACK_NEEDS_OLD_TYPE_EXTENT
      old_type_extent,
#endif
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
XT_MPI_STRP_PRS_DRIVER(const XT_MPI_OFFSET_EXT *v,
                       int count, MPI_Datatype old_type,
                       MPI_Comm comm)
{
  size_t count_ = (size_t)0;
  for (int i=0; i<count; ++i)
    count_ += (size_t)(v[i].size > 0);
  if (count_ < 1) return MPI_DATATYPE_NULL;
  XT_MPI_OFFSET_EXT *v_comp;
  if ((size_t)count != count_) {
    v_comp = xmalloc(count_ * sizeof (*v_comp));
    for (size_t i=0, j=0; i<(size_t)count; ++i) {
      v_comp[j] = v[i];
      j+= v[i].size > 0;
    }
  } else
    v_comp = (XT_MPI_OFFSET_EXT *)v;
  MPI_Datatype dt = XT_MPI_STRP_PRS_ENTRY(v_comp, count_, old_type, comm);
  if ((size_t)count != count_)
    free(v_comp);
  return dt;
}


#undef XT_MPI_STRP_PRS_ENTRY
#undef XT_MPI_STRP_PRS_GEN_FALLBACK_TYPE
#undef XT_MPI_STRP_PRS_MATCH_CONTIGUOUS
#undef XT_MPI_STRP_PRS_MATCH_SIMPLE_VEC
#undef XT_MPI_STRP_PRS_MATCH_INDEXED
#undef XT_MPI_STRP_PRS_MATCH_BLOCK_VEC
#undef XT_MPI_OFFSET_EXT

#undef XT_TOKEN_PASTE2
#undef XT_TOKEN_PASTE2_
#undef XT_TOKEN_PASTE3
#undef XT_TOKEN_PASTE3_

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
