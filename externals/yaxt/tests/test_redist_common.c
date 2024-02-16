/**
 * @file test_redist_common.c
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <unistd.h>

#include "yaxt.h"
#include "tests.h"
#include "test_redist_common.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"

/*
 * build xmap for destination list containing all odd elements of
 * source list dimensioned 1 to src_num_indices
 */
Xt_xmap
build_odd_selection_xmap(int src_num_indices, MPI_Comm comm)
{
  enum {
    selection_stride = 2,
  };
  if (src_num_indices < 0)
    PUT_ERR("error: src_num_indices < 0");
  Xt_int *index_list = xmalloc((size_t)src_num_indices
                               * sizeof (index_list[0]));
  for (int i = 0; i < src_num_indices; ++i)
    index_list[i] = (Xt_int)(i + 1);
  Xt_idxlist src_idxlist = xt_idxvec_new(index_list, src_num_indices);
  int dst_num_indices
    = (int)((src_num_indices + selection_stride - 1) / selection_stride);
  for (int i = 0; i < dst_num_indices; ++i)
    index_list[i] = (Xt_int)(i * selection_stride + 1);
  Xt_idxlist dst_idxlist = xt_idxvec_new(index_list, dst_num_indices);
  free(index_list);
  Xt_xmap xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, comm);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);
  return xmap;
}

int communicators_are_congruent(MPI_Comm comm1, MPI_Comm comm2) {

  int result;

  xt_mpi_call(MPI_Comm_compare(comm1, comm2, &result), MPI_COMM_WORLD);

  return ((result == MPI_IDENT) || (result == MPI_CONGRUENT));
}

void
check_redist_(Xt_redist redist,
              int sync_mode,
              int num_redists,
              const void *src[num_redists],
              size_t dst_num_elems, void *dst[num_redists],
              void *dst_buf_base,
              prepare_dst dst_prep,
              const void *dst_prep_info,
              const void *ref_dst_data,
              MPI_Datatype dst_data_dt,
              MPI_Datatype ref_dst_data_dt,
              const char *file, int line)
{
  MPI_Comm comm = xt_redist_get_MPI_Comm(redist);
  size_t dt_extent;
  {
    MPI_Aint dt_lb, dt_extent_;
    xt_mpi_call(MPI_Type_get_extent(dst_data_dt, &dt_lb, &dt_extent_), comm);
    dt_extent = (size_t)dt_extent_;
  }
  size_t dst_size = dst_num_elems * dt_extent;
  int start_txmode = sync_mode == sync_mode_test_a ? 1 : 0,
    end_txmode = sync_mode == sync_mode_test_s ? 1 : 2;
  for (int txmode = start_txmode; txmode < end_txmode; ++txmode) {
    dst_prep(dst_buf_base, dst_prep_info, dst_num_elems);
    if (txmode == 0) {
      xt_redist_s_exchange(redist, num_redists, src, dst);
    } else {
      wrap_a_exchange(redist, num_redists, src, dst);
    }
    bool compare_failed = false;
    if (dst_data_dt == ref_dst_data_dt) {
      compare_failed = memcmp(dst_buf_base, ref_dst_data, dst_size);
    } else if (dst_data_dt == MPI_DOUBLE && ref_dst_data_dt == MPI_INT) {
      const double *dst_cmp = dst_buf_base;
      const int *ref_dst_cmp = ref_dst_data;
      for (size_t i = 0; i < dst_num_elems; ++i)
        compare_failed |= (dst_cmp[i] != ref_dst_cmp[i]);
    } else if (dst_data_dt == MPI_DOUBLE && ref_dst_data_dt == MPI_LONG) {
      const double *dst_cmp = dst_buf_base;
      const long *ref_dst_cmp = ref_dst_data;
      for (size_t i = 0; i < dst_num_elems; ++i)
        compare_failed |= (dst_cmp[i] != ref_dst_cmp[i]);
    } else if (dst_data_dt == MPI_DOUBLE && ref_dst_data_dt == MPI_SHORT) {
      const double *dst_cmp = dst_buf_base;
      const short *ref_dst_cmp = ref_dst_data;
      for (size_t i = 0; i < dst_num_elems; ++i)
        compare_failed |= (dst_cmp[i] != ref_dst_cmp[i]);
    } else if (dst_data_dt == MPI_DOUBLE && ref_dst_data_dt == MPI_LONG_LONG) {
      const double *dst_cmp = dst_buf_base;
      const long long *ref_dst_cmp = ref_dst_data;
      for (size_t i = 0; i < dst_num_elems; ++i)
        compare_failed |= (dst_cmp[i] != ref_dst_cmp[i]);
    } else if (dst_data_dt == MPI_DOUBLE
               && ref_dst_data_dt == MPI_DATATYPE_NULL) {
      const double *dst_cmp = dst_buf_base;
      for (size_t i = 0; i < dst_num_elems; ++i)
        compare_failed |= (dst_cmp[i] != (double)i);
    } else
      Xt_abort(comm, "internal error: unhandled test case!", file, line);
    if (compare_failed)
      PUT_ERR("error in xt_redist_s/a_exchange, called from %s, line %d\n",
              file, line);
  }
}

void
fill_array_double(void *dst, const void *dst_prep_info, size_t dst_num_elems)
{
  (void)dst_prep_info;
  double *restrict a = dst;
  for (size_t i = 0; i < dst_num_elems; ++i)
    a[i] = -1;
}

void
fill_array_float(void *dst, const void *dst_prep_info, size_t dst_num_elems)
{
  (void)dst_prep_info;
  float *restrict a = dst;
  for (size_t i = 0; i < dst_num_elems; ++i)
    a[i] = -1;
}

void
fill_array_long(void *dst, const void *dst_prep_info, size_t dst_num_elems)
{
  (void)dst_prep_info;
  long *restrict a = dst;
  for (size_t i = 0; i < dst_num_elems; ++i)
    a[i] = -1L;
}

void
fill_array_int(void *dst, const void *dst_prep_info, size_t dst_num_elems)
{
  (void)dst_prep_info;
  int *restrict a = dst;
  for (size_t i = 0; i < dst_num_elems; ++i)
    a[i] = -1;
}

void
fill_array_short(void *dst, const void *dst_prep_info, size_t dst_num_elems)
{
  (void)dst_prep_info;
  short *restrict a = dst;
  for (size_t i = 0; i < dst_num_elems; ++i)
    a[i] = (short)-1;
}

void
fill_array_long_long(void *dst, const void *dst_prep_info, size_t dst_num_elems)
{
  (void)dst_prep_info;
  long long *restrict a = dst;
  for (size_t i = 0; i < dst_num_elems; ++i)
    a[i] = -1;
}

void
fill_array_xt_int(void *dst, const void *dst_prep_info, size_t dst_num_elems)
{
  (void)dst_prep_info;
  Xt_int *restrict a = dst;
  for (size_t i = 0; i < dst_num_elems; ++i)
    a[i] = (Xt_int)-1;
}


void
check_wait_request_(Xt_request *request, const char *file, int line)
{
#ifndef VERBOSE
  (void)file;
  (void)line;
#endif
  if (*request == XT_REQUEST_NULL)
    PUT_ERR("request ==  XT_REQUEST_NULL before xt_request_wait: %s, line %d\n",
            file, line);
  xt_request_wait(request);
  if (*request != XT_REQUEST_NULL)
    PUT_ERR("request !=  XT_REQUEST_NULL after xt_request_wait: %s, line %d\n",
            file, line);
}

void
wrap_a_exchange(Xt_redist redist, int num_data_p, const void *src_data_p[],
                void *dst_data_p[])
{
  Xt_request request;
  xt_redist_a_exchange(redist, num_data_p, src_data_p, dst_data_p, &request);
  check_wait_request(&request);
}

void
wrap_a_exchange1(Xt_redist redist, const void *src_data_p, void *dst_data_p)
{
  Xt_request request;
  xt_redist_a_exchange1(redist, src_data_p, dst_data_p, &request);
  check_wait_request(&request);
}

Xt_config
redist_exchanger_option(int *argc, char ***argv)
{
  Xt_config config = xt_config_new();
  int opt;
  while ((opt = getopt(*argc, *argv, "m:")) != -1)
    switch (opt) {
    case 'm':
      {
        int exchanger_id = xt_exchanger_id_by_name(optarg);
        if (exchanger_id != -1)
          xt_config_set_exchange_method(config, exchanger_id);
        else {
          fprintf(stderr, "error: unexpected command-line argument for "
                  "option -m: %s\n", optarg);
          xt_config_delete(config);
          return NULL;
        }
      }
    case '?':
      break;
    }
  return config;
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
