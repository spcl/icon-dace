/**
 * @file interpolation_utils.c
 *
 * @copyright Copyright  (C)  2019 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
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

#include <string.h>

#include "utils.h"
#include "interpolation.h"
#include "interpolation_utils.h"

static size_t xt_redist_get_buffer_size(
  Xt_redist redist, MPI_Datatype(*xt_redist_get_MPI_Datatype)(Xt_redist, int)) {

  if (redist == NULL) return 0;

  MPI_Comm comm = xt_redist_get_MPI_Comm(redist);
  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t max_size = 0;

  for (int i = 0; i < comm_size; ++i) {

    MPI_Datatype dt = xt_redist_get_MPI_Datatype(redist, i);
    if (dt == MPI_DATATYPE_NULL) continue;
    MPI_Aint lb, extent;
    yac_mpi_call(MPI_Type_get_extent(dt, &lb, &extent), comm);
    size_t curr_size = (size_t)extent + (size_t)lb;
    if (curr_size > max_size) max_size = curr_size;
    yac_mpi_call(MPI_Type_free(&dt), comm);
  }
  return max_size;
}

static size_t xt_redist_get_send_buffer_size(Xt_redist redist) {

  return xt_redist_get_buffer_size(redist, xt_redist_get_send_MPI_Datatype);
}

static size_t xt_redist_get_recv_buffer_size(Xt_redist redist) {

  return xt_redist_get_buffer_size(redist, xt_redist_get_recv_MPI_Datatype);
}

static size_t * get_buffer_sizes(
  Xt_redist * redists, size_t num_fields,
  enum yac_interpolation_buffer_type type) {

  YAC_ASSERT(
    (type == SEND_BUFFER) || (type == RECV_BUFFER),
    "ERROR(get_buffer_sizes): invalid buffer type");

  size_t (*xt_redist_get_buffer_size)(Xt_redist) =
    (type == SEND_BUFFER)?
      xt_redist_get_send_buffer_size:xt_redist_get_recv_buffer_size;

  size_t * buffer_sizes;
  if (redists == NULL) {
    buffer_sizes = xcalloc(num_fields, sizeof(*buffer_sizes));
  } else {
    buffer_sizes = xmalloc(num_fields * sizeof(*buffer_sizes));
    for (size_t i = 0; i < num_fields; ++i)
      buffer_sizes[i] = xt_redist_get_buffer_size(redists[i]);
  }
  return buffer_sizes;
}

static double ** allocate_buffer(
  size_t * buffer_sizes, size_t num_fields, size_t collection_size) {

  size_t total_buffer_size = 0;
  for (int i = 0; i < num_fields; ++i)
    total_buffer_size += buffer_sizes[i];
  double ** buffer_data =
    xmalloc(collection_size * num_fields * sizeof(*buffer_data));
  buffer_data[0] = xmalloc(collection_size * total_buffer_size);
  for (size_t i = 0, offset = 0; i < collection_size; ++i) {
    for (size_t j = 0; j < num_fields; offset += buffer_sizes[j++]) {
      buffer_data[i * num_fields + j] =
        (double*)((char*)(buffer_data[0]) + offset);
    }
  }
  return buffer_data;
}

struct yac_interpolation_buffer yac_interpolation_buffer_init(
  Xt_redist * redists, size_t num_fields, size_t collection_size,
  enum yac_interpolation_buffer_type type) {

  size_t * buffer_sizes = get_buffer_sizes(redists, num_fields, type);

  return
    (struct yac_interpolation_buffer) {
      .buffer_sizes = buffer_sizes,
      .buffer = allocate_buffer(buffer_sizes, num_fields, collection_size)};
}

struct yac_interpolation_buffer yac_interpolation_buffer_init_2(
  Xt_redist * redists, size_t * min_buffer_sizes, size_t num_fields,
  size_t collection_size, enum yac_interpolation_buffer_type type) {


  size_t * buffer_sizes = get_buffer_sizes(redists, num_fields, type);
  for (size_t i = 0; i < num_fields; ++i)
    if (min_buffer_sizes[i] > buffer_sizes[i])
      buffer_sizes[i] = min_buffer_sizes[i];

  return
    (struct yac_interpolation_buffer) {
      .buffer_sizes = buffer_sizes,
      .buffer = allocate_buffer(buffer_sizes, num_fields, collection_size)};
}

struct yac_interpolation_buffer yac_interpolation_buffer_copy(
  struct yac_interpolation_buffer src, size_t num_fields,
  size_t collection_size) {

  return
    (struct yac_interpolation_buffer) {
      .buffer_sizes = COPY_DATA(src.buffer_sizes, num_fields),
      .buffer =
        allocate_buffer(src.buffer_sizes, num_fields, collection_size)};
}

void yac_interpolation_buffer_free(struct yac_interpolation_buffer * buffer) {

  free(buffer->buffer[0]);
  free(buffer->buffer);
  free(buffer->buffer_sizes);
}
