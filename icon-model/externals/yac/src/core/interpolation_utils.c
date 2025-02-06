// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>

#include "utils_core.h"
#include "interpolation.h"
#include "interpolation_utils.h"
#include "yac_mpi_internal.h"

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
  for (size_t i = 0; i < num_fields; ++i)
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
