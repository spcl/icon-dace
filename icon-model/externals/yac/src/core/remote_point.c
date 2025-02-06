// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include "yac_mpi_common.h"
#include "yac_mpi_internal.h"
#include "remote_point.h"

MPI_Datatype yac_get_remote_point_info_mpi_datatype(MPI_Comm comm) {

  struct remote_point_info dummy;
  MPI_Datatype remote_point_info_dt;
  int array_of_blocklengths[] = {1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.rank) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.orig_pos) -
     (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] =
    {MPI_INT, MPI_UINT64_T};
  yac_mpi_call(
    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           array_of_types, &remote_point_info_dt), comm);
  return yac_create_resized(remote_point_info_dt, sizeof(dummy), comm);
}

int yac_remote_point_infos_get_pack_size(
  struct remote_point_infos const * infos, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  int pack_size_count,
      pack_size_data;

  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &pack_size_count), comm);
  yac_mpi_call(
    MPI_Pack_size(infos->count, point_info_dt, comm, &pack_size_data), comm);

  return pack_size_count + pack_size_data;
}

int yac_remote_point_get_pack_size(
  struct remote_point * point, MPI_Datatype point_info_dt, MPI_Comm comm) {

  int pack_size_id,
      pack_size_remote_point_infos;

  yac_mpi_call(MPI_Pack_size(1, yac_int_dt, comm, &pack_size_id), comm);
  pack_size_remote_point_infos =
    yac_remote_point_infos_get_pack_size(&(point->data), point_info_dt, comm);

  return pack_size_id + pack_size_remote_point_infos;
}

void yac_remote_point_infos_pack(
  struct remote_point_infos const * infos, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  int count = infos->count;

  struct remote_point_info const * info =
    (count == 1)?(&(infos->data.single)):(infos->data.multi);

  yac_mpi_call(
      MPI_Pack(&count, 1, MPI_INT, buffer, buffer_size, position, comm), comm);
  yac_mpi_call(
      MPI_Pack(info, count, point_info_dt, buffer, buffer_size, position, comm),
      comm);
}

void yac_remote_point_pack(
  struct remote_point * point, void * buffer, int buffer_size, int * position,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  yac_mpi_call(
      MPI_Pack(&(point->global_id), 1, yac_int_dt, buffer,
               buffer_size, position, comm), comm);

  yac_remote_point_infos_pack(
    &(point->data), buffer, buffer_size, position, point_info_dt, comm);
}

void yac_remote_point_infos_unpack(
  void * buffer, int buffer_size, int * position,
  struct remote_point_infos * infos, MPI_Datatype point_info_dt, MPI_Comm comm) {

  int count;
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, &count, 1, MPI_INT, comm), comm);

  infos->count = count;

  YAC_ASSERT(
    count > 0, "ERROR(yac_remote_point_infos_unpack): invalid count")

  struct remote_point_info * point_infos;
  if (count == 1)
    point_infos = &(infos->data.single);
  else
    point_infos =
      (infos->data.multi = xmalloc((size_t)count * sizeof(*point_infos)));

  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, point_infos, count,
                point_info_dt, comm), comm);
}

void yac_remote_point_unpack(
  void * buffer, int buffer_size, int * position, struct remote_point * point,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &(point->global_id), 1, yac_int_dt, comm),
    comm);
  yac_remote_point_infos_unpack(
    buffer, buffer_size, position, &(point->data), point_info_dt, comm);
}

int yac_remote_points_get_pack_size(
  struct remote_points * points, MPI_Datatype point_info_dt, MPI_Comm comm) {

  size_t count = points->count;
  struct remote_point * points_data = points->data;

  int count_pack_size,
      remote_points_pack_size;

  yac_mpi_call(MPI_Pack_size(2, MPI_UINT64_T, comm, &count_pack_size), comm);
  remote_points_pack_size = 0;
  for (size_t i = 0; i < count; ++i)
    remote_points_pack_size +=
      yac_remote_point_get_pack_size(points_data + i, point_info_dt, comm);

  return count_pack_size + remote_points_pack_size;
}

void yac_remote_points_pack(
  struct remote_points * points, void * buffer, int buffer_size, int * position,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  size_t count = points->count;
  struct remote_point * point_data = points->data;
  uint64_t counts[2] = {(uint64_t)count, 0};
  for (size_t i = 0; i < count; ++i)
    if (point_data[i].data.count > 1)
      counts[1] += (uint64_t)(point_data[i].data.count);

  yac_mpi_call(
      MPI_Pack(counts, 2, MPI_UINT64_T, buffer,
               buffer_size, position, comm), comm);
  for (size_t i = 0; i < count; ++i)
    yac_remote_point_pack(
      point_data + i, buffer, buffer_size, position, point_info_dt, comm);
}

static void yac_remote_point_infos_unpack_info_buffer(
  void * buffer, int buffer_size, int * position,
  struct remote_point_info * info_buffer, size_t * info_buffer_position,
  struct remote_point_infos * infos, MPI_Datatype point_info_dt, MPI_Comm comm) {

  int count;
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, &count, 1, MPI_INT, comm), comm);

  infos->count = count;

  struct remote_point_info * point_infos;
  if (count == 1)
    point_infos = &(infos->data.single);
  else {
    point_infos =
      (infos->data.multi = info_buffer + *info_buffer_position);
    *info_buffer_position += count;
  }

  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, point_infos, count,
               point_info_dt, comm), comm);
}

void yac_remote_point_unpack_info_buffer(
  void * buffer, int buffer_size, int * position,
  struct remote_point_info * info_buffer, size_t * info_buffer_position,
  struct remote_point * point, MPI_Datatype point_info_dt, MPI_Comm comm) {

  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &(point->global_id), 1, yac_int_dt, comm),
    comm);
  yac_remote_point_infos_unpack_info_buffer(
    buffer, buffer_size, position, info_buffer, info_buffer_position,
    &(point->data), point_info_dt, comm);
}

void yac_remote_points_unpack(
  void * buffer, int buffer_size, int * position,
  struct remote_points ** points, MPI_Datatype point_info_dt, MPI_Comm comm) {

  uint64_t counts[2];

  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, counts, 2, MPI_UINT64_T, comm), comm);

  *points = xmalloc(((size_t)counts[1]) * sizeof(struct remote_point_infos) +
                    sizeof(**points));

  size_t count = ((*points)->count = (size_t)(counts[0]));
  struct remote_point * point_data =
    ((*points)->data = xmalloc(count * sizeof(*((*points)->data))));
  struct remote_point_info * remote_point_info_buffer =
    &((*points)->buffer[0]);

  for (size_t i = 0, offset = 0; i < count; ++i) {

    yac_remote_point_unpack_info_buffer(
      buffer, buffer_size, position, remote_point_info_buffer, &offset,
      point_data + i, point_info_dt, comm);
  }
}
