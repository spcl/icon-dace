// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef REMOTE_POINT_H
#define REMOTE_POINT_H

#include <stdint.h>
#include <mpi.h>

#include "yac_types.h"

//! single location information of a point
struct remote_point_info {
  int rank;          // MPI rank
  uint64_t orig_pos; // local id on the remote process
};

//! location information about a point that is located on one or
//  more processes
struct remote_point_infos {
  int count; // number of processes
  union {
    struct remote_point_info single; // valid if count == 1;
    struct remote_point_info * multi; // valid if count != -1;
  } data; // remote process information
};

//! information (global id and location) about a point that
//  is located on one or more processes
struct remote_point {
  yac_int global_id;
  struct remote_point_infos data;
};

//! structure containing the information (global id and location)
//  about a list of points that are located on one or more processes
struct remote_points {
  struct remote_point * data;
  size_t count;
  struct remote_point_info buffer[]; // buffer required to
                                     // store location information
};

/**
 * generates an MPI Datatype for struct remote_point_info
 * @param[in] comm communicator
 * @return MPI Datatype for struct remote_point_info
 * @remark the user has to free the returned MPI Datatype using MPI_Type_free
 */
MPI_Datatype yac_get_remote_point_info_mpi_datatype(MPI_Comm comm);

/**
 * computes the maximum size required by MPI to pack the provided point
 * information of type struct remote_point_infos
 * @param[in] infos         point information for which the pack size is
 *                          to be determined
 * @param[in] point_info_dt MPI Datatype for packing struct point_info
 * @param[in] comm          communicator
 * @return maximum packing size
 */
int yac_remote_point_infos_get_pack_size(
  struct remote_point_infos const * infos, MPI_Datatype point_info_dt,
  MPI_Comm comm);

/**
 * packs a provided remote_point_infos; this is simlar to MPI_Pack
 * @param[in]    infos         remote_point to be packed
 * @param[inout] buffer        packing buffer
 * @param[in]    buffer_size   size of packing buffer
 * @param[inout] position      packing position
 * @param[in]    point_info_dt MPI Datatype for packing struct point_info
 * @param[in]    comm          communicator
 */
void yac_remote_point_infos_pack(
  struct remote_point_infos const * infos, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * unpacks and allocates a remote_point_infos from a buffer; this is
 * similar to MPI_Unpack
 * @param[in]    buffer        packing buffer
 * @param[in]    buffer_size   size of packing buffer
 * @param[inout] position      unpacking position
 * @param[out]   infos         unpacked point information
 * @param[in]    point_info_dt MPI Datatype for unpacking struct point_info
 * @param[in]    comm          communicator
 */
void yac_remote_point_infos_unpack(
  void * buffer, int buffer_size, int * position,
  struct remote_point_infos * infos, MPI_Datatype point_info_dt,
  MPI_Comm comm);

/**
 * computes the maximum size required by MPI to pack the provided point
 * of type struct remote_point
 * @param[in] point         point for which the pack size is to be determined
 * @param[in] point_info_dt MPI Datatype for packing struct point_info
 * @param[in] comm          communicator
 * @return maximum packing size
 */
int yac_remote_point_get_pack_size(
  struct remote_point * point, MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * packs a provided remote_point; this is simlar to MPI_Pack
 * @param[in]    point         remote_point to be packed
 * @param[inout] buffer        packing buffer
 * @param[in]    buffer_size   size of packing buffer
 * @param[inout] position      packing position
 * @param[in]    point_info_dt MPI Datatype for packing struct point_info
 * @param[in]    comm          communicator
 */
void yac_remote_point_pack(
  struct remote_point * point, void * buffer, int buffer_size, int * position,
  MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * unpacks and allocates a remote_point from a buffer; this is similar to
 * MPI_Unpack
 * @param[in]    buffer        packing buffer
 * @param[in]    buffer_size   size of packing buffer
 * @param[inout] position      unpacking position
 * @param[out]   point         unpacked point
 * @param[in]    point_info_dt MPI Datatype for unpacking struct point_info
 * @param[in]    comm          communicator
 */
void yac_remote_point_unpack(
  void * buffer, int buffer_size, int * position, struct remote_point * point,
  MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * computes the maximum size required by MPI to pack the provided points
 * of type struct remote_points
 * @param[in] points        points for which the pack size is to be determined
 * @param[in] point_info_dt MPI Datatype for unpacking struct point_info
 * @param[in] comm          communicator
 * @return maximum packing size
 */
int yac_remote_points_get_pack_size(
  struct remote_points * points, MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * packs provided remote_points; this is simlar to MPI_Pack
 * @param[in]    points        remote_points to be packed
 * @param[inout] buffer        packing buffer
 * @param[in]    buffer_size   size of packing buffer
 * @param[inout] position      packing position
 * @param[in]    point_info_dt MPI Datatype for packing struct point_info
 * @param[in]    comm          communicator
 */
void yac_remote_points_pack(
  struct remote_points * points, void * buffer, int buffer_size, int * position,
  MPI_Datatype point_info_dt, MPI_Comm comm);

/**
 * unpacks and allocates remote_points from a buffer; this is similar to
 * MPI_Unpack
 * @param[in]    buffer        packing buffer
 * @param[in]    buffer_size   size of packing buffer
 * @param[inout] position      unpacking position
 * @param[out]   points        unpacked points
 * @param[in]    point_info_dt MPI Datatype for unpacking struct point_info
 * @param[in]    comm          communicator
 */
void yac_remote_points_unpack(
  void * buffer, int buffer_size, int * position,
  struct remote_points ** points, MPI_Datatype point_info_dt, MPI_Comm comm);

#endif // REMOTE_POINT_H
