// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "read_icon_grid.h"
#include "utils_common.h"
#include "io_utils.h"
#include "geometry.h"
#include "read_grid.h"

#ifdef YAC_NETCDF_ENABLED
#include <netcdf.h>

static int * get_icon_cell_mask( int ncid, size_t nbr_cells );
static void get_icon_connect(
  int ncid, size_t nbr_cells, int ** vertex_of_cell, size_t * nv);

void yac_read_part_icon_grid_information(
  const char * filename, int * nbr_vertices, int * nbr_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells, int ** global_cell_id,
  int ** cell_mask, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size) {

 // read the global grid
  yac_read_icon_grid_information(
    filename, nbr_vertices, nbr_cells, num_vertices_per_cell, cell_to_vertex,
    x_vertices, y_vertices, x_cells, y_cells, cell_mask);

  // determine local range
  int local_start =
    ((unsigned long)(*nbr_cells) * (unsigned long)rank) /
    (unsigned long)size;
  int num_local_cells =
    ((unsigned long)(*nbr_cells) * ((unsigned long)rank+1)) /
    (unsigned long)size - (unsigned long)local_start;

  // mask for required vertices and cells
  int * required_vertices = xcalloc(*nbr_vertices, sizeof(*required_vertices));
  int * required_cells = xcalloc(*nbr_cells, sizeof(*required_cells));

  int offset = 0;
  for (int i = 0; i < local_start; ++i)
    offset += (*num_vertices_per_cell)[i];

  // mark all local cells and their vertices as required
  for (int i = local_start; i < local_start + num_local_cells; ++i) {

    required_cells[i] = 2;
    for (int j = 0; j < (*num_vertices_per_cell)[i]; ++j)
      required_vertices[(*cell_to_vertex)[offset+j]] = 2;

    offset += (*num_vertices_per_cell)[i];
  }

  // mark all halo cells as required
  offset = 0;
  for (int i = 0; i < *nbr_cells; ++i) {

    if (!required_cells[i]) {

      for (int j = 0; j < (*num_vertices_per_cell)[i]; ++j) {

        if (required_vertices[(*cell_to_vertex)[offset+j]]) {

          required_cells[i] = 1;
          break;
        }
      }

    }

    offset += (*num_vertices_per_cell)[i];
  }

  // mark all halo vertices as required
  offset = 0;
  for (int i = 0; i < *nbr_cells; ++i) {

    if (required_cells[i] == 1) {

      for (int j = 0; j < (*num_vertices_per_cell)[i]; ++j) {

        if (!required_vertices[(*cell_to_vertex)[offset+j]]) {

          required_vertices[(*cell_to_vertex)[offset+j]] = 1;
        }
      }
    }

    offset += (*num_vertices_per_cell)[i];
  }

  // count the number of cells and vertices
  int part_num_vertices = 0;
  int part_num_cells = 0;
  for (int i = 0; i < *nbr_vertices; ++i)
    if (required_vertices[i])
      part_num_vertices++;
  for (int i = 0; i < *nbr_cells; ++i)
    if(required_cells[i])
      part_num_cells++;

  *global_cell_id = xmalloc(part_num_cells * sizeof(**global_cell_id));
  *cell_core_mask = xmalloc(part_num_cells * sizeof(**cell_core_mask));
  *global_corner_id = xmalloc(part_num_vertices * sizeof(**global_corner_id));
  *corner_core_mask = xmalloc(part_num_vertices * sizeof(**corner_core_mask));

  // generate final vertex data
  part_num_vertices = 0;
  int * global_to_local_vertex =
    xmalloc(*nbr_vertices * sizeof(*global_to_local_vertex));
  for (int i = 0; i < *nbr_vertices; ++i) {

    if (required_vertices[i]) {

      (*global_corner_id)[part_num_vertices] = i;
      (*corner_core_mask)[part_num_vertices] = required_vertices[i] == 2;
      (*x_vertices)[part_num_vertices] = (*x_vertices)[i];
      (*y_vertices)[part_num_vertices] = (*y_vertices)[i];
      global_to_local_vertex[i] = part_num_vertices;
      part_num_vertices++;
    }
  }

  *x_vertices = xrealloc(*x_vertices, part_num_vertices * sizeof(**x_vertices));
  *y_vertices = xrealloc(*y_vertices, part_num_vertices * sizeof(**y_vertices));
  *nbr_vertices = part_num_vertices;
  free(required_vertices);

  // generate final cell data
  int num_cell_vertex_dependencies = 0;
  part_num_cells = 0;
  offset = 0;
  for (int i = 0; i < *nbr_cells; ++i) {

    if (required_cells[i]) {

      (*global_cell_id)[part_num_cells] = i;
      (*cell_core_mask)[part_num_cells] = required_cells[i] == 2;
      (*x_cells)[part_num_cells] = (*x_cells)[i];
      (*y_cells)[part_num_cells] = (*y_cells)[i];
      (*cell_mask)[part_num_cells] = (*cell_mask)[i];

      for (int j = 0; j < (*num_vertices_per_cell)[i]; ++j)
        (*cell_to_vertex)[num_cell_vertex_dependencies++] =
          global_to_local_vertex[(*cell_to_vertex)[offset+j]];

      (*num_vertices_per_cell)[part_num_cells] = (*num_vertices_per_cell)[i];

      part_num_cells++;
    }

    offset += (*num_vertices_per_cell)[i];
  }

  *x_cells   = xrealloc(*x_cells, part_num_cells * sizeof(**x_cells));
  *y_cells   = xrealloc(*y_cells, part_num_cells * sizeof(**y_cells));
  *cell_mask = xrealloc(*cell_mask, part_num_cells * sizeof(**cell_mask));

  *num_vertices_per_cell = xrealloc(*num_vertices_per_cell, part_num_cells *
                                    sizeof(**num_vertices_per_cell));
  *cell_to_vertex = xrealloc(*cell_to_vertex, num_cell_vertex_dependencies *
                             sizeof(**cell_to_vertex));
  *nbr_cells = part_num_cells;
  free(required_cells);
  free(global_to_local_vertex);
}

void yac_read_icon_grid_information(const char * filename, int * nbr_vertices,
                                    int * nbr_cells, int ** num_vertices_per_cell,
                                    int ** cell_to_vertex, double ** x_vertices,
                                    double ** y_vertices, double ** x_cells,
                                    double ** y_cells, int ** cell_mask) {

   /* Open file */
   int ncid;
   yac_nc_open(filename, NC_NOWRITE, &ncid);

   /* Get vertex longitudes and latitudes of cells */
   size_t nbr_vertices_;
   yac_read_coords(
     ncid, "vlon", "vlat", x_vertices, y_vertices, &nbr_vertices_);
   *nbr_vertices = (int)nbr_vertices_;

   /* Get cell center longitudes and latitudes of cells */
   size_t nbr_cells_;
   yac_read_coords(ncid, "clon", "clat", x_cells, y_cells, &nbr_cells_);
   *nbr_cells = (int)nbr_cells_;

   /* Get mask of cells */
   *cell_mask = get_icon_cell_mask ( ncid, *nbr_cells );

   /* Get relations between vertices and cells */
   int * vertex_of_cell;
   size_t nv;
   get_icon_connect(ncid, *nbr_cells, &vertex_of_cell, &nv);

   /* Close file */
   YAC_HANDLE_ERROR(nc_close(ncid));

   //-------------------------------------------------------------------------//

   *num_vertices_per_cell = xmalloc(*nbr_cells * sizeof(**num_vertices_per_cell));
   for (int i = 0; i < *nbr_cells; (*num_vertices_per_cell)[i++] = (int)nv);

   *cell_to_vertex = xmalloc(*nbr_cells * nv * sizeof(**cell_to_vertex));

   // Unfortunately the data is only available in Fortran order
   for (int i = 0; i < *nbr_cells; ++i) {

     for (size_t j = 0; j < nv; ++j)
       (*cell_to_vertex)[nv * i + j] =
         vertex_of_cell[i + j * (*nbr_cells)] - 1;
   }

   free(vertex_of_cell);
}

// taken from scales-ppm library
// https://www.dkrz.de/redmine/projects/scales-ppm
static inline int
partition_idx_from_element_idx(unsigned element_idx, unsigned num_elements,
                               int num_partitions) {

  return (int)((((unsigned long)element_idx) * ((unsigned long)num_partitions) +
                (unsigned long)num_partitions - 1) /
               ((unsigned long)num_elements));
}

static void convert_to_rad(int ncid, int varid, double * array, size_t count) {

  if (yac_check_coord_units(ncid, varid))
    for (size_t i = 0; i < count; ++i) array[i] *= YAC_RAD;
}

void yac_read_icon_grid_information_parallel(
  const char * filename, MPI_Comm comm, int * nbr_vertices, int * nbr_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex, int ** global_cell_ids,
  int ** cell_owner, int ** global_vertex_ids, int ** vertex_owner,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells, int ** cell_msk) {

  int comm_rank, comm_size;

  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  int local_is_io, * io_ranks, num_io_ranks;
  yac_get_io_ranks(comm, &local_is_io, &io_ranks, &num_io_ranks);

  size_t ncells, nvertices, nv, ne;

  size_t read_local_start_cell = 0;
  size_t read_num_local_cells = 0;
  size_t read_local_start_vertex = 0;
  size_t read_num_local_vertices = 0;
  double * read_cell_lon = NULL;
  double * read_cell_lat = NULL;
  int * read_cell_mask = NULL;
  double * read_vertex_lon = NULL;
  double * read_vertex_lat = NULL;
  int * read_dist_vertex_of_cell = NULL;
  int * read_dist_cells_of_vertex = NULL;

  if (local_is_io) {

    unsigned long io_proc_idx = ULONG_MAX;
    for (int i = 0; (i < num_io_ranks) && (io_proc_idx == ULONG_MAX); ++i)
      if (io_ranks[i] == comm_rank)
        io_proc_idx = (unsigned long)i;

    // open file
    int ncid;
    yac_nc_open(filename, NC_NOWRITE, &ncid);

    // get number of cells and vertices
    int dim_id;
    yac_nc_inq_dimid(ncid, "cell", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &ncells));
    yac_nc_inq_dimid(ncid, "vertex", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &nvertices));
    yac_nc_inq_dimid(ncid, "nv", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &nv));
    yac_nc_inq_dimid(ncid, "ne", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &ne));

    // determine local range for cell and vertex data
    read_local_start_cell =
      ((unsigned long)ncells * io_proc_idx) / (unsigned long)num_io_ranks;
    read_num_local_cells =
      ((unsigned long)ncells * (io_proc_idx+1)) / (unsigned long)num_io_ranks -
      (unsigned long)read_local_start_cell;
    read_local_start_vertex =
      ((unsigned long)nvertices * io_proc_idx) / (unsigned long)num_io_ranks;
    read_num_local_vertices =
      ((unsigned long)nvertices * (io_proc_idx+1)) / (unsigned long)num_io_ranks -
      (unsigned long)read_local_start_vertex;

    // read basic grid data (each process its individual part)
    read_cell_lon = xmalloc(read_num_local_cells * sizeof(*read_cell_lon));
    read_cell_lat = xmalloc(read_num_local_cells * sizeof(*read_cell_lat));
    read_cell_mask = xmalloc(read_num_local_cells * sizeof(*read_cell_mask));
    read_vertex_lon = xmalloc(read_num_local_vertices * sizeof(*read_vertex_lon));
    read_vertex_lat = xmalloc(read_num_local_vertices * sizeof(*read_vertex_lat));
    read_dist_vertex_of_cell =
      xmalloc(read_num_local_cells * nv * sizeof(*read_dist_vertex_of_cell));
    read_dist_cells_of_vertex =
      xmalloc(read_num_local_vertices * ne * sizeof(*read_dist_cells_of_vertex));
    int varid;
    yac_nc_inq_varid(ncid, "clon", &varid);
    YAC_HANDLE_ERROR(nc_get_vara_double(ncid, varid, &read_local_start_cell,
                                    &read_num_local_cells, read_cell_lon));
    convert_to_rad(ncid, varid, read_cell_lon, read_num_local_cells);
    yac_nc_inq_varid(ncid, "clat", &varid);
    YAC_HANDLE_ERROR(nc_get_vara_double(ncid, varid, &read_local_start_cell,
                                    &read_num_local_cells, read_cell_lat));
    convert_to_rad(ncid, varid, read_cell_lat, read_num_local_cells);
    yac_nc_inq_varid(ncid, "cell_sea_land_mask", &varid);
    YAC_HANDLE_ERROR(nc_get_vara_int(ncid, varid, &read_local_start_cell,
                                 &read_num_local_cells, read_cell_mask));
    yac_nc_inq_varid(ncid, "vlon", &varid);
    YAC_HANDLE_ERROR(nc_get_vara_double(ncid, varid, &read_local_start_vertex,
                                    &read_num_local_vertices, read_vertex_lon));
    convert_to_rad(ncid, varid, read_vertex_lon, read_num_local_vertices);
    yac_nc_inq_varid(ncid, "vlat", &varid);
    YAC_HANDLE_ERROR(nc_get_vara_double(ncid, varid, &read_local_start_vertex,
                                    &read_num_local_vertices, read_vertex_lat));
    convert_to_rad(ncid, varid, read_vertex_lat, read_num_local_vertices);
    {
      int * buffer = xmalloc(MAX(nv*read_num_local_cells,
                                 ne*read_num_local_vertices) * sizeof(*buffer));
      {
        size_t tmp_start[2] = {0, read_local_start_cell};
        size_t tmp_count[2] = {nv, read_num_local_cells};
        yac_nc_inq_varid(ncid, "vertex_of_cell", &varid);
        YAC_HANDLE_ERROR(nc_get_vara_int(ncid, varid, tmp_start, tmp_count, buffer));
        for (size_t i = 0; i < read_num_local_cells; ++i)
          for (size_t j = 0; j < nv; ++j)
            read_dist_vertex_of_cell[i * nv + j] =
              buffer[i + j * read_num_local_cells];
      }
      {
        size_t tmp_start[2] = {0, read_local_start_vertex};
        size_t tmp_count[2] = {ne, read_num_local_vertices};
        yac_nc_inq_varid(ncid, "cells_of_vertex", &varid);
        YAC_HANDLE_ERROR(nc_get_vara_int(ncid, varid, tmp_start, tmp_count, buffer));
        for (size_t i = 0; i < read_num_local_vertices; ++i)
          for (size_t j = 0; j < ne; ++j)
            read_dist_cells_of_vertex[i * ne + j] =
              buffer[i + j * read_num_local_vertices];
      }
      free(buffer);

      // adjust for c indices
      for (size_t i = 0; i < read_num_local_cells * nv; ++i)
        if (read_dist_vertex_of_cell[i] > 0) read_dist_vertex_of_cell[i]--;
      for (size_t i = 0; i < read_num_local_vertices * ne; ++i)
        read_dist_cells_of_vertex[i]--;
    }

    YAC_HANDLE_ERROR(nc_close(ncid));
  } else {
    read_cell_lon = xmalloc(1 * sizeof(*read_cell_lon));
    read_cell_lat = xmalloc(1 * sizeof(*read_cell_lat));
    read_cell_mask = xmalloc(1 * sizeof(*read_cell_mask));
    read_vertex_lon = xmalloc(1 * sizeof(*read_vertex_lon));
    read_vertex_lat = xmalloc(1 * sizeof(*read_vertex_lat));
    read_dist_vertex_of_cell = xmalloc(1 * sizeof(*read_dist_vertex_of_cell));
    read_dist_cells_of_vertex = xmalloc(1 * sizeof(*read_dist_cells_of_vertex));
  }

  free(io_ranks);

  {
    int tmp;
    if (comm_rank == 0) tmp = (int)ncells;
    MPI_Bcast(&tmp, 1, MPI_INT, 0, comm);
    ncells = (size_t)tmp;
    if (comm_rank == 0) tmp = (int)nvertices;
    MPI_Bcast(&tmp, 1, MPI_INT, 0, comm);
    nvertices = (size_t)tmp;
    if (comm_rank == 0) tmp = (int)nv;
    MPI_Bcast(&tmp, 1, MPI_INT, 0, comm);
    nv = (size_t)tmp;
    if (comm_rank == 0) tmp = (int)ne;
    MPI_Bcast(&tmp, 1, MPI_INT, 0, comm);
    ne = (size_t)tmp;
  }

  // determine local range for cell and vertex data
  size_t local_start_cell =
    ((unsigned long)ncells * (unsigned long)comm_rank) /
    (unsigned long)comm_size;
  size_t num_local_cells =
    ((unsigned long)ncells * ((unsigned long)comm_rank+1)) /
    (unsigned long)comm_size - (unsigned long)local_start_cell;
  size_t local_start_vertex =
    ((unsigned long)nvertices * (unsigned long)comm_rank) /
    (unsigned long)comm_size;
  size_t num_local_vertices =
    ((unsigned long)nvertices * ((unsigned long)comm_rank+1)) /
    (unsigned long)comm_size - (unsigned long)local_start_vertex;

  // redistribute basic cell data (from io decomposition)
  double * cell_lon = xmalloc(num_local_cells * sizeof(*cell_lon));
  double * cell_lat = xmalloc(num_local_cells * sizeof(*cell_lat));
  int * cell_mask = xmalloc(num_local_cells * sizeof(*cell_mask));
  int * dist_vertex_of_cell =
    xmalloc(num_local_cells * nv * sizeof(*dist_vertex_of_cell));
  {
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (size_t i = 0; i < read_num_local_cells; ++i)
      send_count[
        partition_idx_from_element_idx(
          read_local_start_cell + i, ncells, comm_size)]++;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    MPI_Alltoallv(read_cell_lon, send_count, send_displ, MPI_DOUBLE,
                  cell_lon, recv_count, recv_displ, MPI_DOUBLE, comm);
    MPI_Alltoallv(read_cell_lat, send_count, send_displ, MPI_DOUBLE,
                  cell_lat, recv_count, recv_displ, MPI_DOUBLE, comm);
    MPI_Alltoallv(read_cell_mask, send_count, send_displ, MPI_INT,
                  cell_mask, recv_count, recv_displ, MPI_INT, comm);

    for (int i = 0; i < comm_size; ++i) {
      send_count[i] *= nv;
      send_displ[i] *= nv;
      recv_count[i] *= nv;
      recv_displ[i] *= nv;
    }

    MPI_Alltoallv(read_dist_vertex_of_cell, send_count, send_displ, MPI_INT,
                  dist_vertex_of_cell, recv_count, recv_displ, MPI_INT, comm);

    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
    free(read_cell_lon);
    free(read_cell_lat);
    free(read_cell_mask);
    free(read_dist_vertex_of_cell);
  }

  // redistribute basic vertex data (from io decomposition)
  double * vertex_lon = xmalloc(num_local_vertices * sizeof(*vertex_lon));
  double * vertex_lat = xmalloc(num_local_vertices * sizeof(*vertex_lat));
  int * dist_cells_of_vertex =
    xmalloc(num_local_vertices * ne * sizeof(*dist_cells_of_vertex));
  {
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (size_t i = 0; i < read_num_local_vertices; ++i)
      send_count[
        partition_idx_from_element_idx(
          read_local_start_vertex + i, nvertices, comm_size)]++;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    MPI_Alltoallv(read_vertex_lon, send_count, send_displ, MPI_DOUBLE,
                  vertex_lon, recv_count, recv_displ, MPI_DOUBLE, comm);
    MPI_Alltoallv(read_vertex_lat, send_count, send_displ, MPI_DOUBLE,
                  vertex_lat, recv_count, recv_displ, MPI_DOUBLE, comm);

    for (int i = 0; i < comm_size; ++i) {
      send_count[i] *= ne;
      send_displ[i] *= ne;
      recv_count[i] *= ne;
      recv_displ[i] *= ne;
    }

    MPI_Alltoallv(read_dist_cells_of_vertex, send_count, send_displ, MPI_INT,
                  dist_cells_of_vertex, recv_count, recv_displ, MPI_INT, comm);

    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
    free(read_vertex_lon);
    free(read_vertex_lat);
    free(read_dist_cells_of_vertex);
  }

  // determine required vertices for core cells
  size_t num_core_vertices = num_local_cells * nv;
  int * core_vertices = xmalloc(num_core_vertices * sizeof(*core_vertices));
  {
    memcpy(core_vertices, dist_vertex_of_cell,
           num_core_vertices * sizeof(*core_vertices));
    yac_quicksort_index(core_vertices, num_core_vertices, NULL);
    yac_remove_duplicates_int(core_vertices, &num_core_vertices);
    core_vertices =
      xrealloc(core_vertices, num_core_vertices * sizeof(*core_vertices));
  }

  // get cells_of_vertex for core vertices and compute dist_vertex_owner
  int * cells_of_vertex_core =
    xmalloc(num_core_vertices * ne * sizeof(*cells_of_vertex_core));
  int * dist_vertex_owner =
    xmalloc(num_local_vertices * sizeof(*dist_vertex_owner));
  {
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (size_t i = 0; i < num_core_vertices; ++i)
      send_count[((unsigned long)(core_vertices[i]) *
                  (unsigned long)comm_size + (unsigned long)comm_size - 1) /
                 (unsigned long)nvertices]++;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    int num_core_vertices_remote = 0;
    for (int i = 0; i < comm_size; ++i)
      num_core_vertices_remote += recv_count[i];

    int * remote_vertex_buffer =
      xmalloc(num_core_vertices_remote * sizeof(*remote_vertex_buffer));

    MPI_Alltoallv(core_vertices, send_count, send_displ, MPI_INT,
                  remote_vertex_buffer, recv_count, recv_displ, MPI_INT, comm);

    for (size_t i = 0; i < num_local_vertices; ++i) dist_vertex_owner[i] = -1;

    for (int i = 0, j = 0; i < comm_size; ++i)
      for (int k = 0; k < recv_count[i]; ++k, ++j)
        dist_vertex_owner[remote_vertex_buffer[j] - local_start_vertex] = i;

    int * send_cell_of_vertex =
      xmalloc(num_core_vertices_remote * ne * sizeof(*send_cell_of_vertex));

    for (int i = 0, l = 0, m = 0; i < comm_size; ++i) {
      for (int j = 0; j < recv_count[i]; ++j, ++l) {
        int idx = remote_vertex_buffer[l] - local_start_vertex;
        for (size_t k = 0; k < ne; ++k, ++m)
          send_cell_of_vertex[m] = dist_cells_of_vertex[idx * ne + k];
      }
    }
    free(dist_cells_of_vertex);

    for (int i = 0; i < comm_size; ++i) {
      send_count[i] *= ne;
      send_displ[i] *= ne;
      recv_count[i] *= ne;
      recv_displ[i] *= ne;
    }

    MPI_Alltoallv(send_cell_of_vertex, recv_count, recv_displ, MPI_INT,
                  cells_of_vertex_core, send_count, send_displ, MPI_INT, comm);

    free(send_cell_of_vertex);
    free(remote_vertex_buffer);
    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
  }

  // determine halo cells
  int * halo_cells;
  size_t num_halo_cells = 0;
  {
    int * tmp_cells = cells_of_vertex_core;
    size_t num_tmp_cells = num_core_vertices * ne;

    yac_quicksort_index(tmp_cells, num_tmp_cells, NULL);
    yac_remove_duplicates_int(tmp_cells, &num_tmp_cells);

    if ((num_tmp_cells > 0) && (tmp_cells[0] == -1)) {
      num_tmp_cells--;
      tmp_cells++;
    }

    halo_cells = xmalloc((num_tmp_cells - num_local_cells) * sizeof(*halo_cells));

    size_t i = 0;
    for (; i < num_tmp_cells && tmp_cells[i] < (int)local_start_cell; ++i)
      halo_cells[num_halo_cells++] = tmp_cells[i];
    i += num_local_cells;
    for (; i < num_tmp_cells; ++i) halo_cells[num_halo_cells++] = tmp_cells[i];

    assert(num_halo_cells == num_tmp_cells - num_local_cells);

    free(cells_of_vertex_core);
  }

  // determine all vertices and get coordinates of halo cells
  size_t num_all_local_vertices = num_halo_cells * nv + num_core_vertices;
  int * all_cell_to_vertex;
  int * all_local_vertices = xrealloc(core_vertices, num_all_local_vertices *
                                      sizeof(*all_local_vertices));
  cell_lon = xrealloc(cell_lon, (num_local_cells + num_halo_cells) *
                      sizeof(*cell_lon));
  cell_lat = xrealloc(cell_lat, (num_local_cells + num_halo_cells) *
                      sizeof(*cell_lat));
  cell_mask = xrealloc(cell_mask, (num_local_cells + num_halo_cells) *
                       sizeof(*cell_mask));
  {
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (size_t i = 0; i < num_halo_cells; ++i)
      send_count[((unsigned long)(halo_cells[i]) *
                  (unsigned long)comm_size + (unsigned long)comm_size - 1) /
                 (unsigned long)ncells]++;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    int num_halo_cells_remote = 0;
    for (int i = 0; i < comm_size; ++i)
      num_halo_cells_remote += recv_count[i];

    int * remote_halo_cell_buffer =
      xmalloc(num_halo_cells_remote * sizeof(*remote_halo_cell_buffer));

    MPI_Alltoallv(halo_cells, send_count, send_displ, MPI_INT,
                  remote_halo_cell_buffer, recv_count, recv_displ, MPI_INT,
                  comm);

    int * send_halo_cell_vertices =
      xmalloc(num_halo_cells_remote * nv * sizeof(*send_halo_cell_vertices));
    double * send_cell_lon =
      xmalloc(num_halo_cells_remote * sizeof(*send_cell_lon));
    double * send_cell_lat =
      xmalloc(num_halo_cells_remote * sizeof(*send_cell_lat));
    int * send_cell_mask =
      xmalloc(num_halo_cells_remote * sizeof(*send_cell_mask));

    for (int i = 0, l = 0, m = 0; i < comm_size; ++i) {
      for (int j = 0; j < recv_count[i]; ++j, ++l) {
        int idx = remote_halo_cell_buffer[l] - local_start_cell;
        for (size_t k = 0; k < nv; ++k, ++m)
          send_halo_cell_vertices[m] = dist_vertex_of_cell[idx * nv + k];
        send_cell_lon[l] = cell_lon[idx];
        send_cell_lat[l] = cell_lat[idx];
        send_cell_mask[l] = cell_mask[idx];
      }
    }

    MPI_Alltoallv(send_cell_lon, recv_count, recv_displ, MPI_DOUBLE,
                  cell_lon + num_local_cells, send_count, send_displ,
                  MPI_DOUBLE, comm);
    MPI_Alltoallv(send_cell_lat, recv_count, recv_displ, MPI_DOUBLE,
                  cell_lat + num_local_cells, send_count, send_displ,
                  MPI_DOUBLE, comm);
    MPI_Alltoallv(send_cell_mask, recv_count, recv_displ, MPI_INT,
                  cell_mask + num_local_cells, send_count, send_displ,
                  MPI_INT, comm);

    for (int i = 0; i < comm_size; ++i) {
      send_count[i] *= nv;
      send_displ[i] *= nv;
      recv_count[i] *= nv;
      recv_displ[i] *= nv;
    }

    all_cell_to_vertex =
      xrealloc(dist_vertex_of_cell, (num_local_cells + num_halo_cells) * nv *
               sizeof(*all_cell_to_vertex));

    MPI_Alltoallv(send_halo_cell_vertices, recv_count, recv_displ, MPI_INT,
                  all_cell_to_vertex + num_local_cells * nv, send_count,
                  send_displ, MPI_INT, comm);

    memcpy(all_local_vertices + num_core_vertices,
           all_cell_to_vertex + num_local_cells * nv,
           num_halo_cells * nv * sizeof(*all_local_vertices));

    free(send_cell_mask);
    free(send_cell_lat);
    free(send_cell_lon);
    free(send_halo_cell_vertices);
    free(remote_halo_cell_buffer);
    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);

    yac_quicksort_index(
      all_local_vertices, num_all_local_vertices, NULL);
    yac_remove_duplicates_int(all_local_vertices, &num_all_local_vertices);
  }

  // determine owner and coordinates for all vertices
  double * all_vertex_lon =
    xmalloc(num_all_local_vertices * sizeof(*all_vertex_lon));
  double * all_vertex_lat =
    xmalloc(num_all_local_vertices * sizeof(*all_vertex_lat));
  int * all_local_vertices_owner =
    xmalloc(num_all_local_vertices * sizeof(*all_local_vertices_owner));
  {
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (size_t i = 0; i < num_all_local_vertices; ++i)
      send_count[((unsigned long)(all_local_vertices[i]) *
                  (unsigned long)comm_size + (unsigned long)comm_size - 1) /
                 (unsigned long)nvertices]++;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    int num_all_local_vertices_remote = 0;
    for (int i = 0; i < comm_size; ++i)
      num_all_local_vertices_remote += recv_count[i];

    int * remote_vertex_buffer =
      xmalloc(num_all_local_vertices_remote * sizeof(*remote_vertex_buffer));

    MPI_Alltoallv(all_local_vertices, send_count, send_displ, MPI_INT,
                  remote_vertex_buffer, recv_count, recv_displ, MPI_INT, comm);

    int * send_vertex_owner = remote_vertex_buffer;
    double * send_vertex_lon =
      xmalloc(num_all_local_vertices_remote * sizeof(*send_vertex_lon));
    double * send_vertex_lat =
      xmalloc(num_all_local_vertices_remote * sizeof(*send_vertex_lat));

    for (int i = 0, l = 0; i < comm_size; ++i) {
      for (int j = 0; j < recv_count[i]; ++j, ++l) {
        int idx = remote_vertex_buffer[l] - local_start_vertex;
        send_vertex_owner[l] = dist_vertex_owner[idx];
        send_vertex_lon[l] = vertex_lon[idx];
        send_vertex_lat[l] = vertex_lat[idx];
      }
    }

    free(dist_vertex_owner);
    free(vertex_lon);
    free(vertex_lat);

    MPI_Alltoallv(send_vertex_owner, recv_count, recv_displ, MPI_INT,
                  all_local_vertices_owner, send_count, send_displ, MPI_INT,
                  comm);
    MPI_Alltoallv(send_vertex_lon, recv_count, recv_displ, MPI_DOUBLE,
                  all_vertex_lon, send_count, send_displ, MPI_DOUBLE, comm);
    MPI_Alltoallv(send_vertex_lat, recv_count, recv_displ, MPI_DOUBLE,
                  all_vertex_lat, send_count, send_displ, MPI_DOUBLE, comm);

    free(send_vertex_lat);
    free(send_vertex_lon);
    free(send_vertex_owner);
    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
  }

  // convert global ids within all_cell_to_vertex into local ids
  {
    size_t vertex_count = nv * (num_local_cells + num_halo_cells);
    int * permutation = xmalloc(vertex_count * sizeof(*permutation));
    for (size_t i = 0; i < vertex_count; ++i) permutation[i] = (int)i;

    yac_quicksort_index(all_cell_to_vertex, vertex_count, permutation);

    for (size_t i = 0, j = 0; i < vertex_count; ++i) {
      while (all_local_vertices[j] != all_cell_to_vertex[i]) ++j;
      all_cell_to_vertex[i] = (int)j;
    }

    yac_quicksort_index(permutation, vertex_count, all_cell_to_vertex);
    free(permutation);
  }

  *nbr_vertices = num_all_local_vertices;
  *nbr_cells = num_local_cells + num_halo_cells;

  *num_vertices_per_cell =
    xmalloc((num_local_cells + num_halo_cells) * sizeof(**num_vertices_per_cell));
  for (size_t i = 0; i < num_local_cells + num_halo_cells; ++i)
    (*num_vertices_per_cell)[i] = nv;

  *cell_to_vertex = all_cell_to_vertex;

  *global_cell_ids =
    xmalloc((num_local_cells + num_halo_cells) * sizeof(**global_cell_ids));
  for (size_t i = 0; i < num_local_cells; ++i)
    (*global_cell_ids)[i] = local_start_cell + i;
  memcpy(*global_cell_ids + num_local_cells, halo_cells,
         num_halo_cells * sizeof(*halo_cells));
  *global_vertex_ids = xrealloc(all_local_vertices, num_all_local_vertices *
                                sizeof(*all_local_vertices));

  *cell_owner =
    xmalloc((num_local_cells + num_halo_cells) * sizeof(**cell_owner));
  for (size_t i = 0; i < num_local_cells; ++i)
    (*cell_owner)[i] = -1;
  for (size_t i = 0; i < num_halo_cells; ++i)
    (*cell_owner)[num_local_cells + i] =
      ((unsigned long)(halo_cells[i]) * (unsigned long)comm_size +
      (unsigned long)comm_size - 1) / (unsigned long)ncells;
  free(halo_cells);

  for (size_t i = 0; i < num_all_local_vertices; ++i)
    if (all_local_vertices_owner[i] == comm_rank)
      all_local_vertices_owner[i] = -1;
  *vertex_owner = all_local_vertices_owner;

  *x_vertices = all_vertex_lon;
  *y_vertices = all_vertex_lat;
  *x_cells = cell_lon;
  *y_cells = cell_lat;
  *cell_msk = cell_mask;
}

void yac_read_icon_grid_information_parallel_f2c(
  const char * filename, MPI_Fint comm, int * nbr_vertices, int * nbr_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex, int ** global_cell_ids,
  int ** cell_owner, int ** global_vertex_ids, int ** vertex_owner,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells, int ** cell_msk) {

  yac_read_icon_grid_information_parallel(
    filename, MPI_Comm_f2c(comm), nbr_vertices, nbr_cells,
    num_vertices_per_cell, cell_to_vertex, global_cell_ids,
    cell_owner, global_vertex_ids, vertex_owner,
    x_vertices, y_vertices, x_cells, y_cells, cell_msk);
}

static int * generate_simple_core_mask(size_t N) {
  int * mask = xmalloc(N * sizeof(*mask));
  for (size_t i = 0; i < N; ++i) mask[i] = 1;
  return mask;
}

static size_t * generate_offsets(size_t N, int * counts) {

  size_t * offsets = xmalloc(N * sizeof(*offsets));
  for (size_t i = 0, accu = 0; i < N; ++i) {
    offsets[i] = accu;
    accu += (size_t)(counts[i]);
  }
  return offsets;
}

struct yac_basic_grid_data yac_read_icon_basic_grid_data_parallel(
  const char * filename, MPI_Comm comm) {

  int comm_rank, comm_size;

  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  int local_is_io, * io_ranks, num_io_ranks;
  yac_get_io_ranks(comm, &local_is_io, &io_ranks, &num_io_ranks);

  size_t ncells, nvertices, nedges, nv, ne;

  size_t read_local_start_cell = 0;
  size_t read_num_local_cells = 0;
  size_t read_local_start_vertex = 0;
  size_t read_num_local_vertices = 0;
  size_t read_local_start_edge = 0;
  size_t read_num_local_edges = 0;
  double * read_vertex_lon = NULL;
  double * read_vertex_lat = NULL;
  int * read_dist_vertex_of_cell = NULL;
  int * read_dist_edge_of_cell = NULL;
  int * read_dist_edge_vertices = NULL;

  if (local_is_io) {

    unsigned long io_proc_idx = ULONG_MAX;
    for (int i = 0; (i < num_io_ranks) && (io_proc_idx == ULONG_MAX); ++i)
      if (io_ranks[i] == comm_rank)
        io_proc_idx = (unsigned long)i;

    // open file
    int ncid;
    yac_nc_open(filename, NC_NOWRITE, &ncid);

    // get number of cells and vertices
    int dim_id;
    yac_nc_inq_dimid(ncid, "cell", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &ncells));
    yac_nc_inq_dimid(ncid, "vertex", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &nvertices));
    yac_nc_inq_dimid(ncid, "edge", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &nedges));
    yac_nc_inq_dimid(ncid, "nv", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &nv));
    yac_nc_inq_dimid(ncid, "ne", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &ne));

    // determine local range for cell and vertex data
    read_local_start_cell =
      ((unsigned long)ncells * io_proc_idx) / (unsigned long)num_io_ranks;
    read_num_local_cells =
      ((unsigned long)ncells * (io_proc_idx+1)) / (unsigned long)num_io_ranks -
      (unsigned long)read_local_start_cell;
    read_local_start_vertex =
      ((unsigned long)nvertices * io_proc_idx) / (unsigned long)num_io_ranks;
    read_num_local_vertices =
      ((unsigned long)nvertices * (io_proc_idx+1)) / (unsigned long)num_io_ranks -
      (unsigned long)read_local_start_vertex;
    read_local_start_edge =
      ((unsigned long)nedges * io_proc_idx) / (unsigned long)num_io_ranks;
    read_num_local_edges =
      ((unsigned long)nedges * (io_proc_idx+1)) / (unsigned long)num_io_ranks -
      (unsigned long)read_local_start_edge;

    // read basic grid data (each process its individual part)
    read_vertex_lon = xmalloc(read_num_local_vertices * sizeof(*read_vertex_lon));
    read_vertex_lat = xmalloc(read_num_local_vertices * sizeof(*read_vertex_lat));
    read_dist_vertex_of_cell =
      xmalloc(read_num_local_cells * nv * sizeof(*read_dist_vertex_of_cell));
    read_dist_edge_of_cell =
      xmalloc(read_num_local_cells * nv * sizeof(*read_dist_edge_of_cell));
    read_dist_edge_vertices =
      xmalloc(read_num_local_edges * 2 * sizeof(*read_dist_edge_vertices));
    int varid;
    yac_nc_inq_varid(ncid, "vlon", &varid);
    YAC_HANDLE_ERROR(nc_get_vara_double(ncid, varid, &read_local_start_vertex,
                                    &read_num_local_vertices, read_vertex_lon));
    convert_to_rad(ncid, varid, read_vertex_lon, read_num_local_vertices);
    yac_nc_inq_varid(ncid, "vlat", &varid);
    YAC_HANDLE_ERROR(nc_get_vara_double(ncid, varid, &read_local_start_vertex,
                                    &read_num_local_vertices, read_vertex_lat));
    convert_to_rad(ncid, varid, read_vertex_lat, read_num_local_vertices);
    {
      int * buffer =
        xmalloc(
          MAX(nv*read_num_local_cells, 2*read_num_local_edges) *
          sizeof(*buffer));
      {
        size_t tmp_start[2] = {0, read_local_start_cell};
        size_t tmp_count[2] = {nv, read_num_local_cells};
        yac_nc_inq_varid(ncid, "vertex_of_cell", &varid);
        YAC_HANDLE_ERROR(nc_get_vara_int(ncid, varid, tmp_start, tmp_count, buffer));
        for (size_t i = 0; i < read_num_local_cells; ++i)
          for (size_t j = 0; j < nv; ++j)
            read_dist_vertex_of_cell[i * nv + j] =
              buffer[i + j * read_num_local_cells];
        yac_nc_inq_varid(ncid, "edge_of_cell", &varid);
        YAC_HANDLE_ERROR(nc_get_vara_int(ncid, varid, tmp_start, tmp_count, buffer));
        for (size_t i = 0; i < read_num_local_cells; ++i)
          for (size_t j = 0; j < nv; ++j)
            read_dist_edge_of_cell[i * nv + j] =
              buffer[i + j * read_num_local_cells];
      }
      {
        size_t tmp_start[2] = {0, read_local_start_edge};
        size_t tmp_count[2] = {2, read_num_local_edges};
        yac_nc_inq_varid(ncid, "edge_vertices", &varid);
        YAC_HANDLE_ERROR(nc_get_vara_int(ncid, varid, tmp_start, tmp_count, buffer));
        for (size_t i = 0; i < read_num_local_edges; ++i)
          for (int j = 0; j < 2; ++j)
            read_dist_edge_vertices[i * 2 + j] =
              buffer[i + j * read_num_local_edges];
      }
      free(buffer);

      // adjust for c indices
      for (size_t i = 0; i < read_num_local_cells * nv; ++i)
        if (read_dist_vertex_of_cell[i] > 0) read_dist_vertex_of_cell[i]--;
      for (size_t i = 0; i < read_num_local_cells * nv; ++i)
        if (read_dist_edge_of_cell[i] > 0) read_dist_edge_of_cell[i]--;
      for (size_t i = 0; i < read_num_local_edges * 2; ++i)
        if (read_dist_edge_vertices[i] > 0) read_dist_edge_vertices[i]--;
    }

    YAC_HANDLE_ERROR(nc_close(ncid));
  } else {
    read_vertex_lon = xmalloc(1 * sizeof(*read_vertex_lon));
    read_vertex_lat = xmalloc(1 * sizeof(*read_vertex_lat));
    read_dist_vertex_of_cell = xmalloc(1 * sizeof(*read_dist_vertex_of_cell));
    read_dist_edge_of_cell = xmalloc(1 * sizeof(*read_dist_edge_of_cell));
    read_dist_edge_vertices = xmalloc(1 * sizeof(*read_dist_edge_vertices));
  }

  free(io_ranks);

  {
    int tmp;
    if (comm_rank == 0) tmp = (int)ncells;
    MPI_Bcast(&tmp, 1, MPI_INT, 0, comm);
    ncells = (size_t)tmp;
    if (comm_rank == 0) tmp = (int)nvertices;
    MPI_Bcast(&tmp, 1, MPI_INT, 0, comm);
    nvertices = (size_t)tmp;
    if (comm_rank == 0) tmp = (int)nedges;
    MPI_Bcast(&tmp, 1, MPI_INT, 0, comm);
    nedges = (size_t)tmp;
    if (comm_rank == 0) tmp = (int)nv;
    MPI_Bcast(&tmp, 1, MPI_INT, 0, comm);
    nv = (size_t)tmp;
    if (comm_rank == 0) tmp = (int)ne;
    MPI_Bcast(&tmp, 1, MPI_INT, 0, comm);
    ne = (size_t)tmp;
  }

  // determine local range for cell and vertex data
  size_t local_start_cell =
    ((unsigned long)ncells * (unsigned long)comm_rank) /
    (unsigned long)comm_size;
  size_t num_local_cells =
    ((unsigned long)ncells * ((unsigned long)comm_rank+1)) /
    (unsigned long)comm_size - (unsigned long)local_start_cell;
  size_t local_start_vertex =
    ((unsigned long)nvertices * (unsigned long)comm_rank) /
    (unsigned long)comm_size;
  size_t num_local_vertices =
    ((unsigned long)nvertices * ((unsigned long)comm_rank+1)) /
    (unsigned long)comm_size - (unsigned long)local_start_vertex;
  size_t local_start_edge =
    ((unsigned long)nedges * (unsigned long)comm_rank) /
    (unsigned long)comm_size;
  size_t num_local_edges =
    ((unsigned long)nedges * ((unsigned long)comm_rank+1)) /
    (unsigned long)comm_size - (unsigned long)local_start_edge;

  // redistribute basic cell data (from io decomposition)
  int * dist_vertex_of_cell =
    xmalloc(num_local_cells * nv * sizeof(*dist_vertex_of_cell));
  int * dist_edge_of_cell =
    xmalloc(num_local_cells * nv * sizeof(*dist_edge_of_cell));
  {
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (unsigned i = 0; i < read_num_local_cells; ++i)
      send_count[
        partition_idx_from_element_idx(
          read_local_start_cell + i, ncells, comm_size)] += nv;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    MPI_Alltoallv(read_dist_vertex_of_cell, send_count, send_displ, MPI_INT,
                  dist_vertex_of_cell, recv_count, recv_displ, MPI_INT, comm);
    MPI_Alltoallv(read_dist_edge_of_cell, send_count, send_displ, MPI_INT,
                  dist_edge_of_cell, recv_count, recv_displ, MPI_INT, comm);

    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
    free(read_dist_vertex_of_cell);
    free(read_dist_edge_of_cell);
  }

  // redistribute basic vertex data (from io decomposition)
  double * vertex_lon = xmalloc(num_local_vertices * sizeof(*vertex_lon));
  double * vertex_lat = xmalloc(num_local_vertices * sizeof(*vertex_lat));
  {
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (unsigned i = 0; i < read_num_local_vertices; ++i)
      send_count[
        partition_idx_from_element_idx(
          read_local_start_vertex + i, nvertices, comm_size)]++;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    MPI_Alltoallv(read_vertex_lon, send_count, send_displ, MPI_DOUBLE,
                  vertex_lon, recv_count, recv_displ, MPI_DOUBLE, comm);
    MPI_Alltoallv(read_vertex_lat, send_count, send_displ, MPI_DOUBLE,
                  vertex_lat, recv_count, recv_displ, MPI_DOUBLE, comm);

    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
    free(read_vertex_lon);
    free(read_vertex_lat);
  }

  // redistribute basic edge data (from io decomposition)
  int * dist_edge_vertices =
    xmalloc(num_local_edges * 2 * sizeof(*dist_edge_vertices));
  {
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (unsigned i = 0; i < read_num_local_edges; ++i)
      send_count[
        partition_idx_from_element_idx(
          read_local_start_edge + i, nedges, comm_size)] += 2;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    MPI_Alltoallv(read_dist_edge_vertices, send_count, send_displ, MPI_INT,
                  dist_edge_vertices, recv_count, recv_displ, MPI_INT, comm);

    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
    free(read_dist_edge_vertices);
  }

  // determine required vertices for core cells
  // in additional compute num_cells_per_vertex, vertex_to_cell,
  // and cell_to_vertex
  size_t num_core_vertices;
  yac_int * core_vertices;
  int * num_cells_per_vertex;
  size_t * vertex_to_cell;
  size_t * cell_to_vertex;
  {
    size_t N = num_local_cells * nv;
    core_vertices = xmalloc(N * sizeof(*core_vertices));
    num_cells_per_vertex = xmalloc(N * sizeof(*num_cells_per_vertex));
    vertex_to_cell = xmalloc(N * sizeof(*vertex_to_cell));
    cell_to_vertex = xmalloc(N * sizeof(*cell_to_vertex));
    for (size_t i = 0; i < N; ++i)
      core_vertices[i] = (yac_int)dist_vertex_of_cell[i];
    size_t * permutation = vertex_to_cell;
    for (size_t i = 0; i < N; ++i) permutation[i] = i;
    yac_quicksort_index_yac_int_size_t(core_vertices, N, permutation);
    // remove duplicated core vertices and count number of cells per vertex
    yac_int prev_vertex_id = XT_INT_MAX;
    num_core_vertices = 0;
    for (size_t i = 0; i < N; ++i) {
      yac_int curr_vertex_id = core_vertices[i];
      if (prev_vertex_id == curr_vertex_id) {
        num_cells_per_vertex[num_core_vertices-1]++;
      } else {
        num_cells_per_vertex[num_core_vertices] = 1;
        core_vertices[num_core_vertices] = (prev_vertex_id = curr_vertex_id);
        ++num_core_vertices;
      }
      cell_to_vertex[permutation[i]] = num_core_vertices-1;
      permutation[i] /= nv;
    }
    core_vertices =
      xrealloc(core_vertices, num_core_vertices * sizeof(*core_vertices));
    num_cells_per_vertex =
      xrealloc(num_cells_per_vertex,
               num_core_vertices * sizeof(*num_cells_per_vertex));
    free(dist_vertex_of_cell);
  }

  // determine required edges for core cells
  // in additional compute edge_to_cell and cell_to_edge
  size_t num_core_edges;
  yac_int * core_edges;
  size_t * cell_to_edge;
  {
    size_t N = num_local_cells * nv;
    core_edges = xmalloc(N * sizeof(*core_edges));
    size_t * permutation = xmalloc(N * sizeof(*permutation));
    cell_to_edge = xmalloc(N * sizeof(*cell_to_edge));
    for (size_t i = 0; i < N; ++i)
      core_edges[i] = (yac_int)dist_edge_of_cell[i];
    for (size_t i = 0; i < N; ++i) permutation[i] = i;
    yac_quicksort_index_yac_int_size_t(core_edges, N, permutation);
    // remove duplicated core edges and count number of cells per edge
    yac_int prev_edge_id = XT_INT_MAX;
    num_core_edges = 0;
    for (size_t i = 0; i < N; ++i) {
      yac_int curr_edge_id = core_edges[i];
      if (prev_edge_id != curr_edge_id)
        core_edges[num_core_edges++] = (prev_edge_id = curr_edge_id);
      cell_to_edge[permutation[i]] = num_core_edges-1;
      permutation[i] /= nv;
    }
    free(permutation);
    core_edges =
      xrealloc(core_edges, num_core_edges * sizeof(*core_edges));
    free(dist_edge_of_cell);
  }

  // generate vertex coordinate data
  yac_coordinate_pointer vertex_coordinates =
    xmalloc(num_core_vertices * sizeof(*vertex_coordinates));
  {
    double * local_vertex_lon =
      xmalloc(num_core_vertices * sizeof(*local_vertex_lon));
    double * local_vertex_lat =
      xmalloc(num_core_vertices * sizeof(*local_vertex_lat));
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (size_t i = 0; i < num_core_vertices; ++i)
      send_count[((unsigned long)(core_vertices[i]) *
                  (unsigned long)comm_size + (unsigned long)comm_size - 1) /
                 (unsigned long)nvertices]++;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    int num_all_local_vertices_remote = 0;
    for (int i = 0; i < comm_size; ++i)
      num_all_local_vertices_remote += recv_count[i];

    yac_int * remote_vertex_buffer =
      xmalloc(num_all_local_vertices_remote * sizeof(*remote_vertex_buffer));

    MPI_Alltoallv(
      core_vertices, send_count, send_displ, yac_int_dt,
      remote_vertex_buffer, recv_count, recv_displ, yac_int_dt, comm);

    double * send_vertex_lon =
      xmalloc(num_all_local_vertices_remote * sizeof(*send_vertex_lon));
    double * send_vertex_lat =
      xmalloc(num_all_local_vertices_remote * sizeof(*send_vertex_lat));

    for (int i = 0, l = 0; i < comm_size; ++i) {
      for (int j = 0; j < recv_count[i]; ++j, ++l) {
        size_t idx = (size_t)(remote_vertex_buffer[l]) - local_start_vertex;
        send_vertex_lon[l] = vertex_lon[idx];
        send_vertex_lat[l] = vertex_lat[idx];
      }
    }

    free(remote_vertex_buffer);
    free(vertex_lon);
    free(vertex_lat);

    MPI_Alltoallv(send_vertex_lon, recv_count, recv_displ, MPI_DOUBLE,
                  local_vertex_lon, send_count, send_displ, MPI_DOUBLE, comm);
    MPI_Alltoallv(send_vertex_lat, recv_count, recv_displ, MPI_DOUBLE,
                  local_vertex_lat, send_count, send_displ, MPI_DOUBLE, comm);

    for (unsigned i = 0; i < num_core_vertices; ++i)
      LLtoXYZ(
        local_vertex_lon[i], local_vertex_lat[i], &(vertex_coordinates[i][0]));

    free(send_vertex_lat);
    free(send_vertex_lon);
    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
    free(local_vertex_lon);
    free(local_vertex_lat);
  }

  // generate edge vertex data
  size_t * edge_to_vertex =
    xmalloc(2 * num_core_edges * sizeof(*edge_to_vertex));
  {
    int * local_edge_to_vertex =
      xmalloc(2 * num_core_edges * sizeof(*local_edge_to_vertex));
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (size_t i = 0; i < num_core_edges; ++i)
      send_count[((unsigned long)(core_edges[i]) *
                  (unsigned long)comm_size + (unsigned long)comm_size - 1) /
                 (unsigned long)nedges]++;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    int num_all_local_edges_remote = 0;
    for (int i = 0; i < comm_size; ++i)
      num_all_local_edges_remote += recv_count[i];

    yac_int * remote_edge_buffer =
      xmalloc(num_all_local_edges_remote * sizeof(*remote_edge_buffer));

    MPI_Alltoallv(
      core_edges, send_count, send_displ, yac_int_dt,
      remote_edge_buffer, recv_count, recv_displ, yac_int_dt, comm);

    int * send_edge_vertices =
      xmalloc(2 * num_all_local_edges_remote * sizeof(*send_edge_vertices));

    for (int i = 0, l = 0; i < comm_size; ++i) {
      for (int j = 0; j < recv_count[i]; ++j, ++l) {
        size_t idx = (size_t)(remote_edge_buffer[l]) - local_start_edge;
        send_edge_vertices[2*l+0] = dist_edge_vertices[2*idx+0];
        send_edge_vertices[2*l+1] = dist_edge_vertices[2*idx+1];

      }
      send_count[i] *= 2;
      recv_count[i] *= 2;
      send_displ[i] *= 2;
      recv_displ[i] *= 2;
    }

    free(remote_edge_buffer);
    free(dist_edge_vertices);

    MPI_Alltoallv(send_edge_vertices, recv_count, recv_displ, MPI_INT,
                  local_edge_to_vertex, send_count, send_displ, MPI_INT, comm);

    size_t * permutation = xmalloc(2 * num_core_edges * sizeof(*permutation));
    for (size_t i = 0; i < 2 * num_core_edges; ++i) permutation[i] = i;

    yac_quicksort_index_int_size_t(
      local_edge_to_vertex, 2 * num_core_edges, permutation);

    for (size_t i = 0, j = 0; i < 2 * num_core_edges; ++i) {
      yac_int curr_vertex = (yac_int)(local_edge_to_vertex[i]);
      while ((j < nvertices) && (core_vertices[j] < curr_vertex)) ++j;
      YAC_ASSERT(
        (j < nvertices) && (core_vertices[j] == curr_vertex),
        "ERROR(yac_read_icon_basic_grid_data_parallel): vertex id missmatch")
      edge_to_vertex[permutation[i]] = j;
    }

    free(permutation);
    free(send_edge_vertices);
    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
    free(local_edge_to_vertex);
  }

  // generate cell ids for local partition
  yac_int * cell_ids = xmalloc(num_local_cells * sizeof(*cell_ids));
  for (size_t i = 0; i < num_local_cells; ++i)
    cell_ids[i] = (yac_int)(local_start_cell + i);

  // generate num_vertices_per_cell
  int * num_vertices_per_cell =
    xmalloc(num_local_cells * sizeof(*num_vertices_per_cell));
  for (size_t i = 0; i < num_local_cells; ++i) num_vertices_per_cell[i] = nv;

  // generate edge_type (for icon grids this is always YAC_GREAT_CIRCLE_EDGE)
  enum yac_edge_type * edge_type = xmalloc(num_core_edges * sizeof(*edge_type));
  for (size_t i = 0; i < num_core_edges; ++i) edge_type[i] = YAC_GREAT_CIRCLE_EDGE;

  struct yac_basic_grid_data grid_data;
  grid_data.vertex_coordinates      = vertex_coordinates;
  grid_data.cell_ids                = cell_ids;
  grid_data.vertex_ids              = core_vertices;
  grid_data.edge_ids                = core_edges;
  grid_data.num_cells               = num_local_cells;
  grid_data.num_vertices            = num_core_vertices;
  grid_data.num_edges               = num_core_edges;
  grid_data.core_cell_mask          = generate_simple_core_mask(num_local_cells);
  grid_data.core_vertex_mask        = generate_simple_core_mask(num_core_vertices);
  grid_data.core_edge_mask          = generate_simple_core_mask(num_core_edges);
  grid_data.num_vertices_per_cell   = num_vertices_per_cell;
  grid_data.num_cells_per_vertex    = num_cells_per_vertex;
  grid_data.cell_to_vertex          = cell_to_vertex;
  grid_data.cell_to_vertex_offsets  = generate_offsets(num_local_cells, num_vertices_per_cell);
  grid_data.cell_to_edge            = cell_to_edge;
  grid_data.cell_to_edge_offsets    = grid_data.cell_to_vertex_offsets;
  grid_data.vertex_to_cell          = vertex_to_cell;
  grid_data.vertex_to_cell_offsets  = generate_offsets(num_core_vertices, num_cells_per_vertex);
  grid_data.edge_to_vertex          = (yac_size_t_2_pointer)&(edge_to_vertex[0]);
  grid_data.edge_type               = edge_type;
  grid_data.num_total_cells         = num_local_cells;
  grid_data.num_total_vertices      = num_core_vertices;
  grid_data.num_total_edges         = num_core_edges;

  return grid_data;
}

struct yac_basic_grid * yac_read_icon_basic_grid_parallel(
  char const * filename, char const * gridname, MPI_Comm comm) {

  return
    yac_basic_grid_new(
      gridname,
      yac_read_icon_basic_grid_data_parallel(filename, comm));
}

struct yac_basic_grid * yac_read_icon_basic_grid_parallel_f2c(
  char const * filename, char const * gridname, MPI_Fint comm) {

  return
    yac_read_icon_basic_grid_parallel(
      filename, gridname, MPI_Comm_f2c(comm));
}

struct yac_basic_grid_data yac_read_icon_basic_grid_data (
  char const * filename) {

  int nbr_vertices;
  int nbr_cells;
  int * num_vertices_per_cell = NULL;
  int * cell_to_vertex = NULL;
  int * cell_mask = NULL;

  double * x_vertices = NULL;
  double * y_vertices = NULL;
  double * x_cells = NULL;
  double * y_cells = NULL;

  yac_read_icon_grid_information(filename, &nbr_vertices, &nbr_cells,
                                 &num_vertices_per_cell, &cell_to_vertex,
                                 &x_vertices, &y_vertices,
                                 &x_cells, &y_cells,
                                 &cell_mask);

  free(x_cells);
  free(y_cells);
  free(cell_mask);

  struct yac_basic_grid_data grid =
    yac_generate_basic_grid_data_unstruct(
      (size_t)nbr_vertices, (int)nbr_cells, num_vertices_per_cell,
      x_vertices, y_vertices, cell_to_vertex);
  free(num_vertices_per_cell);
  free(x_vertices);
  free(y_vertices);
  free(cell_to_vertex);

  return grid;
}

struct yac_basic_grid * yac_read_icon_basic_grid(
  char const * filename, char const * gridname) {

  return
    yac_basic_grid_new(
      gridname, yac_read_icon_basic_grid_data(filename));
}

/* ---------------------------------------------------------------- */

static int * get_icon_cell_mask ( int ncid, size_t nbr_cells ) {

  // get variable id
  int mask_id;
  yac_nc_inq_varid (ncid, "cell_sea_land_mask", &mask_id);

  // check number of dimension (has to be 1)
  int ndims;
  YAC_HANDLE_ERROR(nc_inq_varndims(ncid, mask_id, &ndims));
  YAC_ASSERT(
    ndims == 1,
    "ERROR(get_icon_cell_mask): mask array has more than one dimension")

  // get id of dimension
  int dimid;
  YAC_HANDLE_ERROR(nc_inq_vardimid(ncid, mask_id, &dimid));

  // check size of mask (has to be equal to nbr_cells)
  size_t dimlen;
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &dimlen));
  YAC_ASSERT(
    dimlen == nbr_cells,
    "ERROR(get_icon_cell_mask): invalid size of mask array")

  // check mask type (has to be NC_INT)
  nc_type mask_type;
  YAC_HANDLE_ERROR(nc_inq_vartype(ncid, mask_id, &mask_type));
  YAC_ASSERT(
    mask_type == NC_INT, "ERROR(get_icon_cell_mask): invalid mask type")

  // get and return mask
  int * cell_mask = xmalloc(nbr_cells * sizeof(*cell_mask));
  YAC_HANDLE_ERROR(nc_get_var_int (ncid, mask_id, cell_mask));
  return cell_mask;
}

/* ---------------------------------------------------------------- */

static void get_icon_connect(
  int ncid, size_t nbr_cells, int ** vertex_of_cell, size_t * nv) {

  // get variable id
  int conn_id;
  yac_nc_inq_varid(ncid, "vertex_of_cell", &conn_id);

  // check number of dimension (has to be 1)
  int ndims;
  YAC_HANDLE_ERROR(nc_inq_varndims(ncid, conn_id, &ndims));
  YAC_ASSERT(
    ndims == 2,
    "ERROR(get_icon_connect): "
    "connectivity array has invalid number of dimensions")

  // get ids of dimensions
  int dimids[2];
  YAC_HANDLE_ERROR(nc_inq_vardimid(ncid, conn_id, dimids));

  // check size of dimensions
  size_t dimlen;
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimids[0], nv));
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimids[1], &dimlen));
  YAC_ASSERT(
    dimlen == nbr_cells,
    "ERROR(get_icon_connect): invalid size of second dimension of "
    "connectivity array (has to be nbr_cells)")

  // get and return connectivity array
  *vertex_of_cell = xmalloc(*nv * nbr_cells * sizeof(**vertex_of_cell));
  YAC_HANDLE_ERROR(nc_get_var_int (ncid, conn_id, *vertex_of_cell));
}

void yac_delete_icon_grid_data( int ** cell_mask,
                                int ** global_cell_id,
                                int ** global_cell_id_rank,
                                int ** num_vertices_per_cell,
                                int ** global_corner_id,
                                int ** global_corner_id_rank,
                                int ** cell_to_vertex,
                                double ** x_cells,
                                double ** y_cells,
                                double ** x_vertices,
                                double ** y_vertices) {

    free (*cell_mask);
    free (*global_cell_id);
    free (*global_cell_id_rank);
    free (*num_vertices_per_cell);
    free (*global_corner_id);
    free (*global_corner_id_rank);
    free (*cell_to_vertex);
    free (*x_cells);
    free (*y_cells);
    free (*x_vertices);
    free (*y_vertices);
}

#else

struct yac_basic_grid_data yac_read_icon_basic_grid_data(
  char const * filename) {

   UNUSED(filename);
   die(
     "ERROR(yac_basic_grid_data): "
     "YAC is built without the NetCDF support");

   return
    yac_generate_basic_grid_data_reg_2d(
      (size_t[]){0,0}, (int[]){0,0}, NULL, NULL);
}

struct yac_basic_grid * yac_read_icon_basic_grid(
  char const * filename, char const * gridname) {

   UNUSED(filename);
   UNUSED(gridname);
   die(
     "ERROR(yac_read_icon_basic_grid): "
     "YAC is built without the NetCDF support");

   return NULL;
}

void yac_read_icon_grid_information(const char * filename, int * nbr_vertices,
                                    int * nbr_cells, int ** num_vertices_per_cell,
                                    int ** cell_to_vertex, double ** x_vertices,
                                    double ** y_vertices, double ** x_cells,
                                    double ** y_cells, int ** cell_mask) {

   UNUSED(filename);
   UNUSED(nbr_vertices);
   UNUSED(nbr_cells);
   UNUSED(num_vertices_per_cell);
   UNUSED(cell_to_vertex);
   UNUSED(x_vertices);
   UNUSED(y_vertices);
   UNUSED(x_cells);
   UNUSED(y_cells);
   UNUSED(cell_mask);
   die(
     "ERROR(yac_read_icon_grid_information): "
     "YAC is built without the NetCDF support");
}

void yac_read_part_icon_grid_information(
  const char * filename, int * num_vertices, int * num_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells, int ** global_cell_id,
  int ** cell_mask, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size) {

   UNUSED(filename);
   UNUSED(num_vertices);
   UNUSED(num_cells);
   UNUSED(num_vertices_per_cell);
   UNUSED(cell_to_vertex);
   UNUSED(x_vertices);
   UNUSED(y_vertices);
   UNUSED(x_cells);
   UNUSED(y_cells);
   UNUSED(global_cell_id);
   UNUSED(cell_mask);
   UNUSED(cell_core_mask);
   UNUSED(global_corner_id);
   UNUSED(corner_core_mask);
   UNUSED(rank);
   UNUSED(size);
   die(
     "ERROR(yac_read_part_icon_grid_information): "
     "YAC is built without the NetCDF support");
}

void yac_read_icon_grid_information_parallel(
  const char * filename, MPI_Comm comm, int * nbr_vertices, int * nbr_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex, int ** global_cell_ids,
  int ** cell_owner, int ** global_vertex_ids, int ** vertex_owner,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells, int ** cell_msk) {

   UNUSED(filename);
   UNUSED(comm);
   UNUSED(nbr_vertices);
   UNUSED(nbr_cells);
   UNUSED(num_vertices_per_cell);
   UNUSED(cell_to_vertex);
   UNUSED(global_cell_ids);
   UNUSED(cell_owner);
   UNUSED(global_vertex_ids);
   UNUSED(vertex_owner);
   UNUSED(x_vertices);
   UNUSED(y_vertices);
   UNUSED(x_cells);
   UNUSED(y_cells);
   UNUSED(cell_msk);
   die(
     "ERROR(yac_read_icon_grid_information_parallel): "
     "YAC is built without the NetCDF support");
}

struct yac_basic_grid_data yac_read_icon_basic_grid_data_parallel(
  const char * filename, MPI_Comm comm) {

   UNUSED(filename);
   UNUSED(comm);
   die(
     "ERROR(yac_read_icon_basic_grid_data_parallel): "
     "YAC is built without the NetCDF support");

   return
    yac_generate_basic_grid_data_reg_2d(
      (size_t[]){0,0}, (int[]){0,0}, NULL, NULL);
}

struct yac_basic_grid * yac_read_icon_basic_grid_parallel(
  char const * filename, char const * gridname, MPI_Comm comm) {

   UNUSED(filename);
   UNUSED(gridname);
   UNUSED(comm);
   die(
     "ERROR(yac_read_icon_basic_grid_parallel): "
     "YAC is built without the NetCDF support");

   return NULL;
}

void yac_delete_icon_grid_data( int ** cell_mask,
                                int ** global_cell_id,
                                int ** cell_core_mask,
                                int ** num_vertices_per_cell,
                                int ** global_corner_id,
                                int ** corner_core_mask,
                                int ** cell_to_vertex,
                                double ** x_cells,
                                double ** y_cells,
                                double ** x_vertices,
                                double ** y_vertices) {

   UNUSED(cell_mask);
   UNUSED(global_cell_id);
   UNUSED(cell_core_mask);
   UNUSED(num_vertices_per_cell);
   UNUSED(global_corner_id);
   UNUSED(corner_core_mask);
   UNUSED(cell_to_vertex);
   UNUSED(x_vertices);
   UNUSED(y_vertices);
   UNUSED(x_cells);
   UNUSED(y_cells);
   die(
     "ERROR(yac_delete_icon_grid_data): "
     "YAC is built without the NetCDF support");
}

#endif // YAC_NETCDF_ENABLED
