// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

#include "test_common.h"
#include "geometry.h"
#include "tests.h"

static struct yac_grid_cell generate_cell_func(
  double * lon, double * lat, enum yac_edge_type * edge_type,
  size_t num_corners, void (*fun_LLtoXYZ)(double,double,double[])) {

  struct yac_grid_cell cell;

  cell.coordinates_xyz = xmalloc(num_corners * sizeof(*(cell.coordinates_xyz)));
  cell.edge_type = xmalloc(num_corners * sizeof(*(cell.edge_type)));
  cell.num_corners = num_corners;
  cell.array_size = num_corners;
  for (size_t i = 0; i < num_corners; ++i)
    (*fun_LLtoXYZ)(lon[i], lat[i], cell.coordinates_xyz[i]);
  memcpy(cell.edge_type, edge_type, num_corners * sizeof(*edge_type));

  return cell;
}

struct yac_grid_cell generate_cell_deg(
  double * lon, double * lat, enum yac_edge_type * edge_type,
  size_t num_corners) {
  return generate_cell_func(lon, lat, edge_type, num_corners, LLtoXYZ_deg);
}

struct yac_grid_cell generate_cell_rad(
  double * lon, double * lat, enum yac_edge_type * edge_type,
  size_t num_corners) {
  return generate_cell_func(lon, lat, edge_type, num_corners, LLtoXYZ);
}

struct yac_grid_cell generate_cell_3d(
  yac_coordinate_pointer coords, enum yac_edge_type * edge_type,
  size_t num_corners) {

  struct yac_grid_cell cell;

  cell.coordinates_xyz = xmalloc(num_corners * sizeof(*(cell.coordinates_xyz)));
  cell.edge_type = xmalloc(num_corners * sizeof(*(cell.edge_type)));
  cell.num_corners = num_corners;
  cell.array_size = num_corners;
  memcpy(cell.coordinates_xyz, coords, num_corners * sizeof(*coords));
  memcpy(cell.edge_type, edge_type, num_corners * sizeof(*edge_type));

  return cell;
}

int intersect(enum yac_edge_type edge_type_a,
              double lon_a, double lat_a, double lon_b, double lat_b,
              enum yac_edge_type edge_type_b,
              double lon_c, double lat_c, double lon_d, double lat_d,
              double * intersection) {

  double a[3], b[3], c[3], d[3], p[3], q[3];

  LLtoXYZ_deg(lon_a, lat_a, a);
  LLtoXYZ_deg(lon_b, lat_b, b);
  LLtoXYZ_deg(lon_c, lat_c, c);
  LLtoXYZ_deg(lon_d, lat_d, d);

  int ret = yac_intersect_vec(edge_type_a, a, b, edge_type_b, c, d, p, q);

  switch (ret) {
    case ((1 << 0)            | (1 << 2)):
    case ((1 << 0) | (1 << 1) | (1 << 2)):
    case ((1 << 0)            | (1 << 2) | (1 << 3)):
    case ((1 << 0) | (1 << 1) | (1 << 2) | (1 << 3)):
      intersection[0] = p[0];
      intersection[1] = p[1];
      intersection[2] = p[2];
      return 1;
    case (           (1 << 1)            | (1 << 3)):
    case ((1 << 0) | (1 << 1)            | (1 << 3)):
    case (           (1 << 1) | (1 << 2) | (1 << 3)):
      intersection[0] = q[0];
      intersection[1] = q[1];
      intersection[2] = q[2];
      return 1;
    default:
      return 0;
  }
}

void * to_pointer(void * data, size_t size_data) {

  void * ret_value = xmalloc(size_data);
  memcpy(ret_value, data, size_data);
  return ret_value;
}

int double_are_equal(double a, double b) {

  return (a > b) == (a < b);
}

int double_are_unequal(double a, double b) {

  return (a > b) != (a < b);
}

void set_even_io_rank_list(MPI_Comm comm) {

  int comm_size;
  MPI_Comm_size(comm, &comm_size);

  char * io_rank_list = xmalloc(16 * comm_size);
  io_rank_list[0] = '\0';
  for (int i = 0; i < comm_size; i += 2) {
    char rank[16];
    snprintf(rank, sizeof(rank), "%d,", i);
    strcat(io_rank_list, rank);
  }
  char comm_size_str[16];
  snprintf(comm_size_str, sizeof(comm_size_str), "%d", comm_size);
  setenv("YAC_IO_RANK_LIST", io_rank_list, 1);
  setenv("YAC_IO_MAX_NUM_RANKS_PER_NODE", comm_size_str, 1);
  free(io_rank_list);
}

void clear_yac_io_env() {
  unsetenv("YAC_IO_RANK_LIST");
  unsetenv("YAC_IO_MAX_NUM_RANKS");
  unsetenv("YAC_IO_RANK_EXCLUDE_LIST");
  unsetenv("YAC_IO_MAX_NUM_RANKS_PER_NODE");
}

static int check_indices(
  size_t const * indices_a, size_t const * indices_b, size_t count) {

  if (count == 0) return 1;

  if (count == 1) return *indices_a == *indices_b;

  size_t i;
  for (i = 0; i < count; ++i) if (indices_a[i] == indices_b[0]) break;
  if (i == count) return 0;

  if (indices_a[(i+1)%count] == indices_b[1]) {

    for (size_t j = 2; j < count; ++j)
      if (indices_a[(i+j)%count] !=indices_b[j])
        return 0;

  } else if (indices_a[(i+1)%count] == indices_b[count - 1]) {

    for (size_t j = 2; j < count; ++j)
      if (indices_a[(i+j)%count] != indices_b[count - j])
        return 0;

  } else return 0;

  return 1;
}

void check_basic_grid_data(
  struct yac_basic_grid_data grid_a, struct yac_basic_grid_data grid_b,
  char const * grid_name) {

  printf("testing grid: %s\n", grid_name);

  if ((grid_a.num_cells != grid_b.num_cells) ||
      (grid_a.num_total_cells != grid_b.num_total_cells))
    PUT_ERR("error in grid.num_cells or grid.num_total_cells\n")
  if ((grid_a.num_vertices != grid_b.num_vertices) ||
      (grid_a.num_total_vertices != grid_b.num_total_vertices))
    PUT_ERR("error in grid.num_vertices or grid.num_total_vertices\n")
  if ((grid_a.num_edges != grid_b.num_edges) ||
      (grid_a.num_total_edges != grid_b.num_total_edges))
    PUT_ERR("error in grid.num_edges or grid.num_total_edges\n")

  for (size_t i = 0; i < grid_a.num_cells; ++i)
    if (grid_a.num_vertices_per_cell[i] != grid_b.num_vertices_per_cell[i])
      PUT_ERR("error in grid.num_vertices_per_cell\n")

  for (size_t i = 0; i < grid_a.num_vertices; ++i) {
    if (grid_a.num_cells_per_vertex[i] != grid_b.num_cells_per_vertex[i])
      PUT_ERR("error in grid.num_cells_per_vertex\n")

    if (get_vector_angle(
          grid_a.vertex_coordinates[i],
          grid_b.vertex_coordinates[i]) > yac_angle_tol)
      PUT_ERR("error in grid.coordinates_xyz\n")
  }

  for (size_t i = 0, offset = 0; i < grid_a.num_cells; ++i) {

    size_t num_vertices = grid_a.num_vertices_per_cell[i];

    if ((grid_a.cell_to_vertex_offsets[i] !=
         grid_b.cell_to_vertex_offsets[i]) ||
        (grid_a.cell_to_vertex_offsets[i] != offset))
      PUT_ERR("error in grid.cell_to_vertex_offsets\n")

    if (!check_indices(grid_a.cell_to_vertex + offset,
                       grid_b.cell_to_vertex + offset, num_vertices))
      PUT_ERR("error in grid.cell_to_vertex\n")

    if ((grid_a.cell_to_edge_offsets[i] !=
         grid_b.cell_to_edge_offsets[i]) ||
        (grid_a.cell_to_edge_offsets[i] != offset))
      PUT_ERR("error in grid.cell_to_edge_offsets\n")

    if (!check_indices(grid_a.cell_to_edge + offset,
                       grid_b.cell_to_edge + offset, num_vertices))
      PUT_ERR("error in grid.cell_to_edge\n")

    offset += num_vertices;
  }

  for (size_t i = 0, offset = 0; i < grid_a.num_vertices; ++i) {

    size_t num_cells = grid_a.num_cells_per_vertex[i];

    if ((grid_a.vertex_to_cell_offsets[i] !=
         grid_b.vertex_to_cell_offsets[i]) ||
        (grid_a.vertex_to_cell_offsets[i] != offset))
      PUT_ERR("error in grid.vertex_to_cell_offsets\n")

    if (!check_indices(grid_a.vertex_to_cell + offset,
                       grid_b.vertex_to_cell + offset, num_cells))
      PUT_ERR("error in grid.vertex_to_cell\n")

    offset += num_cells;
  }

  for (size_t i = 0; i < grid_a.num_edges; ++i) {
    if (!check_indices(grid_a.edge_to_vertex[i], grid_b.edge_to_vertex[i], 2))
      PUT_ERR("error in grid.vertex_to_cell\n")
    if (grid_a.edge_type[i] != grid_b.edge_type[i])
      PUT_ERR("error in grid.edge_type");
  }
}
