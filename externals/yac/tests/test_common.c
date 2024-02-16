/**
 * @file test_common.c
 *
 * @copyright Copyright  (C)  2018 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler  <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

#include "geometry.h"
#include "test_common.h"
#include "utils.h"

static struct grid_cell generate_cell_func(
  double * lon, double * lat, enum yac_edge_type * edge_type,
  size_t num_corners, void (*fun_LLtoXYZ)(double,double,double[])) {

  struct grid_cell cell;

  cell.coordinates_xyz = xmalloc(num_corners * sizeof(*(cell.coordinates_xyz)));
  cell.edge_type = xmalloc(num_corners * sizeof(*(cell.edge_type)));
  cell.num_corners = num_corners;
  cell.array_size = num_corners;
  for (size_t i = 0; i < num_corners; ++i)
    (*fun_LLtoXYZ)(lon[i], lat[i], cell.coordinates_xyz[i]);
  memcpy(cell.edge_type, edge_type, num_corners * sizeof(*edge_type));

  return cell;
}

struct grid_cell generate_cell_deg(
  double * lon, double * lat, enum yac_edge_type * edge_type,
  size_t num_corners) {
  return generate_cell_func(lon, lat, edge_type, num_corners, LLtoXYZ_deg);
}

struct grid_cell generate_cell_rad(
  double * lon, double * lat, enum yac_edge_type * edge_type,
  size_t num_corners) {
  return generate_cell_func(lon, lat, edge_type, num_corners, LLtoXYZ);
}

struct grid_cell generate_cell_3d(
  coordinate_pointer coords, enum yac_edge_type * edge_type,
  size_t num_corners) {

  struct grid_cell cell;

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
