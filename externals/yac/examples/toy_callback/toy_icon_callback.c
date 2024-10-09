// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "yac.h"
#include "yac_utils.h"

/* ------------------------------------------------- */

const char fieldName[] = "icon_to_cube";

struct grid_info {
  yac_int * global_cell_id;
  yac_int * global_corner_id;
  int * cell_to_vertex;
  int * cell_to_vertex_offset;
  int * cell_core_mask;
  int * num_vertices_per_cell;
};

static void compute_weights_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data);
static inline void LLtoXYZ(double lon, double lat, double p_out[]);
static inline void XYZtoLL (double const p_in[], double * lon, double * lat);

#define STR_USAGE "Usage: %s -c configFilename -g gridFilename\n"
#define YAC_ASSERT_ARGS(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(
  int argc, char ** argv,
  char const ** configFilename, char const ** gridFilename);

int main (int argc, char *argv[]) {

  // Initialisation of MPI

  MPI_Init(0, NULL);
  xt_initialize(MPI_COMM_WORLD);

  char const * configFilename = "toy_callback.yaml"; // default configuration file
  char const * gridFilename = "icon_grid_0043_R02B04_G.nc"; // default grid file
  parse_arguments(argc, argv, &configFilename, &gridFilename);
  yac_cinit ();
  yac_cread_config_yaml(configFilename);

  int comp_id;
  char * comp_name = "ICON";

  yac_cdef_comp ( comp_name, &comp_id );

  MPI_Comm local_comm;
  yac_cget_comp_comm(comp_id, &local_comm);

  int rank, size;
  MPI_Comm_rank(local_comm,&rank);
  MPI_Comm_size(local_comm,&size);

  MPI_Comm_free(&local_comm);


  int cell_point_id;
  int corner_point_id;

  int field_id;
  int grid_id;

  /* Grid definition for the first component (ICON) */

  int num_vertices;
  int num_cells;
  int * num_vertices_per_cell;
  int * cell_to_vertex;
  int * cell_to_vertex_offset;
  double * x_vertices;
  double * y_vertices;
  double * x_cells;
  double * y_cells;

  int * cell_mask;
  int * global_cell_id;
  int * cell_core_mask;
  int * global_corner_id;
  int * corner_core_mask;

  yac_read_part_icon_grid_information(gridFilename, &num_vertices, &num_cells,
                                      &num_vertices_per_cell, &cell_to_vertex,
                                      &x_vertices, &y_vertices, &x_cells,
                                      &y_cells, &global_cell_id,
                                      &cell_mask,
                                      &cell_core_mask, &global_corner_id,
                                      &corner_core_mask, rank, size);

  cell_to_vertex_offset = malloc(num_cells * sizeof(*cell_to_vertex_offset));
  for (int i = 0, accu = 0; i < num_cells; ++i) {
    cell_to_vertex_offset[i] = accu;
    accu += num_vertices_per_cell[i];
  }

  double * x_points, * y_points;

  x_points = x_vertices;
  y_points = y_vertices;

  double * x_center, * y_center;

  x_center = x_cells;
  y_center = y_cells;

  yac_cdef_grid_unstruct(
    "icon_grid", num_vertices, num_cells, num_vertices_per_cell,
    x_vertices, y_vertices, cell_to_vertex, &grid_id);

  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);
  yac_cset_global_index(global_corner_id, YAC_LOCATION_CORNER, grid_id);
  yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, grid_id);

  yac_cdef_points_unstruct(
    grid_id, num_cells, YAC_LOCATION_CELL, x_center, y_center, &cell_point_id );
  yac_cdef_points_unstruct(
    grid_id, num_vertices, YAC_LOCATION_CORNER, x_points, y_points, &corner_point_id );

  /* Field definition for the first component (ICON-atmosphere) */
  int point_ids[2] = {cell_point_id, corner_point_id};
  int num_points = 2;
  yac_cdef_field(
    fieldName, comp_id, point_ids, num_points, 1, "2", YAC_TIME_UNIT_SECOND,
    &field_id);

  // register callback routine
  struct grid_info grid_info;
  grid_info.global_cell_id = global_cell_id;
  grid_info.global_corner_id = global_corner_id;
  grid_info.cell_to_vertex = cell_to_vertex;
  grid_info.cell_to_vertex_offset = cell_to_vertex_offset;
  grid_info.cell_core_mask = cell_core_mask;
  grid_info.num_vertices_per_cell = num_vertices_per_cell;
  yac_cadd_compute_weights_callback(
    compute_weights_callback, (void*)&grid_info, "compute_weights_callback");

  /* Search. */
  yac_cenddef( );

  double * cell_out = malloc(num_cells * sizeof(*cell_out));
  double * corner_out = malloc(num_vertices * sizeof(*corner_out));

  int err, info;

  for (int i = 0; i < num_cells; ++i) {

    double middle_point[3] = {0, 0, 0};

    for (int j = 0; j < num_vertices_per_cell[i]; ++j) {

      double curr_point[3];

      LLtoXYZ(x_vertices[cell_to_vertex[cell_to_vertex_offset[i] + j]],
              y_vertices[cell_to_vertex[cell_to_vertex_offset[i] + j]],
              curr_point);

      middle_point[0] += curr_point[0];
      middle_point[1] += curr_point[1];
      middle_point[2] += curr_point[2];
    }

    double scale = 1.0 / sqrt(middle_point[0] * middle_point[0] + 
                              middle_point[1] * middle_point[1] + 
                              middle_point[2] * middle_point[2]);

    middle_point[0] *= scale;
    middle_point[1] *= scale;
    middle_point[2] *= scale;

    double lon, lat;

    XYZtoLL(middle_point, &lon, &lat);

    cell_out[i] = yac_test_func(lon, lat);
  }
  for (int i = 0; i < num_vertices; ++i)
    corner_out[i] = yac_test_func(x_vertices[i], y_vertices[i]);

  {
    double *point_set_data[2];
    double **collection_data[1] = {point_set_data};

    point_set_data[0] = cell_out;
    point_set_data[1] = corner_out;
    yac_cput(field_id, 1, collection_data, &info, &err);
  }

#ifdef VTK_OUTPUT
  //----------------------------------------------------------
  // write field to vtk output file
  //----------------------------------------------------------

  char vtk_filename[32];

  sprintf(vtk_filename, "toy_icon_callback_%d.vtk", rank);

  YAC_VTK_FILE *vtk_file = yac_vtk_open(vtk_filename, "unstruct_out");
  yac_vtk_write_point_data_ll(
    vtk_file, x_vertices, y_vertices, num_vertices);
  yac_vtk_write_cell_data(vtk_file, (unsigned *)cell_to_vertex,
                          (unsigned*)num_vertices_per_cell, num_cells);
  yac_vtk_write_point_scalars_int(
    vtk_file, corner_core_mask, num_vertices, "corner_core_mask");
  yac_vtk_write_point_scalars_int(
    vtk_file, global_corner_id, num_vertices, "global_corner_id");
  yac_vtk_write_cell_scalars_int(
    vtk_file, cell_core_mask, num_cells, "cell_core_mask");
  yac_vtk_write_cell_scalars_int(
    vtk_file, global_cell_id, num_cells, "global_cell_id");

  yac_vtk_write_cell_scalars_double(vtk_file, cell_out, num_cells, "cell_out");
  yac_vtk_write_point_scalars_double(
    vtk_file, corner_out, num_vertices, "corner_out");

  yac_vtk_close(vtk_file);

#endif // VTK_OUTPUT

  yac_cfinalize();

  xt_finalize();

  MPI_Finalize();

  free(num_vertices_per_cell);
  free(cell_to_vertex);
  free(cell_to_vertex_offset);
  free(x_cells);
  free(y_cells);
  free(x_vertices);
  free(y_vertices);
  free(cell_mask);
  free(global_cell_id);
  free(cell_core_mask);
  free(global_corner_id);
  free(corner_core_mask);
  free(cell_out);
  free(corner_out);

  return EXIT_SUCCESS;
}

static void compute_weights_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data) {

  (void)(tgt_coords);

  struct grid_info * info = (struct grid_info*)user_data;

  for (size_t i = 0; i < 2; ++i) {
    global_results_points[i] = NULL;
    result_weights[i] = NULL;
    result_count[i] = 0;
  }

  // consistency check
  if ((info->global_cell_id[src_cell_idx] != src_cell_id) ||
      (!info->cell_core_mask[src_cell_idx])) {
    fputs("ERROR(compute_weights_callback): inconsistent data\n", stderr);
    exit(EXIT_FAILURE);
  }

  static double cell_weight[1] = {0.5};
  static double vertex_weights[4] = {0.125, 0.125, 0.125, 0.125};

  static int cell_result_points[1];
  static int vertex_result_points[4];

  for (int i = 0; i < info->num_vertices_per_cell[src_cell_idx]; ++i)
    vertex_result_points[i] =
      (int)(info->global_corner_id[
              info->cell_to_vertex[
                info->cell_to_vertex_offset[src_cell_idx] + i]]);
  cell_result_points[0] = (int)src_cell_id;

  global_results_points[0] = cell_result_points;
  global_results_points[1] = vertex_result_points;
  result_weights[0] = cell_weight;
  result_weights[1] = vertex_weights;
  result_count[0] = 1;
  result_count[1] = 4;
}

static void parse_arguments(
  int argc, char ** argv,
  char const ** configFilename, char const ** gridFilename) {

  int opt;
  while ((opt = getopt(argc, argv, "c:g:")) != -1) {
    YAC_ASSERT_ARGS((opt == 'c') || (opt == 'g'), "invalid command argument")
    switch (opt) {
      default:
      case 'c':
        *configFilename = optarg;
        break;
      case 'g':
        *gridFilename = optarg;
        break;
    }
  }
}

static inline void LLtoXYZ(double lon, double lat, double p_out[]) {

   while (lon < -M_PI) lon += 2.0 * M_PI;
   while (lon >= M_PI) lon -= 2.0 * M_PI;

   double cos_lat = cos(lat);
   p_out[0] = cos_lat * cos(lon);
   p_out[1] = cos_lat * sin(lon);
   p_out[2] = sin(lat);
}

static inline void XYZtoLL (double const p_in[], double * lon, double * lat) {

   *lon = atan2(p_in[1] , p_in[0]);
   *lat = M_PI_2 - acos(p_in[2]);
}
