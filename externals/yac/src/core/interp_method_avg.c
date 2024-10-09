// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <string.h>

#include "interp_method_internal.h"
#include "grid_cell.h"
#include "interp_method_avg.h"
#include "yac_lapack_interface.h"

#define SWAP(x, y, T) do { T SWAP = x; x = y; y = SWAP; } while (0)

static size_t do_search_avg(struct interp_method * method,
                            struct yac_interp_grid * interp_grid,
                            size_t * tgt_points, size_t count,
                            struct yac_interp_weights * weights);
static void delete_avg(struct interp_method * method);

typedef int (*func_compute_weights)(
  double[3], size_t, yac_const_coordinate_pointer, int*, double*,
  enum yac_cell_type cell_type);

static struct interp_method_vtable
  interp_method_avg_vtable = {
    .do_search = do_search_avg,
    .delete = delete_avg};

struct interp_method_avg {

  struct interp_method_vtable * vtable;
  func_compute_weights compute_weights;
};

static int check_src_field_cell(
  double * point, size_t field_cell, size_t cell_size,
  size_t const * vertex_to_cell, size_t const * vertex_to_cell_offsets,
  yac_const_coordinate_pointer field_coordinates,
  struct yac_grid_cell * cell_buffer) {

  if (cell_size == 0) return 0;

  cell_buffer->num_corners = cell_size;

  yac_coordinate_pointer coordinates_xyz;

  if (cell_size > cell_buffer->array_size) {
    cell_buffer->coordinates_xyz =
      ((coordinates_xyz =
          xrealloc(
            cell_buffer->coordinates_xyz,
            cell_size * sizeof(*coordinates_xyz))));
    cell_buffer->edge_type =
      xrealloc(
        cell_buffer->edge_type, cell_size * sizeof(cell_buffer->edge_type));
    for (size_t i = 0; i < cell_size; ++i)
      cell_buffer->edge_type[i] = YAC_GREAT_CIRCLE_EDGE;
    cell_buffer->array_size = cell_size;
  } else {
    coordinates_xyz = cell_buffer->coordinates_xyz;
  }

  size_t const * field_cell_vertices =
    vertex_to_cell + vertex_to_cell_offsets[field_cell];
  for (size_t i = 0; i < cell_size; ++i)
    memcpy(coordinates_xyz[i], field_coordinates[field_cell_vertices[i]],
           sizeof(*coordinates_xyz));

  return yac_point_in_cell(point, *cell_buffer);
}

#define IS_GC(x) (((x) == YAC_GREAT_CIRCLE_EDGE) || ((x) == YAC_LON_CIRCLE_EDGE))
#define IS_LON_LAT(x) (((x) == YAC_LAT_CIRCLE_EDGE) || ((x) == YAC_LON_CIRCLE_EDGE))

static enum yac_cell_type determine_cell_type(
  enum yac_edge_type const * edge_types,
  size_t const * edge_indices, size_t num_edges) {

  if (edge_types == NULL) return YAC_GREAT_CIRCLE_CELL;

  int edge_type_flag = 0;

  if (num_edges == 4) {

    enum yac_edge_type temp_edges[4] =
      {edge_types[edge_indices[0]],
       edge_types[edge_indices[1]],
       edge_types[edge_indices[2]],
       edge_types[edge_indices[3]]};
    if (IS_LON_LAT(temp_edges[0]) &&
        IS_LON_LAT(temp_edges[1]) &&
        (temp_edges[0] == temp_edges[2]) &&
        (temp_edges[1] == temp_edges[3]) &&
        (temp_edges[0] != temp_edges[1]))
      return YAC_LON_LAT_CELL;
  }

  for (size_t i = 0; i < num_edges; ++i)
    edge_type_flag |= (1 << (edge_types[edge_indices[i]] == YAC_LAT_CIRCLE_EDGE));

  enum yac_cell_type cell_types[4] =
    {YAC_GREAT_CIRCLE_CELL,
     YAC_GREAT_CIRCLE_CELL,
     YAC_LAT_CELL,
     YAC_MIXED_CELL};
  return cell_types[edge_type_flag];
}
#undef IS_LON_LAT
#undef IS_GC

static size_t do_search_avg (struct interp_method * method,
                             struct yac_interp_grid * interp_grid,
                             size_t * tgt_points, size_t count,
                             struct yac_interp_weights * weights) {

  struct interp_method_avg * method_avg = (struct interp_method_avg *)method;
  func_compute_weights compute_weights = method_avg->compute_weights;

  YAC_ASSERT(
    yac_interp_grid_get_num_src_fields(interp_grid) == 1,
    "ERROR(do_search_avg): invalid number of source fields")

  enum yac_location src_field_location =
    yac_interp_grid_get_src_field_location(interp_grid, 0);

  YAC_ASSERT(
    (src_field_location == YAC_LOC_CORNER) ||
    (src_field_location == YAC_LOC_CELL),
    "ERROR(do_search_avg): unsupported source field location type")

  // get coordinates of target points
  yac_coordinate_pointer tgt_coords = xmalloc(count * sizeof(*tgt_coords));
  yac_interp_grid_get_tgt_coordinates(
    interp_grid, tgt_points, count, tgt_coords);

  size_t * size_t_buffer = xmalloc(2 * count * sizeof(*size_t_buffer));
  size_t * src_field_cells = size_t_buffer;
  size_t * reorder_idx = size_t_buffer + count;

  // get matching source cells for all target points
  yac_interp_grid_do_points_search(
    interp_grid, tgt_coords, count, src_field_cells);

  size_t const * src_field_cell_to_vertex;
  size_t const * src_field_cell_to_vertex_offsets;
  size_t const * src_field_cell_to_edge;
  size_t const * src_field_cell_to_edge_offsets;
  int const * src_field_num_vertices_per_cell;
  enum yac_edge_type const * src_edge_types;

  // if the source field is location at cell points
  if (src_field_location == YAC_LOC_CELL) {

    // generate auxiliary grid for all search result cells
    yac_interp_grid_get_aux_grid_src(
      interp_grid, src_field_cells, count,
      (size_t**)&src_field_cell_to_vertex,
      (size_t**)&src_field_cell_to_vertex_offsets,
      (int**)&src_field_num_vertices_per_cell);

    struct yac_const_basic_grid_data * src_grid_data =
      yac_interp_grid_get_basic_grid_data_src(interp_grid);
    yac_const_coordinate_pointer src_field_coordinates =
      yac_interp_grid_get_src_field_coords(interp_grid, 0);

    struct yac_grid_cell cell_buffer;
    yac_init_grid_cell(&cell_buffer);

    // for all target point
    for (size_t i = 0; i < count; ++i) {

      size_t curr_src_cell = src_field_cells[i];
      if (curr_src_cell == SIZE_MAX) continue;

      double * curr_tgt_coord = tgt_coords[i];
      size_t const * curr_vertices =
        src_grid_data->cell_to_vertex +
        src_grid_data->cell_to_vertex_offsets[curr_src_cell];
      size_t curr_num_vertices =
        src_grid_data->num_vertices_per_cell[curr_src_cell];

      size_t result_src_field_cell = SIZE_MAX;

      // for all auxiliary of the current source result cell
      for (size_t j = 0; j < curr_num_vertices; ++j) {

        size_t src_field_cell = curr_vertices[j];
        size_t src_field_cell_size =
          (size_t)src_field_num_vertices_per_cell[src_field_cell];

        // check whether the target point is in the current source field cell
        if (check_src_field_cell(
              curr_tgt_coord, src_field_cell, src_field_cell_size,
              src_field_cell_to_vertex, src_field_cell_to_vertex_offsets,
              src_field_coordinates, &cell_buffer)) {
          result_src_field_cell = src_field_cell;
          break;
        }
      }

      src_field_cells[i] = result_src_field_cell;
    }
    yac_free_grid_cell(&cell_buffer);

    src_edge_types = NULL;
    src_field_cell_to_edge = NULL;
    src_field_cell_to_edge_offsets = NULL;
  } else {

    struct yac_const_basic_grid_data * grid_data =
      yac_interp_grid_get_basic_grid_data_src(interp_grid);
    src_field_cell_to_vertex = grid_data->cell_to_vertex;
    src_field_cell_to_vertex_offsets = grid_data->cell_to_vertex_offsets;
    src_field_num_vertices_per_cell = grid_data->num_vertices_per_cell;
    src_edge_types = grid_data->edge_type;
    src_field_cell_to_edge = grid_data->cell_to_edge;
    src_field_cell_to_edge_offsets = grid_data->cell_to_edge_offsets;
  }
  yac_const_coordinate_pointer src_field_coordinates =
    yac_interp_grid_get_src_field_coords(interp_grid, 0);
  const_yac_int_pointer src_global_ids =
    yac_interp_grid_get_src_field_global_ids(interp_grid, 0);

  // sort target points, for which we found a source cell, to the beginning of
  // the array
  for (size_t i = 0; i < count; ++i) reorder_idx[i] = i;
  yac_quicksort_index_size_t_size_t(src_field_cells, count, reorder_idx);

  size_t result_count = 0;
  for (result_count = 0; result_count < count; ++result_count)
    if (src_field_cells[result_count] == SIZE_MAX) break;

  size_t total_num_weights = 0;
  size_t max_num_vertices_per_cell = 0;
  size_t * num_weights_per_tgt =
    xmalloc(result_count * sizeof(*num_weights_per_tgt));

  // get the number of vertices per target cell and the total number of vertices
  // required for the interpolation (may contain duplicated vertices)
  for (size_t i = 0; i < result_count; ++i) {
    size_t curr_num_vertices =
      (size_t)(src_field_num_vertices_per_cell[src_field_cells[i]]);
    num_weights_per_tgt[i] = curr_num_vertices;
    total_num_weights += curr_num_vertices;
    if (curr_num_vertices > max_num_vertices_per_cell)
      max_num_vertices_per_cell = curr_num_vertices;
  }
  if (src_field_location == YAC_LOC_CELL)
    free((void*)src_field_num_vertices_per_cell);

  const_int_pointer src_field_mask =
    yac_interp_grid_get_src_field_mask(interp_grid, 0);

  yac_coordinate_pointer src_coord_buffer =
    xmalloc(max_num_vertices_per_cell * sizeof(*src_coord_buffer));
  int * mask_buffer =
    (src_field_mask != NULL)?
      xmalloc(max_num_vertices_per_cell * sizeof(*mask_buffer)):NULL;
  double * w = xmalloc(total_num_weights * sizeof(*w));
  size_t * src_points = xmalloc(total_num_weights * sizeof(*src_points));

  // for each target point, extract relevant source data and compute the weights
  // based on that
  total_num_weights = 0;
  for (size_t i = 0; i < result_count;) {

    size_t curr_num_vertices = num_weights_per_tgt[i];

    const_size_t_pointer curr_cell_to_vertex =
      src_field_cell_to_vertex +
      src_field_cell_to_vertex_offsets[src_field_cells[i]];
    const_size_t_pointer curr_cell_to_edge =
      (src_edge_types != NULL)?
        (src_field_cell_to_edge +
         src_field_cell_to_edge_offsets[src_field_cells[i]]):NULL;
    double * curr_weights = w + total_num_weights;

    // get the index of the source cell corner with the lowest global id
    // (this is being used, that the order in which the source corners are
    //  processes is decomposition independent)
    size_t lowest_global_id_idx = 0;
    {
      yac_int lowest_global_id = src_global_ids[curr_cell_to_vertex[0]];
      for (size_t j = 1; j < curr_num_vertices; ++j) {
        if (src_global_ids[curr_cell_to_vertex[j]] < lowest_global_id) {
          lowest_global_id = src_global_ids[curr_cell_to_vertex[j]];
          lowest_global_id_idx = j;
        }
      }
    }

    // get the mask and the source coordinates of the source vertices required
    // for the current target point
    for (size_t j = 0, l = lowest_global_id_idx; j < curr_num_vertices;
         ++j, ++total_num_weights, ++l) {

      if (l == curr_num_vertices) l = 0;

      size_t curr_vertex_idx = curr_cell_to_vertex[l];
      src_points[total_num_weights] = curr_vertex_idx;
      for (size_t k = 0; k < 3; ++k)
        src_coord_buffer[j][k] = src_field_coordinates[curr_vertex_idx][k];
      if (src_field_mask != NULL)
        mask_buffer[j] = src_field_mask[curr_vertex_idx];
    }

    enum yac_cell_type cell_type =
      determine_cell_type(
        src_edge_types, curr_cell_to_edge, curr_num_vertices);

    // compute the weights
    if (compute_weights(tgt_coords[reorder_idx[i]], curr_num_vertices,
                        (yac_const_coordinate_pointer)src_coord_buffer,
                        mask_buffer, curr_weights, cell_type)) {

      // if weight computation was successful
      num_weights_per_tgt[i] = curr_num_vertices;
      ++i;
    } else {

      // if no weights could be computed
      total_num_weights -= curr_num_vertices;
      result_count--;
      src_field_cells[i] = src_field_cells[result_count];
      size_t temp_reorder_idx = reorder_idx[i];
      reorder_idx[i] = reorder_idx[result_count];
      reorder_idx[result_count] = temp_reorder_idx;
    }
  }
  if (src_field_location == YAC_LOC_CELL) {
    free((void*)src_field_cell_to_vertex);
    free((void*)src_field_cell_to_vertex_offsets);
  }

  // move the non-interpolated target points to the end
  for (size_t i = 0; i < count; ++i) src_field_cells[reorder_idx[i]] = i;
  yac_quicksort_index_size_t_size_t(src_field_cells, count, tgt_points);

  free(tgt_coords);
  free(mask_buffer);
  free(src_coord_buffer);
  free(size_t_buffer);

  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(
        interp_grid, tgt_points, result_count),
    .count = result_count};
  struct remote_point * srcs =
    yac_interp_grid_get_src_remote_points(
      interp_grid, 0, src_points, total_num_weights);

  // store weights
  yac_interp_weights_add_wsum(
    weights, &tgts, num_weights_per_tgt, srcs, w);

  free(tgts.data);
  free(srcs);
  free(w);
  free(num_weights_per_tgt);
  free(src_points);

  return result_count;
}

static int compute_weights_avg_yes(
  double tgt_coords[3], size_t num_vertices,
  yac_const_coordinate_pointer src_coords, int * src_mask, double * weights,
  enum yac_cell_type cell_type) {

  UNUSED(tgt_coords);
  UNUSED(cell_type);
  UNUSED(src_coords);

  size_t num_unmasked_points = 0;

  // if we have a source mask, check whether any points is masked
  if (src_mask != NULL) {
    for (size_t i = 0; i < num_vertices; ++i)
      if (src_mask[i]) num_unmasked_points++;
    if (num_unmasked_points == 0) return 0;

    double weight = 1.0 / (double)num_unmasked_points;
    for (size_t i = 0; i < num_vertices; ++i)
      weights[i] = (src_mask[i])?weight:0.0;

  } else {
    double weight = 1.0 / (double)num_vertices;
    for (size_t i = 0; i < num_vertices; ++i) weights[i] = weight;
  }
  return 1;
}

static int compute_weights_avg_no(
  double tgt_coords[3], size_t num_vertices,
  yac_const_coordinate_pointer src_coords, int * src_mask, double * weights,
  enum yac_cell_type cell_type) {

  UNUSED(tgt_coords);
  UNUSED(cell_type);
  UNUSED(src_coords);

  // if we have a source mask, check whether any points is masked
  if (src_mask != NULL)
    for (size_t i = 0; i < num_vertices; ++i)
      if (!src_mask[i]) return 0;

  double weight = 1.0 / (double)num_vertices;
  for (size_t i = 0; i < num_vertices; ++i) weights[i] = weight;

  return 1;
}

static int compute_weights_dist_yes(
  double tgt_coords[3], size_t num_vertices,
  yac_const_coordinate_pointer src_coords, int * src_mask, double * weights,
  enum yac_cell_type cell_type) {

  UNUSED(cell_type);

  // if there is a source mask, check if there are unmasked vertices
  if (src_mask != NULL) {
    int has_unmasked = 0;
    for (size_t i = 0; i < num_vertices; ++i) has_unmasked |= src_mask[i];
    if (!has_unmasked) return 0;
  }

  if (src_mask != NULL) {
    for (size_t i = 0; i < num_vertices; ++i) {

      if (src_mask[i]) {

        double distance =
          get_vector_angle(tgt_coords, (double*)(src_coords[i]));

        // if the target and source point are nearly identical
        if (distance < yac_angle_tol) {
          for (size_t j = 0; j < num_vertices; ++j) weights[j] = 0.0;
          weights[i] = 1.0;
          return 1;
        }

        weights[i] = 1.0 / distance;
      } else {
        weights[i] = 0.0;
      }
    }
  } else {
    for (size_t i = 0; i < num_vertices; ++i) {

      double distance =
        get_vector_angle(tgt_coords, (double*)(src_coords[i]));

      // if the target and source point are nearly identical
      if (distance < yac_angle_tol) {
        for (size_t j = 0; j < num_vertices; ++j) weights[j] = 0.0;
        weights[i] = 1.0;
        return 1;
      }

      weights[i] = 1.0 / distance;
    }
  }

  // compute scaling factor for the weights
  double inv_distance_sum = 0.0;
  for (size_t i = 0; i < num_vertices; ++i)
    inv_distance_sum += weights[i];
  double scale = 1.0 / inv_distance_sum;

  for (size_t i = 0; i < num_vertices; ++i) weights[i] *= scale;

  return 1;
}

static int compute_weights_dist_no(
  double tgt_coords[3], size_t num_vertices,
  yac_const_coordinate_pointer src_coords, int * src_mask, double * weights,
  enum yac_cell_type cell_type) {

  UNUSED(cell_type);

  for (size_t i = 0; i < num_vertices; ++i) {

    double distance =
      get_vector_angle(tgt_coords, (double*)(src_coords[i]));

    // if the target and source point are nearly identical
    if (distance < yac_angle_tol) {
      if ((src_mask != NULL) && !src_mask[i]) return 0;
      for (size_t j = 0; j < num_vertices; ++j) weights[j] = 0.0;
      weights[i] = 1.0;
      return 1;
    }

    weights[i] = 1.0 / distance;
  }

  // if there is a source mask, check if there are masked vertices
  // (we do this here because there may be a matching source point
  //  that is not masked)
  if (src_mask != NULL)
    for (size_t i = 0; i < num_vertices; ++i) if(!src_mask[i]) return 0;

  // compute scaling factor for the weights
  double inv_distance_sum = 0.0;
  for (size_t i = 0; i < num_vertices; ++i)
    inv_distance_sum += weights[i];
  double scale = 1.0 / inv_distance_sum;

  for (size_t i = 0; i < num_vertices; ++i) weights[i] *= scale;

  return 1;
}

// compute the spherical barycentric coordinates of the given vertices with
// respect to the given triangle
static void compute_barycentric_coords(
  double * barycentric_coords, size_t triangle_indices[3],
  yac_const_coordinate_pointer grid_coords) {

  double A[3][3];
  lapack_int n = 3, nrhs = 1, lda = n, ldx = n, ipiv[3];
  for (int i = 0; i < 3; ++i)
    memcpy(A[i], grid_coords[triangle_indices[i]], sizeof(*grid_coords));

  // for a vertex v the spherical barycentric coordinates b are defined as
  // follows
  // A * b = v
  // where: A is the matrix consisting of the vertex coordinates of the three
  // corners of the triangle
  // we compute b by solving this linear system using LAPACK
  YAC_ASSERT(
    !LAPACKE_dgesv(
      LAPACK_COL_MAJOR, n, nrhs, &A[0][0], lda,
      ipiv, barycentric_coords, ldx),
    "ERROR: internal error (could not solve linear 3x3 system)")
}

// returns 0 if edge_direction_vec points south; 1 otherwise
static inline int get_lat_edge_ordering(double const * a, double const * b) {

  return (a[0] * b[1] - a[1] * b[0]) < 0.0;
}

static void get_cell_lon_lat_bounds(
  yac_const_coordinate_pointer coords,
  double lon[2], double lat[2], int reorder[4]) {

  // determine whether the lat neighbour to vertex 0 is at position 1 or 3
  int lat_neigh_offset =
    (fabs(coords[0][2] - coords[1][2]) <
     fabs(coords[0][2] - coords[3][2]))?1:3;
  // get the index of the coord which is closes to the pole
  // (the one with the lowest absolut z coordinate)
  int closest_to_equator_idx =
    (fabs(coords[0][2]) < fabs(coords[1][2]))?0:2;
  int upper_edge_idx = ((coords[0][2]) > coords[2][2])?0:2;
  int lower_edge_idx = upper_edge_idx^2;
  int upper_edge_is_closer_to_pole =
    closest_to_equator_idx == upper_edge_idx;
  int closest_to_equator_edge_ordering =
    get_lat_edge_ordering(
      coords[closest_to_equator_idx],
      coords[(closest_to_equator_idx+lat_neigh_offset)%4]);
  int cell_ordering =
    closest_to_equator_edge_ordering ^ upper_edge_is_closer_to_pole;

  XYZtoLL(
    coords[closest_to_equator_idx], &lon[closest_to_equator_edge_ordering],
    &lat[upper_edge_is_closer_to_pole]);
  XYZtoLL(
    coords[(closest_to_equator_idx + lat_neigh_offset)%4],
    &lon[closest_to_equator_edge_ordering^1],
    &lat[upper_edge_is_closer_to_pole]);
  lat[!upper_edge_is_closer_to_pole] =
    M_PI_2 - acos(coords[closest_to_equator_idx ^ 2][2]);

  reorder[cell_ordering] = lower_edge_idx;
  reorder[cell_ordering^1] = (lower_edge_idx+lat_neigh_offset)%4;
  reorder[2+cell_ordering] = upper_edge_idx;
  reorder[2+(cell_ordering^1)] = (upper_edge_idx+lat_neigh_offset)%4;

  if (lon[1] < lon[0]) lon[1] += 2.0 * M_PI;
}

static void get_point_lon_lat(
  double point_coord[3], double cell_lon[2], double cell_lat[2],
  double * point_lon, double * point_lat) {

  UNUSED(cell_lat);

  double lon, lat;
  XYZtoLL(point_coord, &lon, &lat);

  // adjust target point lon to source cell lon
  if (lon < cell_lon[0]) {
    while (fabs(cell_lon[0] - lon) > M_PI) lon += 2.0 * M_PI;
  } else {
    while (fabs(lon - cell_lon[1]) > M_PI) lon -= 2.0 * M_PI;
  }

  *point_lon = lon;
  *point_lat = lat;
}

static int determine_triangle_idx(
  double point_lon, double point_lat, double cell_lon[2], double cell_lat[2]) {

  return ((cell_lat[1] - cell_lat[0]) * (point_lat - cell_lat[0]) -
          (cell_lat[1] - cell_lat[0]) * (point_lon - cell_lon[0])) > 0.0;
}

static void compute_weights_bary_reg(
  double tgt_coords[3],
  yac_const_coordinate_pointer src_coords, double * weights) {

  // get lon lat bounds of the cell and the target point
  double src_lon[2], src_lat[2], tgt_lon, tgt_lat;
  int src_reorder[4];
  get_cell_lon_lat_bounds(src_coords, src_lon, src_lat, src_reorder);
  get_point_lon_lat(tgt_coords, src_lon, src_lat, &tgt_lon, &tgt_lat);

  int triangle_idx =
    determine_triangle_idx(tgt_lon, tgt_lat, src_lon, src_lat);

  if (triangle_idx) {

    double w[2];
    w[0] = (tgt_lat - src_lat[1]) / (src_lat[0] - src_lat[1]);
    w[1] = ((src_lat[1] - src_lat[0]) * (tgt_lon - src_lon[0])) /
           ((src_lon[0] - src_lon[1]) * (src_lat[0] - src_lat[1]));

    weights[src_reorder[0]] = w[0];
    weights[src_reorder[2]] = w[1];
    weights[src_reorder[3]] = 1.0 - (w[0] + w[1]);
    weights[src_reorder[1]] = 0.0;
  } else {

    double w[2];
    w[0] = (tgt_lon - src_lon[1]) / (src_lon[0] - src_lon[1]);
    w[1] = ((src_lat[1] - src_lat[0]) * (tgt_lon - src_lon[1]) +
            (src_lon[0] - src_lon[1]) * (tgt_lat - src_lat[1])) /
           ((src_lat[0] - src_lat[1]) * (src_lon[0] - src_lon[1]));

    weights[src_reorder[0]] = w[0];
    weights[src_reorder[1]] = w[1];
    weights[src_reorder[2]] = 1.0 - (w[0] + w[1]);
    weights[src_reorder[3]] = 0.0;
  }
}

static int compute_weights_bary_yes(
  double tgt_coords[3], size_t num_vertices,
  yac_const_coordinate_pointer src_coords, int * src_mask, double * weights,
  enum yac_cell_type cell_type) {

  // if all source points are masked
  if (src_mask != NULL) {
    size_t i;
    for (i = 0; i < num_vertices; ++i) if (src_mask[i] != 0) break;
    if (i == num_vertices) return 0;
  }

  double const tol = 1e-9;

  YAC_ASSERT(
    (cell_type == YAC_LON_LAT_CELL) ||
    (cell_type == YAC_GREAT_CIRCLE_CELL),
    "ERROR(compute_weights_bary_yes): barycentric coordinates "
    "are only supported for great circle edge cells and lon lat cells")

  if (cell_type == YAC_LON_LAT_CELL) {

    // compute weights
    compute_weights_bary_reg(tgt_coords, src_coords, weights);

    // apply mask handling
    if (src_mask != NULL) {
      double weight_sum = 1.0;

      // check if any of the source points is masked
      for (int i = 0; i < 4; ++i) if (!src_mask[i]) {
        weight_sum -= weights[i];
        weights[i] = 0.0;
      }

      // if all relevant points are masked
      if (weight_sum < tol) return 0;

      // if only some of the relevant points are masked
      if (weight_sum < 1.0) {
        double scale = 1.0 / weight_sum;
        for (int i = 0; i < 4; ++i) weights[i] *= scale;
      }
    }

    return 1;

  } else {

    size_t corner_indices[num_vertices];
    size_t triangle_indices[num_vertices-2][3];

    for (size_t i = 0; i < num_vertices; ++i) corner_indices[i] = i;

    // triangulate source polygon
    if (num_vertices > 3) {
      yac_triangulate_cell_indices(
        corner_indices, num_vertices, 0, triangle_indices);
    } else {
      for (size_t i = 0; i < 3; ++i) triangle_indices[0][i] = i;
    }

    // find the best matching triangle
    double min_barycentric_coord = -DBL_MAX;
    double barycentric_coords[3];
    size_t match_index = SIZE_MAX;

    for (size_t i = 0; i < num_vertices - 2; ++i) {

      double temp_barycentric_coords[3];
      memcpy(temp_barycentric_coords, tgt_coords, 3 * sizeof(double));

      // compute barycentric coordinates for the current triangle
      compute_barycentric_coords(
        temp_barycentric_coords, triangle_indices[i], src_coords);

      double curr_min_barycentric_coord =
        MIN(temp_barycentric_coords[0],
            MIN(temp_barycentric_coords[1],
                temp_barycentric_coords[2]));

      if (curr_min_barycentric_coord > min_barycentric_coord) {
        min_barycentric_coord = curr_min_barycentric_coord;
        match_index = i;
        for (int j = 0; j < 3; ++j)
          barycentric_coords[j] = MAX(0.0, temp_barycentric_coords[j]);
      }
    }

    // initialise weights
    if (num_vertices > 3)
      for (size_t j = 0; j < num_vertices; ++j) weights[j] = 0.0;

    // apply mask
    if (src_mask != NULL)
      for (int j = 0; j < 3; ++j)
        if (src_mask[triangle_indices[match_index][j]] == 0)
          barycentric_coords[j] = 0.0;

    double weight_sum = 0.0;
    for (int j = 0; j < 3; ++j) {
      if (barycentric_coords[j] < 1e-4) barycentric_coords[j] = 0.0;
      else weight_sum += barycentric_coords[j];
    }

    if (weight_sum < tol) return 0;

    double scale = 1.0 / weight_sum;

    for (int j = 0; j < 3; ++j)
      weights[triangle_indices[match_index][j]] =
        barycentric_coords[j] * scale;

    return 1;
  }
}

static int compute_weights_bary_no(
  double tgt_coords[3], size_t num_vertices,
  yac_const_coordinate_pointer src_coords, int * src_mask, double * weights,
  enum yac_cell_type cell_type) {

  // if all source points are masked
  if (src_mask != NULL) {
    size_t i;
    for (i = 0; i < num_vertices; ++i) if (src_mask[i] != 0) break;
    if (i == num_vertices) return 0;
  }

  double const tol = 1e-9;

  YAC_ASSERT(
    (cell_type == YAC_LON_LAT_CELL) || (cell_type == YAC_GREAT_CIRCLE_CELL),
    "ERROR(compute_weights_bary_no): barycentric coordinates "
    "are only supported for great circle edge cells and lon lat cells")

  if (cell_type == YAC_LON_LAT_CELL) {

    // compute weights
    compute_weights_bary_reg(tgt_coords, src_coords, weights);

    for (int i = 0; i < 4; ++i) if (weights[i] < tol) weights[i] = 0.0;

    // check if any of the source points is masked
    if (src_mask != NULL)
      for (int i = 0; i < 4; ++i)
        if (!src_mask[i] && (weights[i] != 0.0))
          return 0;
    return 1;

  } else {

    size_t corner_indices[num_vertices];
    size_t triangle_indices[num_vertices-2][3];

    for (size_t i = 0; i < num_vertices; ++i) corner_indices[i] = i;

    // triangulate source polygon
    if (num_vertices > 3) {
      yac_triangulate_cell_indices(
        corner_indices, num_vertices, 0, triangle_indices);
    } else {
      for (size_t i = 0; i < 3; ++i) triangle_indices[0][i] = i;
    }

    // find the best matching triangle
    double min_barycentric_coord = -DBL_MAX;
    double barycentric_coords[3];
    size_t match_index = SIZE_MAX;

    for (size_t i = 0; i < num_vertices - 2; ++i) {

      double temp_barycentric_coords[3];
      memcpy(temp_barycentric_coords, tgt_coords, 3 * sizeof(double));

      // compute barycentric coordinates for the current triangle
      compute_barycentric_coords(
        temp_barycentric_coords, triangle_indices[i], src_coords);

      double curr_min_barycentric_coord =
        MIN(temp_barycentric_coords[0],
            MIN(temp_barycentric_coords[1],
                temp_barycentric_coords[2]));

      if (curr_min_barycentric_coord > min_barycentric_coord) {
        min_barycentric_coord = curr_min_barycentric_coord;
        match_index = i;
        for (int j = 0; j < 3; ++j)
          barycentric_coords[j] = MAX(0.0, temp_barycentric_coords[j]);
      }
    }

    // initialise weights
    if (num_vertices > 3)
      for (size_t j = 0; j < num_vertices; ++j) weights[j] = 0.0;

    // apply mask
    if (src_mask != NULL)
      for (int j = 0; j < 3; ++j)
        if (!src_mask[triangle_indices[match_index][j]] &&
            (barycentric_coords[j] >= 1e-4)) return 0;

    double weight_sum = 0.0;
    for (int j = 0; j < 3; ++j) {
      if (barycentric_coords[j] < 1e-4) barycentric_coords[j] = 0.0;
      else weight_sum += barycentric_coords[j];
    }

    if (weight_sum < tol) return 0;

    double scale = 1.0 / weight_sum;

    for (int j = 0; j < 3; ++j)
      weights[triangle_indices[match_index][j]] =
        barycentric_coords[j] * scale;

    return 1;
  }
}

static func_compute_weights select_compute_weight_routine(
  enum yac_interp_avg_weight_type weight_type, int partial_coverage) {

  YAC_ASSERT(
    (weight_type == YAC_INTERP_AVG_ARITHMETIC) ||
    (weight_type == YAC_INTERP_AVG_DIST) ||
    (weight_type == YAC_INTERP_AVG_BARY),
    "ERROR(select_compute_weight_routine): invalid weight type")

  switch(weight_type) {
    case(YAC_INTERP_AVG_ARITHMETIC):
      return (partial_coverage)?compute_weights_avg_yes:compute_weights_avg_no;
    case(YAC_INTERP_AVG_DIST):
      return (partial_coverage)?compute_weights_dist_yes:compute_weights_dist_no;
    case(YAC_INTERP_AVG_BARY):
    default:
      return (partial_coverage)?compute_weights_bary_yes:compute_weights_bary_no;
  };
}

struct interp_method * yac_interp_method_avg_new(
  enum yac_interp_avg_weight_type weight_type,
  int partial_coverage) {

  struct interp_method_avg * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_avg_vtable;
  method->compute_weights =
    select_compute_weight_routine(weight_type, partial_coverage);

  return (struct interp_method*)method;
}

static void delete_avg(struct interp_method * method) {
  free(method);
}
