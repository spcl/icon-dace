// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <string.h>

#include "interp_method_internal.h"
#include "interp_method_ncc.h"

static size_t do_search_ncc(struct interp_method * method,
                            struct yac_interp_grid * interp_grid,
                            size_t * tgt_points, size_t count,
                            struct yac_interp_weights * weights);
static void delete_ncc(struct interp_method * method);

typedef int (*func_compute_weights)(
  double[3], size_t, yac_const_coordinate_pointer, int *, double *);

static struct interp_method_vtable
  interp_method_ncc_vtable = {
    .do_search = do_search_ncc,
    .delete = delete_ncc};

struct interp_method_ncc {

  struct interp_method_vtable * vtable;
  func_compute_weights compute_weights;
};

// determines for a given coordinate which corner of a specific
// cell is the closest one
static size_t get_closest_src_corner(
  struct yac_const_basic_grid_data * grid_data, size_t cell,
  double const coord[3]) {

  size_t const * cell_to_vertex =
    grid_data->cell_to_vertex + grid_data->cell_to_vertex_offsets[cell];
  yac_const_coordinate_pointer corner_coords =
    grid_data->vertex_coordinates;
  size_t num_vertices = grid_data->num_vertices_per_cell[cell];
  yac_int const * vertex_ids = grid_data->ids[YAC_LOC_CORNER];

  double min_distance = DBL_MAX;
  size_t min_vertex = SIZE_MAX;
  yac_int min_vertex_id = XT_INT_MAX;

  for (size_t i = 0; i < num_vertices; ++i) {
    size_t curr_vertex = cell_to_vertex[i];
    double curr_distance =
      get_vector_angle(coord, corner_coords[curr_vertex]);
    if (curr_distance < min_distance) {
      min_distance = curr_distance;
      min_vertex = curr_vertex;
      min_vertex_id = vertex_ids[curr_vertex];
    } else if ((curr_distance == min_distance) &&
               (vertex_ids[curr_vertex] < min_vertex_id)) {

      min_vertex = curr_vertex;
      min_vertex_id = vertex_ids[curr_vertex];
    }
  }
  return min_vertex;
}

static size_t do_search_ncc (struct interp_method * method,
                             struct yac_interp_grid * interp_grid,
                             size_t * tgt_points, size_t count,
                             struct yac_interp_weights * weights) {

  struct interp_method_ncc * method_ncc = (struct interp_method_ncc *)method;

  YAC_ASSERT(
    yac_interp_grid_get_num_src_fields(interp_grid) == 1,
    "ERROR(do_search_ncc): invalid number of source fields")

  YAC_ASSERT(
    yac_interp_grid_get_src_field_location(interp_grid, 0) == YAC_LOC_CELL,
    "ERROR(do_search_ncc): unsupported source field location type "
    "(has to be CELL)")

  // get coordinates of target points
  yac_coordinate_pointer tgt_coords = xmalloc(count * sizeof(*tgt_coords));
  yac_interp_grid_get_tgt_coordinates(
    interp_grid, tgt_points, count, tgt_coords);

  size_t * size_t_buffer = xmalloc(3 * count * sizeof(*size_t_buffer));
  size_t * src_cells = size_t_buffer;
  size_t * src_corners = src_cells;
  size_t * reorder_idx = size_t_buffer + count;
  size_t * interp_flag = size_t_buffer + 2 * count;

  // get matching source cells for all target points
  yac_interp_grid_do_points_search(
    interp_grid, tgt_coords, count, src_cells);

  // sort target points, for which we found a source cell, to the beginning of
  // the array
  for (size_t i = 0; i < count; ++i) reorder_idx[i] = i;
  yac_quicksort_index_size_t_size_t(src_cells, count, reorder_idx);

  size_t result_count = 0;
  for (result_count = 0; result_count < count; ++result_count)
    if (src_cells[result_count] == SIZE_MAX) break;
  for (size_t i = result_count; i < count; ++i)
    interp_flag[reorder_idx[i]] = SIZE_MAX;

  // for each target point find the closest corner of its source cell
  {
    struct yac_const_basic_grid_data * src_grid_data =
      yac_interp_grid_get_basic_grid_data_src(interp_grid);
    for (size_t i = 0; i < result_count; ++i)
      src_corners[i] =
        get_closest_src_corner(
          src_grid_data, src_cells[i], tgt_coords[reorder_idx[i]]);
  }

  // sort source corners
  yac_quicksort_index_size_t_size_t(
    src_corners, result_count, reorder_idx);

  // count number of unique source corners
  size_t num_unique_src_corners = 0;
  for (size_t i = 0, prev_src_corner = SIZE_MAX; i < result_count; ++i) {
    size_t curr_src_corner = src_corners[i];
    if (prev_src_corner != curr_src_corner) {
      prev_src_corner = curr_src_corner;
      ++num_unique_src_corners;
    }
  }

  // count number of target points per unique source corner and
  // remove duplicated corners
  size_t * num_tgt_per_corner =
    xcalloc(num_unique_src_corners, sizeof(*num_tgt_per_corner));
  num_unique_src_corners = 0;
  for (size_t i = 0, prev_src_corner = SIZE_MAX; i < result_count; ++i) {
    size_t curr_src_corner = src_corners[i];
    if (prev_src_corner != curr_src_corner) {
      prev_src_corner = curr_src_corner;
      src_corners[num_unique_src_corners] = curr_src_corner;
      ++num_unique_src_corners;
    }
    num_tgt_per_corner[num_unique_src_corners-1]++;
  }

  // get all source cells associated with each unique source corner;
  // results cells for each corner are sorted by their global id
  size_t * src_corner_cells = NULL;
  size_t * num_cells_per_corner =
    xmalloc(num_unique_src_corners * sizeof(*num_cells_per_corner));
  yac_interp_grid_get_src_corner_cells(
    interp_grid, src_corners, num_unique_src_corners, &src_corner_cells,
    num_cells_per_corner);

  // get maximum number of source cells per target
  size_t max_num_cell_per_corner = 0;
  size_t total_num_weights = 0;
  for (size_t i = 0; i < num_unique_src_corners; ++i) {
    if (max_num_cell_per_corner < num_cells_per_corner[i])
      max_num_cell_per_corner = num_cells_per_corner[i];
    total_num_weights += num_cells_per_corner[i] * num_tgt_per_corner[i];
  }

  yac_const_coordinate_pointer src_field_coordinates =
    yac_interp_grid_get_src_field_coords(interp_grid, 0);
  const_int_pointer src_field_mask =
    yac_interp_grid_get_src_field_mask(interp_grid, 0);

  yac_coordinate_pointer src_coord_buffer =
    xmalloc(max_num_cell_per_corner * sizeof(*src_coord_buffer));
  int * mask_buffer =
    (src_field_mask != NULL)?
      xmalloc(max_num_cell_per_corner * sizeof(*mask_buffer)):NULL;
  double * w = xmalloc(total_num_weights * sizeof(*w));
  size_t * src_points = xmalloc(total_num_weights * sizeof(*src_points));
  size_t * num_weights_per_tgt =
    xmalloc(result_count * sizeof(*num_weights_per_tgt));

  // for each target point, extract relevant source data and compute the weights
  // based on that
  func_compute_weights compute_weights = method_ncc->compute_weights;
  total_num_weights = 0;
  result_count = 0;
  for (size_t i = 0, l = 0, src_corner_cell_offset = 0;
       i < num_unique_src_corners; ++i) {

    size_t curr_num_cells = num_cells_per_corner[i];

    const_size_t_pointer curr_src_corner_cells =
      src_corner_cells + src_corner_cell_offset;
    src_corner_cell_offset += curr_num_cells;

    // get field coordinates and mask values for current source cells
    for (size_t j = 0; j < curr_num_cells; ++j) {
      size_t const curr_src_cell = curr_src_corner_cells[j];
      for (size_t k = 0; k < 3; ++k)
        src_coord_buffer[j][k] = src_field_coordinates[curr_src_cell][k];
      if (src_field_mask != NULL)
        mask_buffer[j] = src_field_mask[curr_src_cell];
    }

    for (size_t j = 0; j < num_tgt_per_corner[i]; ++j, ++l) {

      // compute weights
      if (compute_weights(
            tgt_coords[reorder_idx[l]], curr_num_cells,
            (yac_const_coordinate_pointer)src_coord_buffer,
            mask_buffer, w + total_num_weights)) {

        interp_flag[reorder_idx[l]] = result_count;
        memcpy(
          src_points + total_num_weights, curr_src_corner_cells,
          curr_num_cells * sizeof(*src_points));
        num_weights_per_tgt[result_count] = curr_num_cells;
        total_num_weights += curr_num_cells;
        result_count++;
      } else {
        interp_flag[reorder_idx[l]] = SIZE_MAX;
      }
    }
  }

  free(src_corner_cells);
  free(num_tgt_per_corner);
  free(num_cells_per_corner);

  // move target points for which the interpolation failed to the end of
  // the array while keeping the origin order for the remaining target points
  yac_quicksort_index_size_t_size_t(interp_flag, count, tgt_points);

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
  double tgt_coords[3], size_t num_src,
  yac_const_coordinate_pointer src_coords, int * src_mask, double * weights) {

  UNUSED(tgt_coords);
  UNUSED(src_coords);

  size_t num_unmasked_points = 0;

  // if we have a source mask, check whether any points is masked
  if (src_mask != NULL) {
    for (size_t i = 0; i < num_src; ++i)
      if (src_mask[i]) num_unmasked_points++;
    if (num_unmasked_points == 0) return 0;

    double weight = 1.0 / (double)num_unmasked_points;
    for (size_t i = 0; i < num_src; ++i)
      weights[i] = (src_mask[i])?weight:0.0;

  } else {
    double weight = 1.0 / (double)num_src;
    for (size_t i = 0; i < num_src; ++i) weights[i] = weight;
  }
  return 1;
}

static int compute_weights_avg_no(
  double tgt_coords[3], size_t num_src,
  yac_const_coordinate_pointer src_coords, int * src_mask, double * weights) {

  UNUSED(tgt_coords);
  UNUSED(src_coords);

  // if we have a source mask, check whether any points is masked
  if (src_mask != NULL)
    for (size_t i = 0; i < num_src; ++i)
      if (!src_mask[i]) return 0;

  double weight = 1.0 / (double)num_src;
  for (size_t i = 0; i < num_src; ++i) weights[i] = weight;

  return 1;
}

static int compute_weights_dist_yes(
  double tgt_coords[3], size_t num_src,
  yac_const_coordinate_pointer src_coords, int * src_mask, double * weights) {

  // if there is a source mask, check if there are unmasked vertices
  if (src_mask != NULL) {
    int has_unmasked = 0;
    for (size_t i = 0; i < num_src; ++i) has_unmasked |= src_mask[i];
    if (!has_unmasked) return 0;
  }

  if (src_mask != NULL) {
    for (size_t i = 0; i < num_src; ++i) {

      if (src_mask[i]) {

        double distance =
          get_vector_angle(tgt_coords, (double*)(src_coords[i]));

        // if the target and source point are nearly identical
        if (distance < yac_angle_tol) {
          for (size_t j = 0; j < num_src; ++j) weights[j] = 0.0;
          weights[i] = 1.0;
          return 1;
        }

        weights[i] = 1.0 / distance;
      } else {
        weights[i] = 0.0;
      }
    }
  } else {
    for (size_t i = 0; i < num_src; ++i) {

      double distance =
        get_vector_angle(tgt_coords, (double*)(src_coords[i]));

      // if the target and source point are nearly identical
      if (distance < yac_angle_tol) {
        for (size_t j = 0; j < num_src; ++j) weights[j] = 0.0;
        weights[i] = 1.0;
        return 1;
      }

      weights[i] = 1.0 / distance;
    }
  }

  // compute scaling factor for the weights
  double inv_distance_sum = 0.0;
  for (size_t i = 0; i < num_src; ++i)
    inv_distance_sum += weights[i];
  double scale = 1.0 / inv_distance_sum;

  for (size_t i = 0; i < num_src; ++i) weights[i] *= scale;

  return 1;
}

static int compute_weights_dist_no(
  double tgt_coords[3], size_t num_src,
  yac_const_coordinate_pointer src_coords, int * src_mask, double * weights) {

  for (size_t i = 0; i < num_src; ++i) {

    double distance =
      get_vector_angle(tgt_coords, (double*)(src_coords[i]));

    // if the target and source point are nearly identical
    if (distance < yac_angle_tol) {
      if ((src_mask != NULL) && !src_mask[i]) return 0;
      for (size_t j = 0; j < num_src; ++j) weights[j] = 0.0;
      weights[i] = 1.0;
      return 1;
    }

    weights[i] = 1.0 / distance;
  }

  // if there is a source mask, check if there are masked vertices
  // (we do this here because there may be a matching source point
  //  that is not masked)
  if (src_mask != NULL)
    for (size_t i = 0; i < num_src; ++i) if(!src_mask[i]) return 0;

  // compute scaling factor for the weights
  double inv_distance_sum = 0.0;
  for (size_t i = 0; i < num_src; ++i)
    inv_distance_sum += weights[i];
  double scale = 1.0 / inv_distance_sum;

  for (size_t i = 0; i < num_src; ++i) weights[i] *= scale;

  return 1;
}

static func_compute_weights select_compute_weight_routine(
  enum yac_interp_ncc_weight_type weight_type, int partial_coverage) {

  YAC_ASSERT(
    (weight_type == YAC_INTERP_NCC_AVG) ||
    (weight_type == YAC_INTERP_NCC_DIST),
    "ERROR(select_compute_weight_routine): invalid weight type")

  switch(weight_type) {
    case(YAC_INTERP_NCC_AVG):
      return
        (partial_coverage)?compute_weights_avg_yes:compute_weights_avg_no;
    case(YAC_INTERP_NCC_DIST):
    default:
      return
        (partial_coverage)?compute_weights_dist_yes:compute_weights_dist_no;
  };
}

struct interp_method * yac_interp_method_ncc_new(
  enum yac_interp_ncc_weight_type weight_type,
  int partial_coverage) {

  struct interp_method_ncc * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_ncc_vtable;
  method->compute_weights =
    select_compute_weight_routine(weight_type, partial_coverage);

  return (struct interp_method*)method;
}

static void delete_ncc(struct interp_method * method) {
  free(method);
}
