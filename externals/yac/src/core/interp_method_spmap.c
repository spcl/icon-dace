// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <float.h>
#include <string.h>

#include "interp_method_internal.h"
#include "interp_method_spmap.h"
#include "ensure_array_size.h"
#include "area.h"

static size_t do_search_spmap(struct interp_method * method,
                              struct yac_interp_grid * interp_grid,
                              size_t * tgt_points, size_t count,
                              struct yac_interp_weights * weights);
static void delete_spmap(struct interp_method * method);

static struct interp_method_vtable
  interp_method_spmap_vtable = {
    .do_search = do_search_spmap,
    .delete = delete_spmap
};

struct interp_method_spmap {

  struct interp_method_vtable * vtable;
  double spread_distance;
  double max_search_distance;
  enum yac_interp_spmap_weight_type weight_type;
  enum yac_interp_spmap_scale_type scale_type;
  double src_area_scale;
  double tgt_area_scale;
};

static inline int compare_size_t(const void * a, const void * b) {

  size_t const * a_ = a, * b_ = b;

  return (*a_ > *b_) - (*b_ > *a_);
}

static int lists_overlap(
  size_t * list_a, size_t count_a, size_t * list_b, size_t count_b) {

  if ((count_a == 0) || (count_b == 0)) return 0;

  size_t i = 0, j = 0;
  size_t curr_a = SIZE_MAX, curr_b = list_b[0];

  do {
    while ((i < count_a) && (((curr_a = list_a[i++])) < curr_b));
    if (curr_a == curr_b) return 1;
    while ((j < count_b) && (((curr_b = list_b[j++])) < curr_a));
    if (curr_a == curr_b) return 1;
  } while ((i < count_a) || (j < count_b));

  return 0;
}

static void merge_lists(
  size_t * list, size_t * list_size, size_t * insert, size_t insert_size) {

  if (insert_size == 0) return;

  size_t new_list_size = *list_size;
  size_t old_list_size = *list_size;

  for (size_t i = 0, j = 0; i < insert_size; ++i) {
    size_t curr_insert = insert[i];
    while ((j < old_list_size) && (list[j] < curr_insert)) ++j;
    if ((j >= old_list_size) || (list[j] != curr_insert))
      list[new_list_size++] = curr_insert;
  }

  if (new_list_size != old_list_size) {
    qsort(list, new_list_size, sizeof(*list), compare_size_t);
    *list_size = new_list_size;
  }
}

static void check_spread_distance(
  yac_const_coordinate_pointer tgt_field_coords, double max_distance,
  size_t tgt_start_point, size_t * tgt_points, size_t * count) {

  double const * start_coord = tgt_field_coords[tgt_start_point];
  size_t new_count = 0;
  for (size_t i = 0, old_count = *count; i < old_count; ++i) {
    if (get_vector_angle(
          start_coord, tgt_field_coords[tgt_points[i]]) <= max_distance) {
      if (new_count != i) tgt_points[new_count] = tgt_points[i];
      ++new_count;
    }
  }

  *count = new_count;
}

static void remove_disconnected_points(
  struct yac_interp_grid * interp_grid, size_t tgt_start_point,
  size_t * from_tgt_points, size_t * to_tgt_points,
  size_t * count, int * flag, size_t * temp_cell_edges,
  size_t ** edges_buffer, size_t * edges_buffer_array_size) {

  struct yac_const_basic_grid_data * tgt_grid_data =
    yac_interp_grid_get_basic_grid_data_tgt(interp_grid);
  const_size_t_pointer cell_to_edge = tgt_grid_data->cell_to_edge;
  const_size_t_pointer cell_to_edge_offsets =
    tgt_grid_data->cell_to_edge_offsets;
  const int * num_vertices_per_cell =
    tgt_grid_data->num_vertices_per_cell;

  size_t old_count = *count;
  memset(flag, 0, old_count * sizeof(*flag));

  for (size_t i = 0; i < old_count; ++i) {
    if (from_tgt_points[i] == tgt_start_point) {
      flag[i] = 1;
      break;
    }
  }

  size_t * edges = *edges_buffer;
  size_t num_edges = 0;
  size_t edges_array_size = *edges_buffer_array_size;

  const_size_t_pointer curr_edges =
    cell_to_edge + cell_to_edge_offsets[tgt_start_point];
  size_t curr_num_edges = num_vertices_per_cell[tgt_start_point];

  ENSURE_ARRAY_SIZE(edges, edges_array_size, curr_num_edges);
  num_edges = curr_num_edges;
  memcpy(edges, curr_edges, num_edges * sizeof(*edges));
  qsort(edges, num_edges, sizeof(*edges), compare_size_t);

  int change_flag = 0;
  do {

    change_flag = 0;

    for (size_t i = 0; i < old_count; ++i) {

      if (flag[i]) continue;

      const_size_t_pointer curr_edges =
        tgt_grid_data->cell_to_edge + cell_to_edge_offsets[from_tgt_points[i]];
      curr_num_edges = num_vertices_per_cell[from_tgt_points[i]];
      memcpy(
        temp_cell_edges, curr_edges, curr_num_edges * sizeof(*temp_cell_edges));
      qsort(temp_cell_edges, curr_num_edges, sizeof(*edges), compare_size_t);
      ENSURE_ARRAY_SIZE(edges, edges_array_size, num_edges + curr_num_edges);

      if (lists_overlap(edges, num_edges, temp_cell_edges, curr_num_edges)) {
        merge_lists(edges, &num_edges, temp_cell_edges, curr_num_edges);
        flag[i] = 1;
        change_flag = 1;
      }
    }
  } while (change_flag);

  *edges_buffer = edges;
  *edges_buffer_array_size = edges_array_size;

  size_t new_count = 0;
  for (size_t i = 0; i < old_count; ++i)
    if (flag[i]) to_tgt_points[new_count++] = from_tgt_points[i];

  *count = new_count;
}

static size_t check_tgt_result_points(
  struct yac_interp_grid * interp_grid, double spread_distance,
  size_t num_src_points, size_t const * const tgt_result_points,
  size_t * num_tgt_per_src, size_t * spread_tgt_result_points) {

  size_t max_num_tgt_per_src = 0;
  for (size_t i = 0; i < num_src_points; ++i)
    if (num_tgt_per_src[i] > max_num_tgt_per_src)
      max_num_tgt_per_src = num_tgt_per_src[i];

  int * flag = xmalloc(max_num_tgt_per_src * sizeof(*flag));

  size_t new_offset = 0;
  size_t * cell_edge_buffer = NULL;
  size_t cell_edge_buffer_array_size = 0;
  size_t * edge_buffer = NULL;
  size_t edge_buffer_array_size = 0;
  int max_num_vertice_per_tgt = 0;
  const int * num_vertices_per_tgt =
    yac_interp_grid_get_basic_grid_data_tgt(interp_grid)->
      num_vertices_per_cell;

  yac_const_coordinate_pointer tgt_field_coords =
    yac_interp_grid_get_tgt_field_coords(interp_grid);

  // for all source points
  for (size_t i = 0, old_offset = 0; i < num_src_points; ++i) {

    size_t * old_results = spread_tgt_result_points + old_offset;
    size_t * new_results = spread_tgt_result_points + new_offset;
    old_offset += num_tgt_per_src[i];

    // remove all target points, which exceed the spread distance from
    // the original tgt
    check_spread_distance(
      tgt_field_coords, spread_distance, tgt_result_points[i],
      old_results, num_tgt_per_src + i);

    // check buffer sizes required by routine remove_disconnected_points
    for (size_t j = 0, curr_num_tgt_per_src = num_tgt_per_src[i];
        j < curr_num_tgt_per_src; ++j)
      if (num_vertices_per_tgt[old_results[j]] > max_num_vertice_per_tgt)
        max_num_vertice_per_tgt = num_vertices_per_tgt[old_results[j]];

    ENSURE_ARRAY_SIZE(
      cell_edge_buffer, cell_edge_buffer_array_size,
      max_num_vertice_per_tgt);

    // remove all tgts that are not directly connected to the original tgt
    remove_disconnected_points(
      interp_grid, tgt_result_points[i],
      old_results, new_results, num_tgt_per_src + i, flag,
      cell_edge_buffer, &edge_buffer, &edge_buffer_array_size);

    new_offset += num_tgt_per_src[i];
  }

  free(edge_buffer);
  free(cell_edge_buffer);
  free(flag);

  return new_offset;
}

static double * compute_weights(
  struct yac_interp_grid * interp_grid,
  enum yac_interp_spmap_weight_type weight_type, size_t num_src_points,
  size_t const * const src_points, size_t const * const num_tgt_per_src,
  size_t total_num_tgt, size_t const * const tgt_result_points) {

  double * weights = xmalloc(total_num_tgt * sizeof(*weights));
  YAC_ASSERT(
    (weight_type == YAC_INTERP_SPMAP_AVG) ||
    (weight_type == YAC_INTERP_SPMAP_DIST),
    "ERROR(do_search_spmap): invalid weight_type")
  switch (weight_type) {
    case (YAC_INTERP_SPMAP_AVG): {
      for (size_t i = 0, offset = 0; i < num_src_points; ++i) {
        size_t curr_num_tgt = num_tgt_per_src[i];
        if (curr_num_tgt == 0) continue;
        if (curr_num_tgt > 1) {
          double curr_weight_data = 1.0 / (double)(curr_num_tgt);
          for (size_t j = 0; j < curr_num_tgt; ++j, ++offset) {
            weights[offset] = curr_weight_data;
          }
        } else {
          weights[offset] = 1.0;
          ++offset;
        }
      }
      break;
    }
    default:
    case (YAC_INTERP_SPMAP_DIST): {

      yac_const_coordinate_pointer src_field_coords =
        yac_interp_grid_get_src_field_coords(interp_grid, 0);
      yac_const_coordinate_pointer tgt_field_coords =
        yac_interp_grid_get_tgt_field_coords(interp_grid);

      for (size_t i = 0, offset = 0; i < num_src_points; ++i) {

        size_t curr_num_tgt = num_tgt_per_src[i];

        if (curr_num_tgt == 0) continue;

        double * curr_weights = weights + offset;

        if (curr_num_tgt > 1) {

          size_t const * const  curr_result_points =
            tgt_result_points + offset;
          double const * curr_src_coord = src_field_coords[src_points[offset]];
          offset += curr_num_tgt;

          int match_flag = 0;

          for (size_t j = 0; j < curr_num_tgt; ++j) {

            double distance =
              get_vector_angle(
                (double*)curr_src_coord,
                (double*)tgt_field_coords[curr_result_points[j]]);

            if (distance < yac_angle_tol) {
              for (size_t k = 0; k < curr_num_tgt; ++k) curr_weights[k] = 0.0;
              curr_weights[j] = 1.0;
              match_flag = 1;
              break;
            }
            curr_weights[j] = 1.0 / distance;
          }

          if (!match_flag) {

            // compute scaling factor for the weights
            double inv_distance_sum = 0.0;
            for (size_t j = 0; j < curr_num_tgt; ++j)
              inv_distance_sum += curr_weights[j];
            double scale = 1.0 / inv_distance_sum;

            for (size_t j = 0; j < curr_num_tgt; ++j) curr_weights[j] *= scale;
          }
        } else {
          *curr_weights = 1.0;
          ++offset;
        }
      }
      break;
    }
  };

  return weights;
}

static void scale_weights(
  struct yac_interp_grid * interp_grid,
  enum yac_interp_spmap_scale_type scale_type,
  double src_area_scale, double tgt_area_scale,
  size_t num_src_points, size_t const * src_points,
  size_t const * num_tgt_per_src, size_t const * tgt_points,
  double * weights) {

  YAC_ASSERT(
    (scale_type == YAC_INTERP_SPMAP_NONE) ||
    (scale_type == YAC_INTERP_SPMAP_SRCAREA) ||
    (scale_type == YAC_INTERP_SPMAP_INVTGTAREA) ||
    (scale_type == YAC_INTERP_SPMAP_FRACAREA),
    "ERROR(scale_weights): invalid scale_type")

  // if there is no scaling
  if (scale_type == YAC_INTERP_SPMAP_NONE) return;

  struct yac_grid_cell grid_cell;
  yac_init_grid_cell(&grid_cell);

  struct yac_const_basic_grid_data * src_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_src(interp_grid);
  struct yac_const_basic_grid_data * tgt_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_tgt(interp_grid);

#define COMPUTE_CELL_AREA(PREFIX, IDX) \
  yac_const_basic_grid_data_get_grid_cell( \
    PREFIX ## _basic_grid_data, PREFIX ## _points[IDX], &grid_cell); \
  double PREFIX ## _cell_area = yac_huiliers_area(grid_cell) * PREFIX ## _area_scale; \
  YAC_ASSERT_F( \
    PREFIX ## _cell_area > YAC_AREA_TOL, \
    "ERROR(scale_weights): " \
    "area of %s cell (global id %"XT_INT_FMT") is close to zero (%e)", \
    # PREFIX, \
    PREFIX ## _basic_grid_data->ids[YAC_LOC_CELL][PREFIX ## _points[IDX]], \
    PREFIX ## _cell_area)
#define NO_COMPUTE_CELL_AREA(PREFIX, IDX)

#define SCALE_WEIGHTS( \
  COMPUTE_SRC_CELL_AREA, \
  COMPUTE_TGT_CELL_AREA, \
  SRC_CELL_AREA, TGT_CELL_AREA) \
{ \
  for (size_t i = 0, offset = 0; i < num_src_points; ++i) { \
    COMPUTE_SRC_CELL_AREA(src, offset) \
    size_t curr_num_tgt = num_tgt_per_src[i]; \
    for (size_t j = 0; j < curr_num_tgt; ++j, ++offset) { \
      COMPUTE_TGT_CELL_AREA(tgt, offset) \
      weights[offset] *= SRC_CELL_AREA / TGT_CELL_AREA; \
    } \
  } \
}

  switch (scale_type) {
    case(YAC_INTERP_SPMAP_SRCAREA):
      SCALE_WEIGHTS(
        COMPUTE_CELL_AREA, NO_COMPUTE_CELL_AREA, src_cell_area, 1.0)
      break;
    case(YAC_INTERP_SPMAP_INVTGTAREA):
      SCALE_WEIGHTS(
        NO_COMPUTE_CELL_AREA, COMPUTE_CELL_AREA, 1.0, tgt_cell_area)
      break;
    default:
    case(YAC_INTERP_SPMAP_FRACAREA):
      SCALE_WEIGHTS(
        COMPUTE_CELL_AREA, COMPUTE_CELL_AREA, src_cell_area, tgt_cell_area)
      break;
  }

  yac_free_grid_cell(&grid_cell);
}

static void spread_src_data(
  struct yac_interp_grid * interp_grid, double spread_distance,
  enum yac_interp_spmap_weight_type weight_type,
  enum yac_interp_spmap_scale_type scale_type,
  double src_area_scale, double tgt_area_scale,
  size_t num_src_points, size_t ** src_points_,
  size_t ** tgt_result_points_, double ** weights_,
  size_t * total_num_weights_) {

  // shortcut in case there is only a single "1.0" weight for all source points
  if ((spread_distance <= 0) && (scale_type == YAC_INTERP_SPMAP_NONE)) {

    *weights_ = NULL;
    *total_num_weights_ = num_src_points;
    return;
  }

  size_t * src_points = *src_points_;
  size_t * tgt_result_points = *tgt_result_points_;

  size_t * num_tgt_per_src =
    xmalloc(num_src_points * sizeof(*num_tgt_per_src));
  size_t total_num_weights;

  // search for additional target points if spread distance is bigger than 0.0
  if (spread_distance > 0.0) {

    struct sin_cos_angle inc_angle =
      sin_cos_angle_new(sin(spread_distance), cos(spread_distance));
    yac_const_coordinate_pointer tgt_field_coords =
      yac_interp_grid_get_tgt_field_coords(interp_grid);

    struct bounding_circle * search_bnd_circles =
      xmalloc(num_src_points * sizeof(*search_bnd_circles));
    for (size_t i = 0; i < num_src_points; ++i) {
      memcpy(
        search_bnd_circles[i].base_vector,
        tgt_field_coords[tgt_result_points[i]], sizeof(*tgt_field_coords));
      search_bnd_circles[i].inc_angle = inc_angle;
      search_bnd_circles[i].sq_crd = DBL_MAX;
    }
    size_t * spread_tgt_result_points = NULL;

    // do bounding circle search around found tgt points
    yac_interp_grid_do_bnd_circle_search_tgt(
      interp_grid, search_bnd_circles, num_src_points,
      &spread_tgt_result_points, num_tgt_per_src);
    free(search_bnd_circles);

    // remove target points which exceed the spread distance and only keep
    // targets that are connected to the original target point or other
    // target that have already been selected
    total_num_weights =
      check_tgt_result_points(
        interp_grid, spread_distance, num_src_points, tgt_result_points,
        num_tgt_per_src, spread_tgt_result_points);
    free((void*)tgt_result_points);
    tgt_result_points =
      xrealloc(
        spread_tgt_result_points,
        total_num_weights * sizeof(*spread_tgt_result_points));

    // adjust src_points (one source per target)
    size_t * new_src_points =
      xmalloc(total_num_weights * sizeof(*new_src_points));
    for (size_t i = 0, offset = 0; i < num_src_points; ++i)
      for (size_t j = 0, curr_src_point = src_points[i];
            j < num_tgt_per_src[i]; ++j, ++offset)
        new_src_points[offset] = curr_src_point;
    free((void*)src_points);
    src_points = new_src_points;

  } else {

    for (size_t i = 0; i < num_src_points; ++i) num_tgt_per_src[i] = 1;
    total_num_weights = num_src_points;
  }

  // compute weights
  double * weights =
    compute_weights(
      interp_grid, weight_type, num_src_points, src_points,
      num_tgt_per_src, total_num_weights, tgt_result_points);

  // scale weights
  scale_weights(
    interp_grid, scale_type, src_area_scale, tgt_area_scale, num_src_points,
    src_points, num_tgt_per_src, tgt_result_points, weights);

  free(num_tgt_per_src);

  // set return values
  *tgt_result_points_ = tgt_result_points;
  *src_points_ = src_points;
  *weights_ = weights;
  *total_num_weights_ = total_num_weights;
}

static size_t do_search_spmap (struct interp_method * method,
                               struct yac_interp_grid * interp_grid,
                               size_t * tgt_points, size_t count,
                               struct yac_interp_weights * weights) {

  YAC_ASSERT(
    yac_interp_grid_get_num_src_fields(interp_grid) == 1,
    "ERROR(do_search_spmap): invalid number of source fields")
  YAC_ASSERT(
    yac_interp_grid_get_src_field_location(interp_grid, 0) == YAC_LOC_CELL,
    "ERROR(do_search_spmap): "
    "invalid source field location (has to be YAC_LOC_CELL)")
  YAC_ASSERT(
    yac_interp_grid_get_tgt_field_location(interp_grid) == YAC_LOC_CELL,
    "ERROR(do_search_spmap): "
    "invalid target field location (has to be YAC_LOC_CELL)")

  // get coordinates of all source points
  size_t * src_points;
  size_t num_src_points;
  yac_interp_grid_get_src_points(
    interp_grid, 0, &src_points, &num_src_points);
  yac_coordinate_pointer src_coords = xmalloc(num_src_points * sizeof(*src_coords));
  yac_interp_grid_get_src_coordinates(
    interp_grid, src_points, num_src_points, 0, src_coords);

  // search for matching tgt points
  size_t * tgt_result_points =
    xmalloc(num_src_points * sizeof(*tgt_result_points));
  yac_interp_grid_do_nnn_search_tgt(
    interp_grid, src_coords, num_src_points, 1, tgt_result_points,
    ((struct interp_method_spmap*)method)->max_search_distance);

  free(src_coords);

  // remove source points for which matching target point was found
  {
    size_t new_num_src_points = 0;
    for (size_t i = 0; i < num_src_points; ++i) {
      if (tgt_result_points[i] != SIZE_MAX) {
        if (i != new_num_src_points) {
          src_points[new_num_src_points] = src_points[i];
          tgt_result_points[new_num_src_points] =
            tgt_result_points[i];
        }
        ++new_num_src_points;
      }
    }
    num_src_points = new_num_src_points;
  }

  // spread the data from each source point to multiple target points
  // (weight_data is set to NULL if no spreading was applied)
  double * weight_data;
  size_t total_num_weights;
  spread_src_data(
    interp_grid, ((struct interp_method_spmap*)method)->spread_distance,
    ((struct interp_method_spmap*)method)->weight_type,
    ((struct interp_method_spmap*)method)->scale_type,
    ((struct interp_method_spmap*)method)->src_area_scale,
    ((struct interp_method_spmap*)method)->tgt_area_scale, num_src_points,
    &src_points, &tgt_result_points, &weight_data, &total_num_weights);

  // relocate source-target-point-pairs to dist owners of the respective
  // target points
  size_t result_count = total_num_weights;
  int to_tgt_owner = 1;
  yac_interp_grid_relocate_src_tgt_pairs(
    interp_grid, to_tgt_owner,
    0, &src_points, &tgt_result_points, &weight_data, &result_count);
  total_num_weights = result_count;

  // sort source-target-point-pairs by target points
  yac_quicksort_index_size_t_size_t_double(
    tgt_result_points, result_count, src_points, weight_data);

  // generate num_src_per_tgt and compact tgt_result_points
  size_t * num_src_per_tgt = xmalloc(result_count * sizeof(*num_src_per_tgt));
  size_t num_unique_tgt_result_points = 0;
  for (size_t i = 0; i < result_count;) {
    size_t prev_i = i;
    size_t curr_tgt = tgt_result_points[i];
    while ((i < result_count) && (curr_tgt == tgt_result_points[i])) ++i;
    num_src_per_tgt[num_unique_tgt_result_points] = i - prev_i;
    tgt_result_points[num_unique_tgt_result_points] = curr_tgt;
    ++num_unique_tgt_result_points;
  }
  result_count = num_unique_tgt_result_points;
  num_src_per_tgt =
    xrealloc(num_src_per_tgt, result_count * sizeof(*num_src_per_tgt));
  tgt_result_points =
    xrealloc(tgt_result_points, result_count * sizeof(*tgt_result_points));

  // match tgt_result_points with available target points
  qsort(tgt_points, count, sizeof(*tgt_points), compare_size_t);
  int * reorder_flag = xmalloc(count * sizeof(*tgt_points));
  {
    size_t j = 0;
    for (size_t i = 0; i < result_count; ++i) {
      size_t curr_result_tgt = tgt_result_points[i];
      while ((j < count) && (tgt_points[j] < curr_result_tgt))
        reorder_flag[j++] = 1;
      YAC_ASSERT(
        (j < count) && (curr_result_tgt == tgt_points[j]),
        "ERROR(do_search_spmap): "
        "required target points already in use or not available")
      reorder_flag[j++] = 0;
    }
    for (; j < count; ++j) reorder_flag[j] = 1;
  }

  // sort used target points to the beginning of the array
  yac_quicksort_index_int_size_t(reorder_flag, count, tgt_points);
  free(reorder_flag);

  struct remote_point * srcs =
    yac_interp_grid_get_src_remote_points(
      interp_grid, 0, src_points, total_num_weights);
  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(
        interp_grid, tgt_result_points, result_count),
    .count = result_count};
  free(tgt_result_points);

  // store results
  if (weight_data == NULL)
    yac_interp_weights_add_sum(
      weights, &tgts, num_src_per_tgt, srcs);
  else
    yac_interp_weights_add_wsum(
      weights, &tgts, num_src_per_tgt, srcs, weight_data);

  free(weight_data);
  free(src_points);
  free(tgts.data);
  free(srcs);
  free(num_src_per_tgt);

  return result_count;
}

struct interp_method * yac_interp_method_spmap_new(
  double spread_distance, double max_search_distance,
  enum yac_interp_spmap_weight_type weight_type,
  enum yac_interp_spmap_scale_type scale_type,
  double src_sphere_radius, double tgt_sphere_radius) {

  struct interp_method_spmap * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_spmap_vtable;
  method->spread_distance = spread_distance;
  method->max_search_distance =
    (max_search_distance == 0.0)?M_PI:max_search_distance;
  method->weight_type = weight_type;
  method->scale_type = scale_type;
  method->src_area_scale = src_sphere_radius * src_sphere_radius;
  method->tgt_area_scale = tgt_sphere_radius * tgt_sphere_radius;

  YAC_ASSERT(
    (spread_distance >= 0.0) && (spread_distance <= M_PI_2),
    "ERROR(yac_interp_method_spmap_new): invalid spread_distance "
    "(has to be >= 0 and <= PI/2")

  YAC_ASSERT(
    (max_search_distance >= 0.0) && (max_search_distance <= M_PI),
    "ERROR(yac_interp_method_spmap_new): invalid max_search_distance "
    "(has to be >= 0 and <= PI")

  YAC_ASSERT(
    src_sphere_radius > 0.0,
    "ERROR(yac_interp_method_spmap_new): invalid src_sphere_radius "
    "(has to be >= 0.0")

  YAC_ASSERT(
    tgt_sphere_radius > 0.0,
    "ERROR(yac_interp_method_spmap_new): invalid tgt_sphere_radius "
    "(has to be >= 0.0")

  return (struct interp_method*)method;
}

static void delete_spmap(struct interp_method * method) {
  free(method);
}
