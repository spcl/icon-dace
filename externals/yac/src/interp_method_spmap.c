/**
 * @file interp_method_spmap.c
 *
 * @copyright Copyright  (C)  2020 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <float.h>

#include "string.h"
#include "interp_method_spmap.h"
#include "ensure_array_size.h"

static size_t do_search_spmap(struct interp_method * method,
                              struct interp_grid * interp_grid,
                              size_t * tgt_points, size_t count,
                              struct interp_weights * weights);
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

static void check_distances(
  const_coordinate_pointer tgt_field_coords, double max_distance,
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
  struct interp_grid * interp_grid, size_t tgt_start_point,
  size_t * from_tgt_points, size_t * to_tgt_points,
  size_t * count, int * flag, size_t * temp_cell_edges,
  size_t ** edges_buffer, size_t * edges_buffer_array_size) {

  struct const_basic_grid_data * tgt_grid_data =
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

static size_t do_search_spmap (struct interp_method * method,
                               struct interp_grid * interp_grid,
                               size_t * tgt_points, size_t count,
                               struct interp_weights * weights) {

  YAC_ASSERT(
    yac_interp_grid_get_num_src_fields(interp_grid) == 1,
    "ERROR(do_search_spmap): invalid number of source fields")
  YAC_ASSERT(
    yac_interp_grid_get_src_field_location(interp_grid, 0) == CELL,
    "ERROR(do_search_spmap): invalid source field location (has to be CELL)")
  YAC_ASSERT(
    yac_interp_grid_get_tgt_field_location(interp_grid) == CELL,
    "ERROR(do_search_spmap): invalid target field location (has to be CELL)")

  // get coordinates of all source points
  size_t * src_points;
  size_t num_src_points;
  yac_interp_grid_get_src_points(
    interp_grid, 0, &src_points, &num_src_points);
  coordinate_pointer src_coords = xmalloc(num_src_points * sizeof(*src_coords));
  yac_interp_grid_get_src_coordinates(
    interp_grid, src_points, num_src_points, 0, src_coords);

  // search for matching tgt points
  size_t * tgt_result_points =
    xmalloc(num_src_points * sizeof(*tgt_result_points));
  yac_interp_grid_do_nnn_search_tgt(
    interp_grid, src_coords, num_src_points, 1, tgt_result_points);

  // check that we found a target for each source point
  for (size_t i = 0; i < num_src_points; ++i)
    YAC_ASSERT(
      tgt_result_points[i] != SIZE_MAX,
      "ERROR(do_search_spmap): could not find matching target point")

  // ensure that all results are within the maximum search distance
  double max_search_distance =
    ((struct interp_method_spmap*)method)->max_search_distance;
  if (max_search_distance > 0.0) {
    const_coordinate_pointer tgt_field_coords =
      yac_interp_grid_get_tgt_field_coords(interp_grid);
    size_t new_num_src_points = 0;
    for (size_t i = 0; i < num_src_points; ++i) {
      if (get_vector_angle(
            src_coords[i], tgt_field_coords[tgt_result_points[i]]) <=
          max_search_distance) {
        if (i != new_num_src_points) {
          src_points[new_num_src_points] = src_points[i];
          tgt_result_points[new_num_src_points] = tgt_result_points[i];
        }
        new_num_src_points++;
      }
    }
    num_src_points = new_num_src_points;
  }

  free(src_coords);

  double * weight_data = NULL;

  double spread_distance =
    ((struct interp_method_spmap*)method)->spread_distance;
  if (spread_distance > 0.0) {

    struct sin_cos_angle inc_angle =
      sin_cos_angle_new(sin(spread_distance), cos(spread_distance));
    const_coordinate_pointer tgt_field_coords =
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
    size_t * temp_tgt_result_points = NULL;

    // do bounding circle search around found tgt points
    size_t * num_tgt_per_src =
      xmalloc(num_src_points * sizeof(*num_tgt_per_src));
    yac_interp_grid_do_bnd_circle_search_tgt(
      interp_grid, search_bnd_circles, num_src_points,
      &temp_tgt_result_points, num_tgt_per_src);
    free(search_bnd_circles);

    size_t max_num_tgt_per_src = 0;
    for (size_t i = 0; i < num_src_points; ++i)
      if (num_tgt_per_src[i] > max_num_tgt_per_src)
        max_num_tgt_per_src = num_tgt_per_src[i];

    int * flag = xmalloc(max_num_tgt_per_src * sizeof(*flag));

    // remove tgt points not directly connected to original tgt and whose
    // distance exceed the spread distance (search results at this point
    // are based on bounding circles of the target cells)
    size_t new_offset = 0;
    size_t * cell_edge_buffer = NULL;
    size_t cell_edge_buffer_array_size = 0;
    size_t * edge_buffer = NULL;
    size_t edge_buffer_array_size = 0;
    size_t max_num_vertice_per_tgt = 0;
    const int * num_vertices_per_tgt =
      yac_interp_grid_get_basic_grid_data_tgt(interp_grid)->
        num_vertices_per_cell;

    {
      // the call to yac_interp_grid_do_bnd_circle_search_tgt might have
      // changed the target field coordinate pointer array
      const_coordinate_pointer tgt_field_coords =
        yac_interp_grid_get_tgt_field_coords(interp_grid);

      // for all source points
      for (size_t i = 0, old_offset = 0; i < num_src_points; ++i) {
        size_t * old_results = temp_tgt_result_points + old_offset;
        size_t * new_results = temp_tgt_result_points + new_offset;
        old_offset += num_tgt_per_src[i];

        // remove all tgts, which exceed the spread distance from
        // the original tgt
        check_distances(
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
    }
    free(edge_buffer);
    free(cell_edge_buffer);
    free(flag);


    // adjust src_points
    size_t * new_src_points = xmalloc(new_offset * sizeof(*new_src_points));
    for (size_t i = 0, offset = 0; i < num_src_points; ++i)
      for (size_t j = 0, curr_src_point = src_points[i];
           j < num_tgt_per_src[i]; ++j, ++offset)
        new_src_points[offset] = curr_src_point;

    // compute weight data
    weight_data = xmalloc(new_offset * sizeof(*weight_data));
    YAC_ASSERT(
      (((struct interp_method_spmap*)method)->weight_type == SPMAP_AVG) ||
      (((struct interp_method_spmap*)method)->weight_type == SPMAP_DIST),
      "ERROR(do_search_spmap): invalid weight_type")
    switch (((struct interp_method_spmap*)method)->weight_type) {
      case (SPMAP_AVG): {
        for (size_t i = 0, offset = 0; i < num_src_points; ++i) {
          double curr_weight_data = 1.0 / (double)(num_tgt_per_src[i]);
          for (size_t j = 0; j < num_tgt_per_src[i]; ++j, ++offset) {
            weight_data[offset] = curr_weight_data;
          }
        }
        break;
      }
      default:
      case (SPMAP_DIST): {

        const_coordinate_pointer src_field_coords =
          yac_interp_grid_get_src_field_coords(interp_grid, 0);
        const_coordinate_pointer tgt_field_coords =
          yac_interp_grid_get_tgt_field_coords(interp_grid);

        for (size_t i = 0, offset = 0; i < num_src_points; ++i) {

          size_t curr_num_tgt = num_tgt_per_src[i];
          size_t * curr_result_points = temp_tgt_result_points + offset;
          double const * curr_src_coord = src_field_coords[src_points[i]];
          double * curr_weights = weight_data + offset;
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
        }
        break;
      }
    };

    free(num_tgt_per_src);
    free(tgt_result_points);
    tgt_result_points =
      xrealloc(
        temp_tgt_result_points, new_offset * sizeof(*temp_tgt_result_points));
    num_src_points = new_offset;
    free(src_points);
    src_points = new_src_points;
  }

  // relocate source-target-point-pairs to dist owners of the respective
  // target points
  size_t result_count = num_src_points;
  int to_tgt_owner = 1;
  yac_interp_grid_relocate_src_tgt_pairs(
    interp_grid, to_tgt_owner,
    0, &src_points, &tgt_result_points, &weight_data, &result_count);
  num_src_points = result_count;

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
      interp_grid, 0, src_points, num_src_points);
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
  enum yac_interp_spmap_weight_type weight_type) {

  struct interp_method_spmap * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_spmap_vtable;
  method->spread_distance = spread_distance;
  method->max_search_distance = max_search_distance;
  method->weight_type = weight_type;

  YAC_ASSERT(
    (spread_distance >= 0.0) && (spread_distance <= M_PI_2),
    "ERROR(yac_interp_method_spmap_new): invalid spread_distance "
    "(has to be >= 0 and <= PI/2")

  YAC_ASSERT(
    (max_search_distance >= 0.0) && (max_search_distance <= M_PI),
    "ERROR(yac_interp_method_spmap_new): invalid max_search_distance "
    "(has to be >= 0 and <= PI")

  return (struct interp_method*)method;
}

static void delete_spmap(struct interp_method * method) {
  free(method);
}
