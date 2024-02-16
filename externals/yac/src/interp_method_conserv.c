/**
 * @file interp_method_conserv.c
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

#include "string.h"

#include "interp_method_conserv.h"
#include "clipping.h"
#include "area.h"
#include "float.h"

#define AREA_TOL_FACTOR (1e-6)

static size_t do_search_conserv_1st_order(struct interp_method * method,
                                          struct interp_grid * interp_grid,
                                          size_t * tgt_points, size_t count,
                                          struct interp_weights * weights);
static size_t do_search_conserv_2nd_order(struct interp_method * method,
                                          struct interp_grid * interp_grid,
                                          size_t * tgt_points, size_t count,
                                          struct interp_weights * weights);
static void delete_conserv(struct interp_method * method);

static struct interp_method_vtable
  interp_method_conserv_1st_order_vtable = {
    .do_search = do_search_conserv_1st_order,
    .delete = delete_conserv};

static struct interp_method_vtable
  interp_method_conserv_2nd_order_vtable = {
    .do_search = do_search_conserv_2nd_order,
    .delete = delete_conserv};

struct interp_method_conserv {

  struct interp_method_vtable * vtable;
  int partial_coverage;
  enum yac_interp_method_conserv_normalisation normalisation;
  int enforced_conserv;
};

struct weight_vector_data {
  double weight;
  size_t local_id;
  yac_int global_id;
};

struct weight_vector_data_3d {
  double weight[3];
  size_t local_id;
  yac_int global_id;
};

struct weight_vector_3d {
  struct weight_vector_data_3d * data;
  size_t n;
};

struct supermesh_cell {
  struct {
    size_t local_id;
    yac_int global_id;
  } src, tgt;
  double norm_area;
  double area;
  double barycenter[3];
  struct weight_vector_3d * src_cell_gradient;
};

static int get_max_num_vertices_per_cell(
  struct const_basic_grid_data * basic_grid_data) {

  int max_num_vertices_per_cell = 0;
  for (size_t i = 0; i < basic_grid_data->count[CELL]; ++i)
    if (basic_grid_data->num_vertices_per_cell[i] >
        max_num_vertices_per_cell)
      max_num_vertices_per_cell =
        basic_grid_data->num_vertices_per_cell[i];
  return max_num_vertices_per_cell;
}

static void get_cell_buffers(
  struct interp_grid * interp_grid, size_t max_num_src_per_tgt,
  struct grid_cell * tgt_grid_cell, struct grid_cell ** src_grid_cells) {

  struct const_basic_grid_data * src_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_src(interp_grid);
  struct const_basic_grid_data * tgt_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_tgt(interp_grid);

  *src_grid_cells = xmalloc(max_num_src_per_tgt * sizeof(**src_grid_cells));
  enum yac_edge_type * edge_type_buffer;
  double (*coordinates_xyz_buffer)[3];

  // prepare source grid cell buffer
  {
    int max_num_vertices_per_cell =
      MAX(get_max_num_vertices_per_cell(src_basic_grid_data),
          get_max_num_vertices_per_cell(tgt_basic_grid_data));

    edge_type_buffer =
      xmalloc((max_num_src_per_tgt + 1) * (size_t)max_num_vertices_per_cell *
              sizeof(*edge_type_buffer));
    coordinates_xyz_buffer =
      xmalloc((max_num_src_per_tgt + 1) * (size_t)max_num_vertices_per_cell *
              sizeof(*coordinates_xyz_buffer));

    tgt_grid_cell->coordinates_xyz = coordinates_xyz_buffer;
    tgt_grid_cell->edge_type = edge_type_buffer;
    tgt_grid_cell->array_size = max_num_vertices_per_cell;
    for (size_t i = 0; i < max_num_src_per_tgt; ++i) {
      (*src_grid_cells)[i].coordinates_xyz =
        coordinates_xyz_buffer + (i + 1) * max_num_vertices_per_cell;
      (*src_grid_cells)[i].edge_type =
        edge_type_buffer + (i + 1) * max_num_vertices_per_cell;
      (*src_grid_cells)[i].array_size = max_num_vertices_per_cell;
    }
  }
}

static void get_cell_buffers_(
  struct interp_grid * interp_grid,
  struct grid_cell * tgt_grid_cell, struct grid_cell * src_grid_cell) {

  struct const_basic_grid_data * src_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_src(interp_grid);
  struct const_basic_grid_data * tgt_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_tgt(interp_grid);

  int max_num_vertices_per_cell =
    MAX(get_max_num_vertices_per_cell(src_basic_grid_data),
        get_max_num_vertices_per_cell(tgt_basic_grid_data));

  enum yac_edge_type * edge_type_buffer =
    xmalloc(2 * (size_t)max_num_vertices_per_cell *
            sizeof(*edge_type_buffer));
  coordinate_pointer coordinates_xyz_buffer =
    xmalloc(2 * (size_t)max_num_vertices_per_cell *
            sizeof(*coordinates_xyz_buffer));

  tgt_grid_cell->coordinates_xyz = coordinates_xyz_buffer;
  tgt_grid_cell->edge_type = edge_type_buffer;
  tgt_grid_cell->array_size = max_num_vertices_per_cell;

  src_grid_cell->coordinates_xyz =
    coordinates_xyz_buffer + max_num_vertices_per_cell;
  src_grid_cell->edge_type =
    edge_type_buffer + max_num_vertices_per_cell;
  src_grid_cell->array_size = max_num_vertices_per_cell;
}

static int compute_1st_order_weights(
  struct const_basic_grid_data * tgt_basic_grid_data, size_t tgt_cell,
  struct const_basic_grid_data * src_basic_grid_data, size_t src_count,
  size_t * src_cells, struct grid_cell tgt_grid_cell_buffer,
  struct grid_cell * src_grid_cell_buffer,
  double * weights, size_t * num_weights, int partial_coverage,
  enum yac_interp_method_conserv_normalisation normalisation,
  int enforced_conserv) {

  yac_const_basic_grid_data_get_grid_cell(
    tgt_basic_grid_data, tgt_cell, &tgt_grid_cell_buffer);
  for (size_t i = 0; i < src_count; ++i)
    yac_const_basic_grid_data_get_grid_cell(
      src_basic_grid_data, src_cells[i], src_grid_cell_buffer + i);

  double * area = weights;
  yac_compute_overlap_areas(
    src_count, src_grid_cell_buffer, tgt_grid_cell_buffer, area);

  size_t num_valid_weights = 0;
  for (size_t i = 0; i < src_count; ++i) {

    if (area[i] > 0.0) {
      if (i != num_valid_weights) {
        area[num_valid_weights] = area[i];
        src_cells[num_valid_weights] = src_cells[i];
      }
      ++num_valid_weights;
    }
  }
  *num_weights = num_valid_weights;
  if (num_valid_weights == 0) return 0;

  double tgt_cell_area = yac_huiliers_area(tgt_grid_cell_buffer);
  double norm_factor;

  YAC_ASSERT(
    (normalisation == CONSERV_DESTAREA) || (normalisation == CONSERV_FRACAREA),
    "ERROR(compute_weights_order_first_conserv_no_partial): "
    "invalid normalisation option in conservative remapping")
  switch(normalisation) {
    case(CONSERV_DESTAREA):
      norm_factor = 1.0 / tgt_cell_area;
      break;
    default:
    case(CONSERV_FRACAREA): {
      double fracarea = 0.0;
      for (size_t i = 0; i < num_valid_weights; ++i) fracarea += area[i];
      norm_factor = 1.0 / fracarea;
      break;
    }
  };

  if (partial_coverage) {
    for (size_t i = 0; i < num_valid_weights; ++i) weights[i] *= norm_factor;
    return 1;
  } else {
    double tgt_cell_area_diff = tgt_cell_area;
    double area_tol = tgt_cell_area * AREA_TOL_FACTOR;
    for (size_t i = 0; i < num_valid_weights; ++i) {
      double curr_area = area[i];
      tgt_cell_area_diff -= curr_area;
      weights[i] = curr_area * norm_factor;
    }
    int successful = fabs(tgt_cell_area_diff) <= area_tol;
    if (successful && enforced_conserv)
      yac_correct_weights(num_valid_weights, weights);
    return successful;
  }
}

static size_t do_search_conserv_1st_order (struct interp_method * method,
                                           struct interp_grid * interp_grid,
                                           size_t * tgt_points, size_t count,
                                           struct interp_weights * weights) {

  struct interp_method_conserv * method_conserv =
    (struct interp_method_conserv *)method;

  YAC_ASSERT(
    yac_interp_grid_get_num_src_fields(interp_grid) == 1,
    "ERROR(do_search_conserv): invalid number of source fields")

  YAC_ASSERT(
    yac_interp_grid_get_src_field_location(interp_grid, 0) == CELL,
    "ERROR(do_search_conserv): unsupported source field location type")

  YAC_ASSERT(
    yac_interp_grid_get_tgt_field_location(interp_grid) == CELL,
    "ERROR(do_search_conserv): unsupported target field location type")


  size_t * src_cells = NULL;
  size_t * num_src_per_tgt = xmalloc(count * sizeof(*num_src_per_tgt));

  // search matching cells
  yac_interp_grid_do_cell_search_src(
    interp_grid, tgt_points, count, &src_cells, num_src_per_tgt);

  // we did a search on the interp_grid, therefore we have to re-get the basic
  // grid data
  struct const_basic_grid_data * tgt_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_tgt(interp_grid);
  struct const_basic_grid_data * src_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_src(interp_grid);

  size_t total_num_weights = 0;
  size_t max_num_src_per_tgt = 0;
  for (size_t i = 0; i < count; ++i) {
    size_t curr_num_src_per_tgt = num_src_per_tgt[i];
    if (curr_num_src_per_tgt > max_num_src_per_tgt)
      max_num_src_per_tgt = curr_num_src_per_tgt;
    total_num_weights += num_src_per_tgt[i];
  }

  // to ensure that the interpolation always procduces the same result, we
  // sort the source cells for each target point by their global ids
  {
    yac_int * temp_src_global_ids =
      xmalloc(max_num_src_per_tgt * sizeof(*temp_src_global_ids));

    for (size_t i = 0, offset = 0; i < count; ++i) {

      size_t curr_num_src_per_tgt = num_src_per_tgt[i];
      size_t * curr_src_cells = src_cells + offset;
      offset += curr_num_src_per_tgt;

      for (size_t j = 0; j < curr_num_src_per_tgt; ++j)
        temp_src_global_ids[j] =
          src_basic_grid_data->ids[CELL][curr_src_cells[j]];

      yac_quicksort_index_yac_int_size_t(
        temp_src_global_ids, curr_num_src_per_tgt, curr_src_cells);
    }

    free(temp_src_global_ids);
  }

  double * w = xmalloc(total_num_weights * sizeof(*w));
  size_t result_count = 0;
  size_t * failed_tgt = xmalloc(count * sizeof(*failed_tgt));
  total_num_weights = 0;

  int partial_coverage = method_conserv->partial_coverage;
  enum yac_interp_method_conserv_normalisation normalisation =
    method_conserv->normalisation;
  int enforced_conserv = method_conserv->enforced_conserv;

  struct grid_cell tgt_grid_cell;
  struct grid_cell * src_grid_cells;
  get_cell_buffers(
    interp_grid, max_num_src_per_tgt, &tgt_grid_cell, &src_grid_cells);

  // compute overlaps
  for (size_t i = 0, offset = 0, result_offset = 0; i < count; ++i) {

    size_t curr_src_count = num_src_per_tgt[i];
    size_t curr_tgt_point = tgt_points[i];
    size_t num_weights;

    // if weight computation was successful
    if (compute_1st_order_weights(
          tgt_basic_grid_data, curr_tgt_point,
          src_basic_grid_data, curr_src_count, src_cells + offset,
          tgt_grid_cell, src_grid_cells, w + result_offset, &num_weights,
          partial_coverage, normalisation, enforced_conserv)) {

      if (offset != result_offset) {

        memmove(
          src_cells + result_offset, src_cells + offset,
          num_weights * sizeof(*src_cells));
      }
      tgt_points[result_count] = curr_tgt_point;
      num_src_per_tgt[result_count] = num_weights;
      result_count++;
      result_offset += num_weights;
      total_num_weights += num_weights;
    } else {
      failed_tgt[i - result_count] = curr_tgt_point;
    }

    offset += curr_src_count;
  }

  free(tgt_grid_cell.edge_type);
  free(tgt_grid_cell.coordinates_xyz);
  free(src_grid_cells);

  if (result_count != count)
    memcpy(tgt_points + result_count, failed_tgt,
           (count - result_count) * sizeof(*tgt_points));
  free(failed_tgt);

  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(
        interp_grid, tgt_points, result_count),
    .count = result_count};
  struct remote_point * srcs =
    yac_interp_grid_get_src_remote_points(
      interp_grid, 0, src_cells, total_num_weights);

  // store weights
  yac_interp_weights_add_wsum(
    weights, &tgts, num_src_per_tgt, srcs, w);

  free(tgts.data);
  free(srcs);
  free(src_cells);
  free(num_src_per_tgt);
  free(w);

  return result_count;
}

static int
compare_supermesh_cell_src_local_ids(const void * a, const void * b) {

  struct supermesh_cell * a_ = (struct supermesh_cell *)a;
  struct supermesh_cell * b_ = (struct supermesh_cell *)b;

  int ret = (a_->src.local_id > b_->src.local_id) -
            (a_->src.local_id < b_->src.local_id);
  if (ret) return ret;
  return (a_->tgt.global_id > b_->tgt.global_id) -
         (a_->tgt.global_id < b_->tgt.global_id);
}

static int
compare_supermesh_cell_tgt_local_ids(const void * a, const void * b) {

  struct supermesh_cell * a_ = (struct supermesh_cell *)a;
  struct supermesh_cell * b_ = (struct supermesh_cell *)b;

  int ret = (a_->tgt.local_id > b_->tgt.local_id) -
            (a_->tgt.local_id < b_->tgt.local_id);
  if (ret) return ret;
  return (a_->src.global_id > b_->src.global_id) -
         (a_->src.global_id < b_->src.global_id);
}

static inline void orthogonalise_weight_vector(
  double M[3][3], struct weight_vector_3d * G_i,
  struct weight_vector_data_3d * buffer) {

  // This routine computes: (I_3 - C_i * C_i^-1) * G_i
  //
  // O(g_i) = g_i - C_i * (C_i^-1 * g_i)
  // where: C_i is the centeroid of a source cell
  //        g_i is the gradient of the source field in C_i
  //        O(g_i) is the projection of g_i into the plane perpendicular to C_i
  // g_i = G_i * f
  // where: G_i is the weight matrix to compute g_i
  //        f is the source field vector
  // => O(g_i) = (I_3 - C_i * C_i^-1) * G_i * f
  // where: I_3 is the identity matrix of size 3 x 3
  //
  // M = I_3 - C_i * C_i^-1

  struct weight_vector_data_3d * G_i_data = G_i->data;

  size_t N = G_i->n;
  for (size_t i = 0; i < N; ++i) {
    buffer[i].local_id = G_i_data[i].local_id;
    buffer[i].global_id = G_i_data[i].global_id;
    for (size_t j = 0; j < 3; ++j) buffer[i].weight[j] = 0.0;
  }

  for (size_t n = 0; n < N; ++n)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        buffer[n].weight[i] += G_i_data[n].weight[j] * M[i][j];

  memcpy(G_i->data, buffer, N * sizeof(*buffer));
}

static int compare_weight_vector_data_weight(
  void const * a, void const * b) {

  struct weight_vector_data const * weight_a =
    (struct weight_vector_data const *)a;
  struct weight_vector_data const * weight_b =
    (struct weight_vector_data const *)b;

  int ret = weight_a->global_id - weight_b->global_id;
  if (ret) return ret;
  double abs_weight_a = fabs(weight_a->weight);
  double abs_weight_b = fabs(weight_b->weight);
  ret = (abs_weight_a > abs_weight_b) - (abs_weight_a < abs_weight_b);
  if (ret) return ret;
  return (weight_a->weight > weight_b->weight) -
         (weight_a->weight < weight_b->weight);
}

static int compare_weight_vector_data(
  void const * a, void const * b) {

  struct weight_vector_data const * weight_a =
    (struct weight_vector_data const *)a;
  struct weight_vector_data const * weight_b =
    (struct weight_vector_data const *)b;

  return weight_a->global_id - weight_b->global_id;
}

static void compact_weight_vector_data(
  struct weight_vector_data * weights, size_t * n) {

  size_t n_ = *n;

  if (n_ <= 1) return;

  // sort weights by global_id then by weight
  qsort(weights, n_, sizeof(*weights), compare_weight_vector_data_weight);

  size_t new_n = 1;
  struct weight_vector_data * prev_weight_data = weights;
  struct weight_vector_data * curr_weight_data = weights + 1;
  for (size_t i = 1; i < n_; ++i, ++curr_weight_data) {

    // if both weights refer to the same source point (by global_id)
    if (!compare_weight_vector_data(prev_weight_data, curr_weight_data)) {
      prev_weight_data->weight += curr_weight_data->weight;
    } else {
      ++new_n;
      ++prev_weight_data;
      *prev_weight_data = *curr_weight_data;
    }
  }

  *n = new_n;
}

static size_t compute_2nd_order_tgt_cell_weights(
  struct supermesh_cell * super_cell, struct weight_vector_data * weights) {

  weights[0].weight       = super_cell->norm_area;
  weights[0].global_id    = super_cell->src.global_id;
  weights[0].local_id     = super_cell->src.local_id;

  struct weight_vector_3d * src_cell_gradient = super_cell->src_cell_gradient;

  size_t N = src_cell_gradient->n;
  // in case we have no gradient for the current supermesh cell,
  // we assume a constant field across the whole associated source cell
  if (N > 0) {
    struct weight_vector_data_3d * gradient_weights = src_cell_gradient->data - 1;
    double * overlap_barycenter = super_cell->barycenter;

    for (size_t n = 1; n <= N; ++n) {
      weights[n].weight       =
        (gradient_weights[n].weight[0] * overlap_barycenter[0] +
         gradient_weights[n].weight[1] * overlap_barycenter[1] +
         gradient_weights[n].weight[2] * overlap_barycenter[2]) *
        super_cell->norm_area;
      weights[n].global_id    = gradient_weights[n].global_id;
      weights[n].local_id     = gradient_weights[n].local_id;
    }
  }

  return 1 + N;
}

static void compute_cell_barycenter(
  struct const_basic_grid_data * grid_data, size_t cell_idx,
  double barycenter[3]) {

  size_t num_vertices = grid_data->num_vertices_per_cell[cell_idx];
  size_t const * vertices =
    grid_data->cell_to_vertex + grid_data->cell_to_vertex_offsets[cell_idx];

  barycenter[0] = 0.0;
  barycenter[1] = 0.0;
  barycenter[2] = 0.0;

  for (size_t i = 0; i < num_vertices; ++i) {
    double const * curr_vertex_coordinate =
      grid_data->vertex_coordinates[vertices[i]];
    barycenter[0] += curr_vertex_coordinate[0];
    barycenter[1] += curr_vertex_coordinate[1];
    barycenter[2] += curr_vertex_coordinate[2];
  }
  normalise_vector(barycenter);
}

static void compute_super_cells(
  struct interp_grid * interp_grid, size_t * tgt_points, size_t count,
  struct supermesh_cell ** super_cells_, size_t * num_super_cells,
  int * interp_fail_flag, size_t ** src_cells, size_t * num_src_cells,
  enum yac_interp_method_conserv_normalisation normalisation,
  int partial_coverage) {

  YAC_ASSERT(
    (normalisation == CONSERV_DESTAREA) || (normalisation == CONSERV_FRACAREA),
    "ERROR(compute_super_cells): "
    "invalid normalisation option in conservative remapping")

  size_t * num_src_per_tgt = xmalloc(count * sizeof(*num_src_per_tgt));

  // search for all source cell overlapping with the the target cells
  yac_interp_grid_do_cell_search_src(
    interp_grid, tgt_points, count, src_cells, num_src_per_tgt);

  // determine the number of unique matching source cells
  size_t total_num_overlaps = 0;
  for (size_t i = 0; i < count; ++i) total_num_overlaps += num_src_per_tgt[i];
  yac_quicksort_index_size_t_int(*src_cells, total_num_overlaps, NULL);
  *num_src_cells = total_num_overlaps;
  yac_remove_duplicates_size_t(*src_cells, num_src_cells);

  size_t * num_tgt_per_src =
    xrealloc(num_src_per_tgt, *num_src_cells * sizeof(*num_tgt_per_src));

  // for some required source cells we may have not all required supermesh cells
  // in that case we have to find the respective target cells and compute the
  // missing supermesh cells from them
  size_t * tgt_cells = NULL;
  yac_interp_grid_do_cell_search_tgt(
    interp_grid, *src_cells, *num_src_cells, &tgt_cells, num_tgt_per_src);

  total_num_overlaps = 0;
  for (size_t i = 0; i < *num_src_cells; ++i)
    total_num_overlaps += num_tgt_per_src[i];

  struct supermesh_cell * super_cells =
    xmalloc(total_num_overlaps * sizeof(*super_cells));

  struct const_basic_grid_data * src_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_src(interp_grid);
  struct const_basic_grid_data * tgt_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_tgt(interp_grid);

  for (size_t i = 0, j = 0; i < *num_src_cells; ++i) {

    size_t curr_num_overlaps = num_tgt_per_src[i];
    size_t curr_src_cell = (*src_cells)[i];
    yac_int curr_src_global_id =
      src_basic_grid_data->ids[CELL][curr_src_cell];

    for (size_t k = 0; k < curr_num_overlaps; ++k, ++j) {

      size_t curr_tgt_cell = tgt_cells[j];
      struct supermesh_cell * curr_super_cell = super_cells + j;
      curr_super_cell->src.local_id = curr_src_cell;
      curr_super_cell->src.global_id = curr_src_global_id;
      curr_super_cell->tgt.local_id = curr_tgt_cell;
      curr_super_cell->tgt.global_id =
        tgt_basic_grid_data->ids[CELL][curr_tgt_cell];
    }
  }
  free(tgt_cells);
  free(num_tgt_per_src);

  // sort supermesh_cell first by local ids of the target cells and
  // second by global cell id
  qsort(super_cells, total_num_overlaps, sizeof(*super_cells),
        compare_supermesh_cell_tgt_local_ids);

  struct grid_cell tgt_grid_cell;
  struct grid_cell src_grid_cell;
  get_cell_buffers_(interp_grid, &tgt_grid_cell, &src_grid_cell);

  src_basic_grid_data = yac_interp_grid_get_basic_grid_data_src(interp_grid);
  tgt_basic_grid_data = yac_interp_grid_get_basic_grid_data_tgt(interp_grid);

  // For all supermesh cells compute the area and normalised area.
  // Additionally, remove all empty supermesh cells.
  size_t new_num_super_cells = 0;
  size_t tgt_idx = 0;
  for (size_t i = 0, j = 0; i < total_num_overlaps;) {

    // get information about the current target cell
    size_t curr_tgt_cell = super_cells[i].tgt.local_id;
    yac_const_basic_grid_data_get_grid_cell(
      tgt_basic_grid_data, curr_tgt_cell, &tgt_grid_cell);
    double curr_tgt_cell_coverage = 0.0;

    // for all supermesh cells overlapping with the current target cell
    for (;(i < total_num_overlaps) &&
           (super_cells[i].tgt.local_id == curr_tgt_cell); ++i)  {

      // get the current source cell
      yac_const_basic_grid_data_get_grid_cell(
        src_basic_grid_data, super_cells[i].src.local_id, &src_grid_cell);

      // compute area of the current supermesh cell
      double super_cell_area;
      double barycenter[3];
      yac_compute_overlap_info(
        1, &src_grid_cell, tgt_grid_cell, &super_cell_area, &barycenter);

      // if there is an overlap between the current source and target cell
      if (super_cell_area > 0.0) {

        super_cells[new_num_super_cells].src = super_cells[i].src;
        super_cells[new_num_super_cells].tgt = super_cells[i].tgt;
        super_cells[new_num_super_cells].area = super_cell_area;
        memcpy(super_cells[new_num_super_cells].barycenter, barycenter,
               3 * sizeof(double));
        super_cells[new_num_super_cells].src_cell_gradient = NULL;
        ++new_num_super_cells;

        curr_tgt_cell_coverage += super_cell_area;
      }
    }

    double curr_tgt_cell_area = yac_huiliers_area(tgt_grid_cell);

    // if there was an overlap
    if (new_num_super_cells != j) {

      double norm_factor;
      YAC_ASSERT(
        (normalisation == CONSERV_DESTAREA) ||
        (normalisation == CONSERV_FRACAREA),
        "ERROR(compute_super_cells): invalid normalisation")
      switch (normalisation) {
        default:
        case (CONSERV_DESTAREA):
          norm_factor = 1.0 / curr_tgt_cell_area;
          break;
        case (CONSERV_FRACAREA):
          norm_factor = 1.0 / curr_tgt_cell_coverage;
          break;
      }
      // compute normalised area
      for (; j < new_num_super_cells; ++j)
        super_cells[j].norm_area = super_cells[j].area * norm_factor;
    }


    // for target cell that do not overlap with any source cell
    while ((tgt_idx < count) && (tgt_points[tgt_idx] < curr_tgt_cell))
      interp_fail_flag[tgt_idx++] = 1;

    if ((tgt_idx < count) && (tgt_points[tgt_idx] == curr_tgt_cell)) {

      double area_tol = curr_tgt_cell_area * AREA_TOL_FACTOR;

      if (partial_coverage) {
        interp_fail_flag[tgt_idx] = curr_tgt_cell_coverage < area_tol;
      } else {
        interp_fail_flag[tgt_idx] =
          fabs(curr_tgt_cell_area - curr_tgt_cell_coverage) > area_tol;
      }
      ++tgt_idx;
    }
  }
  // for all remaining target cell that do not overlap with any source cell
  for (; tgt_idx < count; ++tgt_idx) interp_fail_flag[tgt_idx] = 1;
  *num_super_cells = new_num_super_cells;
  *super_cells_ =
    xrealloc(super_cells, new_num_super_cells * sizeof(*super_cells));
  free(tgt_grid_cell.coordinates_xyz);
  free(tgt_grid_cell.edge_type);
}

static coordinate_pointer compute_src_cell_centroids(
  struct interp_grid * interp_grid,
  size_t * src_cells, int * skip_src_cell, size_t num_src_cells,
  struct supermesh_cell * super_cells, size_t num_super_cells) {

  coordinate_pointer src_cell_centroids =
    xmalloc(num_src_cells * sizeof(*src_cell_centroids));

  // sort supermesh cells first by source local id and second by
  // target global id
  qsort(super_cells, num_super_cells, sizeof(*super_cells),
        compare_supermesh_cell_src_local_ids);

  struct grid_cell src_grid_cell, dummy;
  get_cell_buffers_(interp_grid, &src_grid_cell, &dummy);

  struct const_basic_grid_data * src_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_src(interp_grid);

  // compute centroids of source cells
  // C_i = N(S_(U_k) (A_k*C_k))
  // where: N(C) = C*(C*C)^-0.5 // normalisation
  //        U_k all supermesh cells overlaping with the respective source cell
  //        A_k area of supermesh cell
  //        C_k barycenter of supermesh cell (normalisation of the sum of all
  //                                          vertices of the cell)
  for (size_t i = 0, offset = 0; i < num_src_cells; ++i) {

    if (skip_src_cell[i]) continue;

    size_t curr_src_cell = src_cells[i];
    yac_const_basic_grid_data_get_grid_cell(
      src_basic_grid_data, curr_src_cell, &src_grid_cell);
    double src_cell_area = yac_huiliers_area(src_grid_cell);

    struct supermesh_cell * curr_super_cells = super_cells + offset;
    size_t curr_num_super_cells = offset;
    while ((offset < num_super_cells) &&
           (super_cells[offset].src.local_id == curr_src_cell)) ++offset;
    curr_num_super_cells = offset - curr_num_super_cells;

    double src_cell_area_diff = src_cell_area;
    double src_cell_centroid[3] = {0.0, 0.0, 0.0};

    for (size_t j = 0; j < curr_num_super_cells; ++j) {

      double super_cell_area = curr_super_cells[j].area;
      double * super_cell_barycenter = curr_super_cells[j].barycenter;
      src_cell_centroid[0] += super_cell_area * super_cell_barycenter[0];
      src_cell_centroid[1] += super_cell_area * super_cell_barycenter[1];
      src_cell_centroid[2] += super_cell_area * super_cell_barycenter[2];
      src_cell_area_diff -= super_cell_area;
    }

    if (!(skip_src_cell[i] =
          (fabs(src_cell_area_diff) > src_cell_area * AREA_TOL_FACTOR))) {

      normalise_vector(src_cell_centroid);
      src_cell_centroids[i][0] = src_cell_centroid[0];
      src_cell_centroids[i][1] = src_cell_centroid[1];
      src_cell_centroids[i][2] = src_cell_centroid[2];
    } /* else {
      MoHa: Does it make sense to compute a centroid only based on the vertices
            of the source cell as a fallback?
    }*/

  }
  free(src_grid_cell.coordinates_xyz);
  free(src_grid_cell.edge_type);

  return src_cell_centroids;
}

static struct weight_vector_3d * compute_src_cell_gradients(
  struct interp_grid * interp_grid, size_t * src_cells,
  coordinate_pointer src_cell_centroids, int * skip_src_cell,
  size_t num_src_cells, size_t * src_cell_neighbours) {

  struct weight_vector_3d * src_cell_gradients =
    xmalloc(num_src_cells * sizeof(*src_cell_gradients));

  struct const_basic_grid_data * src_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_src(interp_grid);

  size_t total_num_gradient_weights = 0;
  size_t max_num_neigh_per_src = 0;
  for (size_t i = 0, offset = 0; i < num_src_cells; ++i) {
    src_cell_gradients[i].data = NULL;
    src_cell_gradients[i].n = 0;
    size_t curr_num_neigh =
      src_basic_grid_data->num_vertices_per_cell[src_cells[i]];
    size_t * curr_neighs = src_cell_neighbours + offset;
    offset += curr_num_neigh;
    for (size_t j = 0; (!skip_src_cell[i]) && (j < curr_num_neigh); ++j)
      skip_src_cell[i] |= curr_neighs[j] == SIZE_MAX;
    if (skip_src_cell[i]) continue;
    total_num_gradient_weights += curr_num_neigh + 1;
    if (max_num_neigh_per_src < curr_num_neigh)
      max_num_neigh_per_src = curr_num_neigh;
  }
  struct weight_vector_data_3d * weight_vector_data_buffer =
    (total_num_gradient_weights > 0)?
      (xmalloc(total_num_gradient_weights *
               sizeof(*weight_vector_data_buffer))):NULL;
  for (size_t i = 0; i < total_num_gradient_weights; ++i)
    for (size_t j = 0; j < 3; ++j)
      weight_vector_data_buffer[i].weight[j] = 0.0;

  struct weight_vector_data_3d * orth_buffer =
    xmalloc((max_num_neigh_per_src + 1) * sizeof(*orth_buffer));

  // compute gradient in the centroid for each source cell
  // g_i = O_i(g'_i)
  // where: O_i(g) = g - C_i*(C_i*g) // makes g orthogonal to C_i
  // g'_i = (A_(C_k))^-1 * S_ijk(((f_j + f_k) * 0.5 - f_i) *
  //                         (C_j x C_k) * |C_j x C_k|^-1 *
  //                         asin(|C_j x C_k|))
  // where: g'_i is the estimated centroid gradient
  //        A_(C_k) is the area of the cell generated by connecting the
  //                barycenters of all neighbour cells
  //        S_ijk is the sum of all edges of the previously described cell in
  //              counterclockwise order
  //        f_i, f_j, and f_k are the mean values of the field over the area of
  //                          the respective source cell area
  // asin(|C_j x C_k|) / |C_j x C_k| ~ 1.0
  // => g'_i = S_ijk(((f_j + f_k) * 0.5 - f_i) * (C_j x C_k)) / A_(C_k)
  //
  // g_i = G_i * f
  // G_i = S_ijk(((I_j + I_k) * 0.5 - I_i) * (C_j x C_k)) / A_(C_k)
  // where: G_i is the gradient weight matrix
  //        I_x is a vector of the same size as f, which is all zero except at
  //            position x, where it is one
  //        f is the source field vector
  for (size_t i = 0, offset = 0, weight_vector_data_buffer_offset = 0;
       i < num_src_cells; ++i) {

    size_t curr_num_neigh =
      src_basic_grid_data->num_vertices_per_cell[src_cells[i]];
    size_t * curr_neighs = src_cell_neighbours + offset;
    offset += curr_num_neigh;
    struct weight_vector_3d * G_i = src_cell_gradients + i;
    G_i->data = weight_vector_data_buffer + weight_vector_data_buffer_offset;

    if (skip_src_cell[i]) continue;

    weight_vector_data_buffer_offset += curr_num_neigh + 1;
    G_i->n = curr_num_neigh + 1;
    G_i->data[0].local_id = src_cells[i];
    G_i->data[0].global_id = src_basic_grid_data->ids[CELL][src_cells[i]];
    for (size_t j = 0; j < curr_num_neigh; ++j) {
      G_i->data[j+1].local_id = curr_neighs[j];
      G_i->data[j+1].global_id = src_basic_grid_data->ids[CELL][curr_neighs[j]];
    }

    // area of the polygon that is formed by connecting the barycenters of
    // the neighbouring cells
    double A_C_k = 0.0;

    struct grid_cell centroid_triangle = {
      .coordinates_xyz = (double[3][3]){{0}},
      .edge_type =
        (enum yac_edge_type[]) {GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE,GREAT_CIRCLE_EDGE},
      .num_corners = 3, .array_size = 0};
    centroid_triangle.coordinates_xyz[0][0] = src_cell_centroids[i][0];
    centroid_triangle.coordinates_xyz[0][1] = src_cell_centroids[i][1];
    centroid_triangle.coordinates_xyz[0][2] = src_cell_centroids[i][2];

    // the following matrix is used by the routine
    // orthogonalise_weight_vector
    // M = I_3 - C_i * C_i^-1
    double M[3][3];
    for (size_t k = 0; k < 3; ++k)
      for (size_t l = 0; l < 3; ++l)
        M[k][l] = - centroid_triangle.coordinates_xyz[0][k] *
                    centroid_triangle.coordinates_xyz[0][l];
    for (size_t k = 0; k < 3; ++k)
      M[k][k] += 1.0;

    // We split the cell that is comprised of the barycenters of the edge
    // neigbours into triangles, which has the centroid of the current cell
    // as one corner. The sum of the areas of these triangles is A_C_K.
    for (size_t j = 0; j < curr_num_neigh; ++j) {

      size_t curr_neighbours[2];
      curr_neighbours[0] = curr_neighs[j];
      curr_neighbours[1] = curr_neighs[(j+1)%curr_num_neigh];

      double neigh_cell_barycenters[2][3];
      compute_cell_barycenter(
        src_basic_grid_data, curr_neighbours[0], neigh_cell_barycenters[0]);
      compute_cell_barycenter(
        src_basic_grid_data, curr_neighbours[1], neigh_cell_barycenters[1]);

      centroid_triangle.coordinates_xyz[1][0] = neigh_cell_barycenters[0][0];
      centroid_triangle.coordinates_xyz[1][1] = neigh_cell_barycenters[0][1];
      centroid_triangle.coordinates_xyz[1][2] = neigh_cell_barycenters[0][2];
      centroid_triangle.coordinates_xyz[2][0] = neigh_cell_barycenters[1][0];
      centroid_triangle.coordinates_xyz[2][1] = neigh_cell_barycenters[1][1];
      centroid_triangle.coordinates_xyz[2][2] = neigh_cell_barycenters[1][2];

      A_C_k += yac_huiliers_area(centroid_triangle);

      // C_j x C_k
      double C_j_x_C_k[3];
      crossproduct_kahan(
        neigh_cell_barycenters[0], neigh_cell_barycenters[1], C_j_x_C_k);
      // this factor ensures that C_j_x_C_k points outwards of the centroid
      // cell
      double factor =
        -copysign(1.0, C_j_x_C_k[0] * centroid_triangle.coordinates_xyz[0][0] +
                       C_j_x_C_k[1] * centroid_triangle.coordinates_xyz[0][1] +
                       C_j_x_C_k[2] * centroid_triangle.coordinates_xyz[0][2]);

      for (size_t l = 0; l < 3; ++l) C_j_x_C_k[l] *= factor;

      // -I_i * (C_j x C_k)
      G_i->data[0].weight[0] -= C_j_x_C_k[0];
      G_i->data[0].weight[1] -= C_j_x_C_k[1];
      G_i->data[0].weight[2] -= C_j_x_C_k[2];

      for (size_t l = 0; l < 3; ++l) C_j_x_C_k[l] *= 0.5;

      // 0.5 * I_j * (C_j x C_k)
      G_i->data[j+1].weight[0] += C_j_x_C_k[0];
      G_i->data[j+1].weight[1] += C_j_x_C_k[1];
      G_i->data[j+1].weight[2] += C_j_x_C_k[2];

      // 0.5 * I_k * (C_j x C_k)
      size_t j_ = (j + 1) % curr_num_neigh + 1;
      G_i->data[j_].weight[0] += C_j_x_C_k[0];
      G_i->data[j_].weight[1] += C_j_x_C_k[1];
      G_i->data[j_].weight[2] += C_j_x_C_k[2];
    }

    // ((I_j + I_k) * 0.5 - I_i) * (C_j x C_k) / A_(C_k)
    double inv_A_C_K = 1.0 / A_C_k;
    for (size_t k = 0; k <= curr_num_neigh; ++k)
      for (size_t l = 0; l < 3; ++l)
        G_i->data[k].weight[l] *= inv_A_C_K;

    orthogonalise_weight_vector(M, G_i, orth_buffer);
  }
  free(orth_buffer);

  return src_cell_gradients;
}

static size_t compute_2nd_order_weights(
  struct interp_grid * interp_grid, size_t * tgt_cells, int * interp_fail_flag,
  size_t num_tgt_cells, struct supermesh_cell * super_cells,
  size_t num_super_cells, size_t ** src_per_tgt, double ** weights,
  size_t * num_src_per_tgt) {

  size_t num_interpolated_tgt = 0;

  // sort supermesh cells first by target local id and second by
  // source global id
  qsort(super_cells, num_super_cells, sizeof(*super_cells),
        compare_supermesh_cell_tgt_local_ids);

  size_t max_num_weights_per_tgt = 0;
  size_t max_num_total_weights = 0;

  // count maximum number of weights
  for (size_t i = 0, j = 0; i < num_tgt_cells; ++i) {

    if (interp_fail_flag[i]) continue;

    size_t curr_tgt_cell = tgt_cells[i];
    size_t curr_num_weights = 0;

    // skip supermesh cells not overlapping with current target cell
    while ((j < num_super_cells) &&
           (super_cells[j].tgt.local_id < curr_tgt_cell)) ++j;

    // for all supermesh cells overlapping with the current target cell
    while ((j < num_super_cells) &&
           (super_cells[j].tgt.local_id == curr_tgt_cell))
      curr_num_weights += 1 + super_cells[j++].src_cell_gradient->n;

    max_num_total_weights += curr_num_weights;
    if (max_num_weights_per_tgt < curr_num_weights)
      max_num_weights_per_tgt = curr_num_weights;
    num_interpolated_tgt++;
  }

  struct weight_vector_data * weight_buffer =
    xmalloc(max_num_weights_per_tgt * sizeof(*weight_buffer));
  *weights = xmalloc(max_num_total_weights * sizeof(**weights));
  *src_per_tgt = xmalloc(max_num_total_weights * sizeof(**src_per_tgt));

  // sort all target points that can be interpolated to the beginning of the
  // tgt_cells array
  yac_quicksort_index_int_size_t(interp_fail_flag, num_tgt_cells, tgt_cells);
  yac_quicksort_index_size_t_int(tgt_cells, num_interpolated_tgt, NULL);

  // compute 2nd order weights
  // f_j = SUM(f_k') / A_j
  // f_k' = A_k*(f_i + g_i * C_k)
  //
  // f_k' = A_k * (I_i + C_k * G_i) * f
  // where: I_i is a vector of the same length as f, that is all 0.0 except at
  //        position i, where it is 1.0
  //        f is the source field vector
  // => f_k' / A_j = M_k * f
  // where: M_k = A_k / A_j * (I_i + C_k * G_i)
  // => f_j = SUM(M_k) * f
  //    f_j = M_j * f
  // where: M_j = SUM(M_k)
  size_t w_idx = 0;
  for (size_t i = 0, j = 0; i < num_interpolated_tgt; ++i) {

    size_t curr_tgt_cell = tgt_cells[i];

    // skip supermesh cells not overlapping with current target cell
    while ((j < num_super_cells) &&
           (super_cells[j].tgt.local_id != curr_tgt_cell)) ++j;

    size_t weight_buffer_offset = 0;

    // for all supermesh cells overlapping with the current target cell
    while ((j < num_super_cells) && 
           (super_cells[j].tgt.local_id == curr_tgt_cell)) {

      size_t num_weights =
        compute_2nd_order_tgt_cell_weights(
          super_cells + j, weight_buffer + weight_buffer_offset);
      weight_buffer_offset += num_weights;
      ++j;
    }

    // compact weights and sort them by global_id
    compact_weight_vector_data(weight_buffer, &weight_buffer_offset);

    num_src_per_tgt[i] = weight_buffer_offset;

    for (size_t j = 0; j < weight_buffer_offset; ++j, ++w_idx) {
      (*weights)[w_idx] = weight_buffer[j].weight;
      (*src_per_tgt)[w_idx] = weight_buffer[j].local_id;
    }
  }
  free(weight_buffer);

  return num_interpolated_tgt;
}

static size_t do_search_conserv_2nd_order (struct interp_method * method,
                                           struct interp_grid * interp_grid,
                                           size_t * tgt_points, size_t count,
                                           struct interp_weights * weights) {

  struct interp_method_conserv * method_conserv =
    (struct interp_method_conserv *)method;

  YAC_ASSERT(
    yac_interp_grid_get_num_src_fields(interp_grid) == 1,
    "ERROR(do_search_conserv): invalid number of source fields")

  YAC_ASSERT(
    yac_interp_grid_get_src_field_location(interp_grid, 0) == CELL,
    "ERROR(do_search_conserv): unsupported source field location type")

  YAC_ASSERT(
    yac_interp_grid_get_tgt_field_location(interp_grid) == CELL,
    "ERROR(do_search_conserv): unsupported target field location type")

  // sort target points
  yac_quicksort_index_size_t_int(tgt_points, count, NULL);

  int * interp_fail_flag = xmalloc(count * sizeof(*interp_fail_flag));
  struct supermesh_cell * super_cells = NULL;
  size_t num_super_cells = 0;
  size_t * src_cells = NULL;
  size_t num_src_cells = 0;

  // compute the overlaps between the target cells and the source grid
  compute_super_cells(
    interp_grid, tgt_points, count, &super_cells, &num_super_cells,
    interp_fail_flag, &src_cells, &num_src_cells,
    method_conserv->normalisation, method_conserv->partial_coverage);

  size_t total_num_src_cell_neighbours = 0;
  struct const_basic_grid_data * src_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_src(interp_grid);
  for (size_t i = 0; i < num_src_cells; ++i)
    total_num_src_cell_neighbours +=
      src_basic_grid_data->num_vertices_per_cell[src_cells[i]];
  size_t * src_cell_neighbours =
    xmalloc(total_num_src_cell_neighbours * sizeof(*src_cell_neighbours));

  // get the neighbours for all soure cells overlapping with a target cell
  yac_interp_grid_get_src_cell_neighbours(
    interp_grid, src_cells, num_src_cells, src_cell_neighbours);

  const_int_pointer src_cell_mask =
    yac_interp_grid_get_src_field_mask(interp_grid, 0);
  int * skip_src_cell = xmalloc(num_src_cells * sizeof(*skip_src_cell));

  // check whether any required source cell or any of its neighbours is masked
  if (src_cell_mask != NULL) {
    for (size_t i = 0, offset = 0; i < num_src_cells; ++i) {
      size_t curr_src_cell = src_cells[i];
      if ((skip_src_cell[i] = !src_cell_mask[curr_src_cell])) continue;
      size_t * curr_neighbours = src_cell_neighbours + offset;
      size_t curr_num_neigh =
        src_basic_grid_data->num_vertices_per_cell[curr_src_cell];
      offset += curr_num_neigh;
      for (size_t j = 0; j < curr_num_neigh; ++j)
        skip_src_cell[i] |=
          (curr_neighbours[j] == SIZE_MAX) || (!src_cell_mask[curr_neighbours[j]]);
    }
  } else {
    memset(skip_src_cell, 0, num_src_cells * sizeof(*skip_src_cell));
  }

  // compute the centroids of all required source cells
  coordinate_pointer src_cell_centroids =
    compute_src_cell_centroids(interp_grid, src_cells, skip_src_cell,
                               num_src_cells, super_cells, num_super_cells);

  // compute the gradients weights for the source cells
  struct weight_vector_3d * src_cell_gradients =
    compute_src_cell_gradients(
      interp_grid, src_cells, src_cell_centroids,
      skip_src_cell, num_src_cells, src_cell_neighbours);
  free(src_cell_neighbours);
  free(src_cell_centroids);

  // sort supermesh cells by local src id
  qsort(super_cells, num_super_cells, sizeof(*super_cells),
        compare_supermesh_cell_src_local_ids);
  for (size_t i = 0, j = 0; i < num_src_cells; ++i) {
    size_t curr_src_cell = src_cells[i];
    struct weight_vector_3d * curr_src_cell_gradient = src_cell_gradients + i;
    // there should not be any supercell whose source cell is not in the list of
    // required source cells
    YAC_ASSERT(
      (j >= num_super_cells) ||
      (super_cells[j].src.local_id >= curr_src_cell),
      "ERROR(do_search_conserv_2nd_order): internal error");
    while ((j < num_super_cells) &&
           (super_cells[j].src.local_id == curr_src_cell)) {
      super_cells[j++].src_cell_gradient = curr_src_cell_gradient;
    }
  }
  free(skip_src_cell);
  free(src_cells);

  size_t * src_per_tgt = NULL;
  double * w = NULL;
  size_t * num_src_per_tgt = xmalloc(count * sizeof(*num_src_per_tgt));

  // compute the weights for the target points
  size_t num_interpolated_tgt =
    compute_2nd_order_weights(
      interp_grid, tgt_points, interp_fail_flag, count,
      super_cells, num_super_cells, &src_per_tgt, &w, num_src_per_tgt);
  if (num_src_cells > 0) free(src_cell_gradients->data);
  free(src_cell_gradients);
  free(super_cells);
  free(interp_fail_flag);

  size_t total_num_weights = 0;
  for (size_t i = 0; i < num_interpolated_tgt; ++i)
    total_num_weights += num_src_per_tgt[i];

  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(
        interp_grid, tgt_points, num_interpolated_tgt),
    .count = num_interpolated_tgt};
  struct remote_point * srcs =
    yac_interp_grid_get_src_remote_points(
      interp_grid, 0, src_per_tgt, total_num_weights);

  // store weights
  yac_interp_weights_add_wsum(
    weights, &tgts, num_src_per_tgt, srcs, w);

  free(tgts.data);
  free(srcs);
  free(w);
  free(num_src_per_tgt);
  free(src_per_tgt);

  return num_interpolated_tgt;
}

struct interp_method * yac_interp_method_conserv_new(
  int order, int enforced_conserv, int partial_coverage,
  enum yac_interp_method_conserv_normalisation normalisation) {

  struct interp_method_conserv * method = xmalloc(1 * sizeof(*method));

  YAC_ASSERT(
    (enforced_conserv != 1) || (order == 1),
    "ERROR(yac_interp_method_conserv_new): interp_method_conserv only "
    "supports enforced_conserv with first order conservative remapping")
  YAC_ASSERT(
    (order == 1) || (order == 2),
    "ERROR(yac_interp_method_conserv_new): invalid order")

  method->vtable =
    (order == 1)?
      &interp_method_conserv_1st_order_vtable:
      &interp_method_conserv_2nd_order_vtable;
  method->partial_coverage = partial_coverage;
  method->normalisation = normalisation;
  method->enforced_conserv = enforced_conserv;

  return (struct interp_method*)method;
}

static void delete_conserv(struct interp_method * method) {
  free(method);
}
