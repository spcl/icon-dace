// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <string.h>

#define WEIGHT_TOL (1e-9)

#include "yac_mpi_internal.h"
#include "interp_weights_internal.h"
#include "ensure_array_size.h"
#include "io_utils.h"
#include "utils_core.h"
#include "interp_method_file.h"

enum yac_interp_weight_stencil_type {
  FIXED         = 0,
  DIRECT        = 1,
  SUM           = 2,
  WEIGHT_SUM    = 3,
  DIRECT_MF     = 4,
  SUM_MF        = 5,
  WEIGHT_SUM_MF = 6,
  WEIGHT_STENCIL_TYPE_SIZE = 7, // last entry
};

struct interp_weight_stencil {

  enum yac_interp_weight_stencil_type type;
  union {
    struct {
      double value;
    } fixed;
    struct {
      struct remote_point src; // src id
    } direct;
    struct {
      struct remote_points * srcs; // src ids
    } sum;
    struct {
      struct remote_points * srcs; // src ids
      double * weights;
    } weight_sum;
    struct {
      struct remote_point src; // src id
      size_t field_idx;
    } direct_mf;
    struct {
      struct remote_points * srcs; // src ids
      size_t * field_indices;
    } sum_mf;
    struct {
      struct remote_points * srcs; // src ids
      double * weights;
      size_t * field_indices;
    } weight_sum_mf;
  } data;
  struct remote_point tgt; //tgt id
};

struct interp_weight_stencil_wsum_mf_weight {

  struct remote_point_info src;
  uint64_t src_field_idx;
  double weight;
};

struct interp_weight_stencil_wsum_mf {
  struct remote_point tgt;
  struct interp_weight_stencil_wsum_mf_weight * data;
  size_t count;
};

struct interp_weight_stencils_wsum_mf {
  struct interp_weight_stencil_wsum_mf * data;
  size_t count;
};

struct interp_weight_stencils_wsum_mf_buffer {
  struct interp_weight_stencils_wsum_mf stencils;
  struct interp_weight_stencil_wsum_mf_weight buffer[];
};

struct interp_weight_stencil_fixed {
  double value;
  uint64_t orig_pos;
};

struct interp_weight_stencil_direct {
  struct remote_point_info src;
  uint64_t orig_pos;
};

struct interp_weight_stencil_direct_mf {
  struct remote_point_info src;
  uint64_t src_field_idx;
  uint64_t orig_pos;
};

struct weighted_global_id {
  yac_int global_id;
  double weight;
};

struct yac_interp_weights {

  MPI_Comm comm;

  enum yac_location tgt_location;
  enum yac_location * src_locations;
  size_t num_src_fields;

  struct interp_weight_stencil * stencils;
  size_t stencils_array_size, stencils_size;
};


struct remote_point_info_reorder {
  struct remote_point_info data;
  size_t field_idx;
  size_t reorder_idx;
};

struct yac_interp_weights * yac_interp_weights_new(
  MPI_Comm comm, enum yac_location tgt_location,
  enum yac_location * src_locations, size_t num_src_fields) {

  struct yac_interp_weights * weights = xmalloc(1 * sizeof(*weights));

  weights->comm = comm;
  weights->tgt_location = tgt_location;
  weights->src_locations = xmalloc(num_src_fields * sizeof(*src_locations));
  memcpy(
    weights->src_locations, src_locations,
    num_src_fields * sizeof(*src_locations));
  weights->num_src_fields = num_src_fields;
  weights->stencils = NULL;
  weights->stencils_array_size = 0;
  weights->stencils_size = 0;

  return weights;
}

static inline struct remote_point copy_remote_point(
  struct remote_point point) {

  int count = point.data.count;
  if (count > 1) {
    struct remote_point_info * point_infos =
      xmalloc((size_t)count * sizeof(*point_infos));
    memcpy(point_infos, point.data.data.multi,
           (size_t)count * sizeof(*point_infos));
    point.data.data.multi = point_infos;
  }
  return point;
}

static inline void copy_remote_points_no_alloc(
  struct remote_point * points_to, struct remote_point * points_from,
  size_t count, struct remote_point_info ** point_info_buffer_) {

  struct remote_point_info * point_info_buffer = *point_info_buffer_;

  for (size_t i = 0; i < count; ++i) {
    int curr_count = points_from[i].data.count;
    points_to[i] = points_from[i];
    if (curr_count > 1) {
      points_to[i].data.data.multi = point_info_buffer;
      memcpy(
        point_info_buffer, points_from[i].data.data.multi,
        (size_t)curr_count * sizeof(*point_info_buffer));
      point_info_buffer += curr_count;
    }
  }
  *point_info_buffer_ = point_info_buffer;
}

static inline struct remote_points * copy_remote_points(
  struct remote_point * points, size_t count) {

  size_t point_info_buffer_size = 0;
  for (size_t i = 0; i < count; ++i)
    if (points[i].data.count > 1)
      point_info_buffer_size += (size_t)(points[i].data.count);

  struct remote_points * points_copy =
    xmalloc(point_info_buffer_size * sizeof(struct remote_point_info) +
            sizeof(*points_copy));
  points_copy->data = xmalloc(count * sizeof(*(points_copy->data)));
  points_copy->count = count;
  struct remote_point_info * point_info_buffer = &(points_copy->buffer[0]);

  copy_remote_points_no_alloc(
    points_copy->data, points, count, &point_info_buffer);

  return points_copy;
}

static inline struct remote_points * copy_remote_points_mf(
  struct remote_point ** points, size_t * counts, size_t num_fields) {

  size_t point_info_buffer_size = 0;
  size_t total_count = 0;
  for (size_t i = 0; i < num_fields; ++i) {
    total_count += counts[i];
    for (size_t j = 0; j < counts[i]; ++j) {
      if (points[i][j].data.count > 1)
        point_info_buffer_size += (size_t)(points[i][j].data.count);
    }
  }

  struct remote_points * points_copy =
    xmalloc(point_info_buffer_size * sizeof(struct remote_point_info) +
            sizeof(*points_copy));
  points_copy->data = xmalloc(total_count * sizeof(*(points_copy->data)));
  points_copy->count = total_count;
  struct remote_point_info * point_info_buffer = &(points_copy->buffer[0]);

  for (size_t i = 0, k = 0; i < num_fields; ++i) {
    for (size_t j = 0; j < counts[i]; ++j, ++k) {
      int curr_count = points[i][j].data.count;
      points_copy->data[k] = points[i][j];
      if (curr_count > 1) {
        points_copy->data[k].data.data.multi = point_info_buffer;
        memcpy(
          point_info_buffer, points[i][j].data.data.multi,
          (size_t)curr_count * sizeof(*point_info_buffer));
        point_info_buffer += curr_count;
      }
    }
  }
  return points_copy;
}

void yac_interp_weights_add_fixed(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  double fixed_value) {

  struct interp_weight_stencil * stencils = weights->stencils;
  size_t stencils_array_size = weights->stencils_array_size;
  size_t stencils_size = weights->stencils_size;

  ENSURE_ARRAY_SIZE(stencils, stencils_array_size, stencils_size + tgts->count);

  for (size_t i = 0; i < tgts->count; ++i, ++stencils_size) {

    stencils[stencils_size].type = FIXED;
    stencils[stencils_size].tgt = copy_remote_point(tgts->data[i]);
    stencils[stencils_size].data.fixed.value = fixed_value;
  }

  weights->stencils = stencils;
  weights->stencils_array_size = stencils_array_size;
  weights->stencils_size = stencils_size;
}

void yac_interp_weights_add_wsum(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * num_src_per_tgt, struct remote_point * srcs, double * w) {

  if (tgts->count == 0) return;

  // determine whether there are zero-weights
  int pack_flag = 0;
  for (size_t i = 0, k = 0; (i < tgts->count) && !pack_flag; ++i)
    for (size_t j = 0; (j < num_src_per_tgt[i]) && !pack_flag; ++j, ++k)
      pack_flag = (fabs(w[k]) <= WEIGHT_TOL);

  if (pack_flag) {

    for (size_t i = 0, k = 0, l = 0;
         i < tgts->count; i++) {

      size_t curr_count = num_src_per_tgt[i];

      for (size_t j = 0; j < curr_count; j++, k++) {

        if (fabs(w[k]) < WEIGHT_TOL) {
          num_src_per_tgt[i]--;
        } else {
          if (l != k) {
            srcs[l] = srcs[k];
            w[l] = w[k];
          }
          ++l;
        }
      }

      // if all weights were zero
      if ((curr_count != 0) && (num_src_per_tgt[i] == 0)) {

        if (l != k) {
          srcs[l] = srcs[k - curr_count];
          w[l] = 0.0;
        }
        num_src_per_tgt[i] = 1;
        ++l;
      }
    }
  }

  // check whether all weights are 1.0 and whether the number of source
  // points per target is one for all targets
  int flag_weight_one = 1;
  int flag_count_one = 1;
  for (size_t i = 0, j = 0;
       (i < tgts->count) && (flag_weight_one || flag_count_one); ++i) {

    size_t curr_count = num_src_per_tgt[i];
    flag_count_one &= curr_count == 1;

    for (size_t k = 0; (k < curr_count) && flag_weight_one; ++k, ++j)
      flag_weight_one &= fabs(w[j] - 1.0) < WEIGHT_TOL;
  }

  // if all weights are 1.0 -> use more optimised weight type
  if (flag_weight_one) {

    // if the number of source points for all target points is one
    if (flag_count_one)
      yac_interp_weights_add_direct(weights, tgts, srcs);
    else
      yac_interp_weights_add_sum(weights, tgts, num_src_per_tgt, srcs);

  } else {

    struct interp_weight_stencil * stencils = weights->stencils;
    size_t stencils_array_size = weights->stencils_array_size;
    size_t stencils_size = weights->stencils_size;

    ENSURE_ARRAY_SIZE(stencils, stencils_array_size, stencils_size + tgts->count);

    for (size_t i = 0; i < tgts->count; ++i, ++stencils_size) {

      size_t curr_num_src = num_src_per_tgt[i];

      // remove target for which no weights were provided
      if (curr_num_src == 0) {
        --stencils_size;
        continue;
      }

      double * curr_weights =
        xmalloc(curr_num_src * sizeof(*curr_weights));

      stencils[stencils_size].type = WEIGHT_SUM;
      stencils[stencils_size].tgt = copy_remote_point(tgts->data[i]);
      stencils[stencils_size].data.weight_sum.srcs =
        copy_remote_points(srcs, curr_num_src);
      stencils[stencils_size].data.weight_sum.weights = curr_weights;
      memcpy(curr_weights, w, curr_num_src * sizeof(*curr_weights));

      srcs += curr_num_src;
      w += curr_num_src;
    }

    weights->stencils = stencils;
    weights->stencils_array_size = stencils_array_size;
    weights->stencils_size = stencils_size;
  }
}

void yac_interp_weights_add_sum(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * num_src_per_tgt, struct remote_point * srcs) {

  if (tgts->count == 0) return;

  // check whether the number of source points per target is one
  // for all targets
  int flag_count_one = 1;
  for (size_t i = 0; i < tgts->count; ++i) {
    if (num_src_per_tgt[i] != 1) {
      flag_count_one = 0;
      break;
    }
  }

  if (flag_count_one) {

    yac_interp_weights_add_direct(weights, tgts, srcs);

  } else {

    struct interp_weight_stencil * stencils = weights->stencils;
    size_t stencils_array_size = weights->stencils_array_size;
    size_t stencils_size = weights->stencils_size;

    ENSURE_ARRAY_SIZE(stencils, stencils_array_size, stencils_size + tgts->count);

    for (size_t i = 0; i < tgts->count; ++i, ++stencils_size) {

      size_t curr_num_src = num_src_per_tgt[i];

      stencils[stencils_size].type = SUM;
      stencils[stencils_size].tgt = copy_remote_point(tgts->data[i]);
      stencils[stencils_size].data.weight_sum.srcs =
        copy_remote_points(srcs, curr_num_src);
      stencils[stencils_size].data.weight_sum.weights = NULL;

      srcs += curr_num_src;
    }

    weights->stencils = stencils;
    weights->stencils_array_size = stencils_array_size;
    weights->stencils_size = stencils_size;
  }
}

void yac_interp_weights_add_direct(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  struct remote_point * srcs) {

  if (tgts->count == 0) return;

  struct interp_weight_stencil * stencils = weights->stencils;
  size_t stencils_array_size = weights->stencils_array_size;
  size_t stencils_size = weights->stencils_size;

  ENSURE_ARRAY_SIZE(stencils, stencils_array_size, stencils_size + tgts->count);

  for (size_t i = 0; i < tgts->count; ++i, ++stencils_size) {

    stencils[stencils_size].type = DIRECT;
    stencils[stencils_size].tgt = copy_remote_point(tgts->data[i]);
    stencils[stencils_size].data.direct.src = copy_remote_point(srcs[i]);
  }

  weights->stencils = stencils;
  weights->stencils_array_size = stencils_array_size;
  weights->stencils_size = stencils_size;
}

void yac_interp_weights_add_direct_mf(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * src_field_indices, struct remote_point ** srcs_per_field,
  size_t num_src_fields) {

  if (tgts->count == 0) return;

  if (num_src_fields == 1) {
    yac_interp_weights_add_direct(weights, tgts, srcs_per_field[0]);
    return;
  }

  struct interp_weight_stencil * stencils = weights->stencils;
  size_t stencils_array_size = weights->stencils_array_size;
  size_t stencils_size = weights->stencils_size;

  ENSURE_ARRAY_SIZE(
    stencils, stencils_array_size, stencils_size + tgts->count);
  stencils += stencils_size;

  size_t srcs_offsets[num_src_fields];
  memset(srcs_offsets, 0, num_src_fields * sizeof(srcs_offsets[0]));

  for (size_t i = 0; i < tgts->count; ++i) {

    size_t src_field_idx = src_field_indices[i];
    stencils[i].type = DIRECT_MF;
    stencils[i].tgt = copy_remote_point(tgts->data[i]);
    stencils[i].data.direct_mf.src =
      copy_remote_point(
        srcs_per_field[src_field_idx][srcs_offsets[src_field_idx]++]);
    stencils[i].data.direct_mf.field_idx = src_field_idx;
  }

  weights->stencils = stencils;
  weights->stencils_array_size = stencils_array_size;
  weights->stencils_size += tgts->count;
}

void yac_interp_weights_add_sum_mf(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * num_src_per_field_per_tgt, struct remote_point ** srcs_per_field,
  size_t num_src_fields) {

  if (tgts->count == 0) return;

  if (num_src_fields == 1) {
    yac_interp_weights_add_sum(
      weights, tgts, num_src_per_field_per_tgt, srcs_per_field[0]);
    return;
  }

  // check whether the number of source points per target is one
  // for all targets
  int flag_count_one = 1;
  for (size_t i = 0, k = 0; i < tgts->count; ++i) {
    size_t count = 0;
    for (size_t j = 0; j < num_src_fields; ++j, ++k)
      count += num_src_per_field_per_tgt[k];
    if (count != 1) {
      flag_count_one = 0;
      break;
    }
  }

  if (flag_count_one) {

    size_t * src_field_indices =
      xmalloc(tgts->count * sizeof(*src_field_indices));

      for (size_t i = 0, k = 0; i < tgts->count; ++i)
        for (size_t j = 0; j < num_src_fields; ++j, ++k)
          if (num_src_per_field_per_tgt[k])
            src_field_indices[i] = j;

    yac_interp_weights_add_direct_mf(
      weights, tgts, src_field_indices, srcs_per_field, num_src_fields);

    free(src_field_indices);

  } else {
    struct remote_point * curr_srcs_per_field[num_src_fields];
    memcpy(curr_srcs_per_field, srcs_per_field,
           num_src_fields * sizeof(*srcs_per_field));

    struct interp_weight_stencil * stencils = weights->stencils;
    size_t stencils_array_size = weights->stencils_array_size;
    size_t stencils_size = weights->stencils_size;

    ENSURE_ARRAY_SIZE(
      stencils, stencils_array_size, stencils_size + tgts->count);

    for (size_t i = 0; i < tgts->count; ++i, ++stencils_size) {

      size_t * curr_num_src_per_src_field =
        num_src_per_field_per_tgt + i * num_src_fields;
      size_t curr_num_src = 0;
      for (size_t j = 0; j < num_src_fields; ++j)
        curr_num_src += curr_num_src_per_src_field[j];

      stencils[stencils_size].type = SUM_MF;
      stencils[stencils_size].tgt = copy_remote_point(tgts->data[i]);
      stencils[stencils_size].data.sum_mf.field_indices =
        xmalloc(
          curr_num_src *
          sizeof(*(stencils[stencils_size].data.sum_mf.field_indices)));
      for (size_t j = 0, l = 0; j < num_src_fields; ++j) {
        size_t curr_num_src = curr_num_src_per_src_field[j];
        for (size_t k = 0; k < curr_num_src; ++k, ++l) {
          stencils[stencils_size].data.sum_mf.field_indices[l] = j;
        }
      }
      stencils[stencils_size].data.sum_mf.srcs =
        copy_remote_points_mf(
          curr_srcs_per_field, curr_num_src_per_src_field, num_src_fields);

      for (size_t j = 0; j < num_src_fields; ++j)
        curr_srcs_per_field[j] += curr_num_src_per_src_field[j];
    }

    weights->stencils = stencils;
    weights->stencils_array_size = stencils_array_size;
    weights->stencils_size = stencils_size;
  }
}

void yac_interp_weights_add_wsum_mf(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * num_src_per_field_per_tgt, struct remote_point ** srcs_per_field,
  double * w, size_t num_src_fields) {

  if (tgts->count == 0) return;

  if (num_src_fields == 1) {
    yac_interp_weights_add_wsum(
      weights, tgts, num_src_per_field_per_tgt, srcs_per_field[0], w);
    return;
  }

  // check whether all weights are 1.0 and whether the number of source
  // points per target is one for all targets
  int flag_weight_one = 1;
  for (size_t i = 0, j = 0;
       (i < tgts->count) && flag_weight_one; ++i) {

    for (size_t src_field_idx = 0; src_field_idx < num_src_fields;
         ++src_field_idx) {

      size_t curr_count =
        num_src_per_field_per_tgt[i * num_src_fields + src_field_idx];

      for (size_t k = 0; (k < curr_count) && flag_weight_one; ++k, ++j)
        flag_weight_one &= fabs(w[j] - 1.0) < 1e-12;
    }
  }

  // if all weights are 1.0 -> use more optimised weight type
  if (flag_weight_one) {

    yac_interp_weights_add_sum_mf(
      weights, tgts, num_src_per_field_per_tgt, srcs_per_field, num_src_fields);

  } else {

    struct remote_point * curr_srcs_per_field[num_src_fields];
    memcpy(curr_srcs_per_field, srcs_per_field,
           num_src_fields * sizeof(*srcs_per_field));

    struct interp_weight_stencil * stencils = weights->stencils;
    size_t stencils_array_size = weights->stencils_array_size;
    size_t stencils_size = weights->stencils_size;

    ENSURE_ARRAY_SIZE(
      stencils, stencils_array_size, stencils_size + tgts->count);

    for (size_t i = 0; i < tgts->count; ++i, ++stencils_size) {

      size_t * curr_num_src_per_src_field =
        num_src_per_field_per_tgt + i * num_src_fields;
      size_t curr_num_weights = 0;
      for (size_t j = 0; j < num_src_fields; ++j)
        curr_num_weights += curr_num_src_per_src_field[j];
      double * curr_weights =
        xmalloc(curr_num_weights * sizeof(*curr_weights));
      size_t * field_indices =
        xmalloc(curr_num_weights * sizeof(*field_indices));

      stencils[stencils_size].type = WEIGHT_SUM_MF;
      stencils[stencils_size].tgt = copy_remote_point(tgts->data[i]);
      stencils[stencils_size].data.weight_sum_mf.field_indices = field_indices;
      for (size_t j = 0, l = 0; j < num_src_fields; ++j) {
        size_t curr_num_src = curr_num_src_per_src_field[j];
        for (size_t k = 0; k < curr_num_src; ++k, ++l) field_indices[l] = j;
      }
      stencils[stencils_size].data.weight_sum_mf.srcs =
        copy_remote_points_mf(curr_srcs_per_field, curr_num_src_per_src_field, num_src_fields);
      stencils[stencils_size].data.weight_sum_mf.weights = curr_weights;
      memcpy(curr_weights, w, curr_num_weights * sizeof(*curr_weights));

      for (size_t j = 0; j < num_src_fields; ++j)
        curr_srcs_per_field[j] += curr_num_src_per_src_field[j];
      w += curr_num_weights;
    }

    weights->stencils = stencils;
    weights->stencils_array_size = stencils_array_size;
    weights->stencils_size = stencils_size;
  }
}

static int compare_stencils_fixed(const void * a, const void * b) {

  int ret = (((struct interp_weight_stencil_fixed *)a)->value >
             ((struct interp_weight_stencil_fixed *)b)->value) -
            (((struct interp_weight_stencil_fixed *)a)->value <
             ((struct interp_weight_stencil_fixed *)b)->value);

  if (ret) return ret;

  return (((struct interp_weight_stencil_fixed *)a)->orig_pos >
          ((struct interp_weight_stencil_fixed *)b)->orig_pos) -
         (((struct interp_weight_stencil_fixed *)a)->orig_pos <
          ((struct interp_weight_stencil_fixed *)b)->orig_pos);
}

static MPI_Datatype get_fixed_stencil_mpi_datatype(MPI_Comm comm) {

  struct interp_weight_stencil_fixed dummy;
  MPI_Datatype fixed_stencil_dt;
  int array_of_blocklengths[] = {1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.value) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.orig_pos) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] = {MPI_DOUBLE, MPI_UINT64_T};
  yac_mpi_call(
    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           array_of_types, &fixed_stencil_dt), comm);
  return yac_create_resized(fixed_stencil_dt, sizeof(dummy), comm);
}

static void yac_interp_weights_redist_fixed(
  MPI_Comm comm, uint64_t count,
  struct interp_weight_stencil * fixed_stencils,
  struct yac_interpolation * interp) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * sendcounts, * recvcounts, * sdispls, *rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  for (size_t i = 0; i < count; ++i) {
    int curr_count = fixed_stencils[i].tgt.data.count;
    struct remote_point_info * curr_point_infos =
      (curr_count == 1)?
        (&(fixed_stencils[i].tgt.data.data.single)):
        (fixed_stencils[i].tgt.data.data.multi);
    for (int j = 0; j < curr_count; ++j)
      sendcounts[curr_point_infos[j].rank]++;
  }

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t send_buffer_size =
    sdispls[comm_size] + sendcounts[comm_size - 1];
  size_t recv_buffer_size =
    rdispls[comm_size - 1] + recvcounts[comm_size - 1];

  struct interp_weight_stencil_fixed * buffer =
    xmalloc((send_buffer_size + recv_buffer_size) * sizeof(*buffer));
  struct interp_weight_stencil_fixed * send_buffer = buffer + recv_buffer_size;
  struct interp_weight_stencil_fixed * recv_buffer = buffer;

  // pack fixed stencils
  for (size_t i = 0; i < count; ++i) {
    int curr_count = fixed_stencils[i].tgt.data.count;
    struct remote_point_info * curr_point_infos =
      (curr_count == 1)?
        (&(fixed_stencils[i].tgt.data.data.single)):
        (fixed_stencils[i].tgt.data.data.multi);
    double value = fixed_stencils[i].data.fixed.value;
    for (int j = 0; j < curr_count; ++j) {
      size_t pos = sdispls[curr_point_infos[j].rank + 1]++;
      send_buffer[pos].value = value;
      send_buffer[pos].orig_pos = curr_point_infos[j].orig_pos;
    }
  }

  // create MPI Datatype for exchanging fixed stencils
  MPI_Datatype stencil_fixed_dt = get_fixed_stencil_mpi_datatype(comm);

  // redistribute fixed stencils
  yac_alltoallv_p2p(
    send_buffer, sendcounts, sdispls, recv_buffer, recvcounts, rdispls,
    sizeof(*send_buffer), stencil_fixed_dt, comm);

  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  yac_mpi_call(MPI_Type_free(&stencil_fixed_dt), comm);

  if (recv_buffer_size == 0) {
    free(buffer);
    return;
  }

  // sort stencils first by fixed value and second by orig_pos
  qsort(recv_buffer, recv_buffer_size, sizeof(*recv_buffer),
        compare_stencils_fixed);

  size_t * tgt_pos = xmalloc(recv_buffer_size * sizeof(*tgt_pos));
  for (size_t i = 0; i < recv_buffer_size; ++i)
    tgt_pos[i] = (size_t)(recv_buffer[i].orig_pos);

  size_t offset = 0, i = 0;
  while (offset < recv_buffer_size) {
    double fixed_value = recv_buffer[i].value;
    while ((i < recv_buffer_size) && (fixed_value == recv_buffer[i].value)) ++i;
    size_t curr_count = i - offset;
    yac_interpolation_add_fixed(
      interp, fixed_value, curr_count, tgt_pos + offset);
    offset = i;
  }

  free(buffer);
  free(tgt_pos);
}

static inline struct remote_point_info select_src(
  struct remote_point_infos src) {

  if (src.count == 1) return src.data.single;

  int min_rank = INT_MAX;
  size_t min_rank_idx = SIZE_MAX;
  for (int i = 0; i < src.count; ++i) {
    if (src.data.multi[i].rank < min_rank) {
      min_rank = src.data.multi[i].rank;
      min_rank_idx = i;
    }
  }

  return src.data.multi[min_rank_idx];
}

static void xt_redist_msg_free(
  struct Xt_redist_msg * msgs, size_t count, MPI_Comm comm) {
  for (size_t i = 0; i < count; ++i) {
    MPI_Datatype * dt = &(msgs[i].datatype);
    if (*dt != MPI_DATATYPE_NULL) yac_mpi_call(MPI_Type_free(dt), comm);
  }
  free(msgs);
}

static Xt_redist generate_direct_redist(
  uint64_t * src_orig_poses, size_t * sendcounts,
  struct interp_weight_stencil_direct * tgt_stencils,
  size_t * recvcounts, MPI_Comm comm) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t nsend = 0, nrecv = 0;
  size_t max_buffer_size = 0;
  for (int i = 0; i < comm_size; ++i) {
    nsend += sendcounts[i] > 0;
    nrecv += recvcounts[i] > 0;
    if (max_buffer_size < sendcounts[i]) max_buffer_size = sendcounts[i];
    if (max_buffer_size < recvcounts[i]) max_buffer_size = recvcounts[i];
  }

  size_t total_num_msg = nsend + nrecv;

  struct Xt_redist_msg * msgs_buffer =
    xmalloc(total_num_msg * sizeof(*msgs_buffer));
  struct Xt_redist_msg * send_msgs = msgs_buffer;
  struct Xt_redist_msg * recv_msgs = msgs_buffer + nsend;

  int * pos_buffer = xmalloc((size_t)max_buffer_size * sizeof(*pos_buffer));

  nsend = 0;
  nrecv = 0;
  for (int i = 0; i < comm_size; ++i) {
    if (recvcounts[i] > 0) {
      for (size_t j = 0; j < recvcounts[i]; ++j)
        pos_buffer[j] = (int)tgt_stencils[j].orig_pos;
      tgt_stencils += recvcounts[i];
      recv_msgs[nrecv].rank = i;
      recv_msgs[nrecv].datatype =
        xt_mpi_generate_datatype(pos_buffer, recvcounts[i], MPI_DOUBLE, comm);
      nrecv++;
    }
    if (sendcounts[i] > 0) {
      for (size_t j = 0; j < sendcounts[i]; ++j)
        pos_buffer[j] = (int)src_orig_poses[j];
      src_orig_poses += sendcounts[i];
      send_msgs[nsend].rank = i;
      send_msgs[nsend].datatype =
        xt_mpi_generate_datatype(pos_buffer, sendcounts[i], MPI_DOUBLE, comm);
      nsend++;
    }
  }

  free(pos_buffer);

  Xt_redist redist;
  MPI_Comm split_comm;

  if (total_num_msg > 0) {

    yac_mpi_call(MPI_Comm_split(comm, 1, 0, &split_comm), comm);

    int * rank_buffer =
      xmalloc(2 * total_num_msg * sizeof(*rank_buffer));
    int * orig_ranks = rank_buffer;
    int * split_ranks = rank_buffer + total_num_msg;

    for (size_t i = 0; i < total_num_msg; ++i)
      orig_ranks[i] = msgs_buffer[i].rank;

    MPI_Group orig_group, split_group;
    yac_mpi_call(MPI_Comm_group(comm, &orig_group), comm);
    yac_mpi_call(MPI_Comm_group(split_comm, &split_group), comm);

    yac_mpi_call(
      MPI_Group_translate_ranks(orig_group, total_num_msg, orig_ranks,
                                split_group, split_ranks), split_comm);

    for (size_t i = 0; i < total_num_msg; ++i)
      msgs_buffer[i].rank = split_ranks[i];

    free(rank_buffer);

    yac_mpi_call(MPI_Group_free(&split_group), comm);
    yac_mpi_call(MPI_Group_free(&orig_group), comm);

    // generate redist
    redist =
      xt_redist_single_array_base_new(
        nsend, nrecv, send_msgs, recv_msgs, split_comm);

  } else {
    yac_mpi_call(MPI_Comm_split(comm, 0, 0, &split_comm), comm);
    redist = NULL;
  }

  yac_mpi_call(MPI_Comm_free(&split_comm), comm);
  xt_redist_msg_free(msgs_buffer, total_num_msg, comm);

  return redist;
}

static Xt_redist * generate_direct_mf_redists(
  uint64_t * src_orig_pos, size_t * sendcounts,
  struct interp_weight_stencil_direct_mf * tgt_stencils,
  size_t * recvcounts, size_t num_src_fields, MPI_Comm comm) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t nsends[num_src_fields], nrecvs[num_src_fields];
  size_t max_buffer_size = 0;
  memset(nsends, 0, num_src_fields * sizeof(nsends[0]));
  memset(nrecvs, 0, num_src_fields * sizeof(nrecvs[0]));
  for (int i = 0; i < comm_size; ++i) {
    for (size_t j = 0; j < num_src_fields; ++j) {
      size_t idx = (size_t)i * num_src_fields + j;
      if (sendcounts[idx] > 0) nsends[j]++;
      if (recvcounts[idx] > 0) nrecvs[j]++;
      if (max_buffer_size < sendcounts[idx]) max_buffer_size = sendcounts[idx];
      if (max_buffer_size < recvcounts[idx]) max_buffer_size = recvcounts[idx];
    }
  }

  size_t nsend = 0, nrecv = 0;
  size_t send_offsets[num_src_fields];
  size_t recv_offsets[num_src_fields];
  for (size_t i = 0; i < num_src_fields; ++i) {
    send_offsets[i] = nsend;
    recv_offsets[i] = nrecv;
    nsend += nsends[i];
    nrecv += nrecvs[i];
  }

  size_t total_num_msg = nsend + nrecv;

  struct Xt_redist_msg * msgs_buffer =
    xmalloc(total_num_msg * sizeof(*msgs_buffer));
  struct Xt_redist_msg * send_msgs = msgs_buffer;
  struct Xt_redist_msg * recv_msgs = msgs_buffer + nsend;

  int * pos_buffer = xmalloc(max_buffer_size * sizeof(*pos_buffer));

  for (int i = 0; i < comm_size; ++i) {
    for (size_t src_field_idx = 0; src_field_idx < num_src_fields;
         ++src_field_idx) {
      size_t idx = (size_t)i * num_src_fields + src_field_idx;
      if (recvcounts[idx] > 0) {
        for (size_t j = 0; j < recvcounts[idx]; ++j)
          pos_buffer[j] = (int)tgt_stencils[j].orig_pos;
        tgt_stencils += recvcounts[idx];
        recv_msgs[recv_offsets[src_field_idx]].rank = i;
        recv_msgs[recv_offsets[src_field_idx]].datatype =
          xt_mpi_generate_datatype(
            pos_buffer, recvcounts[idx], MPI_DOUBLE, comm);
        recv_offsets[src_field_idx]++;
      }
      if (sendcounts[idx] > 0) {
        for (size_t j = 0; j < sendcounts[idx]; ++j)
          pos_buffer[j] = (int)src_orig_pos[j];
        src_orig_pos += sendcounts[idx];
        send_msgs[send_offsets[src_field_idx]].rank = i;
        send_msgs[send_offsets[src_field_idx]].datatype =
          xt_mpi_generate_datatype(
            pos_buffer, sendcounts[idx], MPI_DOUBLE, comm);
        send_offsets[src_field_idx]++;
      }
    }
  }

  free(pos_buffer);

  Xt_redist * redists;
  MPI_Comm split_comm;

  if (total_num_msg > 0) {

    yac_mpi_call(MPI_Comm_split(comm, 1, 0, &split_comm), comm);

    int * rank_buffer =
      xmalloc(2 * total_num_msg * sizeof(*rank_buffer));
    int * orig_ranks = rank_buffer;
    int * split_ranks = rank_buffer + total_num_msg;

    for (size_t i = 0; i < total_num_msg; ++i)
      orig_ranks[i] = msgs_buffer[i].rank;

    MPI_Group orig_group, split_group;
    yac_mpi_call(MPI_Comm_group(comm, &orig_group), comm);
    yac_mpi_call(MPI_Comm_group(split_comm, &split_group), comm);

    yac_mpi_call(
      MPI_Group_translate_ranks(orig_group, total_num_msg, orig_ranks,
                                split_group, split_ranks), split_comm);

    for (size_t i = 0; i < total_num_msg; ++i)
      msgs_buffer[i].rank = split_ranks[i];

    free(rank_buffer);

    yac_mpi_call(MPI_Group_free(&split_group), comm);
    yac_mpi_call(MPI_Group_free(&orig_group), comm);

    // generate redists
    redists = xmalloc(num_src_fields * sizeof(*redists));
    for (size_t src_field_idx = 0; src_field_idx < num_src_fields;
         ++src_field_idx) {
      redists[src_field_idx] =
        xt_redist_single_array_base_new(
          nsends[src_field_idx], nrecvs[src_field_idx],
          send_msgs, recv_msgs, split_comm);
      send_msgs += nsends[src_field_idx];
      recv_msgs += nrecvs[src_field_idx];
    }

  } else {
    yac_mpi_call(MPI_Comm_split(comm, 0, 0, &split_comm), comm);
    redists = NULL;
  }

  yac_mpi_call(MPI_Comm_free(&split_comm), comm);
  xt_redist_msg_free(msgs_buffer, total_num_msg, comm);

  return redists;
}

static MPI_Datatype get_direct_stencil_mpi_datatype(MPI_Comm comm) {

  struct interp_weight_stencil_direct dummy;
  MPI_Datatype direct_stencil_dt;
  int array_of_blocklengths[] = {1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.src) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.orig_pos) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  MPI_Datatype array_of_types[] =
    {yac_get_remote_point_info_mpi_datatype(comm), MPI_UINT64_T};
  yac_mpi_call(
    MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements,
                           array_of_types, &direct_stencil_dt), comm);
  yac_mpi_call(MPI_Type_free(&(array_of_types[0])), comm);
  return yac_create_resized(direct_stencil_dt, sizeof(dummy), comm);
}

static int compare_stencils_direct(const void * a, const void * b) {

  int ret = ((struct interp_weight_stencil_direct *)a)->src.rank -
            ((struct interp_weight_stencil_direct *)b)->src.rank;

  if (ret) return ret;

  return (((struct interp_weight_stencil_direct *)a)->orig_pos >
          ((struct interp_weight_stencil_direct *)b)->orig_pos) -
         (((struct interp_weight_stencil_direct *)a)->orig_pos <
          ((struct interp_weight_stencil_direct *)b)->orig_pos);
}

static void yac_interp_weights_redist_direct(
  MPI_Comm comm, uint64_t count,
  struct interp_weight_stencil * direct_stencils,
  struct yac_interpolation * interp) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * sendcounts, * recvcounts, * sdispls, *rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  for (size_t i = 0; i < count; ++i) {
    int curr_count = direct_stencils[i].tgt.data.count;
    struct remote_point_info * curr_point_info =
      (curr_count == 1)?
        (&(direct_stencils[i].tgt.data.data.single)):
        (direct_stencils[i].tgt.data.data.multi);
    for (int j = 0; j < curr_count; ++j)
      sendcounts[curr_point_info[j].rank]++;
  }

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t send_buffer_size =
    sdispls[comm_size] + sendcounts[comm_size - 1];
  size_t recv_buffer_size =
    rdispls[comm_size - 1] + recvcounts[comm_size - 1];
  size_t tgt_count = recv_buffer_size;

  struct interp_weight_stencil_direct * stencil_buffer =
    xmalloc((send_buffer_size + recv_buffer_size) * sizeof(*stencil_buffer));
  struct interp_weight_stencil_direct * send_stencil_buffer =
    stencil_buffer + recv_buffer_size;
  struct interp_weight_stencil_direct * recv_stencil_buffer = stencil_buffer;

  // pack direct stencils
  for (size_t i = 0; i < count; ++i) {
    int curr_count = direct_stencils[i].tgt.data.count;
    struct remote_point_info * curr_point_infos =
      (curr_count == 1)?
        (&(direct_stencils[i].tgt.data.data.single)):
        (direct_stencils[i].tgt.data.data.multi);
    struct remote_point_info src =
      select_src(direct_stencils[i].data.direct.src.data);
    for (int j = 0; j < curr_count; ++j) {
      size_t pos = sdispls[curr_point_infos[j].rank + 1]++;
      send_stencil_buffer[pos].src = src;
      send_stencil_buffer[pos].orig_pos = curr_point_infos[j].orig_pos;
    }
  }

  // create MPI Datatype for exchanging direct stencils
  MPI_Datatype stencil_direct_dt = get_direct_stencil_mpi_datatype(comm);

  // redistribute stencils based on target owners
  yac_alltoallv_p2p(
    send_stencil_buffer, sendcounts, sdispls,
    recv_stencil_buffer, recvcounts, rdispls,
    sizeof(*stencil_buffer), stencil_direct_dt, comm);

  yac_mpi_call(MPI_Type_free(&stencil_direct_dt), comm);

  // sort stencils based on src rank first and by target orig pos second
  qsort(recv_stencil_buffer, tgt_count, sizeof(*recv_stencil_buffer),
        compare_stencils_direct);

  memset(sendcounts, 0, (size_t)comm_size * sizeof(*sendcounts));

  for (size_t i = 0; i < tgt_count; ++i)
    sendcounts[recv_stencil_buffer[i].src.rank]++;

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  send_buffer_size = sdispls[comm_size] + sendcounts[comm_size - 1];
  recv_buffer_size = rdispls[comm_size - 1] + recvcounts[comm_size - 1];

  uint64_t * orig_pos_buffer =
    xmalloc((send_buffer_size + recv_buffer_size) * sizeof(*orig_pos_buffer));
  uint64_t * send_orig_pos_buffer = orig_pos_buffer + recv_buffer_size;
  uint64_t * recv_orig_pos_buffer = orig_pos_buffer;

  for (size_t i = 0; i < tgt_count; ++i)
    send_orig_pos_buffer[sdispls[recv_stencil_buffer[i].src.rank + 1]++] =
      recv_stencil_buffer[i].src.orig_pos;

  // inform source processes about their requested points
  yac_alltoallv_p2p(
    send_orig_pos_buffer, sendcounts, sdispls,
    recv_orig_pos_buffer, recvcounts, rdispls,
    sizeof(*send_orig_pos_buffer), MPI_UINT64_T, comm);

  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  // generate redist
  Xt_redist redist =
    generate_direct_redist(
      recv_orig_pos_buffer, recvcounts, recv_stencil_buffer, sendcounts, comm);
  free(orig_pos_buffer);
  free(stencil_buffer);

  yac_interpolation_add_direct(interp, redist);

  if (redist != NULL) xt_redist_delete(redist);
}

static MPI_Datatype get_direct_mf_stencil_mpi_datatype(MPI_Comm comm) {

  struct interp_weight_stencil_direct_mf dummy;
  MPI_Datatype direct_stencil_mf_dt;
  int array_of_blocklengths[] = {1, 1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.src) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.src_field_idx) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.orig_pos) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  MPI_Datatype array_of_types[] =
    {yac_get_remote_point_info_mpi_datatype(comm), MPI_UINT64_T, MPI_UINT64_T};
  yac_mpi_call(
    MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements,
                           array_of_types, &direct_stencil_mf_dt), comm);
  yac_mpi_call(MPI_Type_free(&(array_of_types[0])), comm);
  return yac_create_resized(direct_stencil_mf_dt, sizeof(dummy), comm);
}

static int compare_stencils_direct_mf(const void * a, const void * b) {

  int ret = ((struct interp_weight_stencil_direct_mf *)a)->src.rank -
            ((struct interp_weight_stencil_direct_mf *)b)->src.rank;

  if (ret) return ret;

  ret = (((struct interp_weight_stencil_direct_mf *)a)->src_field_idx >
         ((struct interp_weight_stencil_direct_mf *)b)->src_field_idx) -
        (((struct interp_weight_stencil_direct_mf *)a)->src_field_idx <
         ((struct interp_weight_stencil_direct_mf *)b)->src_field_idx);

  if (ret) return ret;

  return (((struct interp_weight_stencil_direct_mf *)a)->orig_pos >
          ((struct interp_weight_stencil_direct_mf *)b)->orig_pos) -
         (((struct interp_weight_stencil_direct_mf *)a)->orig_pos <
          ((struct interp_weight_stencil_direct_mf *)b)->orig_pos);
}

static void yac_interp_weights_redist_direct_mf(
  MPI_Comm comm, uint64_t count,
  struct interp_weight_stencil * direct_mf_stencils,
  struct yac_interpolation * interp) {

  // determine the number of source fields
  uint64_t num_src_fields = 0;
  for (size_t i = 0; i < count; ++i) {
    uint64_t src_field_idx =
      (uint64_t)(direct_mf_stencils[i].data.direct_mf.field_idx);
    if (src_field_idx >= num_src_fields) num_src_fields = src_field_idx + 1;
  }
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &num_src_fields, 1, MPI_UINT64_T, MPI_MAX, comm), comm);

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * sendcounts, * recvcounts, * sdispls, *rdispls;
  yac_get_comm_buffers(
    (size_t)num_src_fields, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);
  size_t * size_t_buffer =
    xmalloc(4 * (size_t)comm_size * sizeof(*size_t_buffer));
  size_t * total_sendcounts = size_t_buffer + 0 * comm_size;
  size_t * total_recvcounts = size_t_buffer + 1 * comm_size;
  size_t * total_sdispls =    size_t_buffer + 2 * comm_size;
  size_t * total_rdispls =    size_t_buffer + 3 * comm_size;

  for (size_t i = 0; i < count; ++i) {
    int curr_count = direct_mf_stencils[i].tgt.data.count;
    struct remote_point_info * curr_point_info =
      (curr_count == 1)?
        (&(direct_mf_stencils[i].tgt.data.data.single)):
        (direct_mf_stencils[i].tgt.data.data.multi);
    uint64_t src_field_idx =
      (uint64_t)(direct_mf_stencils[i].data.direct_mf.field_idx);
    for (int j = 0; j < curr_count; ++j)
      sendcounts[
        (uint64_t)(curr_point_info[j].rank) * num_src_fields + src_field_idx]++;
  }

  yac_generate_alltoallv_args(
    (size_t)num_src_fields, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t saccu = 0, raccu = 0;
  for (int i = 0; i < comm_size; ++i) {
    total_sdispls[i] = saccu;
    total_rdispls[i] = raccu;
    total_sendcounts[i] = 0;
    total_recvcounts[i] = 0;
    for (size_t j = 0; j < num_src_fields; ++j) {
      total_sendcounts[i] += sendcounts[num_src_fields * i + j];
      total_recvcounts[i] += recvcounts[num_src_fields * i + j];
    }
    saccu += total_sendcounts[i];
    raccu += total_recvcounts[i];
  }

  size_t send_buffer_size = total_sdispls[comm_size - 1] +
                            total_sendcounts[comm_size - 1];
  size_t recv_buffer_size = total_rdispls[comm_size - 1] +
                            total_recvcounts[comm_size - 1];
  size_t tgt_count = recv_buffer_size;

  struct interp_weight_stencil_direct_mf * stencil_buffer =
    xmalloc((send_buffer_size + recv_buffer_size) * sizeof(*stencil_buffer));
  struct interp_weight_stencil_direct_mf * send_stencil_buffer =
    stencil_buffer + recv_buffer_size;
  struct interp_weight_stencil_direct_mf * recv_stencil_buffer = stencil_buffer;

  // pack direct_mf stencils
  for (size_t i = 0; i < count; ++i) {
    int curr_count = direct_mf_stencils[i].tgt.data.count;
    struct remote_point_info * curr_point_infos =
      (curr_count == 1)?
        (&(direct_mf_stencils[i].tgt.data.data.single)):
        (direct_mf_stencils[i].tgt.data.data.multi);
    struct remote_point_info src =
      select_src(direct_mf_stencils[i].data.direct_mf.src.data);
    uint64_t src_field_idx =
      (uint64_t)(direct_mf_stencils[i].data.direct_mf.field_idx);
    for (int j = 0; j < curr_count; ++j) {
      size_t pos =
        sdispls[(uint64_t)(curr_point_infos[j].rank) * num_src_fields +
                src_field_idx + 1]++;
      send_stencil_buffer[pos].src = src;
      send_stencil_buffer[pos].src_field_idx = src_field_idx;
      send_stencil_buffer[pos].orig_pos = curr_point_infos[j].orig_pos;
    }
  }

  // create MPI Datatype for exchanging direct_mf stencils
  MPI_Datatype stencil_direct_mf_dt = get_direct_mf_stencil_mpi_datatype(comm);

  // redistribute stencils based on target owners
  yac_alltoallv_p2p(
    send_stencil_buffer, total_sendcounts, total_sdispls,
    recv_stencil_buffer, total_recvcounts, total_rdispls,
    sizeof(*stencil_buffer), stencil_direct_mf_dt, comm);

  yac_mpi_call(MPI_Type_free(&stencil_direct_mf_dt), comm);

  // sort stencils based on src rank first, src_field_idx, and
  // by target orig pos second
  qsort(recv_stencil_buffer, tgt_count, sizeof(*recv_stencil_buffer),
        compare_stencils_direct_mf);

  memset(sendcounts, 0,
         (size_t)comm_size * (size_t)num_src_fields * sizeof(*sendcounts));

  for (size_t i = 0; i < tgt_count; ++i)
    sendcounts[(uint64_t)(recv_stencil_buffer[i].src.rank) * num_src_fields +
               recv_stencil_buffer[i].src_field_idx]++;

  yac_generate_alltoallv_args(
    (size_t)num_src_fields, sendcounts, recvcounts, sdispls, rdispls, comm);

  saccu = 0, raccu = 0;
  for (int i = 0; i < comm_size; ++i) {
    total_sdispls[i] = saccu;
    total_rdispls[i] = raccu;
    total_sendcounts[i] = 0;
    total_recvcounts[i] = 0;
    for (size_t j = 0; j < num_src_fields; ++j) {
      total_sendcounts[i] += sendcounts[num_src_fields * i + j];
      total_recvcounts[i] += recvcounts[num_src_fields * i + j];
    }
    saccu += total_sendcounts[i];
    raccu += total_recvcounts[i];
  }

  send_buffer_size = total_sdispls[comm_size - 1] +
                     total_sendcounts[comm_size - 1];
  recv_buffer_size = total_rdispls[comm_size - 1] +
                     total_recvcounts[comm_size - 1];

  uint64_t * orig_pos_buffer =
    xmalloc((send_buffer_size + recv_buffer_size) * sizeof(*orig_pos_buffer));
  uint64_t * send_orig_pos_buffer = orig_pos_buffer + recv_buffer_size;
  uint64_t * recv_orig_pos_buffer = orig_pos_buffer;

  for (size_t i = 0; i < tgt_count; ++i)
    send_orig_pos_buffer[
      sdispls[(uint64_t)(recv_stencil_buffer[i].src.rank) * num_src_fields +
              recv_stencil_buffer[i].src_field_idx + 1]++] =
        recv_stencil_buffer[i].src.orig_pos;

  // inform source processes about their requested points
  yac_alltoallv_p2p(
    send_orig_pos_buffer, total_sendcounts, total_sdispls,
    recv_orig_pos_buffer, total_recvcounts, total_rdispls,
    sizeof(*send_orig_pos_buffer), MPI_UINT64_T, comm);
  free(size_t_buffer);

  // generate redist
  Xt_redist * redists =
    generate_direct_mf_redists(
      recv_orig_pos_buffer, recvcounts, recv_stencil_buffer, sendcounts,
      (size_t)num_src_fields, comm);

  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  free(orig_pos_buffer);
  free(stencil_buffer);

  yac_interpolation_add_direct_mf(interp, redists, (size_t)num_src_fields);

  if (redists != NULL) {
    for (size_t i = 0; i < (size_t)num_src_fields; ++i)
      xt_redist_delete(redists[i]);
    free(redists);
  }
}

static int get_stencil_pack_size_fixed(
  struct interp_weight_stencil * stencil, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  UNUSED(stencil);
  UNUSED(point_info_dt);

  int pack_size_value;

  yac_mpi_call(MPI_Pack_size(1, MPI_DOUBLE, comm, &pack_size_value), comm);

  return pack_size_value;
}

static int get_stencil_pack_size_direct(
  struct interp_weight_stencil * stencil, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  return
    yac_remote_point_get_pack_size(
      &(stencil->data.direct.src), point_info_dt, comm);
}

static int get_stencil_pack_size_sum(
  struct interp_weight_stencil * stencil, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  return
    yac_remote_points_get_pack_size(
      stencil->data.sum.srcs, point_info_dt, comm);
}

static int get_stencil_pack_size_wsum(
  struct interp_weight_stencil * stencil, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  int pack_size_weights;
  yac_mpi_call(
    MPI_Pack_size(
      (int)(stencil->data.weight_sum.srcs->count),
      MPI_DOUBLE, comm, &pack_size_weights), comm);

  return
    yac_remote_points_get_pack_size(
      stencil->data.weight_sum.srcs, point_info_dt, comm) +
    pack_size_weights;
}

static int get_stencil_pack_size_direct_mf(
  struct interp_weight_stencil * stencil, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  int pack_size_src_field_idx;
  yac_mpi_call(
    MPI_Pack_size(
      1, MPI_UINT64_T, comm, &pack_size_src_field_idx), comm);

  return
    yac_remote_point_get_pack_size(
      &(stencil->data.direct_mf.src), point_info_dt, comm) +
    pack_size_src_field_idx;
}

static int get_stencil_pack_size_wsum_mf(
  struct interp_weight_stencil * stencil, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  int pack_size_weights, pack_size_field_indices;
  int count = (int)(stencil->data.weight_sum_mf.srcs->count);
  yac_mpi_call(
    MPI_Pack_size(
      count, MPI_DOUBLE, comm, &pack_size_weights), comm);
  yac_mpi_call(
    MPI_Pack_size(
      count, MPI_UINT64_T, comm, &pack_size_field_indices), comm);

  return
    yac_remote_points_get_pack_size(
      stencil->data.weight_sum_mf.srcs, point_info_dt, comm) +
    pack_size_weights + pack_size_field_indices;
}

static int get_stencil_pack_size_sum_mf(
  struct interp_weight_stencil * stencil, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  int pack_size_field_indices;
  yac_mpi_call(
    MPI_Pack_size(
      (int)(stencil->data.sum_mf.srcs->count),
      MPI_UINT64_T, comm, &pack_size_field_indices), comm);

  return
    yac_remote_points_get_pack_size(
      stencil->data.sum_mf.srcs, point_info_dt, comm) +
    pack_size_field_indices;
}

static struct interp_weight_stencil copy_interp_weight_stencil(
  struct interp_weight_stencil * stencil, struct remote_point point) {

  struct interp_weight_stencil stencil_copy = *stencil;
  stencil_copy.tgt = copy_remote_point(point);

  YAC_ASSERT(
    (stencil->type == FIXED) ||
    (stencil->type == DIRECT) ||
    (stencil->type == SUM) ||
    (stencil->type == WEIGHT_SUM) ||
    (stencil->type == DIRECT_MF) ||
    (stencil->type == SUM_MF) ||
    (stencil->type == WEIGHT_SUM_MF),
    "ERROR(copy_interp_weight_stencil): invalid stencil type")

  switch (stencil->type) {
    case(FIXED):
      // nothing to be done
      break;
    case(DIRECT):
      stencil_copy.data.direct.src =
        copy_remote_point(stencil->data.direct.src);
      break;
    case(SUM):
      stencil_copy.data.weight_sum.weights = NULL;
      stencil_copy.data.sum.srcs =
        copy_remote_points(
          stencil->data.sum.srcs->data, stencil->data.sum.srcs->count);
      break;
    case(WEIGHT_SUM): {
      stencil_copy.data.weight_sum.srcs =
        copy_remote_points(
          stencil->data.weight_sum.srcs->data,
          stencil->data.weight_sum.srcs->count);
      size_t weight_size =
        stencil->data.weight_sum.srcs->count *
        sizeof(*(stencil_copy.data.weight_sum.weights));
      stencil_copy.data.weight_sum.weights = xmalloc(weight_size);
      memcpy(stencil_copy.data.weight_sum.weights,
             stencil->data.weight_sum.weights, weight_size);
      break;
    }
    case(DIRECT_MF):
      stencil_copy.data.direct_mf.src =
        copy_remote_point(stencil->data.direct_mf.src);
      stencil_copy.data.direct_mf.field_idx =
        stencil->data.direct_mf.field_idx;
      break;
    case(SUM_MF): {
      stencil_copy.data.sum_mf.srcs =
        copy_remote_points(
          stencil->data.sum_mf.srcs->data,
          stencil->data.sum_mf.srcs->count);
      size_t field_indices_size =
        stencil->data.sum_mf.srcs->count *
        sizeof(*(stencil_copy.data.sum_mf.field_indices));
      stencil_copy.data.sum_mf.field_indices = xmalloc(field_indices_size);
      memcpy(stencil_copy.data.sum_mf.field_indices,
             stencil->data.sum_mf.field_indices, field_indices_size);
      break;
    }
    default:
    case(WEIGHT_SUM_MF): {
      stencil_copy.data.weight_sum_mf.srcs =
        copy_remote_points(
          stencil->data.weight_sum_mf.srcs->data,
          stencil->data.weight_sum_mf.srcs->count);
      size_t weight_size =
        stencil->data.weight_sum_mf.srcs->count *
        sizeof(*(stencil_copy.data.weight_sum_mf.weights));
      stencil_copy.data.weight_sum_mf.weights = xmalloc(weight_size);
      memcpy(stencil_copy.data.weight_sum_mf.weights,
             stencil->data.weight_sum_mf.weights, weight_size);
      size_t field_indices_size =
        stencil->data.weight_sum_mf.srcs->count *
        sizeof(*(stencil_copy.data.weight_sum_mf.field_indices));
      stencil_copy.data.weight_sum_mf.field_indices =
        xmalloc(field_indices_size);
      memcpy(stencil_copy.data.weight_sum_mf.field_indices,
             stencil->data.weight_sum_mf.field_indices, field_indices_size);
      break;
    }
  };
  return stencil_copy;
}

static struct interp_weight_stencil wcopy_interp_weight_stencil(
  struct interp_weight_stencil * stencil, struct remote_point point,
  double weight) {

  if (weight == 1.0) return copy_interp_weight_stencil(stencil, point);

  struct remote_point * srcs;
  size_t src_count;
  double * weights;

  YAC_ASSERT(
    (stencil->type == FIXED) ||
    (stencil->type == DIRECT) ||
    (stencil->type == SUM) ||
    (stencil->type == WEIGHT_SUM),
    "ERROR(wcopy_interp_weight_stencil): invalid stencil type")

  switch (stencil->type) {
    case (FIXED):
      return
        (struct interp_weight_stencil) {
          .type = FIXED,
          .data.fixed.value = stencil->data.fixed.value * weight,
          .tgt = copy_remote_point(point)};
    case (DIRECT):
      src_count = 1;
      srcs = &(stencil->data.direct.src);
      weights = NULL;
      break;
    case (SUM):
      src_count = stencil->data.sum.srcs->count;
      srcs = stencil->data.sum.srcs->data;
      weights = NULL;
      break;
    default:
    case (WEIGHT_SUM):
      src_count = stencil->data.weight_sum.srcs->count;
      srcs = stencil->data.weight_sum.srcs->data;
      weights = stencil->data.weight_sum.weights;
      break;
  };

  double * new_weights = xmalloc(src_count * sizeof(*new_weights));
  if (weights == NULL)
    for (size_t i = 0; i < src_count; ++i) new_weights[i] = weight;
  else
    for (size_t i = 0; i < src_count; ++i) new_weights[i] = weights[i] * weight;

  struct interp_weight_stencil stencil_wcopy;
  stencil_wcopy.type = WEIGHT_SUM;
  stencil_wcopy.data.weight_sum.srcs = copy_remote_points(srcs, src_count);
  stencil_wcopy.data.weight_sum.weights = new_weights;
  stencil_wcopy.tgt = copy_remote_point(point);

  return stencil_wcopy;
}

static int compare_w_global_id(const void * a, const void * b) {

  int ret = (((struct weighted_global_id *)a)->global_id >
             ((struct weighted_global_id *)b)->global_id) -
            (((struct weighted_global_id *)a)->global_id <
             ((struct weighted_global_id *)b)->global_id);

  if (ret) return ret;

  return (((struct weighted_global_id *)a)->weight >
          ((struct weighted_global_id *)b)->weight) -
         (((struct weighted_global_id *)a)->weight <
          ((struct weighted_global_id *)b)->weight);
}

static int compare_remote_point(const void * a, const void * b) {

  return ((const struct remote_point*)a)->global_id -
         ((const struct remote_point*)b)->global_id;
}

static void compact_srcs_w(
  struct remote_points * srcs, double ** w) {

  struct remote_point * data = srcs->data;
  size_t count = srcs->count;

  struct weighted_global_id * w_global_id =
    xmalloc(count * sizeof(*w_global_id));

  // extract global ids and weights
  for (size_t i = 0; i < count; ++i) {
    w_global_id[i].global_id = data[i].global_id;
    w_global_id[i].weight = (*w)[i];
  }

  // sort by global ids and weights
  qsort(w_global_id, count, sizeof(*w_global_id), compare_w_global_id);

  // sort sources by global ids
  qsort(data, count, sizeof(*data), compare_remote_point);

  size_t new_count = 0;

  // compact sources
  for (size_t i = 0; i < count;) {

    data[new_count] = data[i];

    yac_int curr_global_id = w_global_id[i].global_id;
    double curr_weight = w_global_id[i].weight;

    ++i;

    while((i < count) && (curr_global_id == w_global_id[i].global_id)) {

      curr_weight += w_global_id[i].weight;
      ++i;
    }

    (*w)[new_count] = curr_weight;
    ++new_count;
  }

  free(w_global_id);

  srcs->data = xrealloc(data, new_count * sizeof(*data));
  srcs->count = new_count;
  *w = xrealloc(*w, new_count * sizeof(**w));
}

static struct interp_weight_stencil stencils_merge_wsum(
  struct interp_weight_stencil ** stencils, double * w, size_t num_stencils) {

  size_t src_count = 0;
  size_t point_info_buffer_size = 0;

  for (size_t i = 0; i < num_stencils; ++i) {
    size_t curr_src_count;
    struct remote_point * srcs;
    YAC_ASSERT(
      (stencils[i]->type == DIRECT) ||
      (stencils[i]->type == SUM) ||
      (stencils[i]->type == WEIGHT_SUM),
      "ERROR(stencils_merge_wsum): invalid stencil type")
    switch (stencils[i]->type) {
      case (DIRECT):
        curr_src_count = 1;
        srcs = &(stencils[i]->data.direct.src);
        break;
      case (SUM):
        curr_src_count = stencils[i]->data.sum.srcs->count;
        srcs = stencils[i]->data.sum.srcs->data;
        break;
      default:
      case (WEIGHT_SUM):
        curr_src_count = stencils[i]->data.weight_sum.srcs->count;
        srcs = stencils[i]->data.weight_sum.srcs->data;
        break;
    };
    src_count += curr_src_count;
    for (size_t j = 0, curr_src_data_count; j < curr_src_count; ++j)
      if (((curr_src_data_count = srcs[j].data.count)) > 1)
        point_info_buffer_size += curr_src_data_count;
  }

  struct remote_points * srcs =
    xmalloc(point_info_buffer_size * sizeof(struct remote_point_info) +
            sizeof(*srcs));
  srcs->data = xmalloc(src_count * sizeof(*(srcs->data)));
  srcs->count = src_count;
  struct remote_point_info * point_info_buffer = &(srcs->buffer[0]);
  double * new_w = xmalloc(src_count * sizeof(*new_w));

  for (size_t i = 0, offset = 0; i < num_stencils; ++i) {
    size_t curr_src_count;
    struct remote_point * curr_srcs;
    double * stencil_w;
    YAC_ASSERT(
      (stencils[i]->type == DIRECT) ||
      (stencils[i]->type == SUM) ||
      (stencils[i]->type == WEIGHT_SUM),
      "ERROR(stencils_merge_wsum): invalid stencil type")
    switch (stencils[i]->type) {
      case (DIRECT):
        curr_src_count = 1;
        curr_srcs = &(stencils[i]->data.direct.src);
        stencil_w = NULL;
        break;
      case (SUM):
        curr_src_count = stencils[i]->data.sum.srcs->count;
        curr_srcs = stencils[i]->data.sum.srcs->data;
        stencil_w = NULL;
        break;
      default:
      case (WEIGHT_SUM):
        curr_src_count = stencils[i]->data.weight_sum.srcs->count;
        curr_srcs = stencils[i]->data.weight_sum.srcs->data;
        stencil_w = stencils[i]->data.weight_sum.weights;
        break;
    };
    copy_remote_points_no_alloc(
      srcs->data + offset, curr_srcs, curr_src_count, &point_info_buffer);
    if (stencil_w == NULL)
      for (size_t j = 0; j < curr_src_count; ++j, ++offset)
        new_w[offset] = w[i];
    else
      for (size_t j = 0; j < curr_src_count; ++j, ++offset)
        new_w[offset] = w[i] * stencil_w[j];
  }

  compact_srcs_w(srcs, &new_w);

  struct interp_weight_stencil merge_stencil;
  merge_stencil.type = WEIGHT_SUM;
  merge_stencil.data.weight_sum.srcs = srcs;
  merge_stencil.data.weight_sum.weights = new_w;

  return merge_stencil;
}

static struct interp_weight_stencil stencils_merge_sum(
  struct interp_weight_stencil ** stencils, double * w, size_t num_stencils) {

  for (size_t i = 0; i < num_stencils; ++i)
    if (w[i] != 1.0)
      return stencils_merge_wsum(stencils, w, num_stencils);

  size_t src_count = 0;
  size_t point_info_buffer_size = 0;

  for (size_t i = 0; i < num_stencils; ++i) {
    size_t curr_src_count;
    struct remote_point * srcs;
    YAC_ASSERT(
      (stencils[i]->type == DIRECT) ||
      (stencils[i]->type == SUM),
      "ERROR(stencils_merge_sum): invalid stencil type")
    switch (stencils[i]->type) {
      case (DIRECT):
        curr_src_count = 1;
        srcs = &(stencils[i]->data.direct.src);
        break;
      default:
      case (SUM):
        curr_src_count = stencils[i]->data.sum.srcs->count;
        srcs = stencils[i]->data.sum.srcs->data;
        break;
    };
    src_count += curr_src_count;
    for (size_t j = 0, curr_src_data_count; j < curr_src_count; ++j)
      if (((curr_src_data_count = srcs[j].data.count)) > 1)
        point_info_buffer_size += curr_src_data_count;
  }

  struct remote_points * srcs =
    xmalloc(point_info_buffer_size * sizeof(struct remote_point_info) +
            sizeof(*srcs));
  srcs->data = xmalloc(src_count * sizeof(*(srcs->data)));
  srcs->count = src_count;
  struct remote_point_info * point_info_buffer = &(srcs->buffer[0]);

  for (size_t i = 0, offset = 0; i < num_stencils; ++i) {
    size_t curr_src_count;
    struct remote_point * curr_srcs;
    YAC_ASSERT(
      (stencils[i]->type == DIRECT) ||
      (stencils[i]->type == SUM),
      "ERROR(stencils_merge_sum): invalid stencil type")
    switch (stencils[i]->type) {
      case (DIRECT):
        curr_src_count = 1;
        curr_srcs = &(stencils[i]->data.direct.src);
        break;
      default:
      case (SUM):
        curr_src_count = stencils[i]->data.sum.srcs->count;
        curr_srcs = stencils[i]->data.sum.srcs->data;
        break;
    };
    copy_remote_points_no_alloc(
      srcs->data + offset, curr_srcs, curr_src_count, &point_info_buffer);
    offset += curr_src_count;
  }

  qsort(srcs->data, srcs->count, sizeof(*(srcs->data)), compare_remote_point);

  struct interp_weight_stencil merge_stencil;
  merge_stencil.type = SUM;
  merge_stencil.data.sum.srcs = srcs;

  return merge_stencil;
}

static struct interp_weight_stencil stencils_merge(
  struct interp_weight_stencil ** stencils, double * w, size_t num_stencils,
  struct remote_point point) {

  if (num_stencils == 1)
    return wcopy_interp_weight_stencil(*stencils, point, *w);

  int fixed_count = 0;
  int direct_count = 0;
  int sum_count = 0;
  int wsum_count = 0;
  double fixed_value = 0.0;

  for (size_t i = 0; i < num_stencils; ++i) {
    YAC_ASSERT(
      (stencils[i]->type != DIRECT_MF) &&
      (stencils[i]->type != SUM_MF) &&
      (stencils[i]->type != WEIGHT_SUM_MF),
      "ERROR(stencils_merge): multiple source fields not yet supported")
    YAC_ASSERT(
      (stencils[i]->type == FIXED) ||
      (stencils[i]->type == DIRECT) ||
      (stencils[i]->type == SUM) ||
      (stencils[i]->type == WEIGHT_SUM),
      "ERROR(stencils_merge): unsupported stencil type")
    switch (stencils[i]->type) {
      case (FIXED):
        fixed_value += stencils[i]->data.fixed.value * w[i];
        fixed_count++;
        break;
      case (DIRECT):
        direct_count++;
        break;
      case (SUM):
        sum_count++;
        break;
      default:
      case (WEIGHT_SUM):
        wsum_count++;
        break;
    };
  }

  struct interp_weight_stencil merge_stencil;

  YAC_ASSERT(
    (fixed_count > 0) || (wsum_count > 0) ||
    (sum_count > 0) || (direct_count > 0),
    "ERROR(stencils_merge): unknown error")
  if (fixed_count > 0) {

    YAC_ASSERT(
      (direct_count + sum_count + wsum_count) <= 0,
      "ERROR(stencils_merge): invalid stencil combination")

    merge_stencil = **stencils;
    merge_stencil.data.fixed.value = fixed_value;
  } else if (wsum_count > 0)
    merge_stencil =
      stencils_merge_wsum(stencils, w, num_stencils);
  else if ((sum_count > 0) || (direct_count > 0))
    merge_stencil =
      stencils_merge_sum(stencils, w, num_stencils);

  merge_stencil.tgt = copy_remote_point(point);

  return merge_stencil;
}

static void get_stencils_pack_sizes(
  struct interp_weight_stencil * stencils, size_t count, size_t * pack_order,
  int * pack_sizes, MPI_Datatype point_info_dt, MPI_Comm comm) {

  int pack_size_type;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &pack_size_type), comm);

  for (size_t i = 0; i < count; ++i) {

    struct interp_weight_stencil * curr_stencil = stencils + pack_order[i];
    int (*func_pack_size)(
      struct interp_weight_stencil * stencil, MPI_Datatype point_info_dt,
      MPI_Comm comm);
    YAC_ASSERT(
      (curr_stencil->type == FIXED) ||
      (curr_stencil->type == DIRECT) ||
      (curr_stencil->type == SUM) ||
      (curr_stencil->type == WEIGHT_SUM) ||
      (curr_stencil->type == DIRECT_MF) ||
      (curr_stencil->type == SUM_MF) ||
      (curr_stencil->type == WEIGHT_SUM_MF),
      "ERROR(get_stencils_pack_sizes): invalid stencil type")
    switch (curr_stencil->type) {
      case(FIXED):
        func_pack_size = get_stencil_pack_size_fixed;
        break;
      case(DIRECT):
        func_pack_size = get_stencil_pack_size_direct;
        break;
      case(SUM):
        func_pack_size = get_stencil_pack_size_sum;
        break;
      default:
      case(WEIGHT_SUM):
        func_pack_size = get_stencil_pack_size_wsum;
        break;
      case(DIRECT_MF):
        func_pack_size = get_stencil_pack_size_direct_mf;
        break;
      case(SUM_MF):
        func_pack_size = get_stencil_pack_size_sum_mf;
        break;
      case(WEIGHT_SUM_MF):
        func_pack_size = get_stencil_pack_size_wsum_mf;
        break;
    };
    pack_sizes[i] = pack_size_type +
                    yac_remote_point_get_pack_size(
                      &(curr_stencil->tgt), point_info_dt, comm) +
                    func_pack_size(curr_stencil, point_info_dt, comm);
  }
}

static void pack_stencil_fixed(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  UNUSED(point_info_dt);

  // fixed value
  yac_mpi_call(
    MPI_Pack(&(stencil->data.fixed.value), 1, MPI_DOUBLE, buffer, buffer_size,
             position, comm), comm);
}

static void pack_stencil_direct(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // src
  yac_remote_point_pack(
    &(stencil->data.direct.src), buffer, buffer_size, position, point_info_dt,
    comm);
}

static void pack_stencil_sum(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // srcs
  yac_remote_points_pack(
    stencil->data.sum.srcs, buffer, buffer_size, position, point_info_dt, comm);
}

static void pack_stencil_wsum(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // srcs
  yac_remote_points_pack(
    stencil->data.weight_sum.srcs, buffer, buffer_size, position,
    point_info_dt, comm);
  // weights
  yac_mpi_call(
    MPI_Pack(stencil->data.weight_sum.weights,
             (int)(stencil->data.weight_sum.srcs->count), MPI_DOUBLE,
             buffer, buffer_size, position, comm), comm);
}

static void pack_stencil_direct_mf(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // src
  yac_remote_point_pack(
    &(stencil->data.direct_mf.src), buffer, buffer_size, position,
    point_info_dt, comm);

  // field_idx
  uint64_t temp_field_idx = (uint64_t)(stencil->data.direct_mf.field_idx);
  yac_mpi_call(
    MPI_Pack(&temp_field_idx, 1, MPI_UINT64_T,
             buffer, buffer_size, position, comm), comm);
}

static void pack_stencil_sum_mf(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // srcs
  yac_remote_points_pack(
    stencil->data.sum_mf.srcs, buffer, buffer_size, position,
    point_info_dt, comm);

  size_t count = stencil->data.sum_mf.srcs->count;
  // field_indices
  uint64_t * temp_field_indices = xmalloc(count * sizeof(*temp_field_indices));
  for (size_t i = 0; i < count; ++i)
    temp_field_indices[i] =
      (uint64_t)(stencil->data.sum_mf.field_indices[i]);
  yac_mpi_call(
    MPI_Pack(temp_field_indices, (int)count, MPI_UINT64_T,
             buffer, buffer_size, position, comm), comm);
  free(temp_field_indices);
}

static void pack_stencil_wsum_mf(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // srcs
  yac_remote_points_pack(
    stencil->data.weight_sum_mf.srcs, buffer, buffer_size, position,
    point_info_dt, comm);

  size_t count = stencil->data.weight_sum_mf.srcs->count;
  // weights
  yac_mpi_call(
    MPI_Pack(stencil->data.weight_sum_mf.weights, (int)count, MPI_DOUBLE,
             buffer, buffer_size, position, comm), comm);
  // field_indices
  uint64_t * temp_field_indices = xmalloc(count * sizeof(*temp_field_indices));
  for (size_t i = 0; i < count; ++i)
    temp_field_indices[i] =
      (uint64_t)(stencil->data.weight_sum_mf.field_indices[i]);
  yac_mpi_call(
    MPI_Pack(temp_field_indices, (int)count, MPI_UINT64_T,
             buffer, buffer_size, position, comm), comm);
  free(temp_field_indices);
}

static void pack_stencils(
  struct interp_weight_stencil * stencils, size_t count, size_t * pack_order,
  void ** pack_data, int * pack_sizes, MPI_Datatype point_info_dt,
  MPI_Comm comm) {

  get_stencils_pack_sizes(
    stencils, count, pack_order, pack_sizes, point_info_dt, comm);

  size_t pack_buffer_size = 0;
  for (size_t i = 0; i < count; ++i)
    pack_buffer_size += (size_t)(pack_sizes[i]);

  void * pack_data_ = xmalloc(pack_buffer_size);
  size_t total_pack_size = 0;

  for (size_t i = 0; i < count; ++i) {

    struct interp_weight_stencil * curr_stencil = stencils + pack_order[i];
    void (*func_pack)(
      struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
      int * position, MPI_Datatype point_info_dt, MPI_Comm comm);

    YAC_ASSERT(
      (curr_stencil->type == FIXED) ||
      (curr_stencil->type == DIRECT) ||
      (curr_stencil->type == SUM) ||
      (curr_stencil->type == WEIGHT_SUM) ||
      (curr_stencil->type == DIRECT_MF) ||
      (curr_stencil->type == SUM_MF) ||
      (curr_stencil->type == WEIGHT_SUM_MF),
      "ERROR(pack_stencils): invalid stencil type")
    switch (curr_stencil->type) {
      default:
      case(FIXED):
        func_pack = pack_stencil_fixed;
        break;
      case(DIRECT):
        func_pack = pack_stencil_direct;
        break;
      case(SUM):
        func_pack = pack_stencil_sum;
        break;
      case(WEIGHT_SUM):
        func_pack = pack_stencil_wsum;
        break;
      case(DIRECT_MF):
        func_pack = pack_stencil_direct_mf;
        break;
      case(SUM_MF):
        func_pack = pack_stencil_sum_mf;
        break;
      case(WEIGHT_SUM_MF):
        func_pack = pack_stencil_wsum_mf;
        break;
    };

    int position = 0;
    int type = (int)curr_stencil->type;
    void * buffer = (void*)((char*)pack_data_ + total_pack_size);
    int buffer_size = pack_sizes[i];

    // type
    yac_mpi_call(
      MPI_Pack(&type, 1, MPI_INT, buffer, buffer_size, &position, comm), comm);
    // tgt
    yac_remote_point_pack(&(curr_stencil->tgt), buffer, buffer_size,
                          &position, point_info_dt, comm);
    // stencil data
    func_pack(curr_stencil, buffer, buffer_size, &position, point_info_dt, comm);

    YAC_ASSERT_F(
      pack_sizes[i] >= position,
      "ERROR(pack_stencils): "
      "actual pack size is bigger then computed one (%d > %d)",
      position, pack_sizes[i]);

    pack_sizes[i] = position;
    total_pack_size += (size_t)position;
  }

  *pack_data = xrealloc(pack_data_, total_pack_size);
}

static void unpack_stencil_fixed(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  UNUSED(point_info_dt);

  // fixed value
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, &(stencil->data.fixed.value), 1,
               MPI_DOUBLE, comm), comm);
}

static void unpack_stencil_direct(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // src
  yac_remote_point_unpack(
    buffer, buffer_size, position, &stencil->data.direct.src,
    point_info_dt, comm);
}

static void unpack_stencil_sum(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // srcs
  stencil->data.weight_sum.weights = NULL;
  yac_remote_points_unpack(
    buffer, buffer_size, position, &(stencil->data.sum.srcs), point_info_dt,
    comm);
}

static void unpack_stencil_wsum(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // srcs
  yac_remote_points_unpack(
    buffer, buffer_size, position, &(stencil->data.weight_sum.srcs),
    point_info_dt, comm);

  size_t count = stencil->data.weight_sum.srcs->count;

  stencil->data.weight_sum.weights =
    xmalloc(count * sizeof(*(stencil->data.weight_sum.weights)));

  // weights
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position, stencil->data.weight_sum.weights,
               (int)count, MPI_DOUBLE, comm), comm);
}

static void unpack_stencil_direct_mf(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // src
  yac_remote_point_unpack(
    buffer, buffer_size, position, &stencil->data.direct_mf.src,
    point_info_dt, comm);

  // field_idx
  uint64_t temp_field_idx;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &temp_field_idx,
      1, MPI_UINT64_T, comm), comm);
  stencil->data.direct_mf.field_idx = (size_t)temp_field_idx;
}

static void unpack_stencil_sum_mf(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // srcs
  yac_remote_points_unpack(
    buffer, buffer_size, position, &(stencil->data.sum_mf.srcs),
    point_info_dt, comm);

  size_t count = stencil->data.sum_mf.srcs->count;

  uint64_t * temp_field_indices = xmalloc(count * sizeof(*temp_field_indices));
  stencil->data.sum_mf.field_indices =
    xmalloc(count * sizeof(*(stencil->data.sum_mf.field_indices)));

  // field_indices
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, temp_field_indices,
      (int)count, MPI_UINT64_T, comm), comm);
  for (size_t i = 0; i < count; ++i)
    stencil->data.sum_mf.field_indices[i] = (size_t)(temp_field_indices[i]);
  free(temp_field_indices);
}

static void unpack_stencil_wsum_mf(
  struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype point_info_dt, MPI_Comm comm) {

  // srcs
  yac_remote_points_unpack(
    buffer, buffer_size, position, &(stencil->data.weight_sum_mf.srcs),
    point_info_dt, comm);

  size_t count = stencil->data.weight_sum_mf.srcs->count;

  stencil->data.weight_sum_mf.weights =
    xmalloc(count * sizeof(*(stencil->data.weight_sum_mf.weights)));

  // weights
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, stencil->data.weight_sum_mf.weights,
      (int)count, MPI_DOUBLE, comm), comm);

  uint64_t * temp_field_indices = xmalloc(count * sizeof(*temp_field_indices));
  stencil->data.weight_sum_mf.field_indices =
    xmalloc(count * sizeof(*(stencil->data.weight_sum_mf.field_indices)));

  // field_indices
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, temp_field_indices,
      (int)count, MPI_UINT64_T, comm), comm);
  for (size_t i = 0; i < count; ++i)
    stencil->data.weight_sum_mf.field_indices[i] =
      (size_t)(temp_field_indices[i]);
  free(temp_field_indices);
}

static void unpack_stencils(
  struct interp_weight_stencil * stencils, size_t count,
  void * packed_data, size_t packed_data_size,
  MPI_Datatype point_info_dt, MPI_Comm comm) {

  for (size_t i = 0, offset = 0; i < count; ++i) {

    YAC_ASSERT(
      packed_data_size >= offset,
      "ERROR(unpack_stencils): invalid offset");

    int position = 0;
    void * curr_buffer = (void*)((unsigned char*)packed_data + offset);
    int buffer_size = (int)(MIN(packed_data_size - offset, INT_MAX));
    struct interp_weight_stencil * curr_stencil = stencils + i;

    int type;
    yac_mpi_call(
      MPI_Unpack(
        curr_buffer, buffer_size, &position, &type, 1, MPI_INT, comm), comm);

    void (*func_unpack)(
      struct interp_weight_stencil * stencil, void * buffer, int buffer_size,
      int * position, MPI_Datatype point_info_dt, MPI_Comm comm);

    YAC_ASSERT(
      (type == FIXED) ||
      (type == DIRECT) || (type == SUM) || (type == WEIGHT_SUM) ||
      (type == DIRECT_MF) || (type == SUM_MF) || (type == WEIGHT_SUM_MF),
      "ERROR(unpack_stencils): invalid stencil type")
    switch (type) {
      case(FIXED):
        func_unpack = unpack_stencil_fixed;
        break;
      case(DIRECT):
        func_unpack = unpack_stencil_direct;
        break;
      case(SUM):
        func_unpack = unpack_stencil_sum;
        break;
      default:
      case(WEIGHT_SUM):
        func_unpack = unpack_stencil_wsum;
        break;
      case(DIRECT_MF):
        func_unpack = unpack_stencil_direct_mf;
        break;
      case(SUM_MF):
        func_unpack = unpack_stencil_sum_mf;
        break;
      case(WEIGHT_SUM_MF):
        func_unpack = unpack_stencil_wsum_mf;
        break;
    };

    curr_stencil->type =
      (enum yac_interp_weight_stencil_type)type;
    yac_remote_point_unpack(
      curr_buffer, buffer_size, &position, &(curr_stencil->tgt),
      point_info_dt, comm);
    func_unpack(
      curr_stencil, curr_buffer, buffer_size, &position, point_info_dt, comm);
    offset += (size_t)position;
  }
}

static struct interp_weight_stencil * exchange_stencils(
  MPI_Comm comm, struct interp_weight_stencil * stencils,
  size_t * stencil_indices,
  size_t * stencil_sendcounts, size_t * stencil_recvcounts) {

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  YAC_ASSERT(
    stencil_sendcounts[comm_rank] == stencil_recvcounts[comm_rank],
    "ERROR(exchange_stencils): error in arguments")

  size_t send_count = 0, recv_count = 0;
  size_t local_send_offset = 0;
  size_t local_recv_offset = 0;
  size_t local_count = (size_t)(stencil_sendcounts[comm_rank]);
  for (int i = 0; i < comm_rank; ++i) {
    send_count += stencil_sendcounts[i];
    recv_count += stencil_recvcounts[i];
    local_send_offset += stencil_sendcounts[i];
    local_recv_offset += stencil_recvcounts[i];
  }
  local_send_offset = send_count;
  local_recv_offset = recv_count;
  stencil_sendcounts[comm_rank] = 0;
  stencil_recvcounts[comm_rank] = 0;
  for (int i = comm_rank + 1; i < comm_size; ++i) {
    send_count += stencil_sendcounts[i];
    recv_count += stencil_recvcounts[i];
  }

  struct interp_weight_stencil * new_stencils =
    xmalloc((recv_count + local_count) * sizeof(*new_stencils));
  size_t * local_stencil_indices =
    xmalloc(local_count * sizeof(*local_stencil_indices));
  memcpy(local_stencil_indices, stencil_indices + local_send_offset,
         local_count * sizeof(*local_stencil_indices));

  // remove the local stencil indices
  memmove(
    stencil_indices + local_send_offset,
    stencil_indices + local_send_offset + local_count,
    (send_count - local_send_offset) * sizeof(*stencil_indices));

  // pack the stencils that need to be send to other processes
  void * send_buffer;
  int * pack_sizes = xmalloc(send_count * sizeof(*pack_sizes));
  MPI_Datatype point_info_dt = yac_get_remote_point_info_mpi_datatype(comm);
  pack_stencils(
    stencils, send_count, stencil_indices, &send_buffer, pack_sizes,
    point_info_dt, comm);

  size_t * sendcounts, * recvcounts, * sdispls, *rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  send_count = 0;
  for (int rank = 0; rank < comm_size; ++rank) {
    size_t sendcount = 0;
    int curr_num_stencils = stencil_sendcounts[rank];
    for (int j = 0; j < curr_num_stencils; ++j, ++send_count)
      sendcount += (size_t)(pack_sizes[send_count]);
    sendcounts[rank] = sendcount;
  }
  free(pack_sizes);

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t recv_size = recvcounts[comm_size - 1] + rdispls[comm_size - 1];

  void * recv_buffer = xmalloc(recv_size);

  // exchange stencils
  yac_alltoallv_packed_p2p(
    send_buffer, sendcounts, sdispls+1,
    recv_buffer, recvcounts, rdispls, comm);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  free(send_buffer);

  // unpack stencils
  unpack_stencils(
    new_stencils, recv_count,
    recv_buffer, recv_size, point_info_dt, comm);
  yac_mpi_call(MPI_Type_free(&point_info_dt), comm);
  free(recv_buffer);

  memmove(new_stencils + local_recv_offset + local_count,
          new_stencils + local_recv_offset ,
          (recv_count - local_recv_offset ) * sizeof(*new_stencils));
  for (size_t i = 0; i < local_count; ++i, ++local_recv_offset )
    new_stencils[local_recv_offset] =
      copy_interp_weight_stencil(
        stencils + local_stencil_indices[i],
        stencils[local_stencil_indices[i]].tgt);
  free(local_stencil_indices);

  return new_stencils;
}

static struct interp_weight_stencil *  yac_interp_weights_get_stencils(
  struct yac_interp_weights * weights, size_t * stencil_indices,
  int * stencil_ranks, size_t count) {

  MPI_Comm comm = weights->comm;
  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  YAC_ASSERT(
    count <= INT_MAX,
    "ERROR(yac_interp_weights_get_stencils): count exceeds INT_MAX");

  size_t * reorder_idx = xmalloc(count * sizeof(*reorder_idx));
  for (size_t i = 0; i < count; ++i) reorder_idx[i] = i;

  yac_quicksort_index_int_size_t_size_t(
    stencil_ranks, count, stencil_indices, reorder_idx);

  // exchange requested stencils indices
  size_t * sendcounts, * recvcounts, * sdispls, *rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);
  for (size_t i = 0; i < count; ++i) sendcounts[stencil_ranks[i]]++;
  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);
  size_t recv_count =
    rdispls[comm_size - 1] + recvcounts[comm_size - 1];
  uint64_t * uint64_t_buffer =
    xmalloc((count + recv_count) * sizeof(*uint64_t_buffer));
  uint64_t * send_stencil_indices = uint64_t_buffer;
  uint64_t * recv_stencil_indices = uint64_t_buffer + count;
  for (size_t i = 0; i < count; ++i)
    send_stencil_indices[i] = (uint64_t)(stencil_indices[i]);
  yac_alltoallv_uint64_p2p(
    send_stencil_indices, sendcounts, sdispls+1,
    recv_stencil_indices, recvcounts, rdispls, comm);

  // exchange stencils
  size_t * exchange_stencil_indices =
    xmalloc(recv_count * sizeof(*exchange_stencil_indices));
  for (size_t i = 0; i < recv_count; ++i) {
    YAC_ASSERT(
      (size_t)(recv_stencil_indices[i]) < weights->stencils_size,
      "ERROR(yac_interp_weights_get_stencils): invalid stencil index");
    exchange_stencil_indices[i] = (size_t)(recv_stencil_indices[i]);
  }
  free(uint64_t_buffer);
  struct interp_weight_stencil * stencils =
    exchange_stencils(comm, weights->stencils, exchange_stencil_indices,
                      recvcounts, sendcounts);
  free(exchange_stencil_indices);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  // sort received stencils into original order
  struct interp_weight_stencil * sorted_stencils =
    xmalloc(count * sizeof(*sorted_stencils));
  for (size_t i = 0; i < count; ++i)
    sorted_stencils[reorder_idx[i]] = stencils[i];
  free(stencils);
  free(reorder_idx);

  return sorted_stencils;
}

static void yac_interp_weight_stencils_delete(
  struct interp_weight_stencil * stencils, size_t count);

void yac_interp_weights_wcopy_weights(
  struct yac_interp_weights * weights, struct remote_points * tgts,
  size_t * num_stencils_per_tgt, size_t * stencil_indices,
  int * stencil_ranks, double * w) {

  size_t count = (tgts != NULL)?tgts->count:0;
  MPI_Comm comm = weights->comm;
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // count the number of missing stencils
  size_t total_num_stencils = 0;
  size_t max_num_stencils_per_tgt = 0;
  for (size_t i = 0; i < count; ++i) {
    size_t curr_num_stencils_per_tgt = num_stencils_per_tgt[i];
    if (curr_num_stencils_per_tgt > max_num_stencils_per_tgt)
      max_num_stencils_per_tgt = curr_num_stencils_per_tgt;
    total_num_stencils += num_stencils_per_tgt[i];
  }
  size_t num_missing_stencils = 0;
  for (size_t i = 0; i < total_num_stencils; ++i)
    if (stencil_ranks[i] != comm_rank) num_missing_stencils++;

  // get missing stencils
  size_t * missing_stencil_indices =
    xmalloc(num_missing_stencils * sizeof(*missing_stencil_indices));
  int * missing_stencil_ranks =
    xmalloc(num_missing_stencils * sizeof(*missing_stencil_ranks));
  for (size_t i = 0, j = 0; i < total_num_stencils; ++i) {
    if (stencil_ranks[i] != comm_rank) {
      missing_stencil_indices[j] = stencil_indices[i];
      missing_stencil_ranks[j] = stencil_ranks[i];
      ++j;
    }
  }
  struct interp_weight_stencil * missing_stencils =
    yac_interp_weights_get_stencils(
      weights, missing_stencil_indices, missing_stencil_ranks,
      num_missing_stencils);
  free(missing_stencil_ranks);
  free(missing_stencil_indices);

  // merge stencils to generate new ones
  {
    struct interp_weight_stencil * stencils = weights->stencils;
    size_t stencils_array_size = weights->stencils_array_size;
    size_t stencils_size = weights->stencils_size;

    struct interp_weight_stencil ** stencils_buffer =
      xmalloc(max_num_stencils_per_tgt * sizeof(*stencils_buffer));

    ENSURE_ARRAY_SIZE(
      stencils, stencils_array_size, stencils_size + count);

    for (size_t i = 0, j = 0; i < count;
         ++i, ++stencils_size) {

      size_t curr_num_stencils = num_stencils_per_tgt[i];
      for (size_t k = 0; k < curr_num_stencils; ++k)
        stencils_buffer[k] =
          (stencil_ranks[k] == comm_rank)?
            (stencils + stencil_indices[k]):(missing_stencils + (j++));

      stencils[stencils_size] =
        stencils_merge(stencils_buffer, w, curr_num_stencils, tgts->data[i]);
      w += curr_num_stencils;
      stencil_indices += curr_num_stencils;
      stencil_ranks += curr_num_stencils;
    }

    weights->stencils = stencils;
    weights->stencils_array_size = stencils_array_size;
    weights->stencils_size = stencils_size;

    free(stencils_buffer);
  }

  yac_interp_weight_stencils_delete(missing_stencils, num_missing_stencils);
}

static int compute_owner(int * ranks, size_t count) {

  YAC_ASSERT(count != 0, "ERROR(compute_owner): count == 0")

  yac_quicksort_index(ranks, count, NULL);

  int best_rank = -1;
  size_t best_rank_count = 0;

  size_t curr_rank_count = 1;
  int prev_rank = ranks[0];

  for (size_t i = 1; i < count; ++i, ++curr_rank_count) {
    int curr_rank = ranks[i];
    if (prev_rank != curr_rank) {
      if (curr_rank_count > best_rank_count) {
        best_rank = prev_rank;
        best_rank_count = curr_rank_count;
      }
      prev_rank = curr_rank;
      curr_rank_count = 0;
    }
  }

  return (curr_rank_count > best_rank_count)?prev_rank:best_rank;
}

static struct interp_weight_stencils_wsum_mf * generate_w_sum_mf_stencils(
  struct interp_weight_stencil * stencils, size_t count,
  enum yac_interp_weight_stencil_type stencil_type) {

  // compute total number of links
  size_t total_num_links = 0;

  for (size_t i = 0; i < count; ++i) {
    YAC_ASSERT(
      stencils[i].type == stencil_type,
      "ERROR(generate_w_sum_mf_stencils): wrong stencil type")
    // due to the data layout this works for "sum" and "weight_sum"
    total_num_links += stencils[i].data.weight_sum.srcs->count;
  }

  struct interp_weight_stencils_wsum_mf_buffer * temp =
    xmalloc(sizeof(*temp) + total_num_links * sizeof(temp->buffer[0]));
  struct interp_weight_stencils_wsum_mf * wsum_stencils =
    (struct interp_weight_stencils_wsum_mf *)temp;
  wsum_stencils->data = xmalloc(count * sizeof(*(wsum_stencils->data)));
  wsum_stencils->count = count;

  // extract data from stencils
  for (size_t i = 0, k = 0; i < count; ++i) {
    struct interp_weight_stencil_wsum_mf * curr_wsum_stencil =
      wsum_stencils->data + i;
    struct interp_weight_stencil_wsum_mf_weight * curr_links =
      &(temp->buffer[k]);
    size_t curr_stencil_size = stencils[i].data.weight_sum.srcs->count;
    struct remote_point * curr_srcs = stencils[i].data.weight_sum.srcs->data;
    curr_wsum_stencil->tgt = copy_remote_point(stencils[i].tgt);
    curr_wsum_stencil->count = curr_stencil_size;
    curr_wsum_stencil->data = curr_links;
    for (size_t j = 0; j < curr_stencil_size; ++j) {
      int curr_count = curr_srcs[j].data.count;
      YAC_ASSERT(
        curr_count >= 1,
        "ERROR(generate_w_sum_mf_stencils): global src id no found")
      curr_links[j].src =
        (curr_count == 1)?
          (curr_srcs[j].data.data.single):(curr_srcs[j].data.data.multi[0]);
      YAC_ASSERT(
        (stencil_type == SUM) || (stencil_type == WEIGHT_SUM) ||
        (stencil_type == SUM_MF) || (stencil_type == WEIGHT_SUM_MF),
        "ERROR(generate_w_sum_mf_stencils): unsupported stencil type")
      switch(stencil_type) {
        default:
        case(SUM):
          curr_links[j].weight = 1.0;
          curr_links[j].src_field_idx = 0;
          break;
        case(WEIGHT_SUM):
          curr_links[j].weight = stencils[i].data.weight_sum.weights[j];
          curr_links[j].src_field_idx = 0;
          break;
        case(SUM_MF):
          curr_links[j].weight = 1.0;
          curr_links[j].src_field_idx =
            stencils[i].data.sum_mf.field_indices[j];
          break;
        case(WEIGHT_SUM_MF):
          curr_links[j].weight = stencils[i].data.weight_sum_mf.weights[j];
          curr_links[j].src_field_idx =
            stencils[i].data.weight_sum_mf.field_indices[j];
          break;
      };
    }
    k += curr_stencil_size;
  }

  return wsum_stencils;
}

static MPI_Datatype get_wsum_mf_weight_mpi_datatype(MPI_Comm comm) {

  struct interp_weight_stencil_wsum_mf_weight dummy;
  MPI_Datatype dt;
  int array_of_blocklengths[] = {1, 1, 1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.src.rank) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.src.orig_pos) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.src_field_idx) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.weight) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] =
    {MPI_INT, MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE};
  yac_mpi_call(
    MPI_Type_create_struct(4, array_of_blocklengths, array_of_displacements,
                           array_of_types, &dt), comm);
  return yac_create_resized(dt, sizeof(dummy), comm);
}

static int get_stencil_wsum_mf_pack_size(
  struct interp_weight_stencil_wsum_mf * stencil,
  MPI_Datatype wsum_mf_weight_dt, MPI_Datatype point_info_dt, MPI_Comm comm) {

  int pack_size_count,
      pack_size_weights,
      pack_size_tgt;

  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &pack_size_count), comm);
  yac_mpi_call(
    MPI_Pack_size(
      (int)(stencil->count), wsum_mf_weight_dt, comm, &pack_size_weights), comm);
  pack_size_tgt =
    yac_remote_point_get_pack_size(&(stencil->tgt), point_info_dt, comm);

  return pack_size_count + pack_size_weights + pack_size_tgt;
}

static void pack_stencils_wsum_mf(
  struct interp_weight_stencil_wsum_mf * wsum_stencils, size_t count,
  size_t * pack_order, void ** pack_data, int * pack_sizes,
  int * weight_counts, MPI_Comm comm) {

  MPI_Datatype wsum_mf_weight_dt = get_wsum_mf_weight_mpi_datatype(comm);
  MPI_Datatype point_info_dt = yac_get_remote_point_info_mpi_datatype(comm);

  // get the pack sizes and the upper bound for the pack buffer size
  size_t temp_total_pack_size = 0;
  for (size_t i = 0; i < count; ++i) {
    temp_total_pack_size +=
      (pack_sizes[i] =
        get_stencil_wsum_mf_pack_size(
          wsum_stencils + pack_order[i],
          wsum_mf_weight_dt, point_info_dt, comm));
  }

  void * pack_data_ = xmalloc(temp_total_pack_size);
  size_t total_pack_size = 0;

  // pack the stencils
  for (size_t i = 0; i < count; ++i) {

    size_t idx = pack_order[i];

    int position = 0;
    void * buffer = (void*)((unsigned char*)pack_data_ + total_pack_size);
    int buffer_size = pack_sizes[i];
    int curr_count = wsum_stencils[idx].count;

    // tgt
    yac_remote_point_pack(
      &(wsum_stencils[idx].tgt), buffer, buffer_size, &position,
      point_info_dt, comm);
    // weight count
    yac_mpi_call(
      MPI_Pack(&curr_count, 1, MPI_INT, buffer, buffer_size, &position, comm), comm);
    // weights
    yac_mpi_call(
      MPI_Pack(wsum_stencils[idx].data, curr_count, wsum_mf_weight_dt,
               buffer, buffer_size, &position, comm), comm);

    pack_sizes[i] = position;
    weight_counts[i] = curr_count;
    total_pack_size += (size_t)position;
  }

  yac_mpi_call(MPI_Type_free(&point_info_dt), comm);
  yac_mpi_call(MPI_Type_free(&wsum_mf_weight_dt), comm);

  *pack_data = xrealloc(pack_data_, total_pack_size);
}

static size_t unpack_stencils_wsum_mf(
  struct interp_weight_stencil_wsum_mf * wsum_stencils,
  struct interp_weight_stencil_wsum_mf_weight * weight_buffer, size_t count,
  void * packed_data, size_t packed_data_size, MPI_Comm comm) {

  MPI_Datatype wsum_mf_weight_dt = get_wsum_mf_weight_mpi_datatype(comm);
  MPI_Datatype point_info_dt = yac_get_remote_point_info_mpi_datatype(comm);

  size_t weight_offset = 0;
  for (size_t i = 0, offset = 0; i < count; ++i) {

    int position = 0;
    void * curr_buffer = (void*)((char*)packed_data + offset);
    int buffer_size = (int)(packed_data_size - offset);
    struct interp_weight_stencil_wsum_mf * curr_wsum_stencil =
      wsum_stencils + i;

    struct remote_point tgt;
    struct interp_weight_stencil_wsum_mf_weight * curr_weights =
      weight_buffer + weight_offset;
    int weight_count;
    yac_remote_point_unpack(
      curr_buffer, buffer_size, &position, &tgt, point_info_dt, comm);
    yac_mpi_call(
      MPI_Unpack(curr_buffer, buffer_size, &position,
                 &weight_count, 1, MPI_INT, comm),
      comm);
    yac_mpi_call(
      MPI_Unpack(curr_buffer, buffer_size, &position,
                 curr_weights, weight_count, wsum_mf_weight_dt, comm), comm);

    curr_wsum_stencil->tgt = tgt;
    curr_wsum_stencil->data = curr_weights;
    curr_wsum_stencil->count = (size_t)weight_count;

    weight_offset += (size_t)weight_count;
    offset += (size_t)position;
  }

  yac_mpi_call(MPI_Type_free(&point_info_dt), comm);
  yac_mpi_call(MPI_Type_free(&wsum_mf_weight_dt), comm);

  return weight_offset;
}

static struct interp_weight_stencils_wsum_mf * redist_wsum_mf_stencils(
  MPI_Comm comm, struct interp_weight_stencils_wsum_mf * wsum_stencils_data,
  int * stencil_owner, size_t * reorder_idx, size_t num_owners) {

  struct interp_weight_stencil_wsum_mf * wsum_stencils =
    wsum_stencils_data->data;

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t local_weight_count = 0;
  size_t local_count = 0;
  for (size_t i = 0; i < num_owners; ++i) {
    if (stencil_owner[i] == comm_rank) {
      local_weight_count += wsum_stencils[reorder_idx[i]].count;
      stencil_owner[i] = INT_MAX;
      ++local_count;
    }
  }
  yac_quicksort_index_int_size_t(stencil_owner, num_owners, reorder_idx);

  size_t send_count = num_owners - local_count;

  // pack the stencils that need to be send to other processes
  void * send_buffer;
  int * pack_sizes = xmalloc(2 * send_count * sizeof(*pack_sizes));
  int * weight_counts = pack_sizes + send_count;
  pack_stencils_wsum_mf(wsum_stencils, send_count, reorder_idx, &send_buffer,
                        pack_sizes, weight_counts, comm);

  size_t * sendcounts, * recvcounts, * sdispls, *rdispls;
  yac_get_comm_buffers(
    3, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  for (size_t i = 0; i < send_count; ++i) {
    int curr_rank = stencil_owner[i];
    sendcounts[3 * curr_rank + 0]++;
    sendcounts[3 * curr_rank + 1] += (size_t)(pack_sizes[i]);
    sendcounts[3 * curr_rank + 2] += (size_t)(weight_counts[i]);
  }
  free(pack_sizes);

  // exchange the number of stencils to be exchanged and the total pack sizes
  yac_mpi_call(MPI_Alltoall(sendcounts, 3, YAC_MPI_SIZE_T,
                            recvcounts, 3, YAC_MPI_SIZE_T, comm), comm);

  size_t recv_count = 0;
  size_t recv_size = 0;
  size_t recv_weight_count = 0;
  size_t saccu = 0, raccu = 0;
  for (int i = 0; i < comm_size; ++i) {
    sdispls[i] = saccu;
    rdispls[i] = raccu;
    recv_count += recvcounts[3 * i + 0];
    recv_size += recvcounts[3 * i + 1];
    recv_weight_count += recvcounts[3 * i + 2];
    saccu += sendcounts[3 * i + 1];
    raccu += recvcounts[3 * i + 1];
    sendcounts[i] = sendcounts[3 * i + 1];
    recvcounts[i] = recvcounts[3 * i + 1];
  }

  void * recv_buffer = xmalloc(recv_size);

  // exchange stencils
  yac_alltoallv_packed_p2p(
    send_buffer, sendcounts, sdispls, recv_buffer, recvcounts, rdispls, comm);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  free(send_buffer);

  struct interp_weight_stencils_wsum_mf_buffer * temp =
    xmalloc(sizeof(*temp) +
           (local_weight_count + recv_weight_count) * sizeof(temp->buffer[0]));
  struct interp_weight_stencils_wsum_mf * new_wsum_stencils_data =
    (struct interp_weight_stencils_wsum_mf *)temp;
  struct interp_weight_stencil_wsum_mf * new_wsum_stencils =
    ((new_wsum_stencils_data->data =
      xmalloc((local_count + recv_count) *
              sizeof(*(new_wsum_stencils_data->data)))));
  new_wsum_stencils_data->count = local_count + recv_count;

  // unpack stencils
  size_t weight_offset =
    unpack_stencils_wsum_mf(
      new_wsum_stencils, &(temp->buffer[0]), recv_count,
      recv_buffer, recv_size, comm);
  free(recv_buffer);
  new_wsum_stencils += recv_count;
  struct interp_weight_stencil_wsum_mf_weight * weight_buffer =
    &(temp->buffer[weight_offset]);

  // copy the stencils that stay locally into the new stencil array
  yac_quicksort_index_size_t_size_t(reorder_idx + send_count, local_count, NULL);
  for (size_t i = 0, weight_offset = 0; i < local_count; ++i) {
    struct interp_weight_stencil_wsum_mf * curr_wsum_stencil =
      wsum_stencils + reorder_idx[i + send_count];
    struct interp_weight_stencil_wsum_mf * curr_new_wsum_stencil =
      new_wsum_stencils + i;
    struct interp_weight_stencil_wsum_mf_weight * curr_new_weights =
      weight_buffer + weight_offset;
    size_t curr_stencil_size = curr_wsum_stencil->count;
    curr_new_wsum_stencil->tgt = copy_remote_point(curr_wsum_stencil->tgt);
    curr_new_wsum_stencil->count = curr_stencil_size;
    curr_new_wsum_stencil->data = curr_new_weights;
    memcpy(curr_new_weights, curr_wsum_stencil->data,
           curr_stencil_size * sizeof(*curr_new_weights));
    weight_offset += curr_stencil_size;
  }

  return new_wsum_stencils_data;
}

static struct interp_weight_stencils_wsum_mf * redist_wsum_mf_stencils_src(
  MPI_Comm comm, struct interp_weight_stencils_wsum_mf * wsum_stencils_data) {

  struct interp_weight_stencil_wsum_mf * wsum_stencils =
    wsum_stencils_data->data;
  size_t count = wsum_stencils_data->count;

  // determine maximum stencil size
  size_t max_stencil_size = 0;
  for (size_t i = 0; i < count; ++i) {
    size_t curr_stencil_size = wsum_stencils[i].count;
    if (curr_stencil_size > max_stencil_size)
      max_stencil_size = curr_stencil_size;
  }

  // determine source process for each stencil
  int * rank_buffer =
    xmalloc((count + max_stencil_size) * sizeof(*rank_buffer));
  int * stencil_owner = rank_buffer;
  int * stencil_owners = rank_buffer + count;
  size_t * reorder_idx = xmalloc(count * sizeof(*reorder_idx));
  for (size_t i = 0; i < count; ++i) {
    size_t curr_stencil_size = wsum_stencils[i].count;
    struct interp_weight_stencil_wsum_mf_weight * curr_weights =
      wsum_stencils[i].data;
    for (size_t j = 0; j < curr_stencil_size; ++j)
      stencil_owners[j] = curr_weights[j].src.rank;
    stencil_owner[i] = compute_owner(stencil_owners, curr_stencil_size);
    reorder_idx[i] = i;
  }

  struct interp_weight_stencils_wsum_mf * new_wsum_stencils_data =
    redist_wsum_mf_stencils(
      comm, wsum_stencils_data, stencil_owner, reorder_idx, count);

  free(reorder_idx);
  free(rank_buffer);

  return new_wsum_stencils_data;
}

static int compare_remote_point_info(const void * a, const void * b) {

  int ret = ((struct remote_point_info *)a)->rank -
            ((struct remote_point_info *)b)->rank;

  if (ret) return ret;

  return (((struct remote_point_info *)a)->orig_pos >
          ((struct remote_point_info *)b)->orig_pos) -
         (((struct remote_point_info *)a)->orig_pos <
          ((struct remote_point_info *)b)->orig_pos);
}

static struct interp_weight_stencils_wsum_mf * redist_wsum_mf_stencils_tgt(
  MPI_Comm comm, struct interp_weight_stencils_wsum_mf * wsum_stencils_data) {

  struct interp_weight_stencil_wsum_mf * wsum_stencils =
    wsum_stencils_data->data;
  size_t count = wsum_stencils_data->count;

  // determine total number of stencils to be sent to other processes
  // (a single stencil may be sent to multiple target processes)
  size_t total_owner_count = 0;
  for (size_t i = 0; i < count; ++i) {
    int stencil_size = wsum_stencils[i].tgt.data.count;
    if (stencil_size == 1) {
      total_owner_count++;
    } else {
      struct remote_point_info * tgt_point_infos =
        wsum_stencils[i].tgt.data.data.multi;
      qsort(
        tgt_point_infos, stencil_size, sizeof(*tgt_point_infos),
        compare_remote_point_info);
      int prev_rank = INT_MAX;
      for (int j = 0; j < stencil_size; ++j) {
        int curr_rank = tgt_point_infos[j].rank;
        if (curr_rank != prev_rank) {
          ++total_owner_count;
          prev_rank = curr_rank;
        }
      }
    }
  }

  int * stencil_owner = xmalloc(total_owner_count * sizeof(*stencil_owner));
  size_t * reorder_idx = xmalloc(total_owner_count * sizeof(*reorder_idx));
  for (size_t i = 0, k = 0; i < count; ++i) {
    int stencil_size = wsum_stencils[i].tgt.data.count;
    if (stencil_size == 1) {
      stencil_owner[k] = wsum_stencils[i].tgt.data.data.single.rank;
      reorder_idx[k] = i;
      ++k;
    } else {
      struct remote_point_info * tgt_point_infos =
        wsum_stencils[i].tgt.data.data.multi;
      int prev_rank = INT_MAX;
      for (int j = 0; j < stencil_size; ++j) {
        int curr_rank = tgt_point_infos[j].rank;
        if (curr_rank != prev_rank) {
          stencil_owner[k] = tgt_point_infos[j].rank;
          reorder_idx[k] = i;
          ++k;
          prev_rank = curr_rank;
        }
      }
    }
  }

  struct interp_weight_stencils_wsum_mf * new_wsum_stencils_data =
    redist_wsum_mf_stencils(
      comm, wsum_stencils_data, stencil_owner, reorder_idx, total_owner_count);

  wsum_stencils = new_wsum_stencils_data->data;
  count = new_wsum_stencils_data->count;

  free(reorder_idx);
  free(stencil_owner);

  if (count == 0) return new_wsum_stencils_data;

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  // count total number of local target locations
  size_t total_num_tgt_pos = 0;
  for (size_t i = 0; i < count; ++i) {
    size_t curr_count = wsum_stencils[i].tgt.data.count;
    if (curr_count == 1) {
      ++total_num_tgt_pos;
    } else {
      struct remote_point_info * curr_point_infos =
        wsum_stencils[i].tgt.data.data.multi;
      for (size_t j = 0; j < curr_count; ++j)
        if (curr_point_infos[j].rank == comm_rank)
          ++total_num_tgt_pos;
    }
  }

  if (total_num_tgt_pos != count) {
    new_wsum_stencils_data->data =
      ((wsum_stencils =
          xrealloc(wsum_stencils, total_num_tgt_pos * sizeof(*wsum_stencils))));
    new_wsum_stencils_data->count = total_num_tgt_pos;
  }

  // remove all non local target point information
  for (size_t i = 0, offset = count; i < count; ++i) {
    size_t curr_count = wsum_stencils[i].tgt.data.count;
    if (curr_count > 1) {
      struct remote_point_info * curr_point_infos =
        wsum_stencils[i].tgt.data.data.multi;
      // find first local target point
      size_t j;
      for (j = 0; j < curr_count; ++j) {
        if (curr_point_infos[j].rank == comm_rank) {
          wsum_stencils[i].tgt.data.count = 1;
          wsum_stencils[i].tgt.data.data.single.rank = comm_rank;
          wsum_stencils[i].tgt.data.data.single.orig_pos =
            curr_point_infos[j].orig_pos;
          break;
        }
      }
      // make a copy for the remaining local target positions
      for (j = j + 1; j < curr_count; ++j) {
        if (curr_point_infos[j].rank == comm_rank) {
          wsum_stencils[offset] = wsum_stencils[i];
          wsum_stencils[offset].tgt.data.data.single.orig_pos =
            curr_point_infos[j].orig_pos;
          ++offset;
        }
      }
      free(curr_point_infos);
    }
  }

  return new_wsum_stencils_data;
}

static Xt_redist * generate_halo_redists(
  struct remote_point_info_reorder * halo_points, size_t count,
  size_t num_src_fields, MPI_Comm comm) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * sendcounts, * recvcounts, * sdispls, *rdispls;
  yac_get_comm_buffers(
    num_src_fields, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);
  size_t * size_t_buffer =
    xmalloc(4 * (size_t)comm_size * sizeof(*size_t_buffer));
  size_t * total_sendcounts = size_t_buffer + 0 * comm_size;
  size_t * total_recvcounts = size_t_buffer + 1 * comm_size;
  size_t * total_sdispls =    size_t_buffer + 2 * comm_size;
  size_t * total_rdispls =    size_t_buffer + 3 * comm_size;

  for (size_t i = 0; i < count; ++i)
    sendcounts[halo_points[i].data.rank * num_src_fields +
               halo_points[i].field_idx]++;

  yac_generate_alltoallv_args(
    num_src_fields, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t saccu = 0, raccu = 0;
  for (int i = 0; i < comm_size; ++i) {
    total_sdispls[i] = saccu;
    total_rdispls[i] = raccu;
    total_sendcounts[i] = 0;
    total_recvcounts[i] = 0;
    for (size_t j = 0; j < num_src_fields; ++j) {
      total_sendcounts[i] += sendcounts[num_src_fields * i + j];
      total_recvcounts[i] += recvcounts[num_src_fields * i + j];
    }
    saccu += total_sendcounts[i];
    raccu += total_recvcounts[i];
  }

  size_t recv_count = total_recvcounts[comm_size - 1] +
                      total_rdispls[comm_size - 1];

  int * exchange_buffer =
    xmalloc((2 * count + recv_count) * sizeof(*exchange_buffer));
  int * send_buffer = exchange_buffer;
  int * reorder_idx = exchange_buffer + count;
  int * recv_buffer = exchange_buffer + 2 * count;

  // pack the original positions of the requested points
  size_t num_halo_per_src_field[num_src_fields];
  memset(
    num_halo_per_src_field, 0,
    num_src_fields * sizeof(num_halo_per_src_field[0]));
  for (size_t i = 0; i < count; ++i) {
    size_t curr_src_field_idx = (size_t)(halo_points[i].field_idx);
    size_t pos = sdispls[(size_t)(halo_points[i].data.rank) * num_src_fields +
                         curr_src_field_idx + 1]++;
    uint64_t orig_pos = halo_points[i].data.orig_pos;
    YAC_ASSERT(
      orig_pos <= INT_MAX,
      "ERROR(generate_halo_redists): offset not supported by MPI")
    send_buffer[pos] = (int)orig_pos;
    reorder_idx[pos] = num_halo_per_src_field[curr_src_field_idx]++;
  }

  // exchange original positions of the requested points
  yac_alltoallv_int_p2p(
    send_buffer, total_sendcounts, total_sdispls,
    recv_buffer, total_recvcounts, total_rdispls, comm);

  free(size_t_buffer);

  size_t nsend = 0, nsends[num_src_fields];
  size_t nrecv = 0, nrecvs[num_src_fields];
  memset(nsends, 0, num_src_fields * sizeof(nsends[0]));
  memset(nrecvs, 0, num_src_fields * sizeof(nrecvs[0]));
  for (int i = 0; i < comm_size; ++i) {
    for (size_t field_idx = 0; field_idx < num_src_fields; ++field_idx) {
      if (sendcounts[i * num_src_fields + field_idx] > 0) {
        nrecv++;
        nrecvs[field_idx]++;
      }
      if (recvcounts[i * num_src_fields + field_idx] > 0) {
        nsend++;
        nsends[field_idx]++;
      }
    }
  }

  size_t total_num_msg = nsend + nrecv;

  struct Xt_redist_msg * msgs_buffer =
    xmalloc(total_num_msg * sizeof(*msgs_buffer));
  struct Xt_redist_msg * send_msgs = msgs_buffer;
  struct Xt_redist_msg * recv_msgs = msgs_buffer + nsend;

  for (size_t field_idx = 0, nsend = 0, nrecv = 0;
       field_idx < num_src_fields; ++field_idx) {
    for (int rank = 0; rank < comm_size; ++rank) {
      size_t idx = (size_t)rank * num_src_fields + field_idx;
      if (sendcounts[idx] > 0) {
        recv_msgs[nrecv].rank = rank;
        recv_msgs[nrecv].datatype =
          xt_mpi_generate_datatype(
            reorder_idx + sdispls[idx], sendcounts[idx], MPI_DOUBLE, comm);
        nrecv++;
      }
      if (recvcounts[idx] > 0) {
        send_msgs[nsend].rank = rank;
        send_msgs[nsend].datatype =
          xt_mpi_generate_datatype(
            recv_buffer + rdispls[idx], recvcounts[idx], MPI_DOUBLE, comm);
        nsend++;
      }
    }
  }

  Xt_redist * redist;
  MPI_Comm halo_comm;

  if (total_num_msg > 0) {

    yac_mpi_call(MPI_Comm_split(comm, 1, 0, &halo_comm), comm);

    int * rank_buffer = xmalloc(2 * total_num_msg * sizeof(*rank_buffer));
    int * orig_ranks = rank_buffer;
    int * split_ranks = rank_buffer + total_num_msg;

    for (size_t i = 0; i < total_num_msg; ++i)
      orig_ranks[i] = msgs_buffer[i].rank;

    MPI_Group orig_group, split_group;
    yac_mpi_call(MPI_Comm_group(comm, &orig_group), comm);
    yac_mpi_call(MPI_Comm_group(halo_comm, &split_group), comm);

    yac_mpi_call(
      MPI_Group_translate_ranks(orig_group, (int)total_num_msg, orig_ranks,
                                split_group, split_ranks), halo_comm);

    for (size_t i = 0; i < total_num_msg; ++i)
      msgs_buffer[i].rank = split_ranks[i];

    free(rank_buffer);

    yac_mpi_call(MPI_Group_free(&split_group), comm);
    yac_mpi_call(MPI_Group_free(&orig_group), comm);

    // generate redist
    redist = xmalloc(num_src_fields * sizeof(*redist));
    if (num_src_fields == 1) {
      *redist =
        xt_redist_single_array_base_new(
          nsend, nrecv, send_msgs, recv_msgs, halo_comm);
    } else {
      for (size_t field_idx = 0; field_idx < num_src_fields; ++field_idx) {
        redist[field_idx] =
          xt_redist_single_array_base_new(
            nsends[field_idx], nrecvs[field_idx],
            send_msgs, recv_msgs, halo_comm);
        send_msgs += nsends[field_idx];
        recv_msgs += nrecvs[field_idx];
      }
    }

  } else {
    yac_mpi_call(MPI_Comm_split(comm, 0, 0, &halo_comm), comm);
    redist = NULL;
  }

  yac_mpi_call(MPI_Comm_free(&halo_comm), comm);
  free(exchange_buffer);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  xt_redist_msg_free(msgs_buffer, total_num_msg, comm);

  return redist;
}

static int compare_rank_pos_reorder_field_idx(const void * a, const void * b) {

  int ret = (((struct remote_point_info_reorder *)a)->field_idx >
             ((struct remote_point_info_reorder *)b)->field_idx) -
            (((struct remote_point_info_reorder *)a)->field_idx <
             ((struct remote_point_info_reorder *)b)->field_idx);

  if (ret) return ret;

  ret = ((struct remote_point_info_reorder *)a)->data.rank -
        ((struct remote_point_info_reorder *)b)->data.rank;

  if (ret) return ret;

  return (((struct remote_point_info_reorder *)a)->data.orig_pos >
          ((struct remote_point_info_reorder *)b)->data.orig_pos) -
         (((struct remote_point_info_reorder *)a)->data.orig_pos <
          ((struct remote_point_info_reorder *)b)->data.orig_pos);
}

static int compare_interp_weight_stencil_wsum_mf_src_orig_pos(
  const void * a, const void * b) {

  struct interp_weight_stencil_wsum_mf * a_ =
    (struct interp_weight_stencil_wsum_mf *)a;
  struct interp_weight_stencil_wsum_mf * b_ =
    (struct interp_weight_stencil_wsum_mf *)b;

  size_t count = MIN(a_->count, b_->count);

  for (size_t i = 0; i < count; ++i) {
    int ret = (a_->data[i].src_field_idx > b_->data[i].src_field_idx) -
              (a_->data[i].src_field_idx < b_->data[i].src_field_idx);
    if (ret) return ret;
    ret = (a_->data[i].src.orig_pos > b_->data[i].src.orig_pos) -
          (a_->data[i].src.orig_pos < b_->data[i].src.orig_pos);
    if (ret) return ret;
  }
  return 0;
}

static int compare_interp_weight_stencil_wsum_mf_tgt_orig_pos(
  const void * a, const void * b) {

  struct interp_weight_stencil_wsum_mf * a_ =
    (struct interp_weight_stencil_wsum_mf *)a;
  struct interp_weight_stencil_wsum_mf * b_ =
    (struct interp_weight_stencil_wsum_mf *)b;

  YAC_ASSERT(
    (a_->tgt.data.count == 1) && (b_->tgt.data.count == 1),
    "ERROR(compare_interp_weight_stencil_wsum_tgt_orig_pos): invalid data")

  size_t a_orig_pos = a_->tgt.data.data.single.orig_pos;
  size_t b_orig_pos = b_->tgt.data.data.single.orig_pos;

  return (a_orig_pos > b_orig_pos) - (a_orig_pos < b_orig_pos);
}

static void free_remote_point(struct remote_point point) {

  if (point.data.count > 1) free(point.data.data.multi);
}

static Xt_redist generate_redist_put_double(
  struct remote_point_infos * point_infos, size_t count, MPI_Comm comm) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  for (size_t i = 0; i < count; ++i) {
    int curr_count = point_infos[i].count;
    struct remote_point_info * curr_point_infos =
      (curr_count == 1)?
        (&(point_infos[i].data.single)):(point_infos[i].data.multi);
    YAC_ASSERT(
      curr_count >= 1,
      "ERROR(generate_redist_put_double): no owner found for global id")
    for (int j = 0; j < curr_count; ++j)
      sendcounts[curr_point_infos[j].rank]++;
  }

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  size_t send_count =
    sdispls[comm_size] + sendcounts[comm_size - 1];
  size_t recv_count =
    rdispls[comm_size - 1] + recvcounts[comm_size - 1];

  int * exchange_buffer =
    xmalloc((2 * send_count + recv_count) * sizeof(*exchange_buffer));
  int * send_buffer = exchange_buffer;
  int * reorder_idx = exchange_buffer + send_count;
  int * recv_buffer = exchange_buffer + 2 * send_count;

  // pack the original positions of the points that have to updated
  for (size_t i = 0; i < count; ++i) {
    int curr_count = point_infos[i].count;
    struct remote_point_info * curr_point_infos =
      (curr_count == 1)?
        (&(point_infos[i].data.single)):(point_infos[i].data.multi);
    for (int j = 0; j < curr_count; ++j) {
      size_t pos = sdispls[curr_point_infos[j].rank + 1]++;
      uint64_t orig_pos = curr_point_infos[j].orig_pos;
      YAC_ASSERT(
        orig_pos <= INT_MAX,
        "ERROR(generate_redist_put_double): offset not supported by MPI")
      send_buffer[pos] = (int)orig_pos;
      reorder_idx[pos] = i;
    }
  }

  // exchange original positions of the points that have to updated
  yac_alltoallv_int_p2p(
    send_buffer, sendcounts, sdispls, recv_buffer, recvcounts, rdispls, comm);

  size_t nsend = 0;
  size_t nrecv = 0;
  for (int i = 0; i < comm_size; ++i) {
    if (sendcounts[i] > 0) nsend++;
    if (recvcounts[i] > 0) nrecv++;
  }

  struct Xt_redist_msg * send_msgs = xmalloc(nsend * sizeof(*send_msgs));
  struct Xt_redist_msg * recv_msgs = xmalloc(nrecv * sizeof(*send_msgs));

  for (int i = 0, nsend = 0, nrecv = 0; i < comm_size; ++i) {
    if (sendcounts[i] > 0) {
      send_msgs[nsend].rank = i;
      send_msgs[nsend].datatype =
        xt_mpi_generate_datatype(
          reorder_idx + sdispls[i], sendcounts[i], MPI_DOUBLE, comm);
      nsend++;
    }
    if (recvcounts[i] > 0) {
      recv_msgs[nrecv].rank = i;
      recv_msgs[nrecv].datatype =
        xt_mpi_generate_datatype(
          recv_buffer + rdispls[i], recvcounts[i], MPI_DOUBLE, comm);
      nrecv++;
    }
  }

  // generate redist
  Xt_redist redist =
    xt_redist_single_array_base_new(nsend, nrecv, send_msgs, recv_msgs, comm);

  free(exchange_buffer);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  xt_redist_msg_free(recv_msgs, nrecv, comm);
  xt_redist_msg_free(send_msgs, nsend, comm);

  return redist;
}

static void yac_interp_weights_redist_w_sum_mf(
  MPI_Comm comm, struct interp_weight_stencils_wsum_mf * wsum_mf_stencils_data,
  struct yac_interpolation * interp,
  enum yac_interp_weights_reorder_type reorder,
  enum yac_interp_weight_stencil_type stencil_type) {

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  // redistribute stencils to respective owners
  struct interp_weight_stencils_wsum_mf * (*redist_wsum_mf_stencils)(
    MPI_Comm comm, struct interp_weight_stencils_wsum_mf * wsum_mf_stencils_data);
  YAC_ASSERT(
    (reorder == YAC_MAPPING_ON_SRC) || (reorder == YAC_MAPPING_ON_TGT),
    "ERROR(yac_interp_weights_redist_w_sum_mf): invalid reorder type")
  redist_wsum_mf_stencils =
    (reorder == YAC_MAPPING_ON_SRC)?
      redist_wsum_mf_stencils_src:redist_wsum_mf_stencils_tgt;
  struct interp_weight_stencils_wsum_mf * new_wsum_mf_stencils_data =
    redist_wsum_mf_stencils(comm, wsum_mf_stencils_data);

  size_t wsum_mf_count = new_wsum_mf_stencils_data->count;
  struct interp_weight_stencil_wsum_mf * wsum_mf_stencils =
    new_wsum_mf_stencils_data->data;

  // compute the total number of links
  size_t total_num_links = 0, total_num_remote_weights = 0;
  for (size_t i = 0; i < wsum_mf_count; ++i) {
    size_t curr_stencil_size = wsum_mf_stencils[i].count;
    total_num_links += curr_stencil_size;
    for (size_t j = 0; j < curr_stencil_size; ++j)
      if (wsum_mf_stencils[i].data[j].src.rank != comm_rank)
        ++total_num_remote_weights;
  }

  // gather all remote source points and determine number of source fields
  struct remote_point_info_reorder * remote_src_points =
    xmalloc(total_num_remote_weights * sizeof(*remote_src_points));
  size_t num_src_fields = 0;
  for (size_t i = 0, k = 0; i < wsum_mf_count; ++i) {
    size_t curr_stencil_size = wsum_mf_stencils[i].count;
    struct interp_weight_stencil_wsum_mf_weight * curr_weights =
      wsum_mf_stencils[i].data;
    for (size_t j = 0; j < curr_stencil_size; ++j) {
      size_t curr_src_field_idx = curr_weights[j].src_field_idx;
      if (curr_src_field_idx >= num_src_fields)
        num_src_fields = curr_src_field_idx + 1;
      if (curr_weights[j].src.rank != comm_rank) {
        remote_src_points[k].data = curr_weights[j].src;
        remote_src_points[k].field_idx = curr_src_field_idx;
        remote_src_points[k].reorder_idx = i;
        ++k;
      }
    }
  }
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &num_src_fields, 1, YAC_MPI_SIZE_T, MPI_MAX, comm), comm);

  // sort remote points first by field_idx, second by rank, and
  // then by orig_pos
  qsort(remote_src_points, total_num_remote_weights, sizeof(*remote_src_points),
        compare_rank_pos_reorder_field_idx);

  // update stencils: set owner to -1; set orig_pos to position of respecitve
  // point in halo data
  // remove duplicated remote points
  struct remote_point_info * prev_remote_src_point;
  size_t prev_field_idx;
  size_t halo_size;
  if (total_num_remote_weights > 0) {
    prev_remote_src_point = &(remote_src_points[0].data);
    prev_field_idx = remote_src_points[0].field_idx;
    halo_size = 1;
  } else {
    prev_field_idx = SIZE_MAX;
    halo_size = 0;
  }
  for (size_t i = 0; i < total_num_remote_weights; ++i) {
    struct remote_point_info * curr_remote_src_point =
      &(remote_src_points[i].data);
    size_t curr_field_idx = remote_src_points[i].field_idx;
    if (compare_remote_point_info(
          prev_remote_src_point, curr_remote_src_point) ||
        (prev_field_idx != curr_field_idx)) {
      prev_remote_src_point = curr_remote_src_point;
      prev_field_idx = curr_field_idx;
      remote_src_points[halo_size].data = *curr_remote_src_point;
      remote_src_points[halo_size].field_idx = curr_field_idx;
      ++halo_size;
    }
    struct interp_weight_stencil_wsum_mf * curr_stencil =
      wsum_mf_stencils + remote_src_points[i].reorder_idx;
    size_t curr_stencil_size = curr_stencil->count;
    for (size_t j = 0; j < curr_stencil_size; ++j) {
      if ((!compare_remote_point_info(
             &(curr_stencil->data[j].src), curr_remote_src_point)) &&
          (curr_stencil->data[j].src_field_idx == curr_field_idx)) {
        curr_stencil->data[j].src.rank = -1;
        curr_stencil->data[j].src.orig_pos = halo_size - 1;
        curr_stencil->data[j].src_field_idx = SIZE_MAX;
      }
    }
  }

  // sort stencils by their memory access pattern on the local process
  qsort(wsum_mf_stencils, wsum_mf_count, sizeof(*wsum_mf_stencils),
        (reorder == YAC_MAPPING_ON_SRC)?
          compare_interp_weight_stencil_wsum_mf_src_orig_pos:
          compare_interp_weight_stencil_wsum_mf_tgt_orig_pos);

  size_t * num_src_per_tgt = xmalloc(wsum_mf_count * sizeof(*num_src_per_tgt));
  double * weights = xmalloc(total_num_links * sizeof(*weights));
  size_t * src_idx = xmalloc(total_num_links * sizeof(*src_idx));
  size_t * src_field_idx = xmalloc(total_num_links * sizeof(*src_field_idx));

  // extract data from stencil
  for (size_t i = 0, k = 0; i < wsum_mf_count; ++i) {
    size_t curr_stencil_size = wsum_mf_stencils[i].count;
    struct interp_weight_stencil_wsum_mf_weight * curr_weights =
      wsum_mf_stencils[i].data;
    num_src_per_tgt[i] = curr_stencil_size;
    for (size_t j = 0; j < curr_stencil_size; ++j, ++k){
      weights[k] = curr_weights[j].weight;
      src_idx[k] = curr_weights[j].src.orig_pos;
      src_field_idx[k] = curr_weights[j].src_field_idx;
    }
  }

  // generate halo redists (one per source field
  Xt_redist * halo_redists =
    generate_halo_redists(
      remote_src_points, halo_size, num_src_fields, comm);

  YAC_ASSERT(
    (stencil_type == WEIGHT_SUM) ||
    (stencil_type == SUM) ||
    (stencil_type == WEIGHT_SUM_MF) ||
    (stencil_type == SUM_MF),
    "ERROR(yac_interp_weights_redist_w_sum_mf): unsupported stencil type");

  // add weights to interpolation
  if (reorder == YAC_MAPPING_ON_SRC) {

    // generate result redist
    struct remote_point_infos * tgt_infos =
      xmalloc(wsum_mf_count * sizeof(*tgt_infos));
    for (size_t i = 0; i < wsum_mf_count; ++i)
      tgt_infos[i] = wsum_mf_stencils[i].tgt.data;
    Xt_redist result_redist =
      generate_redist_put_double(
        tgt_infos, wsum_mf_count, comm);
    free(tgt_infos);

    switch(stencil_type) {
      default:
      case (WEIGHT_SUM):
      case (WEIGHT_SUM_MF):
        yac_interpolation_add_weight_sum_mvp_at_src(
          interp, halo_redists, wsum_mf_count, num_src_per_tgt, weights,
          src_field_idx, src_idx, ((stencil_type==WEIGHT_SUM)?1:num_src_fields),
          result_redist);
        break;
      case (SUM):
      case (SUM_MF):
        yac_interpolation_add_sum_at_src(
          interp, halo_redists, wsum_mf_count, num_src_per_tgt,
          src_field_idx, src_idx, ((stencil_type == SUM)?1:num_src_fields),
          result_redist);
        break;
    }

    if (result_redist != NULL) xt_redist_delete(result_redist);

  } else {

    size_t * tgt_orig_pos = xmalloc(wsum_mf_count * sizeof(*tgt_orig_pos));
    for (size_t i = 0; i < wsum_mf_count; ++i) {
      YAC_ASSERT(
        wsum_mf_stencils[i].tgt.data.count == 1,
        "ERROR(yac_interp_weights_redist_w_sum): currently unsupported target "
        "point distribution")
      tgt_orig_pos[i] =
        (size_t)(wsum_mf_stencils[i].tgt.data.data.single.orig_pos);
    }

    switch(stencil_type) {
      default:
      case (WEIGHT_SUM):
      case (WEIGHT_SUM_MF):
      yac_interpolation_add_weight_sum_mvp_at_tgt(
        interp, halo_redists, tgt_orig_pos, wsum_mf_count,
        num_src_per_tgt, weights, src_field_idx, src_idx,
        ((stencil_type == WEIGHT_SUM)?1:num_src_fields));
        break;
      case (SUM):
      case (SUM_MF):
        yac_interpolation_add_sum_at_tgt(
          interp, halo_redists, tgt_orig_pos, wsum_mf_count,
          num_src_per_tgt, src_field_idx, src_idx,
          ((stencil_type == SUM)?1:num_src_fields));
        break;
    }
    free(tgt_orig_pos);
  }

  for (size_t i = 0; i < new_wsum_mf_stencils_data->count; ++i)
    free_remote_point(new_wsum_mf_stencils_data->data[i].tgt);
  free(new_wsum_mf_stencils_data->data);
  free(new_wsum_mf_stencils_data);

  free(remote_src_points);
  free(src_field_idx);
  free(src_idx);
  free(weights);
  free(num_src_per_tgt);
  if (halo_redists != NULL) {
    for (size_t i = 0; i < num_src_fields; ++i)
      xt_redist_delete(halo_redists[i]);
    free(halo_redists);
  }
}

static int compare_stencils(const void * a, const void * b) {

  return (int)(((struct interp_weight_stencil *)a)->type) -
         (int)(((struct interp_weight_stencil *)b)->type);
}

struct yac_interpolation * yac_interp_weights_get_interpolation(
  struct yac_interp_weights * weights,
  enum yac_interp_weights_reorder_type reorder,
  size_t collection_size, double frac_mask_fallback_value,
  double scaling_factor, double scaling_summand) {

  MPI_Comm comm = weights->comm;

  // sort stencils by type
  qsort(weights->stencils, weights->stencils_size, sizeof(*(weights->stencils)),
        compare_stencils);

  uint64_t local_stencil_counts[WEIGHT_STENCIL_TYPE_SIZE];
  size_t stencils_offsets[WEIGHT_STENCIL_TYPE_SIZE];

  // count local number of stencils per type
  memset(&(local_stencil_counts[0]), 0, sizeof(local_stencil_counts));
  for (size_t i = 0; i < weights->stencils_size; ++i)
    local_stencil_counts[(int)(weights->stencils[i].type)]++;

  for (size_t i = 0, accu = 0; i < (size_t)WEIGHT_STENCIL_TYPE_SIZE; ++i) {
    stencils_offsets[i] = accu;
    accu += local_stencil_counts[i];
  }

  uint64_t global_stencil_counts[WEIGHT_STENCIL_TYPE_SIZE];

  // determine global number of stencils per type
  yac_mpi_call(
    MPI_Allreduce(
      local_stencil_counts, global_stencil_counts,
      (int)WEIGHT_STENCIL_TYPE_SIZE, MPI_UINT64_T, MPI_SUM, comm), comm);

  { // check whether the collection_size is consistant across all processes
    uint64_t max_collection_size = (uint64_t)collection_size;
    yac_mpi_call(
      MPI_Allreduce(
        MPI_IN_PLACE, &max_collection_size, 1, MPI_UINT64_T, MPI_MAX, comm),
      comm);
    YAC_ASSERT(
      (size_t)max_collection_size == collection_size,
      "ERROR(yac_interp_weights_get_interpolation): "
      "mismatching collection sizes")
  }

  struct yac_interpolation * interp =
    yac_interpolation_new(
      collection_size, frac_mask_fallback_value,
      scaling_factor, scaling_summand);

  if (global_stencil_counts[FIXED] > 0)
    yac_interp_weights_redist_fixed(
      weights->comm, local_stencil_counts[FIXED],
      weights->stencils + stencils_offsets[FIXED], interp);

  if (global_stencil_counts[DIRECT] > 0)
    yac_interp_weights_redist_direct(
      weights->comm, local_stencil_counts[DIRECT],
      weights->stencils + stencils_offsets[DIRECT], interp);

  if (global_stencil_counts[SUM] > 0) {

    struct interp_weight_stencils_wsum_mf * wsum_stencils =
      generate_w_sum_mf_stencils(
        weights->stencils + stencils_offsets[SUM],
        (size_t)(local_stencil_counts[SUM]), SUM);
    YAC_ASSERT(
      (reorder == YAC_MAPPING_ON_SRC) || (reorder == YAC_MAPPING_ON_TGT),
      "ERROR(yac_interp_weights_get_interpolation): invalid reorder type")
    yac_interp_weights_redist_w_sum_mf(
      weights->comm, wsum_stencils, interp, reorder, SUM);
    for (size_t i = 0; i < wsum_stencils->count; ++i)
      free_remote_point(wsum_stencils->data[i].tgt);
    free(wsum_stencils->data);
    free(wsum_stencils);
  }

  if (global_stencil_counts[WEIGHT_SUM] > 0) {

    struct interp_weight_stencils_wsum_mf * wsum_stencils =
      generate_w_sum_mf_stencils(
        weights->stencils + stencils_offsets[WEIGHT_SUM],
        (size_t)(local_stencil_counts[WEIGHT_SUM]), WEIGHT_SUM);
    YAC_ASSERT(
      (reorder == YAC_MAPPING_ON_SRC) || (reorder == YAC_MAPPING_ON_TGT),
      "ERROR(yac_interp_weights_get_interpolation): invalid reorder type")
    yac_interp_weights_redist_w_sum_mf(
      weights->comm, wsum_stencils, interp, reorder, WEIGHT_SUM);
    for (size_t i = 0; i < wsum_stencils->count; ++i)
      free_remote_point(wsum_stencils->data[i].tgt);
    free(wsum_stencils->data);
    free(wsum_stencils);
  }

  if (global_stencil_counts[DIRECT_MF] > 0)
    yac_interp_weights_redist_direct_mf(
      weights->comm, local_stencil_counts[DIRECT_MF],
      weights->stencils + stencils_offsets[DIRECT_MF], interp);

  if (global_stencil_counts[SUM_MF] > 0) {

    struct interp_weight_stencils_wsum_mf * sum_mf_stencils =
      generate_w_sum_mf_stencils(
        weights->stencils + stencils_offsets[SUM_MF],
        (size_t)(local_stencil_counts[SUM_MF]), SUM_MF);
    YAC_ASSERT(
      (reorder == YAC_MAPPING_ON_SRC) || (reorder == YAC_MAPPING_ON_TGT),
      "ERROR(yac_interp_weights_get_interpolation): invalid reorder type")
    yac_interp_weights_redist_w_sum_mf(
      weights->comm, sum_mf_stencils, interp, reorder, SUM_MF);
    for (size_t i = 0; i < sum_mf_stencils->count; ++i)
      free_remote_point(sum_mf_stencils->data[i].tgt);
    free(sum_mf_stencils->data);
    free(sum_mf_stencils);
  }

  if (global_stencil_counts[WEIGHT_SUM_MF] > 0) {

    struct interp_weight_stencils_wsum_mf * wsum_mf_stencils =
      generate_w_sum_mf_stencils(
        weights->stencils + stencils_offsets[WEIGHT_SUM_MF],
        (size_t)(local_stencil_counts[WEIGHT_SUM_MF]), WEIGHT_SUM_MF);
    YAC_ASSERT(
      (reorder == YAC_MAPPING_ON_SRC) || (reorder == YAC_MAPPING_ON_TGT),
      "ERROR(yac_interp_weights_get_interpolation): invalid reorder type")
    yac_interp_weights_redist_w_sum_mf(
      weights->comm, wsum_mf_stencils, interp, reorder, WEIGHT_SUM_MF);
    for (size_t i = 0; i < wsum_mf_stencils->count; ++i)
      free_remote_point(wsum_mf_stencils->data[i].tgt);
    free(wsum_mf_stencils->data);
    free(wsum_mf_stencils);
  }

  return interp;
}

struct yac_interpolation * yac_interp_weights_get_interpolation_f2c(
  struct yac_interp_weights * weights, int reorder,
  size_t collection_size, double frac_mask_fallback_value,
  double scaling_factor, double scaling_summand) {

  YAC_ASSERT(
    (reorder == YAC_MAPPING_ON_SRC) ||
    (reorder == YAC_MAPPING_ON_TGT),
    "ERROR(yac_interp_weights_get_interpolation_f2c): "
    "reorder type must be of YAC_MAPPING_ON_SRC/YAC_MAPPING_ON_TGT");

  return
    yac_interp_weights_get_interpolation(
      weights, (enum yac_interp_weights_reorder_type)reorder,
      collection_size, frac_mask_fallback_value,
      scaling_factor, scaling_summand);
}

static void free_remote_points(struct remote_points * points) {

  free(points->data);
  free(points);
}

static void yac_interp_weight_stencils_delete(
  struct interp_weight_stencil * stencils, size_t count) {

  for (size_t i = 0 ; i < count; ++i) {

    YAC_ASSERT(
      (stencils[i].type == DIRECT) ||
      (stencils[i].type == SUM) ||
      (stencils[i].type == WEIGHT_SUM) ||
      (stencils[i].type == DIRECT_MF) ||
      (stencils[i].type == SUM_MF) ||
      (stencils[i].type == WEIGHT_SUM_MF) ||
      (stencils[i].type == FIXED),
      "ERROR(yac_interp_weights_delete): invalid stencil type")
    switch(stencils[i].type) {

      case(DIRECT):
        free_remote_point(stencils[i].data.direct.src);
        break;
      case(SUM):
        free_remote_points(stencils[i].data.sum.srcs);
        break;
      case(WEIGHT_SUM):
        free_remote_points(stencils[i].data.weight_sum.srcs);
        free(stencils[i].data.weight_sum.weights);
        break;
      case(DIRECT_MF):
        free_remote_point(stencils[i].data.direct_mf.src);
        break;
      case(SUM_MF):
        free_remote_points(stencils[i].data.sum_mf.srcs);
        free(stencils[i].data.sum_mf.field_indices);
        break;
      case (WEIGHT_SUM_MF):
        free_remote_points(stencils[i].data.weight_sum_mf.srcs);
        free(stencils[i].data.weight_sum_mf.weights);
        free(stencils[i].data.weight_sum_mf.field_indices);
        break;
      default:
      case(FIXED):
        break;
    };
    free_remote_point(stencils[i].tgt);
  }
  free(stencils);
}

#ifdef YAC_NETCDF_ENABLED
static int compare_double(void const * a, void const * b) {

  return (*(double const *)a > *(double const *)b) -
         (*(double const *)a < *(double const *)b);
}

/**
 * create a weight file with basic meta data
 */
static void create_weight_file(
  char const * filename, char const * src_grid_name, char const * tgt_grid_name,
  size_t num_fixed_values, double * fixed_values,
  size_t * num_tgt_per_fixed_value, size_t num_links,
  size_t num_weights_per_link, size_t num_src_fields,
  size_t * num_links_per_src_field,
  enum yac_location * src_locations, enum yac_location tgt_location,
  size_t src_grid_size, size_t tgt_grid_size) {

  int ncid;

  // create file
  yac_nc_create(filename, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);

  int dim_weight_id[8];

  // define dimensions
  if (num_links > 0) {
    YAC_HANDLE_ERROR(nc_def_dim(ncid, "num_links", num_links, &dim_weight_id[0]));
    YAC_ASSERT_F(
      num_weights_per_link > 0,
      "ERROR(create_weight_file): number of links is %zu but number of "
      "weights per link is zero for weight file %s", num_links, filename)
    YAC_HANDLE_ERROR(
      nc_def_dim(ncid, "num_wgts", num_weights_per_link, &dim_weight_id[1]));
  }
  YAC_ASSERT_F(
    num_src_fields > 0,
    "ERROR(create_weight_file): number of source fields is zero for "
    "weight file %s", filename)
  YAC_HANDLE_ERROR(
    nc_def_dim(ncid, "num_src_fields", num_src_fields, &dim_weight_id[2]));
  YAC_HANDLE_ERROR(
    nc_def_dim(
      ncid, "max_loc_str_len", YAC_MAX_LOC_STR_LEN, &dim_weight_id[3]));

  if (num_fixed_values > 0) {
    YAC_HANDLE_ERROR(
      nc_def_dim(
        ncid, "num_fixed_values", num_fixed_values, &dim_weight_id[4]));
    size_t num_fixed_dst = 0;
    for (size_t i = 0; i < num_fixed_values; ++i)
      num_fixed_dst += num_tgt_per_fixed_value[i];
    YAC_ASSERT_F(
      num_fixed_dst > 0,
      "ERROR(create_weight_file): number of fixed values is %zu but number "
      "of fixed destination points is zero for weight file %s",
      num_fixed_dst, filename)
    YAC_HANDLE_ERROR(
      nc_def_dim(ncid, "num_fixed_dst", num_fixed_dst, &dim_weight_id[5]));
  }

  if (src_grid_size > 0)
    YAC_HANDLE_ERROR(
      nc_def_dim(ncid, "src_grid_size", src_grid_size, &dim_weight_id[6]));

  if (tgt_grid_size > 0)
    YAC_HANDLE_ERROR(
      nc_def_dim(ncid, "dst_grid_size", tgt_grid_size, &dim_weight_id[7]));

  int var_src_add_id, var_dst_add_id, var_weight_id, var_num_links_id,
      src_var_locs_id, tgt_var_loc_id, var_fixed_values_id,
      var_num_dst_per_fixed_value_id, var_dst_add_fixed_id;

  // define variables
  if (num_links > 0) {
    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, "src_address", NC_INT, 1, dim_weight_id, &var_src_add_id));
    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, "dst_address", NC_INT, 1, dim_weight_id, &var_dst_add_id));
    YAC_HANDLE_ERROR(
      nc_def_var(
        ncid, "remap_matrix", NC_DOUBLE, 2, dim_weight_id, &var_weight_id));
    YAC_HANDLE_ERROR(
      nc_def_var(ncid, "num_links_per_src_field", NC_INT, 1,
                 &dim_weight_id[2], &var_num_links_id));
  }
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "src_locations", NC_CHAR, 2, &dim_weight_id[2], &src_var_locs_id));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, "dst_location", NC_CHAR, 1, &dim_weight_id[3], &tgt_var_loc_id));
  if (num_fixed_values > 0) {
    YAC_HANDLE_ERROR(
      nc_def_var(ncid, "fixed_values", NC_DOUBLE, 1, &dim_weight_id[4],
                 &var_fixed_values_id));
    YAC_HANDLE_ERROR(
      nc_def_var(ncid, "num_dst_per_fixed_value", NC_INT, 1, &dim_weight_id[4],
                 &var_num_dst_per_fixed_value_id));
    YAC_HANDLE_ERROR(
      nc_def_var(ncid, "dst_address_fixed", NC_INT, 1, &dim_weight_id[5],
                 &var_dst_add_fixed_id));
  }

  // put attributes
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, NC_GLOBAL, "version",
                    strlen(YAC_WEIGHT_FILE_VERSION_STRING),
                    YAC_WEIGHT_FILE_VERSION_STRING));
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, NC_GLOBAL, "src_grid_name",
                    strlen(src_grid_name), src_grid_name));
  YAC_HANDLE_ERROR(
    nc_put_att_text(ncid, NC_GLOBAL, "dst_grid_name",
                    strlen(tgt_grid_name), tgt_grid_name));
  {
    char const * str_logical[2] = {"FALSE", "TRUE"};
    YAC_HANDLE_ERROR(nc_put_att_text(ncid, NC_GLOBAL, "contains_links",
                                 strlen(str_logical[num_links > 0]),
                                 str_logical[num_links > 0]));
    YAC_HANDLE_ERROR(nc_put_att_text(ncid, NC_GLOBAL, "contains_fixed_dst",
                                 strlen(str_logical[num_fixed_values > 0]),
                                 str_logical[num_fixed_values > 0]));
  }

  // end definition
  YAC_HANDLE_ERROR(nc_enddef(ncid));

  // write some basic data

  if (num_links > 0) {
    int * num_links_per_src_field_int =
      xmalloc(num_src_fields * sizeof(*num_links_per_src_field_int));
    for (size_t i = 0; i < num_src_fields; ++i) {
      YAC_ASSERT(
        num_links_per_src_field[i] <= INT_MAX,
        "ERROR(create_weight_file): "
        "number of links per source field too big (not yet supported)")
      num_links_per_src_field_int[i] = (int)num_links_per_src_field[i];
    }
    YAC_HANDLE_ERROR(
      nc_put_var_int(ncid, var_num_links_id, num_links_per_src_field_int));
    free(num_links_per_src_field_int);
  }

  for (size_t i = 0; i < num_src_fields; ++i) {
    char const * loc_str = yac_loc2str(src_locations[i]);
    size_t str_start[2] = {i, 0};
    size_t str_count[2] = {1, strlen(loc_str)};
    YAC_HANDLE_ERROR(
      nc_put_vara_text(ncid, src_var_locs_id, str_start, str_count, loc_str));
  }

  {
    char const * loc_str = yac_loc2str(tgt_location);
    size_t str_start[1] = {0};
    size_t str_count[1] = {strlen(loc_str)};
    YAC_HANDLE_ERROR(
      nc_put_vara_text(ncid, tgt_var_loc_id, str_start, str_count, loc_str));
  }
  if (num_fixed_values > 0) {

    int * num_tgt_per_fixed_value_int =
      xmalloc(num_fixed_values * sizeof(*num_tgt_per_fixed_value_int));
    for (unsigned i = 0; i < num_fixed_values; ++i) {
      YAC_ASSERT(
        num_tgt_per_fixed_value[i] <= INT_MAX,
        "ERROR(create_weight_file): "
        "number of targets per fixed value is too big (not yet supported)")
      num_tgt_per_fixed_value_int[i] = (int)num_tgt_per_fixed_value[i];
    }
    YAC_HANDLE_ERROR(nc_put_var_double(ncid, var_fixed_values_id, fixed_values));
    YAC_HANDLE_ERROR(nc_put_var_int(ncid, var_num_dst_per_fixed_value_id,
                                num_tgt_per_fixed_value_int));
    free(num_tgt_per_fixed_value_int);
  }

  // close file
  YAC_HANDLE_ERROR(nc_close(ncid));
}

static int compare_interp_weight_stencil(const void * a, const void * b) {

  int a_is_fixed = (((struct interp_weight_stencil *)a)->type == FIXED);
  int b_is_fixed = (((struct interp_weight_stencil *)b)->type == FIXED);
  int ret = b_is_fixed - a_is_fixed;

  if (ret) return ret;

  // if both are fixed stencils
  if (a_is_fixed) {

    double fixed_value_a =
      ((struct interp_weight_stencil *)a)->data.fixed.value;
    double fixed_value_b =
      ((struct interp_weight_stencil *)b)->data.fixed.value;
    ret = (fixed_value_a > fixed_value_b) -
          (fixed_value_a < fixed_value_b);
    if (ret) return ret;
  }

  return (((struct interp_weight_stencil *)a)->tgt.global_id >
          ((struct interp_weight_stencil *)b)->tgt.global_id) -
         (((struct interp_weight_stencil *)a)->tgt.global_id <
          ((struct interp_weight_stencil *)b)->tgt.global_id);
}

static void stencil_determine_tgt_global_id_range(
  struct interp_weight_stencil * stencils, size_t stencils_size,
  yac_int * min_tgt_global_id, yac_int * max_tgt_global_id, MPI_Comm comm) {

  yac_int min_max[2] = {XT_INT_MAX, XT_INT_MIN};

  for (size_t i = 0; i < stencils_size; ++i) {

    yac_int curr_id = stencils[i].tgt.global_id;
    if (curr_id < min_max[0]) min_max[0] = curr_id;
    if (curr_id > min_max[1]) min_max[1] = curr_id;
  }

  min_max[0] = XT_INT_MAX - min_max[0];

  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, min_max, 2, yac_int_dt, MPI_MAX, comm), comm);

  *min_tgt_global_id = XT_INT_MAX - min_max[0];
  *max_tgt_global_id = min_max[1];
}

static void determine_stencils_io_owner(
  struct interp_weight_stencil * stencils, size_t stencils_size,
  yac_int min_tgt_global_id, yac_int max_tgt_global_id,
  int num_io_procs_int, int * io_owner) {

  long long num_io_procs = (long long)num_io_procs_int;
  long long id_range =
    MAX((long long)(max_tgt_global_id - min_tgt_global_id),1);

  for (size_t i = 0; i < stencils_size; ++i)
    io_owner[i] =
      ((int)(MIN(((long long)(stencils[i].tgt.global_id - min_tgt_global_id) *
                 num_io_procs) / id_range, num_io_procs - 1)));
}

static void stencil_get_fixed_values(
  struct interp_weight_stencil * stencils, size_t stencil_count,
  double ** fixed_values, size_t * num_fixed_values, MPI_Comm comm) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  double * local_fixed_values =
    xmalloc(stencil_count * sizeof(*local_fixed_values));

  int * int_buffer = xmalloc(2 * (size_t)comm_size * sizeof(*int_buffer));
  int * recvcounts = int_buffer + 0 * comm_size;
  int * rdispls = int_buffer + 1 * comm_size;

  size_t local_num_fixed = 0;

  // get all local fixed values
  for (size_t i = 0; i < stencil_count;
       ++i, ++local_num_fixed) {
    if (stencils[i].type != FIXED) break;
    local_fixed_values[i] = stencils[i].data.fixed.value;
  }
  qsort(local_fixed_values, local_num_fixed, sizeof(*local_fixed_values),
        compare_double);
  yac_remove_duplicates_double(local_fixed_values, &local_num_fixed);

  // get number of fixed values per rank
  int local_num_fixed_int = (int)(local_num_fixed);
  yac_mpi_call(
    MPI_Allgather(
      &local_num_fixed_int, 1, MPI_INT, recvcounts, 1,MPI_INT, comm), comm);
  for (int i = 0, accu = 0; i < comm_size; ++i) {
    rdispls[i] = accu;
    accu += recvcounts[i];
  }

  size_t num_all_fixed_values = 0;
  for (int i = 0; i < comm_size; ++i)
    num_all_fixed_values += (size_t)(recvcounts[i]);

  double * all_fixed_values =
    xmalloc(num_all_fixed_values * sizeof(*all_fixed_values));

  // gather all fixed values
  yac_mpi_call(
    MPI_Allgatherv(
      local_fixed_values, local_num_fixed_int, MPI_DOUBLE,
      all_fixed_values, recvcounts, rdispls, MPI_DOUBLE, comm), comm);
  free(int_buffer);
  free(local_fixed_values);

  qsort(all_fixed_values, num_all_fixed_values, sizeof(*all_fixed_values),
        compare_double);
  yac_remove_duplicates_double(all_fixed_values, &num_all_fixed_values);
  *fixed_values = xrealloc(all_fixed_values,
                           num_all_fixed_values * sizeof(*all_fixed_values));
  *num_fixed_values = num_all_fixed_values;
}

static size_t get_num_weights_per_link(struct interp_weight_stencil * stencil) {

  YAC_ASSERT(
    (stencil->type == FIXED) ||
    (stencil->type == DIRECT) ||
    (stencil->type == SUM) ||
    (stencil->type == WEIGHT_SUM) ||
    (stencil->type == DIRECT_MF) ||
    (stencil->type == SUM_MF) ||
    (stencil->type == WEIGHT_SUM_MF),
    "ERROR(get_num_weights_per_link): invalid stencil type")

  return (stencil->type == FIXED)?0:1;
}

static size_t stencil_get_num_weights_per_tgt(
  struct interp_weight_stencil * stencils, size_t stencil_count,
  MPI_Comm comm) {

  size_t num_weights_per_link = 0;
  for (size_t i = 0; i < stencil_count; ++i)
    num_weights_per_link =
      MAX(num_weights_per_link, get_num_weights_per_link(stencils + i));

  uint64_t num_weights_per_link_64_t = (uint64_t)num_weights_per_link;
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &num_weights_per_link_64_t, 1, MPI_UINT64_T,
      MPI_MAX, comm), comm);
  num_weights_per_link = (size_t)num_weights_per_link_64_t;

  return num_weights_per_link;
}

static size_t get_num_links_per_src_field(
  struct interp_weight_stencil * stencil, size_t src_field_idx) {;

  YAC_ASSERT(
    stencil->type != FIXED,
    "ERROR(get_num_links_per_src_field): "
    "stencil type FIXED not supported by this routine")
  YAC_ASSERT(
    (stencil->type == DIRECT) ||
    (stencil->type == SUM) ||
    (stencil->type == WEIGHT_SUM) ||
    (stencil->type == DIRECT_MF) ||
    (stencil->type == SUM_MF) ||
    (stencil->type == WEIGHT_SUM_MF),
    "ERROR(get_num_links_per_src_field): invalid stencil type")
  switch (stencil->type) {
    default:
    case(DIRECT): return (src_field_idx == 0)?1:0;
    case(SUM): return (src_field_idx == 0)?stencil->data.sum.srcs->count:0;
    case(WEIGHT_SUM):
      return (src_field_idx == 0)?stencil->data.weight_sum.srcs->count:0;
    case(DIRECT_MF): return stencil->data.direct_mf.field_idx == src_field_idx;
    case(SUM_MF): {
      size_t count = 0;
      size_t stencil_size = stencil->data.sum_mf.srcs->count;
      size_t * field_indices = stencil->data.sum_mf.field_indices;
      for (size_t i = 0; i < stencil_size; ++i)
        if (field_indices[i] == src_field_idx) ++count;
      return count;
    }
    case(WEIGHT_SUM_MF): {
      size_t count = 0;
      size_t stencil_size = stencil->data.weight_sum_mf.srcs->count;
      size_t * field_indices = stencil->data.weight_sum_mf.field_indices;
      for (size_t i = 0; i < stencil_size; ++i)
        if (field_indices[i] == src_field_idx) ++count;
      return count;
    }
  };
}

static void stencil_get_counts(
  struct interp_weight_stencil * stencils, size_t stencil_count,
  size_t num_fixed_values, double * fixed_values,
  size_t * num_tgt_per_fixed_value,
  size_t * num_fixed_tgt, size_t num_src_fields,
  size_t * num_links_per_src_field, size_t * num_links) {

  *num_fixed_tgt = 0;
  *num_links = 0;
  for (size_t i = 0; i < num_fixed_values; ++i) num_tgt_per_fixed_value[i] = 0;
  for (size_t i = 0; i < num_src_fields; ++i) num_links_per_src_field[i] = 0;

  for (size_t i = 0; i < stencil_count; ++i) {
    if (stencils[i].type == FIXED) {
      double curr_fixed_value = stencils[i].data.fixed.value;
      for (size_t j = 0; j < num_fixed_values; ++j) {
        if (curr_fixed_value == fixed_values[j]) {
          num_tgt_per_fixed_value[j]++;
          break;
        }
      }
      ++*num_fixed_tgt;
    } else {
      for (size_t j = 0; j < num_src_fields; ++j) {
        num_links_per_src_field[j] +=
          get_num_links_per_src_field(stencils + i, j);
      }
    }
  }
  for (size_t i = 0; i < num_src_fields; ++i)
    *num_links += num_links_per_src_field[i];
}

static void stencil_xscan_offsets(
  size_t num_fixed_values, size_t * num_tgt_per_fixed_value,
  size_t num_src_fields, size_t * num_links_per_src_field,
  size_t * fixed_offsets, size_t * link_offsets, MPI_Comm comm) {

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  size_t count = num_fixed_values + num_src_fields;
  uint64_t * uint64_t_buffer = xmalloc(3 * count * sizeof(*uint64_t_buffer));
  uint64_t * global_counts = uint64_t_buffer + 0 * count;
  uint64_t * local_counts = uint64_t_buffer + 1 * count;
  uint64_t * offsets = uint64_t_buffer + 2 * count;

  for (size_t i = 0; i < num_fixed_values; ++i)
    local_counts[i] = (uint64_t)(num_tgt_per_fixed_value[i]);
  for (size_t i = 0; i < num_src_fields; ++i)
    local_counts[num_fixed_values + i] = (uint64_t)(num_links_per_src_field[i]);

  yac_mpi_call(
    MPI_Allreduce(local_counts, global_counts, (int)count, MPI_UINT64_T,
                  MPI_SUM, comm), comm);
  yac_mpi_call(
    MPI_Exscan(local_counts, offsets, (int)count, MPI_UINT64_T, MPI_SUM, comm),
    comm);
  if (comm_rank == 0) memset(offsets, 0, count * sizeof(*offsets));

  for (size_t i = 0, accu = 0; i < num_fixed_values; ++i) {
    fixed_offsets[i] = (size_t)(offsets[i]) + accu;
    accu += (size_t)(global_counts[i]);
  }
  for (size_t i = 0, accu = 0; i < num_src_fields; ++i) {
    link_offsets[i] = (size_t)(offsets[i+num_fixed_values]) + accu;
    accu += (size_t)(global_counts[i+num_fixed_values]);
  }
  free(uint64_t_buffer);
}

static int global_id_to_address(yac_int global_id) {

  YAC_ASSERT(
    (global_id < INT_MAX) && (global_id != XT_INT_MAX),
    "ERROR(global_id_to_address): "
    "a global id cannot be converted into a address; too big")
  return (int)global_id + 1;
}

static void stencil_get_tgt_address(
  struct interp_weight_stencil * stencils, size_t stencil_count,
  int * tgt_address) {

  for (size_t i = 0; i < stencil_count; ++i)
    tgt_address[i] = global_id_to_address(stencils[i].tgt.global_id);
}

static void stencil_get_link_data(
  struct interp_weight_stencil * stencils, size_t stencil_count,
  size_t * num_links_per_src_field, size_t num_src_fields,
  int * src_address, int * tgt_address, double * weight) {

  size_t * src_field_offsets =
    xmalloc(2 * num_src_fields * sizeof(*src_field_offsets));
  size_t * prev_src_field_offsets = src_field_offsets + num_src_fields;
  for (size_t i = 0, accu = 0; i < num_src_fields; ++i) {
    src_field_offsets[i] = accu;
    accu += num_links_per_src_field[i];
  }

  struct interp_weight_stencil * curr_stencil = stencils;
  for (size_t i = 0; i < stencil_count; ++i, ++curr_stencil) {

    memcpy(prev_src_field_offsets, src_field_offsets,
           num_src_fields * sizeof(*prev_src_field_offsets));

    int curr_tgt_address =  global_id_to_address(curr_stencil->tgt.global_id);
    YAC_ASSERT(
      curr_stencil->type != FIXED,
      "ERROR(stencil_get_link_data): this call is invalid for FIXED stencils")
    YAC_ASSERT(
      (curr_stencil->type == DIRECT) ||
      (curr_stencil->type == SUM) ||
      (curr_stencil->type == WEIGHT_SUM) ||
      (curr_stencil->type == DIRECT_MF) ||
      (curr_stencil->type == SUM_MF) ||
      (curr_stencil->type == WEIGHT_SUM_MF),
      "ERROR(stencil_get_link_data): invalid stencil type")
    size_t src_field_offset;
    switch (curr_stencil->type) {
      default:
      case(DIRECT):
        src_field_offset = src_field_offsets[0]++;
        src_address[src_field_offset] =
          global_id_to_address(curr_stencil->data.direct.src.global_id);
        tgt_address[src_field_offset] = curr_tgt_address;
        weight[src_field_offset] = 1.0;
        break;
      case(SUM): {
        size_t curr_count = curr_stencil->data.sum.srcs->count;
        struct remote_point * srcs = curr_stencil->data.sum.srcs->data;
        for (size_t k = 0; k < curr_count; ++k) {
          src_field_offset = src_field_offsets[0]++;
          src_address[src_field_offset] =
            global_id_to_address(srcs[k].global_id);
          tgt_address[src_field_offset] = curr_tgt_address;
          weight[src_field_offset] = 1.0;
        }
        break;
      }
      case(WEIGHT_SUM): {
        size_t curr_count = curr_stencil->data.weight_sum.srcs->count;
        struct remote_point * srcs = curr_stencil->data.weight_sum.srcs->data;
        double * weights = curr_stencil->data.weight_sum.weights;
        for (size_t k = 0; k < curr_count; ++k) {
          src_field_offset = src_field_offsets[0]++;
          src_address[src_field_offset] =
            global_id_to_address(srcs[k].global_id);
          tgt_address[src_field_offset] = curr_tgt_address;
          weight[src_field_offset] = weights[k];
        }
        break;
      }
      case(DIRECT_MF):
        src_field_offset =
          src_field_offsets[curr_stencil->data.direct_mf.field_idx]++;
        src_address[src_field_offset ] =
          global_id_to_address(curr_stencil->data.direct_mf.src.global_id);
        tgt_address[src_field_offset ] = curr_tgt_address;
        weight[src_field_offset ] = 1.0;
        break;
      case(SUM_MF): {
        size_t curr_count = curr_stencil->data.sum_mf.srcs->count;
        struct remote_point * srcs =
          curr_stencil->data.sum_mf.srcs->data;
        size_t * field_indices = curr_stencil->data.sum_mf.field_indices;
        for (size_t k = 0; k < curr_count; ++k) {
          src_field_offset = src_field_offsets[field_indices[k]]++;
          src_address[src_field_offset] =
            global_id_to_address(srcs[k].global_id);
          tgt_address[src_field_offset] = curr_tgt_address;
          weight[src_field_offset] = 1.0;
        }
        break;
      }
      case(WEIGHT_SUM_MF): {
        size_t curr_count = curr_stencil->data.weight_sum_mf.srcs->count;
        struct remote_point * srcs =
          curr_stencil->data.weight_sum_mf.srcs->data;
        double * weights = curr_stencil->data.weight_sum_mf.weights;
        size_t * field_indices = curr_stencil->data.weight_sum_mf.field_indices;
        for (size_t k = 0; k < curr_count; ++k) {
          src_field_offset = src_field_offsets[field_indices[k]]++;
          src_address[src_field_offset] =
            global_id_to_address(srcs[k].global_id);
          tgt_address[src_field_offset] = curr_tgt_address;
          weight[src_field_offset] = weights[k];
        }
        break;
      }
    };

    for (size_t j = 0; j < num_src_fields; ++j)
      yac_quicksort_index_int_double(
        src_address + prev_src_field_offsets[j],
        src_field_offsets[j] - prev_src_field_offsets[j],
        weight + prev_src_field_offsets[j]);
  }
  free(src_field_offsets);
}

static void yac_interp_weights_redist_stencils(
  MPI_Comm comm, size_t count, struct interp_weight_stencil * stencils,
  int * owner_ranks, size_t * new_count,
  struct interp_weight_stencil ** new_stencils) {

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * sendcounts, * recvcounts, * sdispls, *rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  size_t * stencil_indices = xmalloc(count * sizeof(*stencil_indices));
  for (size_t i = 0; i < count; ++i) {
    stencil_indices[i] = i;
    sendcounts[owner_ranks[i]]++;
  }

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  // sort the stencil indices by owner rank
  yac_quicksort_index_int_size_t(owner_ranks, count, stencil_indices);

  *new_count = recvcounts[comm_size - 1] + rdispls[comm_size - 1];
  *new_stencils =
    exchange_stencils(comm, stencils, stencil_indices, sendcounts, recvcounts);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  free(stencil_indices);
}

#endif // YAC_NETCDF_ENABLED

void yac_interp_weights_write_to_file(
  struct yac_interp_weights * weights, char const * filename,
  char const * src_grid_name, char const * tgt_grid_name,
  size_t src_grid_size, size_t tgt_grid_size) {

#ifndef YAC_NETCDF_ENABLED

  UNUSED(weights);
  UNUSED(filename);
  UNUSED(src_grid_name);
  UNUSED(tgt_grid_name);
  UNUSED(src_grid_size);
  UNUSED(tgt_grid_size);

  die(
    "ERROR(yac_interp_weights_write_to_file): "
    "YAC is built without the NetCDF support");
#else

  MPI_Comm comm = weights->comm;
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // determine processes that will do output
  int io_flag;
  int * io_ranks;
  int num_io_ranks;
  yac_get_io_ranks(comm, &io_flag, &io_ranks, &num_io_ranks);

  // determine range of global ids
  yac_int min_tgt_global_id, max_tgt_global_id;
  stencil_determine_tgt_global_id_range(
    weights->stencils, weights->stencils_size,
    &min_tgt_global_id, &max_tgt_global_id, comm);

  // determine io owners for all stencils
  int * io_owner =
    xmalloc(weights->stencils_size * sizeof(*io_owner));
  determine_stencils_io_owner(
    weights->stencils, weights->stencils_size,
    min_tgt_global_id, max_tgt_global_id,
    num_io_ranks, io_owner);
  for (size_t i = 0; i < weights->stencils_size; ++i)
    io_owner[i] = io_ranks[io_owner[i]];
  free(io_ranks);

  size_t io_stencil_count = 0;
  struct interp_weight_stencil * io_stencils = NULL;

  // redistribute stencils into io decomposition
  yac_interp_weights_redist_stencils(
    comm, weights->stencils_size, weights->stencils, io_owner,
    &io_stencil_count, &io_stencils);
  free(io_owner);

  // distribute global grid sizes
  uint64_t grid_sizes[2] = {(uint64_t)src_grid_size, (uint64_t)tgt_grid_size};
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, grid_sizes, 2, MPI_UINT64_T, MPI_MAX, comm), comm);
  src_grid_size = (size_t)(grid_sizes[0]);
  tgt_grid_size = (size_t)(grid_sizes[1]);

  MPI_Comm io_comm;
  yac_mpi_call(MPI_Comm_split(comm, io_flag, comm_rank, &io_comm), comm);

  // all non-io processes exit here, the remaining ones work using their own
  // communicator
  if (!io_flag) {
    yac_mpi_call(MPI_Comm_free(&io_comm), comm);
    free(io_stencils);
    // ensure that the writing of the weight file is complete
    yac_mpi_call(MPI_Barrier(comm), comm);
    return;
  }

  // sort the stencils first by type (fixed first; fixed stencils are sorted
  // by their fixed value) and second by tgt id
  qsort(io_stencils, io_stencil_count, sizeof(*io_stencils),
        compare_interp_weight_stencil);

  yac_mpi_call(MPI_Comm_rank(io_comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(io_comm, &comm_size), comm);

  double * fixed_values = NULL;
  size_t num_fixed_values = 0;
  stencil_get_fixed_values(
    io_stencils, io_stencil_count, &fixed_values, &num_fixed_values, io_comm);
  size_t num_src_fields = weights->num_src_fields;
  size_t num_weights_per_link =
    stencil_get_num_weights_per_tgt(io_stencils, io_stencil_count, io_comm);

  size_t * size_t_buffer =
    xmalloc(2 * (num_fixed_values + num_src_fields) * sizeof(*size_t_buffer));
  size_t * num_tgt_per_fixed_value = size_t_buffer;
  size_t * num_links_per_src_field = size_t_buffer + num_fixed_values;
  size_t * fixed_offsets = size_t_buffer + num_fixed_values + num_src_fields;
  size_t * link_offsets = size_t_buffer + 2 * num_fixed_values + num_src_fields;

  size_t num_fixed_tgt = 0;
  size_t num_links = 0;
  stencil_get_counts(
    io_stencils, io_stencil_count, num_fixed_values, fixed_values,
    num_tgt_per_fixed_value, &num_fixed_tgt, num_src_fields,
    num_links_per_src_field, &num_links);

  stencil_xscan_offsets(
    num_fixed_values, num_tgt_per_fixed_value,
    num_src_fields, num_links_per_src_field,
    fixed_offsets, link_offsets, io_comm);

  if (comm_rank == comm_size - 1) {

    size_t * total_num_tgt_per_fixed_value =
      xmalloc(num_fixed_values * sizeof(*total_num_tgt_per_fixed_value));
    for (size_t i = 0, accu = 0; i < num_fixed_values; ++i) {
      total_num_tgt_per_fixed_value[i] =
        fixed_offsets[i] + num_tgt_per_fixed_value[i] - accu;
      accu += total_num_tgt_per_fixed_value[i];
    }
    size_t total_num_links = link_offsets[num_src_fields-1] +
                             num_links_per_src_field[num_src_fields-1];

    size_t * total_num_links_per_src_field =
      xmalloc(num_src_fields * sizeof(*total_num_links_per_src_field));
    for (size_t i = 0, accu = 0; i < num_src_fields; ++i) {
      total_num_links_per_src_field[i] =
        link_offsets[i] + num_links_per_src_field[i] - accu;
      accu += total_num_links_per_src_field[i];
    }

    create_weight_file(
      filename, src_grid_name, tgt_grid_name,
      num_fixed_values, fixed_values, total_num_tgt_per_fixed_value,
      total_num_links, num_weights_per_link,
      num_src_fields, total_num_links_per_src_field,
      weights->src_locations, weights->tgt_location,
      src_grid_size, tgt_grid_size);

    free(total_num_links_per_src_field);
    free(total_num_tgt_per_fixed_value);
  }
  free(fixed_values);

  // ensure that the basic weight file has been written
  yac_mpi_call(MPI_Barrier(io_comm), comm);
  yac_mpi_call(MPI_Comm_free(&io_comm), comm);

  int ncid;

  // open weight file
  yac_nc_open(filename, NC_WRITE | NC_SHARE, &ncid);

  if (num_fixed_tgt > 0) {

    int * tgt_address_fixed =
      xmalloc(num_fixed_tgt * sizeof(*tgt_address_fixed));
    stencil_get_tgt_address(io_stencils, num_fixed_tgt, tgt_address_fixed);

    // inquire variable ids
    int var_dst_add_fixed_id;
    yac_nc_inq_varid(ncid, "dst_address_fixed", &var_dst_add_fixed_id);

    // target ids that receive a fixed value to file
    for (size_t i = 0, offset = 0; i < num_fixed_values; ++i) {

      if (num_tgt_per_fixed_value[i] == 0) continue;

      size_t start[1] = {fixed_offsets[i]};
      size_t count[1] = {num_tgt_per_fixed_value[i]};
      YAC_HANDLE_ERROR(
        nc_put_vara_int(
          ncid, var_dst_add_fixed_id, start, count, tgt_address_fixed + offset));
      offset += num_tgt_per_fixed_value[i];
    }

    free(tgt_address_fixed);
  }

  if (num_links > 0) {

    int * src_address_link = xmalloc(num_links * sizeof(*src_address_link));
    int * tgt_address_link = xmalloc(num_links * sizeof(*tgt_address_link));
    double * w = xmalloc(num_links * num_weights_per_link * sizeof(*w));
    stencil_get_link_data(
      io_stencils + num_fixed_tgt, io_stencil_count - num_fixed_tgt,
      num_links_per_src_field, num_src_fields,
      src_address_link, tgt_address_link, w);

    int var_src_add_id, var_dst_add_id, var_weight_id;
    yac_nc_inq_varid(ncid, "src_address", &var_src_add_id);
    yac_nc_inq_varid(ncid, "dst_address", &var_dst_add_id);
    yac_nc_inq_varid(ncid, "remap_matrix", &var_weight_id);

    for (size_t i = 0, offset = 0; i < num_src_fields; ++i) {

      if (num_links_per_src_field[i] == 0) continue;

      size_t start[2] = {link_offsets[i], 0};
      size_t count[2] = {num_links_per_src_field[i], num_weights_per_link};

      YAC_HANDLE_ERROR(
        nc_put_vara_int(
          ncid, var_src_add_id, start, count, src_address_link + offset));
      YAC_HANDLE_ERROR(
        nc_put_vara_int(
          ncid, var_dst_add_id, start, count, tgt_address_link + offset));
      YAC_HANDLE_ERROR(
        nc_put_vara_double(
          ncid, var_weight_id, start, count,
          w + num_weights_per_link * offset));

      offset += num_links_per_src_field[i];
    }

    free(w);
    free(tgt_address_link);
    free(src_address_link);
  }

  // close weight file
  YAC_HANDLE_ERROR(nc_close(ncid));

  // ensure that the writing of the weight file is complete
  yac_mpi_call(MPI_Barrier(comm), comm);

  free(size_t_buffer);
  yac_interp_weight_stencils_delete(io_stencils, io_stencil_count);
#endif
}

size_t yac_interp_weights_get_interp_count(
  struct yac_interp_weights * weights) {

  return weights->stencils_size;
}

yac_int * yac_interp_weights_get_interp_tgt(
  struct yac_interp_weights * weights) {

  struct interp_weight_stencil * stencils = weights->stencils;
  size_t stencils_size = weights->stencils_size;

  yac_int * global_ids = xmalloc(stencils_size * sizeof(*global_ids));

  for (size_t i = 0; i < stencils_size; ++i)
    global_ids[i] = stencils[i].tgt.global_id;

  return global_ids;
}

MPI_Comm yac_interp_weights_get_comm(struct yac_interp_weights * weights) {
  return weights->comm;
}

void yac_interp_weights_delete(struct yac_interp_weights * weights) {

  if (weights  == NULL) return;

  yac_interp_weight_stencils_delete(weights->stencils, weights->stencils_size);
  free(weights->src_locations);
  free(weights);
}
