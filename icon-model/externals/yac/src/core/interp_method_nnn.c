// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <string.h>

#include "interp_method_internal.h"
#include "interp_method_nnn.h"
#include "yac_lapack_interface.h"
#include "yac_mpi_internal.h"

static size_t do_search_nnn(struct interp_method * method,
                            struct yac_interp_grid * interp_grid,
                            size_t * tgt_points, size_t count,
                            struct yac_interp_weights * weights);
static size_t do_search_1nn(struct interp_method * method,
                            struct yac_interp_grid * interp_grid,
                            size_t * tgt_points, size_t count,
                            struct yac_interp_weights * weights);
static size_t do_search_zero(struct interp_method * method,
                            struct yac_interp_grid * interp_grid,
                            size_t * tgt_points, size_t count,
                            struct yac_interp_weights * weights);
static void delete_nnn(struct interp_method * method);

static struct interp_method_vtable
  interp_method_nnn_vtable = {
    .do_search = do_search_nnn,
    .delete = delete_nnn},
  interp_method_1nn_vtable = {
    .do_search = do_search_1nn,
    .delete = delete_nnn},
  interp_method_zero_vtable = {
    .do_search = do_search_zero,
    .delete = delete_nnn};

struct interp_method_nnn {

  struct interp_method_vtable * vtable;
  void (*compute_weights)(
    double[3], yac_coordinate_pointer, size_t const, double *, double const);
  size_t n;
  double max_search_distance;
  double scale;
};

struct interp_method_1nn {

  struct interp_method_vtable * vtable;
  void (*interp_weights_add)(
    struct yac_interp_weights * weights, struct remote_points * tgts,
    struct remote_point * srcs);
  double max_search_distance;
};

struct interp_method_zero {

  struct interp_method_vtable * vtable;
};

static void compute_weights_avg(
  double tgt_coord[3], yac_coordinate_pointer src_coords, size_t const n,
  double * weights, double const dummy) {

  UNUSED(tgt_coord);
  UNUSED(src_coords);
  UNUSED(dummy);

  double weight = 1.0 / (double)n;
  for (size_t i = 0; i < n; ++i) weights[i] = weight;
}

static void compute_weights_dist(
  double tgt_coord[3], yac_coordinate_pointer src_coords, size_t const n,
  double * weights, double const dummy) {

  UNUSED(dummy);

  for (size_t i = 0; i < n; ++i) {

    double distance = get_vector_angle(tgt_coord, (double*)(src_coords[i]));

    // if the target and source point are nearly identical
    if (distance < yac_angle_tol) {
      for (size_t j = 0; j < n; ++j) weights[j] = 0.0;
      weights[i] = 1.0;
      return;
    }

    weights[i] = 1.0 / distance;
  }

  // compute scaling factor for the weights
  double inv_distance_sum = 0.0;
  for (size_t i = 0; i < n; ++i) inv_distance_sum += weights[i];
  double scale = 1.0 / inv_distance_sum;

  for (size_t i = 0; i < n; ++i) weights[i] *= scale;
}

static void compute_weights_gauss(
  double tgt_coord[3], yac_coordinate_pointer src_coords, size_t const n,
  double * weights, double const gauss_scale) {

  for (size_t i = 0; i < n; ++i) {

    double distance = get_vector_angle(tgt_coord, src_coords[i]);

    // if the target and source point are nearly identical
    if (distance < yac_angle_tol) {
      for (size_t j = 0; j < n; ++j) weights[j] = 0.0;
      weights[i] = 1.0;
      return;
    }
    weights[i] = distance * distance;
  }

  // a) compute sum of source point distances
  double src_distances_sum = 0.0;
  double src_distances_count = 0.5 * (double)(n * n - n);
  for (size_t i = 0; i < n-1; ++i)
     for (size_t j = i + 1; j < n; ++j)
        src_distances_sum += get_vector_angle(src_coords[i], src_coords[j]);

  // b) C = -1 / (c * d_mean^2)
  double src_distances_mean = src_distances_sum / src_distances_count;
  double scale = -1.0 / (gauss_scale * src_distances_mean * src_distances_mean);

  // c) calculate weights
  // w_i = e^(-d_i^2/(c*s^2))
  // w_i = e^(C * d_i^2)
  double weights_sum = 0.0;
  for (size_t i = 0; i < n; ++i) {
    weights[i] = exp(scale * weights[i]);
    weights_sum += weights[i];
  }

  // If the sum of the weights is very low, which can happen in case
  // the target point is very far away from the group of source points.
  if (fabs(weights_sum) < 1e-9) {

    // Due to limited accuracy the exact contribution of each source
    // point cannot be computed. Therefore, the normalisation would
    // generate NaN's. Hence we fall back to inverse distance weighted
    // averge for this target point.
    compute_weights_dist(
      tgt_coord, src_coords, n, weights, gauss_scale);
    return;
  }

  // d) scale weights such that SUM(w_i) == 1
  for (size_t i = 0; i < n; ++i) weights[i] /= weights_sum;
}

static void inverse(double * A, size_t n) {

// LAPACKE_dsytrf_work and LAPACKE_dsytri might not be available (see yac_lapack_interface.h).
#ifdef YAC_LAPACK_NO_DSYTR
  lapack_int ipiv[n+1];
  double work[n*n];

  for (size_t i = 0; i < n+1; ++i) ipiv[i] = 0;
  for (size_t i = 0; i < n*n; ++i) work[i] = 0;

  YAC_ASSERT(
    !LAPACKE_dgetrf(
      LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) n,
      A, (lapack_int) n, ipiv), "internal ERROR: dgetrf")

  YAC_ASSERT(
    !LAPACKE_dgetri_work(
      LAPACK_COL_MAJOR, (lapack_int) n, A, (lapack_int) n,
      ipiv, work, (lapack_int) (n * n)), "internal ERROR: dgetri")
#else
  lapack_int ipiv[n];
  double work[n];

  YAC_ASSERT(
    !LAPACKE_dsytrf_work(
      LAPACK_COL_MAJOR, 'L', (lapack_int) n , A,
      (lapack_int) n, ipiv, work, (lapack_int) n), "internal ERROR: dsytrf")

  YAC_ASSERT(
    !LAPACKE_dsytri_work(
      LAPACK_COL_MAJOR, 'L', (lapack_int) n , A,
      (lapack_int) n, ipiv, work), "internal ERROR: dsytri")

  for (size_t i = 0; i < n; ++i)
    for (size_t j = i + 1; j < n; ++j)
      A[j*n+i] = A[i*n+j];
#endif
}

static void compute_weights_rbf(
  double tgt_coord[3], yac_coordinate_pointer src_coords, size_t const n,
  double * weights, double const rbf_scale) {

  double A[n][n];
  double a[n];

  double sum_d = 0.0, scale_d = 1.0;

  // compute distance matrix for all found source points
  for (size_t i = 0; i < n; ++i) A[i][i] = 0.0;
  for (size_t i = 0; i < n - 1; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      double d = get_vector_angle(src_coords[i], src_coords[j]);
      A[i][j] = d;
      A[j][i] = d;
      sum_d += d;
    }
  }

  // compute and apply scale factor for distance matrix
  if (sum_d > 0.0)
    scale_d = ((double)((n - 1) * n)) / (2.0 * sum_d);
  scale_d /= rbf_scale;

  double sq_scale_d = scale_d * scale_d;

  // compute matrix A[n][n]
  // with A = rbf(A)
  //      rbf(a) => Radial basis function
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double d = A[i][j];
      // gauß rbf kernel
      A[i][j] = exp(-1.0 * d * d * sq_scale_d);
    }
  }

  // compute inverse of A
  inverse(&A[0][0], n);

  // compute a[NUM_NEIGH]
  // with a_i = rbf(d(x_i, y))
  //      x => vector containing the coordinates of the
  //           n nearest neighbours of the current
  //           target point
  //      y => coordinates of target point
  //      d(a, b) => great circle distance between point a and b
  //      rbf(a) => Radial basis function
  for (size_t i = 0; i < n; ++i) {
    double d = get_vector_angle(tgt_coord, src_coords[i]);
    // gauß rbf kernel
    a[i] = exp(-1.0 * d * d * sq_scale_d);
  }

  // compute weights
  // with w_i = SUM(A_inv[i][j]*a[j]) // j = 0..n-1
  for (size_t i = 0; i < n; ++i) {
    weights[i] = 0.0;
    for (size_t j = 0; j < n; ++j) weights[i] += A[i][j] * a[j];
  }
}

static size_t do_search_nnn (struct interp_method * method,
                             struct yac_interp_grid * interp_grid,
                             size_t * tgt_points, size_t count,
                             struct yac_interp_weights * weights) {

  struct interp_method_nnn * method_nnn = (struct interp_method_nnn *)method;

  YAC_ASSERT(
    yac_interp_grid_get_num_src_fields(interp_grid) == 1,
    "ERROR(do_search_nnn): invalid number of source fields")

  // get coordinates of target points
  yac_coordinate_pointer tgt_coords = xmalloc(count * sizeof(*tgt_coords));
  yac_interp_grid_get_tgt_coordinates(
    interp_grid, tgt_points, count, tgt_coords);

  size_t n = method_nnn->n;

  size_t * size_t_buffer = xmalloc((2 * n + 1) * count * sizeof(*size_t_buffer));
  size_t * src_points = size_t_buffer;
  size_t * src_points_reorder = size_t_buffer + n * count;
  size_t * reorder_idx = size_t_buffer + 2 * n * count;
  int * fail_flag = xmalloc(count * sizeof(*fail_flag));

  // search for the n nearest source points
  yac_interp_grid_do_nnn_search_src(
    interp_grid, tgt_coords, count, n, src_points,
    method_nnn->max_search_distance);

  size_t result_count = 0;

  // check for which cells the search was successful
  for (size_t i = 0, k = 0; i < count; ++i) {
    int flag = 0;
    for (size_t j = 0; j < n; ++j, ++k) flag |= src_points[k] == SIZE_MAX;
    result_count += !(size_t)(fail_flag[i] = flag);
    reorder_idx[i] = i;
  }

  // sort target points, for which the search was successful to the beginning of
  // the array
  for (size_t i = 0; i < count; ++i) reorder_idx[i] = i;
  yac_quicksort_index_int_size_t(fail_flag, count, reorder_idx);
  free(fail_flag);

  void (*compute_weights)(
    double[3], yac_coordinate_pointer, size_t, double *, double const) =
      method_nnn->compute_weights;
  double scale = method_nnn->scale;
  double * w = xmalloc(n * result_count * sizeof(*w));
  size_t * num_weights_per_tgt =
    xmalloc(result_count * sizeof(*num_weights_per_tgt));
  yac_coordinate_pointer src_coord_buffer = xmalloc(n * sizeof(*src_coord_buffer));
  yac_const_coordinate_pointer src_field_coordinates =
    yac_interp_grid_get_src_field_coords(interp_grid, 0);

  // compute the weights
  for (size_t i = 0; i < result_count; ++i) {
    size_t curr_reorder_idx = reorder_idx[i];
    for (size_t j = 0; j < n; ++j) {
      size_t curr_src_point = src_points[curr_reorder_idx * n + j];
      memcpy(src_coord_buffer[j], src_field_coordinates[curr_src_point],
             3 * sizeof(double));
      src_points_reorder[i * n + j] = curr_src_point;
    }

    compute_weights(
      tgt_coords[curr_reorder_idx], src_coord_buffer, n, w + i * n, scale);
    num_weights_per_tgt[i] = n;
  }

  // move the non-interpolated target points to the end
  for (size_t i = 0; i < count; ++i) src_points[reorder_idx[i]] = i;
  yac_quicksort_index_size_t_size_t(src_points, count, tgt_points);

  free(src_coord_buffer);
  free(tgt_coords);

  struct remote_point * srcs =
    yac_interp_grid_get_src_remote_points(
      interp_grid, 0, src_points_reorder, result_count * n);

  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(
        interp_grid, tgt_points, result_count),
    .count = result_count};

  // store weights
  yac_interp_weights_add_wsum(
    weights, &tgts, num_weights_per_tgt, srcs, w);

  free(size_t_buffer);
  free(tgts.data);
  free(srcs);
  free(w);
  free(num_weights_per_tgt);

  return result_count;
}

static size_t do_search_1nn(struct interp_method * method,
                            struct yac_interp_grid * interp_grid,
                            size_t * tgt_points, size_t count,
                            struct yac_interp_weights * weights) {

  YAC_ASSERT(
    yac_interp_grid_get_num_src_fields(interp_grid) == 1,
    "ERROR(do_search_1nn): invalid number of source fields")

  // get coordinates of target points
  yac_coordinate_pointer tgt_coords = xmalloc(count * sizeof(*tgt_coords));
  yac_interp_grid_get_tgt_coordinates(
    interp_grid, tgt_points, count, tgt_coords);

  size_t * src_points = xmalloc(count * sizeof(*src_points));

  // search for the nearest source points
  yac_interp_grid_do_nnn_search_src(
    interp_grid, tgt_coords, count, 1, src_points,
    ((struct interp_method_1nn *)method)->max_search_distance);

  // move the non-interpolated target points to the end
  // (source points == SIZE_MAX if no source point was found)
  yac_quicksort_index_size_t_size_t(src_points, count, tgt_points);
  free(tgt_coords);

  // count the number of target points for which we found a source point
  size_t result_count = 0;
  for (; (result_count < count) && (src_points[result_count] != SIZE_MAX);
       ++result_count);

  struct remote_point * srcs =
    yac_interp_grid_get_src_remote_points(
      interp_grid, 0, src_points, result_count);

  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(
        interp_grid, tgt_points, result_count),
    .count = result_count};

  // store weights
  ((struct interp_method_1nn *)method)->interp_weights_add(
    weights, &tgts, srcs);

  free(src_points);
  free(tgts.data);
  free(srcs);

  return result_count;
}

static void interp_weights_add_zero(
    struct yac_interp_weights * weights, struct remote_points * tgts,
    struct remote_point src) {

  size_t * num_src_per_tgt = xmalloc(tgts->count * sizeof(*num_src_per_tgt));
  double * w = xmalloc(tgts->count * sizeof(*w));
  struct remote_point * srcs = xmalloc(tgts->count * sizeof(*srcs));

  for (size_t i = 0; i < tgts->count; ++i) {
    num_src_per_tgt[i] = 1;
    w[i] = 0.0;
    srcs[i] = src;
  }

  yac_interp_weights_add_wsum(weights, tgts, num_src_per_tgt, srcs, w);

  free(w);
  free(num_src_per_tgt);
}

static size_t do_search_zero(struct interp_method * method,
                             struct yac_interp_grid * interp_grid,
                             size_t * tgt_points, size_t count,
                             struct yac_interp_weights * weights) {

  YAC_ASSERT(
    yac_interp_grid_get_num_src_fields(interp_grid) == 1,
    "ERROR(do_search_zero): invalid number of source fields")

  UNUSED(method);

  // get global ids of all valid local source points
  size_t * src_points, src_count;
  yac_interp_grid_get_src_points(interp_grid, 0, &src_points, &src_count);
  yac_int * src_global_ids = xmalloc(src_count * sizeof(*src_global_ids));
  if (src_count > 0)
    yac_interp_grid_get_src_global_ids(
      interp_grid, src_points, src_count, 0, src_global_ids);
  free(src_points);

  // determine smallest global id of all unmasked source points in the
  // local data
  yac_int min_src_global_id = XT_INT_MAX;
  for (size_t i = 0; i < src_count; ++i)
    if (src_global_ids[i] < min_src_global_id)
      min_src_global_id = src_global_ids[i];
  free(src_global_ids);

  // determine smallest global id of all unmasked source points in the
  // global data
  MPI_Comm comm = yac_interp_grid_get_MPI_Comm(interp_grid);
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &min_src_global_id, 1, yac_int_dt, MPI_MIN, comm), comm);
  YAC_ASSERT(
    min_src_global_id != XT_INT_MAX,
    "ERROR(do_search_zero): unable to find unmasked source points");

  // get the local id of the source point with the smallest global ids
  // (add it to the local data if it does not yet exist there)
  size_t min_src_local_id;
  yac_interp_grid_src_global_to_local(
    interp_grid, 0, &min_src_global_id, 1, &min_src_local_id);

  struct remote_point * src =
    yac_interp_grid_get_src_remote_points(
      interp_grid, 0, &min_src_local_id, 1);

  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(
        interp_grid, tgt_points, count),
    .count = count};

  // store weights
  interp_weights_add_zero(weights, &tgts, *src);

  free(tgts.data);
  free(src);

  return count;
}

struct interp_method * yac_interp_method_nnn_new(struct yac_nnn_config config) {

  YAC_ASSERT_F(
    (config.max_search_distance >= 0.0) &&
    (config.max_search_distance <= M_PI),
    "ERROR(yac_interp_method_nnn_new): "
    "unsupported value for max_search_distance (%lf), has to be in [0.0,PI]",
    config.max_search_distance)

  if (config.type == YAC_INTERP_NNN_ZERO) {

    YAC_ASSERT(
      config.n == 1,
      "ERROR(yac_interp_method_nnn_new): n != 1 does not make sense for zero")
    YAC_ASSERT(
      config.max_search_distance == 0.0, "ERROR(yac_interp_method_nnn_new): "
      "max_search_distance != 0.0 does not make sense for zero")

    struct interp_method_zero * method = xmalloc(1 * sizeof(*method));
    method->vtable = &interp_method_zero_vtable;
    return (struct interp_method*)method;
  }

  double max_search_distance =
    (config.max_search_distance == 0.0)?M_PI:config.max_search_distance;

  if (config.n == 1) {

    struct interp_method_1nn * method = xmalloc(1 * sizeof(*method));
    method->vtable = &interp_method_1nn_vtable;
    method->interp_weights_add = yac_interp_weights_add_direct;
    method->max_search_distance = max_search_distance;
    return (struct interp_method*)method;
  }

  struct interp_method_nnn * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_nnn_vtable;

  method->max_search_distance = max_search_distance;
  method->n = config.n;

  YAC_ASSERT(
    (config.type == YAC_INTERP_NNN_AVG) ||
    (config.type == YAC_INTERP_NNN_DIST) ||
    (config.type == YAC_INTERP_NNN_GAUSS) ||
    (config.type == YAC_INTERP_NNN_RBF),
    "ERROR(yac_interp_method_nnn_new): invalid NNN type")

  switch (config.type) {
    case(YAC_INTERP_NNN_AVG):
      YAC_ASSERT(
        config.n >= 1,
        "ERROR(yac_interp_method_nnn_new): n needs to be at least 1 for "
        "distance weighted average")
      method->scale = 0.0;
      method->compute_weights = compute_weights_avg;
      break;
    case(YAC_INTERP_NNN_DIST):
      YAC_ASSERT(
        config.n >= 2,
        "ERROR(yac_interp_method_nnn_new): n needs to be at least 2 for "
        "distance weighted average")
      method->scale = 0.0;
      method->compute_weights = compute_weights_dist;
      break;
    case(YAC_INTERP_NNN_GAUSS):
      YAC_ASSERT(
        config.n >= 2,
        "ERROR(yac_interp_method_nnn_new): n needs to be at least 2 for gauss")
      method->scale = config.data.gauss_scale;
      method->compute_weights = compute_weights_gauss;
      break;
    default:
    case(YAC_INTERP_NNN_RBF):
      YAC_ASSERT(
        config.n >= 2,
        "ERROR(yac_interp_method_nnn_new): n needs to be at least 2 for RBF")
      method->scale = config.data.rbf_scale;
      method->compute_weights = compute_weights_rbf;
      break;
  };

  return (struct interp_method*)method;
}

static void delete_nnn(struct interp_method * method) {
  free(method);
}
