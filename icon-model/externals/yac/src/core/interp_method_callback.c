// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <string.h>

#include "interp_method_internal.h"
#include "utils_core.h"
#include "yac_mpi_internal.h"
#include "ensure_array_size.h"
#include "interp_method_callback.h"

struct tgt_request_data {
  uint64_t src_orig_pos;
  yac_int src_global_id;
  double tgt_coord[3];
};

static size_t do_search_callback(struct interp_method * method,
                                 struct yac_interp_grid * interp_grid,
                                 size_t * tgt_points, size_t count,
                                 struct yac_interp_weights * weights);
static void delete_callback(struct interp_method * method);

static struct interp_method_vtable
  interp_method_callback_vtable = {
    .do_search = do_search_callback,
    .delete = delete_callback};

struct interp_method_callback {

  struct interp_method_vtable * vtable;
  yac_func_compute_weights compute_weights_callback;
  void * user_data;
};

static void get_orig_data(
  struct remote_point const * remote_point,
  int * owner_rank, size_t * orig_pos) {

  struct remote_point_info const * point_infos = NULL;

  if (remote_point->data.count == 1) {
    point_infos = &(remote_point->data.data.single);
  } else {
    int min_rank = INT_MAX;
    for (int i = 0; i < remote_point->data.count; ++i) {
      int curr_rank = remote_point->data.data.multi[i].rank;
      if (curr_rank < min_rank) {
        point_infos = remote_point->data.data.multi + i;
        min_rank = curr_rank;
      }
    }
  }

  YAC_ASSERT(point_infos != NULL, "ERROR(get_orig_data): internal error")

  *owner_rank = point_infos->rank;
  *orig_pos = point_infos->orig_pos;
}

static MPI_Datatype yac_get_request_data_mpi_datatype(MPI_Comm comm) {

  struct tgt_request_data dummy;
  MPI_Datatype tgt_request_data_dt;
  int array_of_blocklengths[] = {1, 1, 3};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.src_orig_pos) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.src_global_id) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.tgt_coord[0]) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] =
    {MPI_UINT64_T, yac_int_dt, MPI_DOUBLE};
  yac_mpi_call(
    MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements,
                           array_of_types, &tgt_request_data_dt), comm);
  return yac_create_resized(tgt_request_data_dt, sizeof(dummy), comm);
}

// move the target points, which have a valid source cells
// (invalid src_cell == SIZE_MAX) to the front of the array and count the number
// of valid results
static size_t get_valid_results(
  size_t * src_cells, size_t count, size_t * tgt_points) {

  size_t valid_count = 0;
  size_t end_pos = count;

  while (valid_count < end_pos) {

    // find next invalid result
    while ((valid_count < end_pos) && (src_cells[valid_count] != SIZE_MAX))
      ++valid_count;

    // find valid result from the end of the array
    do {
      --end_pos;
    } while ((valid_count < end_pos) && (src_cells[end_pos] == SIZE_MAX));

    // switch valid result with invalid one
    if (valid_count < end_pos) {
      size_t temp_src_cell = src_cells[valid_count];
      size_t temp_tgt_point = tgt_points[valid_count];
      src_cells[valid_count] = src_cells[end_pos];
      tgt_points[valid_count] = tgt_points[end_pos];
      src_cells[end_pos] = temp_src_cell;
      tgt_points[end_pos] = temp_tgt_point;
      ++valid_count;
    }
  }

  return valid_count;
}

static size_t do_search_callback (struct interp_method * method,
                                  struct yac_interp_grid * interp_grid,
                                  size_t * tgt_points, size_t count,
                                  struct yac_interp_weights * weights) {

  struct interp_method_callback * method_callback =
    (struct interp_method_callback *)method;

  MPI_Comm comm = yac_interp_grid_get_MPI_Comm(interp_grid);
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // get coordinates of target points
  yac_coordinate_pointer tgt_coords = xmalloc(count * sizeof(*tgt_coords));
  yac_interp_grid_get_tgt_coordinates(
    interp_grid, tgt_points, count, tgt_coords);

  // get matching source cells for all target points
  size_t * src_cells = xmalloc(count * sizeof(*src_cells));
  yac_interp_grid_do_points_search(
    interp_grid, tgt_coords, count, src_cells);
  free(tgt_coords);

  // move the target points, which have a valid source cells
  // (invalid src_cell == SIZE_MAX) to the front of the array
  size_t temp_result_count = get_valid_results(src_cells, count, tgt_points);

  // get all unique source result cells
  yac_quicksort_index_size_t_size_t(src_cells, temp_result_count, tgt_points);
  size_t num_unique_src_cells = 0;
  size_t * src_to_unique_src =
    xmalloc(temp_result_count * sizeof(*src_to_unique_src));
  for (size_t i = 0, prev_src_cell = SIZE_MAX; i < temp_result_count; ++i) {
    size_t curr_src_cell = src_cells[i];
    if (curr_src_cell != prev_src_cell) {
      src_cells[num_unique_src_cells++] = curr_src_cell;
      prev_src_cell = curr_src_cell;
    }
    src_to_unique_src[i] = num_unique_src_cells - 1;
  }

  // get remote point information for all unique source result cell
  struct remote_point * src_remote_points =
    yac_interp_grid_get_src_remote_points2(
      interp_grid, YAC_LOC_CELL, src_cells, num_unique_src_cells);
  free(src_cells);

  // get original owners of all unique source result cells
  int * orig_src_cell_ranks =
    xmalloc(num_unique_src_cells * sizeof(*orig_src_cell_ranks));
  size_t * orig_src_cell_pos =
    xmalloc(num_unique_src_cells * sizeof(*orig_src_cell_pos));
  for (size_t i = 0; i < num_unique_src_cells; ++i)
    get_orig_data(
      src_remote_points + i, orig_src_cell_ranks + i, orig_src_cell_pos + i);

  // set up communication buffers
  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);
  // count number of target requests per rank
  for (size_t i = 0; i < temp_result_count; ++i)
    sendcounts[orig_src_cell_ranks[src_to_unique_src[i]]]++;
  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);
  size_t request_count = recvcounts[comm_size - 1] + rdispls[comm_size - 1];

  // pack target request data
  yac_const_coordinate_pointer tgt_field_coords =
    yac_interp_grid_get_tgt_field_coords(interp_grid);
  struct tgt_request_data * request_buffer =
    xmalloc((temp_result_count + request_count) * sizeof(*request_buffer));
  struct tgt_request_data * request_send_buffer = request_buffer;
  struct tgt_request_data * request_recv_buffer =
    request_buffer + temp_result_count;
  size_t * new_tgt_points =
    xmalloc(temp_result_count * sizeof(*new_tgt_points));
  for (size_t i = 0; i < temp_result_count; ++i) {
    size_t unique_src_cell_idx = src_to_unique_src[i];
    size_t pos = sdispls[orig_src_cell_ranks[unique_src_cell_idx] + 1]++;
    request_send_buffer[pos].src_orig_pos =
      orig_src_cell_pos[unique_src_cell_idx];
    request_send_buffer[pos].src_global_id =
      src_remote_points[unique_src_cell_idx].global_id;
    memcpy(request_send_buffer[pos].tgt_coord, tgt_field_coords[tgt_points[i]],
           3 * sizeof(double));
    new_tgt_points[pos] = tgt_points[i];
  }
  free(src_remote_points);
  free(orig_src_cell_pos);
  free(src_to_unique_src);
  free(orig_src_cell_ranks);

  // bring target points into the order in which we will receive the results
  memcpy(tgt_points, new_tgt_points, temp_result_count * sizeof(*tgt_points));
  free(new_tgt_points);

  // transfer tgt coords, orig src pos and src id
  MPI_Datatype request_data_dt = yac_get_request_data_mpi_datatype(comm);
  yac_alltoallv_p2p(
    request_send_buffer, sendcounts, sdispls,
    request_recv_buffer, recvcounts, rdispls,
    sizeof(*request_recv_buffer), request_data_dt, comm);
  yac_mpi_call(MPI_Type_free(&request_data_dt), comm);

  yac_func_compute_weights compute_weights =
    method_callback->compute_weights_callback;
  void * user_data = method_callback->user_data;
  size_t num_src_fields = yac_interp_grid_get_num_src_fields(interp_grid);
  uint64_t * uint64_t_buffer =
    xmalloc(num_src_fields * (request_count + temp_result_count) *
            sizeof(*uint64_t_buffer));
  uint64_t * temp_num_results_per_src_field_per_tgt = uint64_t_buffer;
  uint64_t * num_results_per_src_field_per_tgt_uint64 =
    uint64_t_buffer + num_src_fields * request_count;
  yac_int * temp_global_ids = NULL;
  size_t temp_global_ids_array_size = 0;
  double * temp_w = NULL;
  size_t temp_w_array_size = 0;
  size_t temp_weights_count = 0;

  // the weight computation function should be available on
  // all source
  YAC_ASSERT(
    (request_count == 0) || (compute_weights != NULL),
    "ERROR(do_search_callback): "
    "no callback routine defined on source process")

  // compute weights and store results
  for (size_t i = 0, k = 0; i < request_count; ++i) {

    // get weights for the current target point from the user
    int const * curr_global_result_points[num_src_fields];
    double * curr_result_weights[num_src_fields];
    size_t curr_result_counts[num_src_fields];
    compute_weights(
      (double const *)(request_recv_buffer[i].tgt_coord),
      (int)(request_recv_buffer[i].src_global_id),
      (size_t)(request_recv_buffer[i].src_orig_pos),
      curr_global_result_points, curr_result_weights, curr_result_counts,
      user_data);

    // copy current results
    for (size_t j = 0; j < num_src_fields; ++j, ++k) {
      size_t curr_count = curr_result_counts[j];
      temp_num_results_per_src_field_per_tgt[k] = (uint64_t)curr_count;
      ENSURE_ARRAY_SIZE(
        temp_global_ids, temp_global_ids_array_size,
        temp_weights_count + curr_count);
      ENSURE_ARRAY_SIZE(
        temp_w, temp_w_array_size,
        temp_weights_count + curr_count);
      for (size_t l = 0; l < curr_count; ++l)
        temp_global_ids[temp_weights_count + l] =
          (yac_int)(curr_global_result_points[j][l]);
      memcpy(temp_w + temp_weights_count, curr_result_weights[j],
             curr_count * sizeof(*curr_result_weights));
      temp_weights_count += curr_count;
    }
  }
  free(request_buffer);

  // return number of results per source field per target point
  for (int i = 0; i < comm_size; ++i) {
    sendcounts[i] *= num_src_fields;
    recvcounts[i] *= num_src_fields;
    sdispls[i] *= num_src_fields;
    rdispls[i] *= num_src_fields;
  }
  yac_alltoallv_uint64_p2p(
    temp_num_results_per_src_field_per_tgt, recvcounts, rdispls,
    num_results_per_src_field_per_tgt_uint64, sendcounts, sdispls,  comm);

  // set up comm buffers for exchanging of the interpolation results and
  // count the total number of weights per source field
  size_t num_weights = 0;
  size_t  * total_num_results_per_src_field =
    xcalloc(num_src_fields, sizeof(*total_num_results_per_src_field));
  {
    uint64_t * curr_num_send_results = temp_num_results_per_src_field_per_tgt;
    uint64_t * curr_num_recv_results = num_results_per_src_field_per_tgt_uint64;
    size_t saccu = 0, raccu = 0;
    for (int i = 0; i < comm_size; ++i) {
      size_t num_send_results = recvcounts[i] / num_src_fields;
      size_t num_recv_results = sendcounts[i] / num_src_fields;
      size_t num_send_weights = 0;
      size_t num_recv_weights = 0;
      for (size_t j = 0; j < num_send_results; ++j)
        for (size_t k = 0; k < num_src_fields; ++k, ++curr_num_send_results)
          num_send_weights += (size_t)*curr_num_send_results;
      for (size_t j = 0; j < num_recv_results; ++j) {
        for (size_t k = 0; k < num_src_fields; ++k, ++curr_num_recv_results) {
          size_t curr_count = (size_t)*curr_num_recv_results;
          num_recv_weights += curr_count;
          total_num_results_per_src_field[k] += curr_count;
        }
      }
      sdispls[i] = saccu;
      rdispls[i] = raccu;
      sendcounts[i] = num_send_weights;
      recvcounts[i] = num_recv_weights;
      saccu += num_send_weights;
      raccu += num_recv_weights;
      num_weights += num_recv_weights;
    }
  }

  double * w = xmalloc(num_weights * sizeof(*w));
  yac_int * global_ids = xmalloc(num_weights * sizeof(*global_ids));

  // return interpolation results
  yac_alltoallv_dble_p2p(
    temp_w, sendcounts, sdispls, w, recvcounts, rdispls, comm);
  yac_alltoallv_yac_int_p2p(
    temp_global_ids, sendcounts, sdispls,
    global_ids, recvcounts, rdispls, comm);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  free(temp_w);
  free(temp_global_ids);

  // check which target points have results
  size_t * interpolated_flag =
    xmalloc(temp_result_count * sizeof(*interpolated_flag));
  size_t result_count = 0;
  for (size_t i = 0, k = 0; i < temp_result_count; ++i) {
    int flag = 0;
    for (size_t j = 0; j < num_src_fields; ++j, ++k)
      flag |= (num_results_per_src_field_per_tgt_uint64[k] > 0);
    if (flag) {
      if (result_count != i)
        memmove(
          num_results_per_src_field_per_tgt_uint64 + result_count * num_src_fields,
          num_results_per_src_field_per_tgt_uint64 + i * num_src_fields,
          num_src_fields * sizeof(*num_results_per_src_field_per_tgt_uint64));
      interpolated_flag[i] = result_count++;
    } else {
      interpolated_flag[i] = SIZE_MAX;
    }
  }

  // sort the target points that can be interpolated to the beginning
  // of the array
  yac_quicksort_index_size_t_size_t(
    interpolated_flag, temp_result_count, tgt_points);
  free(interpolated_flag);

  // sort global result ids into a per tgt per src_field order
  size_t * global_id_reorder_idx =
    xmalloc((num_weights + num_src_fields) * sizeof(*global_id_reorder_idx));
  size_t * src_field_displ = global_id_reorder_idx + num_weights;
  size_t max_num_results_per_src_field = 0;
  size_t * num_results_per_src_field_per_tgt =
    xmalloc(result_count * num_src_fields *
            sizeof(*num_results_per_src_field_per_tgt));
  for (size_t i = 0, accu = 0; i < num_src_fields; ++i) {
    src_field_displ[i] = accu;
    accu += total_num_results_per_src_field[i];
    if (max_num_results_per_src_field < total_num_results_per_src_field[i])
      max_num_results_per_src_field = total_num_results_per_src_field[i];
  }
  for (size_t i = 0, k = 0, l = 0; i < result_count; ++i) {
    for (size_t j = 0; j < num_src_fields; ++j, ++k) {
      num_results_per_src_field_per_tgt[k] =
        (size_t)(num_results_per_src_field_per_tgt_uint64[k]);
      size_t curr_count =
        (size_t)(num_results_per_src_field_per_tgt_uint64[k]);
      for (size_t m = 0; m < curr_count; ++m, ++l)
        global_id_reorder_idx[l] = src_field_displ[j]++;
    }
  }
  yac_quicksort_index_size_t_yac_int(
    global_id_reorder_idx, num_weights, global_ids);
  free(uint64_t_buffer);
  free(global_id_reorder_idx);

  // get remote data for all required source points
  struct remote_point * srcs_per_field[num_src_fields];
  size_t * result_point_buffer =
    xmalloc(max_num_results_per_src_field * sizeof(*result_point_buffer));
  for (size_t i = 0, offset = 0; i < num_src_fields; ++i) {
    size_t curr_count = total_num_results_per_src_field[i];
    yac_interp_grid_src_global_to_local(
      interp_grid, i, global_ids + offset, curr_count, result_point_buffer);
    srcs_per_field[i] =
      yac_interp_grid_get_src_remote_points(
        interp_grid, i, result_point_buffer, curr_count);
    offset += curr_count;
  }
  free(result_point_buffer);
  free(global_ids);
  free(total_num_results_per_src_field);

  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(
        interp_grid, tgt_points, result_count),
    .count = result_count};

  // generate weights
  yac_interp_weights_add_wsum_mf(
    weights, &tgts, num_results_per_src_field_per_tgt, srcs_per_field, w,
    num_src_fields);

  free(tgts.data);
  for (size_t i = 0; i < num_src_fields; ++i) free(srcs_per_field[i]);
  free(num_results_per_src_field_per_tgt);
  free(w);

  return result_count;
}

struct interp_method * yac_interp_method_callback_new(
  yac_func_compute_weights compute_weights_callback, void * user_data) {

  struct interp_method_callback * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_callback_vtable;
  method->compute_weights_callback = compute_weights_callback;
  method->user_data = user_data;

  return (struct interp_method*)method;
}

static void delete_callback(struct interp_method * method) {
  free(method);
}

typedef void (*func_dummy)(void);
static struct {
  yac_func_compute_weights callback;
  void * user_data;
  char * key;
} * callback_lookup_table = NULL;
static size_t callback_lookup_table_array_size = 0;
static size_t callback_lookup_table_size = 0;

void yac_interp_method_callback_add_compute_weights_callback(
  yac_func_compute_weights compute_weights_callback,
  void * user_data, char const * key) {

  YAC_ASSERT(
    key != NULL,
    "ERROR(yac_interp_method_callback_add_compute_weights_callback): "
    "key is NULL")

  for (size_t i = 0; i < callback_lookup_table_size; ++i) {

    // if the key as already been set
    if (!strcmp(callback_lookup_table[i].key, key)) {

      YAC_ASSERT_F(
        (callback_lookup_table[i].callback == compute_weights_callback) &&
        (callback_lookup_table[i].user_data == user_data),
        "ERROR(interp_method_callback_add_do_search_callback): "
        "identical key has been set before with different callback "
        "or user data (key = \"%s\")", key)

      // nothing needs to be added
      return;
    }
  }

  ENSURE_ARRAY_SIZE(callback_lookup_table, callback_lookup_table_array_size,
                    callback_lookup_table_size + 1);

  callback_lookup_table[callback_lookup_table_size].key =
    xmalloc((strlen(key)+1) * sizeof(*key));
  strcpy(callback_lookup_table[callback_lookup_table_size].key, key);
  callback_lookup_table[callback_lookup_table_size].callback =
    compute_weights_callback;
  callback_lookup_table[callback_lookup_table_size].user_data = user_data;
  callback_lookup_table_size++;
}

void yac_interp_method_callback_get_compute_weights_callback(
  char const * key, yac_func_compute_weights * compute_weights_callback,
  void ** user_data) {

  YAC_ASSERT(
    key != NULL,
    "ERROR(yac_interp_method_callback_get_compute_weights_callback): "
    "key is NULL")

  for (size_t i = 0; i < callback_lookup_table_size; ++i) {
    if (!strcmp(callback_lookup_table[i].key, key)) {
      *compute_weights_callback = callback_lookup_table[i].callback;
      *user_data = callback_lookup_table[i].user_data;
      return;
    }
  }
  *compute_weights_callback = NULL;
  *user_data = NULL;
  return;
}
