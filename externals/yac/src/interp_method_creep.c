/**
 * @file interp_method_creep.c
 *
 * @copyright Copyright  (C)  2021 Moritz Hanke <hanke@dkrz.de>
 * @brief contains creep fill interpolation method (inspired by\n
 *        Kara, A. Birol, Alan J. Wallcraft, and Harley E. Hurlburt.
 *        " A Correction for Land Contamination of Atmospheric Variables
 *        near Land–Sea Boundaries". Journal of Physical Oceanography 37.4
 *        (2007): 803-818. < https://doi.org/10.1175/JPO2984.1>.
 *        Web. 2 Jul. 2021. )
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

#include "interp_method_creep.h"
#include "yac_mpi.h"
#include "ensure_array_size.h"

static size_t do_search_creep(struct interp_method * method,
                              struct interp_grid * interp_grid,
                              size_t * tgt_points, size_t count,
                              struct interp_weights * weights);
static void delete_creep(struct interp_method * method);

static struct interp_method_vtable
  interp_method_creep_vtable = {
    .do_search = do_search_creep,
    .delete = delete_creep};

struct interp_method_creep {

  struct interp_method_vtable * vtable;
  int max_creep_distance;
};

struct stencil_info {
  int rank;          // rank of the owner of the stencil
  uint64_t idx;      // index of the stencil that is supposed to be used
  double weight;     // weight of the stencil
};

struct result_stencil {
  yac_int global_id; // global id of target point for which
                     // the stencil was originally generated
  size_t count;      // number of contributing stencils
  union {
    struct stencil_info single,
                        * multi;
  } data;
};

struct result_stencils {
  struct result_stencil * stencils; // result stencils
  size_t count;                     // number of result stencils
  struct stencil_info stencil_info_buffer[];
};

struct interp_result {
  size_t local_id; // local id of target points that is supposed to
                   // be interpolated by the stencil
  yac_int global_id; // global id of target points that is supposed to
                     // be interpolated by the stencil
  size_t idx;
  struct result_stencil stencil;
};

struct tgt_request {
  yac_int global_id; // global id of requested point
  int rank;          // rank of process which requested
  struct result_stencil stencil;
};

struct comm_stuff {
  int rank, size;
  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  MPI_Datatype stencil_info_dt;
  MPI_Comm comm;
};

static MPI_Datatype yac_get_stencil_info_mpi_datatype(MPI_Comm comm) {

  struct stencil_info dummy;
  MPI_Datatype stencil_info_dt;
  int array_of_blocklengths[] = {1, 1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.rank) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.idx) -
       (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.weight) -
       (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] =
    {MPI_INT, MPI_UINT64_T, MPI_DOUBLE};
  yac_mpi_call(
    MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements,
                           array_of_types, &stencil_info_dt), comm);
  return yac_create_resized(stencil_info_dt, sizeof(dummy), comm);
}

static void get_result_stencil_pack_sizes(
  void * results, size_t result_count, size_t result_size,
  struct result_stencil*(*get_stencil)(void*), size_t * pack_order,
  int * pack_sizes, MPI_Datatype stencil_info_dt, MPI_Comm comm) {

  int pack_size_global_id, pack_size_count;
  yac_mpi_call(MPI_Pack_size(1, yac_int_dt, comm, &pack_size_global_id), comm);
  yac_mpi_call(MPI_Pack_size(1, MPI_UINT64_T, comm, &pack_size_count), comm);

  for (size_t i = 0; i < result_count; ++i) {

    struct result_stencil * stencil =
      (*get_stencil)(
        (void*)((unsigned char*)results + pack_order[i] * result_size));
    size_t curr_count = stencil->count;
    int pack_size_stencils;
    yac_mpi_call(
      MPI_Pack_size(
        (int)curr_count, stencil_info_dt, comm, &pack_size_stencils), comm);
    pack_sizes[i] = pack_size_global_id +
                    pack_size_count +
                    pack_size_stencils;
  }
}

static void pack_result_stencil(
  struct result_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype stencil_info_dt, MPI_Comm comm) {

  // global_id
  yac_mpi_call(
    MPI_Pack(&(stencil->global_id), 1, yac_int_dt,
            buffer, buffer_size, position, comm), comm);

  // count
  uint64_t count_uint64_t = (uint64_t)(stencil->count);
  yac_mpi_call(
    MPI_Pack(&count_uint64_t, 1, MPI_UINT64_T,
             buffer, buffer_size, position, comm), comm);

  // stencils
  struct stencil_info * stencils =
    (count_uint64_t == 1)?(&(stencil->data.single)):(stencil->data.multi);
  yac_mpi_call(
    MPI_Pack(stencils, (int)count_uint64_t, stencil_info_dt,
             buffer, buffer_size, position, comm), comm);
}

static void pack_result_stencils(
  void * results, size_t result_count, size_t result_size,
  struct result_stencil*(*get_stencil)(void*), size_t * pack_order,
  void ** pack_data, int * pack_sizes, MPI_Datatype stencil_info_dt,
  MPI_Comm comm) {

  get_result_stencil_pack_sizes(
    results, result_count, result_size, get_stencil, pack_order,
    pack_sizes, stencil_info_dt, comm);

  size_t temp_total_pack_size = 0;
  for (size_t i = 0; i < result_count; ++i)
    temp_total_pack_size += (size_t)(pack_sizes[i]);

  void * pack_data_ = xmalloc(temp_total_pack_size);
  size_t total_pack_size = 0;

  for (size_t i = 0; i < result_count; ++i) {

    struct result_stencil * curr_stencil =
      (*get_stencil)(
        (void*)((unsigned char*)results + pack_order[i] * result_size));

    int position = 0;
    void * buffer = (void*)((char*)pack_data_ + total_pack_size);
    int buffer_size = (int)(temp_total_pack_size - total_pack_size);

    // stencil
    pack_result_stencil(
      curr_stencil, buffer, buffer_size, &position, stencil_info_dt, comm);

    pack_sizes[i] = position;
    total_pack_size += (size_t)position;
  }

  if (total_pack_size != temp_total_pack_size)
    pack_data_ = xrealloc(pack_data_, total_pack_size);

  *pack_data = pack_data_;
}

static void unpack_result_stencil(
  struct result_stencil * stencil, void * buffer, int buffer_size,
  int * position, MPI_Datatype stencil_info_dt, MPI_Comm comm,
  struct stencil_info ** stencil_info_buffer,
  size_t * stencil_info_buffer_array_size,
  size_t * stencil_info_buffer_size) {

  // global_id
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position,
               &(stencil->global_id), 1, yac_int_dt, comm), comm);

  // count
  uint64_t count_uint64_t;
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position,
               &count_uint64_t, 1, MPI_UINT64_T, comm), comm);
  size_t count = ((stencil->count = (size_t)count_uint64_t));

  // stencils
  struct stencil_info * stencils;
  if (count_uint64_t == 1) {
    stencils = &(stencil->data.single);
  } else {
    ENSURE_ARRAY_SIZE(
      *stencil_info_buffer, *stencil_info_buffer_array_size,
      *stencil_info_buffer_size + count);
    stencils =
      ((stencil->data.multi =
          *stencil_info_buffer + *stencil_info_buffer_size));
    *stencil_info_buffer_size += count;
  }
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, position,
               stencils, (int)count, stencil_info_dt, comm), comm);
}

static struct result_stencils * unpack_result_stencils(
  size_t count, void * packed_data, size_t packed_data_size,
  MPI_Datatype stencil_info_dt, MPI_Comm comm) {

  struct stencil_info * stencil_info_buffer = NULL;
  size_t stencil_info_buffer_array_size = 0;
  size_t stencil_info_buffer_size = 0;

  struct result_stencil * result_stencils =
    xmalloc(count * sizeof(*result_stencils));

  for (size_t i = 0, offset = 0; i < count; ++i) {

    int position = 0;
    void * curr_buffer = (void*)((char*)packed_data + offset);
    int buffer_size = (int)(packed_data_size - offset);
    struct result_stencil * curr_stencil = result_stencils + i;

    unpack_result_stencil(
      curr_stencil, curr_buffer, buffer_size, &position, stencil_info_dt, comm,
      &stencil_info_buffer, &stencil_info_buffer_array_size,
      &stencil_info_buffer_size);
    offset += (size_t)position;
  }

  struct result_stencils * results =
    xmalloc(sizeof(*results) +
            stencil_info_buffer_size * sizeof(results->stencil_info_buffer[0]));
  results->stencils = result_stencils;
  results->count = count;
  memcpy(&(results->stencil_info_buffer[0]), stencil_info_buffer,
         stencil_info_buffer_size * sizeof(*stencil_info_buffer));
  free(stencil_info_buffer);
  for (size_t i = 0, offset = 0; i < count; ++i) {
    size_t curr_count = result_stencils[i].count;
    if (curr_count > 1) {
      result_stencils[i].data.multi = &(results->stencil_info_buffer[offset]);
      offset += curr_count;
    }
  }

  return results;
}

static struct result_stencils * exchange_interp_results(
  void * results, size_t result_count, size_t result_size,
  struct result_stencil*(*result_get_stencil)(void*), size_t * pack_order,
  int * ranks, struct comm_stuff comm) {

  // mark and count local results
  size_t local_count = 0;
  for (size_t i = 0; i < result_count; ++i) {
    if (ranks[i] == comm.rank) {
      ranks[i] = INT_MAX;
      ++local_count;
    }
  }
  size_t send_count = result_count - local_count;

  // sort results by rank (local results go to the end of the array)
  yac_quicksort_index_int_size_t(ranks, result_count, pack_order);

  // pack the result stencils that need to be send to other processes
  void * send_buffer;
  int * pack_sizes = xmalloc(send_count * sizeof(*pack_sizes));
  pack_result_stencils(
    results, send_count, result_size,
    result_get_stencil, pack_order, &send_buffer, pack_sizes,
    comm.stencil_info_dt, comm.comm);

  // set up comm buffers
  size_t * num_results_per_rank =
    xmalloc((size_t)comm.size * sizeof(*num_results_per_rank));
  memset(comm.sendcounts, 0,
         (size_t)comm.size * sizeof(*comm.sendcounts));
  size_t j = 0;
  for (int rank = 0; rank < comm.size; ++rank) {
    size_t curr_num_results = 0;
    size_t curr_sendcount = 0;
    while ((j < send_count) && (ranks[j] == rank)) {
      curr_sendcount += (size_t)(pack_sizes[j++]);
      curr_num_results++;
    }
    num_results_per_rank[rank] = curr_num_results;
    YAC_ASSERT(
      curr_sendcount <= INT_MAX,
      "ERROR(exchange_interp_results): pack size to big")
    comm.sendcounts[rank] = curr_sendcount;
  }
  free(pack_sizes);
  yac_generate_alltoallv_args(
    1, comm.sendcounts, comm.recvcounts,
    comm.sdispls, comm.rdispls, comm.comm);
  yac_mpi_call(
    MPI_Alltoall(MPI_IN_PLACE, 1, YAC_MPI_SIZE_T,
                 num_results_per_rank, 1, YAC_MPI_SIZE_T, comm.comm),
    comm.comm);
  size_t recv_count = 0;
  for (int i = 0; i < comm.size; ++i)
    recv_count += num_results_per_rank[i];
  free(num_results_per_rank);

  size_t recv_size = comm.recvcounts[comm.size - 1] +
                     comm.rdispls[comm.size - 1];
  void * recv_buffer = xmalloc(recv_size);

  // exchange result stencils
  yac_alltoallv_packed_p2p(
    send_buffer, comm.sendcounts, comm.sdispls+1,
    recv_buffer, comm.recvcounts, comm.rdispls, comm.comm);
  free(send_buffer);

  // unpack stencils
  struct result_stencils * result_stencils =
    unpack_result_stencils(
      recv_count, recv_buffer, recv_size, comm.stencil_info_dt, comm.comm);
  free(recv_buffer);

  result_stencils->count += local_count;
  result_stencils->stencils =
    xrealloc(
      result_stencils->stencils,
      (recv_count + local_count) * sizeof(*(result_stencils->stencils)));

  // add the local result stenils
  struct result_stencil * local_stencils =
    result_stencils->stencils + recv_count;
  pack_order += send_count;
  for (size_t i = 0; i < local_count; ++i)
    local_stencils[i] =
      *(*result_get_stencil)(
        (void*)((unsigned char*)results + pack_order[i] * result_size));

  return result_stencils;
}

static struct result_stencil * interp_result_get_stencil(
  void * interp_result) {

  return &(((struct interp_result*)interp_result)->stencil);
}

static struct result_stencils * relocate_interp_results(
  struct interp_grid * interp_grid, struct comm_stuff comm,
  struct interp_result * interp_results, size_t result_count) {

  // get the the newly interpolated target points
  size_t * tgt_points = xmalloc(result_count * sizeof(*tgt_points));
  for (size_t i = 0; i < result_count; ++i)
    tgt_points[i] = interp_results[i].local_id;

  // get the distributed owners for all targets
  int * tgt_points_dist_owner =
    xmalloc(result_count * sizeof(*tgt_points_dist_owner));
  yac_interp_grid_determine_dist_tgt_owners(
    interp_grid, tgt_points, result_count, tgt_points_dist_owner);
  size_t * pack_order = tgt_points;
  for (size_t i = 0; i < result_count; ++i) pack_order[i] = i;

  struct result_stencils * relocated_results =
    exchange_interp_results(
      interp_results, result_count, sizeof(*interp_results),
      interp_result_get_stencil, pack_order, tgt_points_dist_owner,
      comm);
  free(pack_order);
  free(tgt_points_dist_owner);

  return relocated_results;
}

static void get_initial_results(
  struct interp_grid * interp_grid, struct comm_stuff comm,
  struct interp_weights * interp_weights,
  struct interp_result ** interp_results, size_t * result_count) {

  // get list of all already interpolated target points
  size_t num_interpolated_tgt =
    yac_interp_weights_get_interp_count(interp_weights);
  yac_int * interpolated_tgts_global_ids =
    yac_interp_weights_get_interp_tgt(interp_weights);
  size_t * interpolated_tgts_local_ids =
    xmalloc(num_interpolated_tgt * sizeof(*interpolated_tgts_local_ids));
  yac_interp_grid_tgt_global_to_local(
    interp_grid, interpolated_tgts_global_ids, num_interpolated_tgt,
    interpolated_tgts_local_ids);

  // initialise initial interpolation results
  struct interp_result * initial_interp_results =
    (num_interpolated_tgt > 0)?
      xmalloc(num_interpolated_tgt * sizeof(*initial_interp_results)):NULL;
  for (size_t i = 0; i < num_interpolated_tgt; ++i) {
    initial_interp_results[i].local_id = interpolated_tgts_local_ids[i];
    initial_interp_results[i].global_id = interpolated_tgts_global_ids[i];
    initial_interp_results[i].idx = SIZE_MAX;
    initial_interp_results[i].stencil.global_id =
      interpolated_tgts_global_ids[i];
    initial_interp_results[i].stencil.count = 1;
    initial_interp_results[i].stencil.data.single.rank = comm.rank;
    initial_interp_results[i].stencil.data.single.idx = (uint64_t)i;
    initial_interp_results[i].stencil.data.single.weight = 1.0;
  }
  free(interpolated_tgts_local_ids);
  free(interpolated_tgts_global_ids);

  *interp_results = initial_interp_results;
  *result_count = num_interpolated_tgt;
}

static int compare_interp_result_stencil(
  const void * a, const void * b) {

  int ret = (((const struct interp_result *)a)->stencil.count >
             ((const struct interp_result *)b)->stencil.count) -
            (((const struct interp_result *)a)->stencil.count <
             ((const struct interp_result *)b)->stencil.count);

  if (ret) return ret;

  return (((const struct interp_result *)a)->global_id >
          ((const struct interp_result *)b)->global_id) -
          (((const struct interp_result *)a)->global_id <
          ((const struct interp_result *)b)->global_id);
}

static int compare_interp_result_stencil_local_id(
  const void * a, const void * b) {

  return (((const struct interp_result *)a)->local_id >
          ((const struct interp_result *)b)->local_id) -
         (((const struct interp_result *)a)->local_id <
          ((const struct interp_result *)b)->local_id);
}

static int compare_interp_result_global_id(const void * a, const void * b) {

  return (((const struct interp_result *)a)->global_id >
          ((const struct interp_result *)b)->global_id) -
         (((const struct interp_result *)a)->global_id <
          ((const struct interp_result *)b)->global_id);
}

static int compare_tgt_request_stencil_count(
  const void * a, const void * b) {

  return (((const struct tgt_request *)a)->stencil.count >
          ((const struct tgt_request *)b)->stencil.count) -
         (((const struct tgt_request *)a)->stencil.count <
          ((const struct tgt_request *)b)->stencil.count);
}

static int compare_tgt_request_global_id(
  const void * a, const void * b) {

  return (((const struct tgt_request *)a)->global_id >
          ((const struct tgt_request *)b)->global_id) -
         (((const struct tgt_request *)a)->global_id <
          ((const struct tgt_request *)b)->global_id);
}

static int compare_result_stencil_global_id(
  const void * a, const void * b) {

  return (((const struct result_stencil *)a)->global_id >
          ((const struct result_stencil *)b)->global_id) -
         (((const struct result_stencil *)a)->global_id <
          ((const struct result_stencil *)b)->global_id);
}

static inline int compare_yac_int(const void * a, const void * b) {

  return ((*(const yac_int *)a) >
          (*(const yac_int *)b)) -
         ((*(const yac_int *)a) <
          (*(const yac_int *)b));
}

static void extract_interp_info(
  struct interp_grid * interp_grid, struct comm_stuff comm,
  struct interp_result * interp_results, size_t result_count,
  struct remote_points * interp_tgt_remote_points,
  size_t ** num_stencils_per_tgt_, size_t ** stencil_indices_,
  int ** stencil_ranks_, double ** w_) {

  // sort interpolation results by target local id
  qsort(interp_results, result_count, sizeof(*interp_results),
        compare_interp_result_stencil_local_id);

  size_t total_num_stencils = 0;
  size_t * interpolated_tgts_local_ids =
    xmalloc(result_count * sizeof(*interpolated_tgts_local_ids));
  for (size_t i = 0; i < result_count; ++i) {
    interpolated_tgts_local_ids[i] = interp_results[i].local_id;
    total_num_stencils += interp_results[i].stencil.count;
  }
  interp_tgt_remote_points->data =
    yac_interp_grid_get_tgt_remote_points(
      interp_grid, interpolated_tgts_local_ids, result_count);
  interp_tgt_remote_points->count = result_count;
  free(interpolated_tgts_local_ids);

  size_t * num_stencils_per_tgt =
    xmalloc(result_count * sizeof(*num_stencils_per_tgt));
  size_t * stencil_indices =
    xmalloc(total_num_stencils * sizeof(*stencil_indices));
  int * stencil_ranks =
    xmalloc(total_num_stencils * sizeof(*stencil_ranks));
  double * w = xmalloc(total_num_stencils * sizeof(*w));

  for (size_t i = 0, j = 0; i < result_count; ++i) {

    size_t curr_count = interp_results[i].stencil.count;

    num_stencils_per_tgt[i] = curr_count;

    if (curr_count == 1) {
      stencil_indices[j] = interp_results[i].stencil.data.single.idx;
      stencil_ranks[j] = interp_results[i].stencil.data.single.rank;
      w[j] = interp_results[i].stencil.data.single.weight;
      ++j;
    } else {
      for (size_t k = 0; k < curr_count; ++k, ++j) {
        stencil_indices[j] = interp_results[i].stencil.data.multi[k].idx;
        stencil_ranks[j] = interp_results[i].stencil.data.multi[k].rank;
        w[j] = interp_results[i].stencil.data.multi[k].weight;
      }
    }
  }

  *num_stencils_per_tgt_ = num_stencils_per_tgt;
  *stencil_indices_ = stencil_indices;
  *stencil_ranks_ = stencil_ranks;
  *w_ = w;
}

static struct result_stencil * tgt_request_get_stencil(
  void * tgt_request) {

  return &(((struct tgt_request*)tgt_request)->stencil);
}

static struct result_stencils * update_neigh_requests(
  struct interp_grid * interp_grid, struct comm_stuff comm,
  struct tgt_request * neigh_requests, size_t * request_count_,
  struct result_stencils * interp_stencils) {

  size_t request_count = *request_count_;
  size_t match_count = 0;

  struct result_stencil * stencils = interp_stencils->stencils;
  size_t stencil_count = interp_stencils->count;

  // match current result points with neighbour requests
  qsort(
    stencils, stencil_count, sizeof(*stencils),
    compare_result_stencil_global_id);
  for (size_t i = 0, j = 0; i < stencil_count; ++i) {
    yac_int curr_global_id = stencils[i].global_id;

    while ((j < request_count) &&
           (neigh_requests[j].global_id < curr_global_id)) ++j;
    while ((j < request_count) &&
           (neigh_requests[j].global_id == curr_global_id)) {

      ++match_count;
      neigh_requests[j].stencil = stencils[i];
      ++j;
    }
  }

  // sort neighbour request by stencil count
  // (open requests have a stencil count of 0)
  qsort(neigh_requests, request_count, sizeof(*neigh_requests),
        compare_tgt_request_stencil_count);

  request_count -= match_count;
  *request_count_ = request_count;
  struct tgt_request * neigh_request_matches = neigh_requests + request_count;

  // sort open neighbour requests by global id
  qsort(neigh_requests, request_count, sizeof(*neigh_requests),
        compare_tgt_request_global_id);

  size_t * pack_order = xmalloc(match_count * sizeof(*pack_order));
  int * origin_rank = xmalloc(match_count * sizeof(*origin_rank));
  for (size_t i = 0; i < match_count; ++i) {
    pack_order[i] = i;
    origin_rank[i] = neigh_request_matches[i].rank;
  }

  struct result_stencils * relocated_results =
    exchange_interp_results(
      neigh_request_matches, match_count, sizeof(*neigh_request_matches),
      tgt_request_get_stencil, pack_order, origin_rank, comm);
  free(origin_rank);
  free(pack_order);

  return relocated_results;
}

static void get_tgt_neigh_info_cell(
  struct interp_grid * interp_grid, size_t * tgt_local_ids,
  yac_int * tgt_global_ids, size_t count,
  size_t ** neigh_local_ids_, yac_int ** neigh_to_tgt_global_id_,
  size_t * total_num_neighbours_) {

  // get neighbour cells for all target cells
  size_t total_num_neighbours = 0;
  struct const_basic_grid_data * tgt_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_tgt(interp_grid);
  for (size_t i = 0; i < count; ++i)
    total_num_neighbours +=
      tgt_basic_grid_data->num_vertices_per_cell[tgt_local_ids[i]];
  size_t * neigh_local_ids =
    xmalloc(total_num_neighbours * sizeof(*neigh_local_ids));
  yac_interp_grid_get_tgt_cell_neighbours(
    interp_grid, tgt_local_ids, count, neigh_local_ids);

  // generate mapping between neighbour global ids and target points
  yac_int * neigh_to_tgt_global_id =
    xmalloc(total_num_neighbours * sizeof(*neigh_to_tgt_global_id));
  for (size_t i = 0, j = 0; i < count; ++i) {
    int curr_num_neigh =
      tgt_basic_grid_data->num_vertices_per_cell[tgt_local_ids[i]];
    yac_int curr_tgt_global_id = tgt_global_ids[i];
    for (int k = 0; k < curr_num_neigh; ++k, ++j)
      neigh_to_tgt_global_id[j] = curr_tgt_global_id;
  }

  *neigh_local_ids_ = neigh_local_ids;
  *neigh_to_tgt_global_id_ = neigh_to_tgt_global_id;
  *total_num_neighbours_ = total_num_neighbours;
}

static void get_tgt_neigh_info_vertex(
  struct interp_grid * interp_grid, size_t * tgt_local_ids,
  yac_int * tgt_global_ids, size_t count,
  size_t ** neigh_local_ids_, yac_int ** neigh_to_tgt_global_id_,
  size_t * total_num_neighbours_) {

  int * num_neighs_per_vertex = xmalloc(count * sizeof(*num_neighs_per_vertex));
  size_t * neigh_vertices;

  yac_interp_grid_get_tgt_vertex_neighbours(
    interp_grid, tgt_local_ids, count,
    &neigh_vertices, num_neighs_per_vertex);

  size_t total_num_neighbours = 0;
  for (size_t i = 0; i < count; ++i)
    total_num_neighbours += (size_t)(num_neighs_per_vertex[i]);

  // generate mapping between neighbour global ids and target points
  yac_int * neigh_to_tgt_global_id =
    xmalloc(total_num_neighbours * sizeof(*neigh_to_tgt_global_id));
  for (size_t i = 0, j = 0; i < count; ++i) {
    int curr_num_neigh = num_neighs_per_vertex[i];
    yac_int curr_tgt_global_id = tgt_global_ids[i];
    for (int k = 0; k < curr_num_neigh; ++k, ++j)
      neigh_to_tgt_global_id[j] = curr_tgt_global_id;
  }
  free(num_neighs_per_vertex);

  *neigh_local_ids_ = neigh_vertices;
  *neigh_to_tgt_global_id_ = neigh_to_tgt_global_id;
  *total_num_neighbours_ = total_num_neighbours;
}

// extracts basic information about the neighbours of the target points
static void get_tgt_neigh_info(
  struct interp_grid * interp_grid,
  size_t * tgt_local_ids, yac_int * tgt_global_ids, size_t count,
  size_t ** neigh_local_ids_, yac_int ** neigh_global_ids_,
  yac_int ** neigh_to_tgt_global_id_, size_t * total_num_neighbours_) {

  // get basic neighbour information
  size_t * neigh_local_ids;
  yac_int * neigh_to_tgt_global_id;
  size_t total_num_neighbours;
  enum yac_location tgt_field_location =
    yac_interp_grid_get_tgt_field_location(interp_grid);
  YAC_ASSERT(
    (tgt_field_location == CELL) || (tgt_field_location == CORNER),
    "ERROR(get_tgt_neigh_info): unsupported target field location")
  if (tgt_field_location == CELL)
    get_tgt_neigh_info_cell(
      interp_grid, tgt_local_ids, tgt_global_ids, count,
      &neigh_local_ids, &neigh_to_tgt_global_id, &total_num_neighbours);
  else
    get_tgt_neigh_info_vertex(
      interp_grid, tgt_local_ids, tgt_global_ids, count,
      &neigh_local_ids, &neigh_to_tgt_global_id, &total_num_neighbours);

  // remove invalid neighbour indices
  yac_quicksort_index_size_t_yac_int(
    neigh_local_ids, total_num_neighbours, neigh_to_tgt_global_id);
  while((total_num_neighbours > 0) &&
        (neigh_local_ids[total_num_neighbours-1] == SIZE_MAX))
    --total_num_neighbours;

  // get global ids for all neighbours
  yac_int * neigh_global_ids =
    xmalloc(total_num_neighbours * sizeof(*neigh_global_ids));
  yac_interp_grid_get_tgt_global_ids(
    interp_grid, neigh_local_ids, total_num_neighbours, neigh_global_ids);

  *neigh_local_ids_ =
    xrealloc(neigh_local_ids, total_num_neighbours * sizeof(*neigh_local_ids));
  *neigh_global_ids_ = neigh_global_ids;
  *neigh_to_tgt_global_id_ =
    xrealloc(neigh_to_tgt_global_id,
             total_num_neighbours * sizeof(*neigh_to_tgt_global_id));
  *total_num_neighbours_ = total_num_neighbours;
}

// sends request for target points neighbours to the respective
// distributed owners
static void send_neigh_request(
  struct interp_grid * interp_grid, struct comm_stuff comm,
  size_t * neigh_local_ids, yac_int * neigh_global_ids,
  size_t num_neighbours, struct tgt_request ** neigh_requests_,
  size_t * request_count_) {

  // get the distributed owner of the neighbour target points
  int * neigh_dist_owner =
    xmalloc(num_neighbours * sizeof(*neigh_dist_owner));
  yac_interp_grid_determine_dist_tgt_owners(
    interp_grid, neigh_local_ids, num_neighbours,
    neigh_dist_owner);

  yac_int * send_neigh_global_ids =
    xmalloc(num_neighbours * sizeof(*send_neigh_global_ids));
  memcpy(send_neigh_global_ids, neigh_global_ids,
        num_neighbours * sizeof(*send_neigh_global_ids));

  // sort send buffer by rank
  yac_quicksort_index_int_yac_int(
    neigh_dist_owner, num_neighbours, send_neigh_global_ids);

  // remove duplicated global ids
  size_t to = 0, new_to = 0, from = 0, new_from = 0;
  for (int rank = 0; rank < comm.size; ++rank) {
    while ((new_from  < num_neighbours) &&
           (neigh_dist_owner[new_from] == rank)) new_from++;
    size_t curr_count = new_from - from;
    qsort(send_neigh_global_ids + from, curr_count,
          sizeof(*send_neigh_global_ids), compare_yac_int);
    yac_int prev_global_id =
      (curr_count > 0)?(send_neigh_global_ids[from]-1):0;
    for (; from < new_from; ++from) {
      yac_int curr_global_id = send_neigh_global_ids[from];
      if (prev_global_id != curr_global_id) {
        send_neigh_global_ids[new_to++] = curr_global_id;
        prev_global_id = curr_global_id;
      }
    }
    curr_count = new_to - to;
    to = new_to;
    comm.sendcounts[rank] = curr_count;
  }
  free(neigh_dist_owner);

  // send request for all neighbours
  yac_generate_alltoallv_args(
    1, comm.sendcounts, comm.recvcounts, comm.sdispls, comm.rdispls, comm.comm);
  size_t request_count = comm.recvcounts[comm.size-1] +
                         comm.rdispls[comm.size-1];
  yac_int * recv_neigh_global_ids =
    xmalloc(request_count * sizeof(*recv_neigh_global_ids));
  yac_alltoallv_yac_int_p2p(
    send_neigh_global_ids, comm.sendcounts, comm.sdispls+1,
    recv_neigh_global_ids, comm.recvcounts, comm.rdispls, comm.comm);
  struct tgt_request * neigh_requests =
    xmalloc(request_count * sizeof(*neigh_requests));
  for (int i = 0, k = 0; i < comm.size; ++i) {
    for (size_t j = 0; j < comm.recvcounts[i]; ++j, ++k) {
      neigh_requests[k].global_id = recv_neigh_global_ids[k];
      neigh_requests[k].rank = i;
      neigh_requests[k].stencil.count = 0;
    }
  }
  qsort(neigh_requests, request_count, sizeof(*neigh_requests),
        compare_tgt_request_global_id);
  free(send_neigh_global_ids);
  free(recv_neigh_global_ids);

  *neigh_requests_ = neigh_requests;
  *request_count_ = request_count;
}

static struct interp_result * init_interp_results(
  size_t * tgt_local_ids, yac_int * tgt_global_ids, size_t count) {

  struct interp_result * interp_results =
    xmalloc(count * sizeof(*interp_results));
  for (size_t i = 0; i < count; ++i) {
    interp_results[i].local_id = tgt_local_ids[i];
    interp_results[i].global_id = tgt_global_ids[i];
    interp_results[i].idx = i;
    interp_results[i].stencil.count = 0;
  }
  qsort(interp_results, count, sizeof(*interp_results),
        compare_interp_result_global_id);
  return interp_results;
}

static int compare_stencil_info(
  const void * a, const void * b) {

  int ret = ((struct stencil_info *)a)->rank -
            ((struct stencil_info *)b)->rank;
  if (ret) return ret;

  return (((struct stencil_info *)a)->idx >
          ((struct stencil_info *)b)->idx) -
         (((struct stencil_info *)a)->idx <
          ((struct stencil_info *)b)->idx);
}

static struct result_stencil copy_result_stencil_multi(
  struct result_stencil * neigh_stencils, size_t * stencil_indices,
  size_t count, yac_int global_id, double weight) {

  size_t stencil_info_count = 0;

  for (size_t i = 0; i < count; ++i)
    stencil_info_count += neigh_stencils[stencil_indices[i]].count;

  struct result_stencil stencil;

  struct stencil_info * stencil_infos;
  if (stencil_info_count > 1) {
    stencil_infos =
      ((stencil.data.multi =
          xmalloc(stencil_info_count * sizeof(*stencil_infos))));
  } else {
    stencil_infos = &(stencil.data.single);
  }

  stencil.global_id = global_id;
  for (size_t i = 0, j = 0; i < count; ++i) {
    struct result_stencil * curr_stencil =
      neigh_stencils + stencil_indices[i];
    size_t curr_stencil_info_count = curr_stencil->count;
    struct stencil_info * curr_stencil_infos =
      (curr_stencil_info_count == 1)?
        (&(curr_stencil->data.single)):(curr_stencil->data.multi);
    memcpy(stencil_infos + j, curr_stencil_infos,
          curr_stencil_info_count * sizeof(*stencil_infos));
    for (size_t k = 0; k < curr_stencil_info_count; ++k, ++j)
      stencil_infos[j].weight *= weight;
  }

  // in case we have multiple stencil infos, there may be duplicated
  // entries from different source
  // here we remove these duplicated entries
  if (stencil_info_count > 1) {

    // sort the stencils
    qsort(stencil_infos, stencil_info_count, sizeof(*stencil_infos),
          compare_stencil_info);

    // remote duplicated stencils
    struct stencil_info * prev_stencil_info = stencil_infos,
                        * curr_stencil_info = stencil_infos + 1;
    size_t new_stencil_info_count = 1;
    for (size_t i = 1; i < stencil_info_count; ++i, ++curr_stencil_info) {
      if (compare_stencil_info(
            curr_stencil_info, prev_stencil_info)) {
        if (new_stencil_info_count != i)
          stencil_infos[new_stencil_info_count] =
            *curr_stencil_info;
        ++new_stencil_info_count;
        prev_stencil_info = curr_stencil_info;
      } else {
        stencil_infos[new_stencil_info_count-1].weight +=
          curr_stencil_info->weight;
      }
    }
    if (new_stencil_info_count != stencil_info_count) {
      stencil_info_count = new_stencil_info_count;
      if (new_stencil_info_count == 1) {
        stencil.data.single = *stencil_infos;
        free(stencil_infos);
      } else {
        stencil.data.multi =
          xrealloc(
            stencil_infos, stencil_info_count * sizeof(*stencil_infos));
      }
    }
  }
  stencil.count = stencil_info_count;

  return stencil;
}

static size_t match_neigh_answers_with_tgts(
  struct result_stencils * neigh_answer,
  yac_int * neigh_global_ids, yac_int * neigh_to_tgt_global_id,
  size_t * stencil_indices, size_t * num_neighbours_,
  struct interp_result * interp_results, size_t * num_open_tgt_) {

  size_t num_neighbours = *num_neighbours_;
  size_t num_open_tgt = *num_open_tgt_;

  struct result_stencil * neigh_stencils = neigh_answer->stencils;
  size_t answer_count = neigh_answer->count;

  // match received neigh request answers with tgt neighbours
  // (remove duplicated results)
  qsort(
    neigh_stencils, answer_count, sizeof(*neigh_stencils),
    compare_result_stencil_global_id);
  size_t match_count = 0;
  for (size_t i = 0, j = 0; i < answer_count; ++i) {
    yac_int curr_global_id = neigh_stencils[i].global_id;

    while ((j < num_neighbours) &&
            (neigh_global_ids[j] < curr_global_id)) ++j;

    while ((j < num_neighbours) &&
            (neigh_global_ids[j] == curr_global_id)) {
      neigh_global_ids[j] = XT_INT_MAX;
      stencil_indices[j] = i;
      ++match_count;
      ++j;
    }
  }

  // move fullfilled neighbour requests to the end of the array
  yac_quicksort_index_yac_int_yac_int_size_t(
    neigh_global_ids, num_neighbours, neigh_to_tgt_global_id, stencil_indices);
  num_neighbours -= match_count;
  *num_neighbours_ = num_neighbours;

  neigh_to_tgt_global_id += num_neighbours;
  stencil_indices += num_neighbours;

  // sort matches by target global ids
  yac_quicksort_index_yac_int_size_t(
    neigh_to_tgt_global_id, match_count, stencil_indices);

  // set stencils for target matches
  for (size_t i = 0, k = 0; i < match_count;) {

    size_t prev_i = i;

    // count the number of stencils for the current target
    yac_int curr_tgt_global_id = neigh_to_tgt_global_id[i++];
    while ((i < match_count) &&
           (neigh_to_tgt_global_id[i] == curr_tgt_global_id)) ++i;
    size_t curr_stencil_count = i - prev_i;

    while ((k < num_open_tgt) &&
            (interp_results[k].global_id < curr_tgt_global_id)) ++k;

    while ((k < num_open_tgt) &&
            (interp_results[k].global_id == curr_tgt_global_id)) {

      interp_results[k].stencil =
        copy_result_stencil_multi(
          neigh_stencils, stencil_indices + prev_i, curr_stencil_count,
          curr_tgt_global_id, 1.0 / (double)curr_stencil_count);
      ++k;
    }
  }

  // move successfully interpolated target points to the end of
  // the array
  qsort(interp_results, num_open_tgt, sizeof(*interp_results),
        compare_interp_result_stencil);
  size_t new_num_open_tgt = 0;
  while ((new_num_open_tgt < num_open_tgt) &&
         (interp_results[new_num_open_tgt].stencil.count == 0))
    new_num_open_tgt++;
  *num_open_tgt_ = new_num_open_tgt;
  return num_open_tgt - new_num_open_tgt;
}

static struct comm_stuff init_comm_stuff(struct interp_grid * interp_grid) {

  struct comm_stuff comm;

  comm.comm = yac_interp_grid_get_MPI_Comm(interp_grid);
  yac_get_comm_buffers(
    1, &comm.sendcounts, &comm.recvcounts, &comm.sdispls, &comm.rdispls,
    comm.comm);
  yac_mpi_call(MPI_Comm_rank(comm.comm, &comm.rank), comm.comm);
  yac_mpi_call(MPI_Comm_size(comm.comm, &comm.size), comm.comm);
  comm.stencil_info_dt = yac_get_stencil_info_mpi_datatype(comm.comm);

  return comm;
}

static void free_comm_stuff(struct comm_stuff comm) {

  yac_free_comm_buffers(
    comm.sendcounts, comm.recvcounts, comm.sdispls, comm.rdispls);
  yac_mpi_call(MPI_Type_free(&comm.stencil_info_dt), comm.comm);
}

// the local process is the distributed owner for all target points
// passed to this routine
static void do_search_creep_2(
  struct interp_grid * interp_grid, int const max_creep_distance,
  size_t * tgt_points, yac_int * tgt_global_ids, size_t count,
  int * interp_flag, struct interp_weights * weights) {

  struct comm_stuff comm = init_comm_stuff(interp_grid);

  // get information about the neighbours of the target points
  size_t * neigh_local_ids;
  yac_int * neigh_global_ids;
  yac_int * neigh_to_tgt_global_id;
  size_t total_num_neighbours;
  get_tgt_neigh_info(
    interp_grid, tgt_points, tgt_global_ids, count,
    &neigh_local_ids, &neigh_global_ids, &neigh_to_tgt_global_id,
    &total_num_neighbours);
  size_t * stencil_indices =
    xmalloc(total_num_neighbours * sizeof(*stencil_indices));

  // send request for target points neighbours to the respective
  // distributed owners
  struct tgt_request * neigh_requests;
  size_t request_count;
  send_neigh_request(
    interp_grid, comm, neigh_local_ids, neigh_global_ids, total_num_neighbours,
    &neigh_requests, &request_count);
  free(neigh_local_ids);

  // sort global ids of target point neighbours
  yac_quicksort_index_yac_int_yac_int(
    neigh_global_ids, total_num_neighbours, neigh_to_tgt_global_id);

  // initialise interpolation results
  size_t num_open_tgt = count;
  struct interp_result * interp_results =
    init_interp_results(tgt_points, tgt_global_ids, count);

  // get already existing results and relocate them to their respective
  // distributed owners
  struct interp_result * initial_interp_results;
  size_t result_count;
  get_initial_results(
    interp_grid, comm, weights, &initial_interp_results, &result_count);

  for (int creep_distance = 0;
       creep_distance < max_creep_distance; ++creep_distance) {

    // check whether there are initial results on any process
    int result_flag = result_count > 0;
    yac_mpi_call(
      MPI_Allreduce(
        MPI_IN_PLACE, &result_flag, 1, MPI_INT, MPI_MAX, comm.comm), comm.comm);
    if (result_flag == 0) break;

    // relocate interpolation results to distributed owners of associated
    // target points
    struct result_stencils * interp_stencils =
      relocate_interp_results(
        interp_grid, comm,
        (creep_distance == 0)?
          initial_interp_results:(interp_results + num_open_tgt),
        result_count);
    if (creep_distance == 0) free(initial_interp_results);

    // match interpolation results with distributed neighbour requests and
    // inform origins of requests about matches
    struct result_stencils * neigh_matches =
      update_neigh_requests(
        interp_grid, comm, neigh_requests, &request_count,
        interp_stencils);

    // match received neigh request answers with tgts
    result_count =
      match_neigh_answers_with_tgts(
        neigh_matches, neigh_global_ids, neigh_to_tgt_global_id,
        stencil_indices, &total_num_neighbours, interp_results,
        &num_open_tgt);
    free(interp_stencils->stencils);
    free(interp_stencils);
    free(neigh_matches->stencils);
    free(neigh_matches);
  } // creep_distance < max_creep_distance

  free(neigh_requests);
  free(stencil_indices);
  free(neigh_global_ids);
  free(neigh_to_tgt_global_id);

  for (size_t i = 0; i < count; ++i)
    interp_flag[interp_results[i].idx] =
      interp_results[i].stencil.count > 0;

  // copy stencils
  struct remote_points interp_tgt_remote_points;
  size_t * num_stencils_per_tgt;
  int * stencil_ranks;
  double * w;
  extract_interp_info(
    interp_grid, comm, interp_results + num_open_tgt, count - num_open_tgt,
    &interp_tgt_remote_points, &num_stencils_per_tgt, &stencil_indices,
    &stencil_ranks, &w);
  yac_interp_weights_wcopy_weights(
    weights, &interp_tgt_remote_points, num_stencils_per_tgt,
    stencil_indices, stencil_ranks, w);
  free(interp_tgt_remote_points.data);
  free(w);
  free(num_stencils_per_tgt);
  free(stencil_ranks);
  free(stencil_indices);
  for (size_t i = num_open_tgt; i < count; ++i)
    if (interp_results[i].stencil.count > 1)
      free(interp_results[i].stencil.data.multi);
  free(interp_results);

  free_comm_stuff(comm);
}

static size_t do_search_creep (struct interp_method * method,
                               struct interp_grid * interp_grid,
                               size_t * tgt_points, size_t count,
                               struct interp_weights * weights) {

  struct interp_method_creep * creep_method =
    (struct interp_method_creep *)method;

  int const max_creep_distance = creep_method->max_creep_distance;

  if (max_creep_distance == 0) return 0;

  YAC_ASSERT(
    (yac_interp_grid_get_tgt_field_location(interp_grid) == CELL) ||
    (yac_interp_grid_get_tgt_field_location(interp_grid) == CORNER),
    "ERROR(do_search_creep): unsupported target field location "
    "(has to be CELL or CORNER)")

  MPI_Comm comm = yac_interp_grid_get_MPI_Comm(interp_grid);
  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);
  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // get the distributed owners for all tgts
  int * tgt_points_dist_owner =
    xmalloc(count * sizeof(*tgt_points_dist_owner));
  yac_interp_grid_determine_dist_tgt_owners(
    interp_grid, tgt_points, count, tgt_points_dist_owner);
  yac_int * tgt_points_global_ids =
    xmalloc(count * sizeof(*tgt_points_global_ids));
  yac_interp_grid_get_tgt_global_ids(
    interp_grid, tgt_points, count, tgt_points_global_ids);

  // determine which of these points are local
  size_t local_count = 0;
  for (size_t i = 0; i < count; ++i) {
    if (tgt_points_dist_owner[i] == comm_rank) {
      local_count++;
      tgt_points_dist_owner[i] = INT_MAX;
    }
  }
  size_t send_count = count - local_count;

  // sort target points (remote points to the start of the array)
  yac_quicksort_index_int_size_t_yac_int(
    tgt_points_dist_owner, count, tgt_points, tgt_points_global_ids);

  // relocate tgt_points to distributed owners
  for (size_t i = 0; i < send_count; ++i)
    sendcounts[tgt_points_dist_owner[i]]++;
  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);
  size_t recv_count = recvcounts[comm_size-1] + rdispls[comm_size-1];
  yac_int * global_ids_buffer =
    xrealloc(
      tgt_points_global_ids,
      (send_count + local_count + recv_count) * sizeof(*global_ids_buffer));
  yac_int * send_global_ids = global_ids_buffer;
  yac_int * recv_global_ids = global_ids_buffer + send_count;
  yac_alltoallv_yac_int_p2p(
    send_global_ids, sendcounts, sdispls+1,
    recv_global_ids + local_count, recvcounts, rdispls, comm);

  // get local ids for received global ids
  size_t * temp_tgt_points =
    xmalloc((local_count + recv_count) * sizeof(*temp_tgt_points));
  memcpy(temp_tgt_points, tgt_points + send_count,
         local_count * sizeof(*tgt_points));
  yac_interp_grid_tgt_global_to_local(
    interp_grid, recv_global_ids + local_count, recv_count,
    temp_tgt_points + local_count);

  // do the interpolation for the redistributed target points
  int * interp_flag_buffer =
    xmalloc((send_count + local_count + recv_count) *
            sizeof(*interp_flag_buffer));
  int * temp_interp_flag = interp_flag_buffer + send_count;
  int * interp_flag = interp_flag_buffer;
  memset(temp_interp_flag, 0, (local_count + recv_count) * sizeof(*temp_interp_flag));
  do_search_creep_2(
    interp_grid, max_creep_distance, temp_tgt_points, global_ids_buffer + send_count,
    local_count + recv_count, temp_interp_flag, weights);
  free(global_ids_buffer);
  free(temp_tgt_points);

  // relocate interp_flag
  yac_alltoallv_int_p2p(
    temp_interp_flag + local_count, recvcounts, rdispls,
    interp_flag, sendcounts, sdispls+1, comm);

  // count number of points that can be interpolated and reorder
  // tgt_points accordingly (interpolated first)
  size_t num_interpolated_tgt = 0;
  for (size_t i = 0; i < count; ++i) {
    if (interp_flag[i]) {
      interp_flag[i] = 0;
      ++num_interpolated_tgt;
    } else {
      interp_flag[i] = 1;
    }
  }
  yac_quicksort_index_int_size_t(interp_flag, count, tgt_points);
  free(interp_flag_buffer);
  free(tgt_points_dist_owner);

  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  return num_interpolated_tgt;
}

struct interp_method * yac_interp_method_creep_new(int creep_distance) {

  struct interp_method_creep * method_creep =
    xmalloc(1 * sizeof(*method_creep));

  method_creep->vtable = &interp_method_creep_vtable;
  method_creep->max_creep_distance =
    (creep_distance >= 0)?creep_distance:INT_MAX;

  return (struct interp_method*)method_creep;
}

static void delete_creep(struct interp_method * method) {
  free(method);
}
