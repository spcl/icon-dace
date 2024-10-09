// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>
#include "yac_mpi_internal.h"
#include "io_utils.h"
#include "math.h"

#define IO_RANK_LIST_STR "YAC_IO_RANK_LIST"
#define IO_MAX_NUM_RANKS_STR "YAC_IO_MAX_NUM_RANKS"
#define IO_RANK_EXCLUDE_LIST_STR "YAC_IO_RANK_EXCLUDE_LIST"
#define IO_MAX_NUM_RANKS_PER_NODE "YAC_IO_MAX_NUM_RANKS_PER_NODE"
#define DEFAULT_MAX_NUM_IO_RANK_PER_NODE (1);

static inline int compare_int(const void * a, const void * b) {

  int const * a_ = a, * b_ = b;

  return (*a_ > *b_) - (*b_ > *a_);
}

static void read_rank_list(
  char const * env_name, MPI_Comm comm, size_t * num_ranks_, int ** ranks_) {

  int rank;
  yac_mpi_call(MPI_Comm_rank(comm, &rank), comm);

  size_t num_ranks = 0;
  int * ranks = NULL;

  // environment is only checked on rank 0, results are broadcasted
  // to other processes
  if (rank == 0) {

    // check whether the user provided a list of ranks for IO
    char * rank_list_str = getenv(env_name);
    if ((rank_list_str != NULL) && (rank_list_str[0] != '\0')) {

      int temp_num_ranks = 1;
      for (char * curr_char = rank_list_str + 1; *curr_char != '\0'; ++curr_char)
        if (*curr_char == ',') ++temp_num_ranks;

      char * rank_list_copy = strdup(rank_list_str);
      num_ranks = 0;
      int * world_ranks = xmalloc(temp_num_ranks * sizeof(*world_ranks));

      // parse rank list
      char * rank_str = strtok(rank_list_copy, ",");
      while (rank_str != NULL) {
        int curr_rank = atoi(rank_str);
        YAC_ASSERT_F(
          curr_rank >= 0, "ERROR(read_rank_list): \"%s\" is not a valid rank",
          rank_str);
        world_ranks[num_ranks++] = curr_rank;
        rank_str = strtok(NULL, ",");
      }
      free(rank_list_copy);

      // sort ranks
      qsort(world_ranks, num_ranks, sizeof(*world_ranks), compare_int);

      // remove duplicated ranks
      yac_remove_duplicates_int(world_ranks, &num_ranks);

      // translate rank list from MPI_COMM_WORLD to comm
      ranks = xmalloc(num_ranks * sizeof(ranks));
      MPI_Group world_group, comm_group;
      yac_mpi_call(
        MPI_Comm_group(MPI_COMM_WORLD, &world_group), MPI_COMM_WORLD);
      yac_mpi_call(MPI_Comm_group(comm, &comm_group), comm);
      yac_mpi_call(
        MPI_Group_translate_ranks(
          world_group, (int)num_ranks, world_ranks,
          comm_group, ranks), MPI_COMM_WORLD);
      yac_mpi_call(MPI_Group_free(&comm_group), comm);
      yac_mpi_call(MPI_Group_free(&world_group), MPI_COMM_WORLD);
      free(world_ranks);

      // remove ranks, which are not avaible in comm
      size_t new_num_ranks = 0;
      for (size_t i = 0; i < num_ranks; ++i) {
        if (ranks[i] != MPI_UNDEFINED) {
          if (i != new_num_ranks)
            ranks[new_num_ranks] = ranks[i];
          ++new_num_ranks;
        }
      }
      num_ranks = new_num_ranks;

      // sort ranks
      qsort(ranks, num_ranks, sizeof(*ranks), compare_int);
    }

    // broadcast ranks
    yac_mpi_call(
      MPI_Bcast(&num_ranks, 1, YAC_MPI_SIZE_T, 0, comm), comm);
    if (num_ranks > 0) {
      ranks = xrealloc(ranks, num_ranks * sizeof(*ranks));
      yac_mpi_call(
        MPI_Bcast(ranks, (int)num_ranks, MPI_INT, 0, comm), comm);
    } else {
      free(ranks);
    }
  } else {

    // receive ranks from root
    yac_mpi_call(
      MPI_Bcast(&num_ranks, 1, YAC_MPI_SIZE_T, 0, comm), comm);
    if (num_ranks > 0) {
      ranks = xmalloc(num_ranks * sizeof(*ranks));
      yac_mpi_call(
        MPI_Bcast(ranks, (int)num_ranks, MPI_INT, 0, comm), comm);
    }
  }
  *num_ranks_ = num_ranks;
  *ranks_ = ranks;
}

static void read_io_rank_list(
  MPI_Comm comm, size_t * num_io_ranks, int ** io_ranks) {

  read_rank_list(IO_RANK_LIST_STR, comm, num_io_ranks, io_ranks);
}

static void check_io_max_num_ranks_per_node(
  MPI_Comm comm, size_t * num_io_ranks_, int ** io_ranks_) {

  int rank, size;
  yac_mpi_call(MPI_Comm_rank(comm, &rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &size), comm);

  int max_num_io_rank_per_node;

  // environment is only checked on rank 0, results are broadcasted
  // to other processes
  if (rank == 0) {

    // check whether the user provided a maximum number of ranks per node
    char * max_num_io_rank_per_node_str = getenv(IO_MAX_NUM_RANKS_PER_NODE);
    if ((max_num_io_rank_per_node_str != NULL) &&
        (max_num_io_rank_per_node_str[0] != '\0')) {

      max_num_io_rank_per_node = atoi(max_num_io_rank_per_node_str);
      YAC_ASSERT_F(
        (max_num_io_rank_per_node > 0) || (max_num_io_rank_per_node == -1),
        "ERROR(check_io_max_num_ranks_per_node): "
        "\"%s\" is not a valid value for the maximum number of io ranks "
        "per node", max_num_io_rank_per_node_str);
    } else {
      max_num_io_rank_per_node = DEFAULT_MAX_NUM_IO_RANK_PER_NODE;
    }
    yac_mpi_call(
      MPI_Bcast(&max_num_io_rank_per_node, 1, MPI_INT, 0, comm), comm);
  } else {
    yac_mpi_call(
      MPI_Bcast(&max_num_io_rank_per_node, 1, MPI_INT, 0, comm), comm);
  }

  // if there is no limit on the number of io ranks per node
  if (max_num_io_rank_per_node == -1) return;

  // create one communicator per node
  MPI_Comm node_comm;
  yac_mpi_call(
    MPI_Comm_split_type(
      comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node_comm), comm);

  // determine rank within node comm
  int node_rank;
  yac_mpi_call(MPI_Comm_rank(node_comm, &node_rank), node_comm);

  // translate io ranks to node ranks
  int * node_ranks = xmalloc(*num_io_ranks_ * sizeof(*node_ranks));
  MPI_Group io_group, node_group;
  yac_mpi_call(MPI_Comm_group(comm, &io_group), comm);
  yac_mpi_call(MPI_Comm_group(node_comm, &node_group), node_comm);
  yac_mpi_call(
    MPI_Group_translate_ranks(
      io_group, (int)*num_io_ranks_, *io_ranks_,
      node_group, node_ranks), comm);
  yac_mpi_call(MPI_Group_free(&node_group), node_comm);
  yac_mpi_call(MPI_Group_free(&io_group), comm);
  yac_mpi_call(MPI_Comm_free(&node_comm), comm);

  // remove MPI_UNDEFINED entries
  size_t num_node_io_ranks = 0;
  for (size_t i = 0; i < *num_io_ranks_; ++i) {
    if (node_ranks[i] != MPI_UNDEFINED) {
      if (i != num_node_io_ranks)
        node_ranks[num_node_io_ranks] = node_ranks[i];
      ++num_node_io_ranks;
    }
  }

  // sort ranks
  qsort(node_ranks, num_node_io_ranks, sizeof(*node_ranks), compare_int);

  // check whether local process is within the first n io ranks on its node
  int local_is_io_rank = 0;
  for (size_t i = 0;
       (i < num_node_io_ranks) &&
       (i < (size_t)max_num_io_rank_per_node) && !local_is_io_rank; ++i)
    if (node_ranks[i] == node_rank) local_is_io_rank = 1;
  free(node_ranks);

  // determine io ranks on all nodes
  int * is_io_rank = xmalloc((size_t)size * sizeof(*is_io_rank));
  yac_mpi_call(
    MPI_Allgather(
      &local_is_io_rank, 1, MPI_INT, is_io_rank, 1, MPI_INT, comm), comm);

  // compress is_io_rank into io_ranks
  size_t num_io_ranks = 0;
  for (int i = 0; i < size; ++i)
    if (is_io_rank[i]) is_io_rank[num_io_ranks++] = i;

  free(*io_ranks_);
  *num_io_ranks_ = num_io_ranks;
  *io_ranks_ = xrealloc(is_io_rank, num_io_ranks * sizeof(**io_ranks_));
}

static void check_io_max_num_ranks(
  MPI_Comm comm, size_t * num_io_ranks, int ** io_ranks) {

  int rank;
  yac_mpi_call(MPI_Comm_rank(comm, &rank), comm);

  int max_num_ranks = -1;

  // environment is only checked on rank 0, results are broadcasted
  // to other processes
  if (rank == 0) {

    // check whether the user provided a maximum number of io ranks
    char * max_num_ranks_str = getenv(IO_MAX_NUM_RANKS_STR);
    if ((max_num_ranks_str != NULL) && (max_num_ranks_str[0] != '\0')) {
      max_num_ranks = atoi(max_num_ranks_str);
      YAC_ASSERT_F(
        (max_num_ranks > 0) || (max_num_ranks == -1),
        "ERROR(check_io_max_num_ranks): "
        "\"%s\" is not a valid value for the maximum number of io ranks",
        max_num_ranks_str);
    }

    yac_mpi_call(
      MPI_Bcast(&max_num_ranks, 1, MPI_INT, 0, comm), comm);
  } else {
    yac_mpi_call(
      MPI_Bcast(&max_num_ranks, 1, MPI_INT, 0, comm), comm);
  }

  // if there is no limit on the number of io ranks
  if (max_num_ranks == -1) return;

  if ((max_num_ranks > 0) && ((size_t)max_num_ranks < *num_io_ranks)) {
    *io_ranks =
      xrealloc(*io_ranks, (size_t)max_num_ranks * sizeof(**io_ranks));
    *num_io_ranks = (size_t)max_num_ranks;
  }
}

static void check_io_rank_exclude_list(
  MPI_Comm comm, size_t * num_io_ranks, int ** io_ranks) {

  // read in the list of ranks which are not to be used for io
  size_t num_io_ranks_excluded;
  int * io_ranks_excluded;
  read_rank_list(
    IO_RANK_EXCLUDE_LIST_STR, comm, &num_io_ranks_excluded, &io_ranks_excluded);

  // if there are ranks to be excluded
  if (num_io_ranks_excluded > 0) {

    // sort ranks
    qsort(*io_ranks, *num_io_ranks, sizeof(**io_ranks), compare_int);

    // match exclude list with io rank list
    size_t new_num_io_ranks = 0;
    for (size_t i = 0, j = 0; i < *num_io_ranks; ++i) {

      while ((j < num_io_ranks_excluded) &&
             (io_ranks_excluded[j] < (*io_ranks)[i])) ++j;

      if ((j >= num_io_ranks_excluded) ||
          ((*io_ranks)[i] != io_ranks_excluded[j])) {

        if (i != new_num_io_ranks)
          (*io_ranks)[new_num_io_ranks] = (*io_ranks)[i];
        ++new_num_io_ranks;
      }
    }

    if (new_num_io_ranks != *num_io_ranks) {
      *io_ranks = xrealloc(*io_ranks, new_num_io_ranks * sizeof(**io_ranks));
      *num_io_ranks = new_num_io_ranks;
    }
  }

  free(io_ranks_excluded);
}

void yac_get_io_ranks(
  MPI_Comm comm, int * local_is_io_, int ** io_ranks_, int * num_io_ranks_) {

  int rank, size;
  yac_mpi_call(MPI_Comm_rank(comm, &rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &size), comm);

  size_t num_io_ranks = 0;
  int * io_ranks = NULL;

  // check environment for io rank list
  read_io_rank_list(comm, &num_io_ranks, &io_ranks);

  // if no rank list was provided -> generate default rank list (all processes)
  if (num_io_ranks == 0) {
    num_io_ranks = (size_t)size;
    io_ranks = xmalloc(num_io_ranks * sizeof(*io_ranks));
    for (int i = 0; i < size; ++i) io_ranks[i] = i;
  }

  // check whether we have to exclude some ranks
  check_io_rank_exclude_list(comm, &num_io_ranks, &io_ranks);

  // check for the maximum number of io ranks per node
  check_io_max_num_ranks_per_node(comm, &num_io_ranks, &io_ranks);

  // check maximum number of io ranks
  check_io_max_num_ranks(comm, &num_io_ranks, &io_ranks);

  YAC_ASSERT(
    num_io_ranks > 0, "ERROR(yac_get_io_ranks): could not determine io ranks");

  int local_is_io = 0;
  for (size_t i = 0; (i < num_io_ranks) && !local_is_io; ++i)
    if (io_ranks[i] == rank) local_is_io = 1;

  *local_is_io_ = local_is_io;
  *io_ranks_ = io_ranks;
  *num_io_ranks_ = (int)num_io_ranks;
}

void yac_nc_open(const char * path, int omode, int * ncidp) {

#ifndef YAC_NETCDF_ENABLED

  UNUSED(path);
  UNUSED(omode);
  UNUSED(ncidp);
  die("ERROR(yac_nc_open): YAC is built without the NetCDF support");
#else

  YAC_ASSERT_F(
    yac_file_exists(path),
    "ERROR(yac_nc_open): file \"%s\" does not exist", path);
  YAC_HANDLE_ERROR(nc_open(path, omode, ncidp));
#endif
}

void yac_nc_create(const char * path, int cmode, int * ncidp) {

#ifndef YAC_NETCDF_ENABLED

  UNUSED(path);
  UNUSED(cmode);
  UNUSED(ncidp);
  die("ERROR(yac_nc_create): YAC is built without the NetCDF support");
#else

  int status = nc_create(path, cmode, ncidp);
  YAC_ASSERT_F(
    status == NC_NOERR,
    "ERROR(yac_nc_create): failed to create file \"%s\" "
    "(NetCDF error message: \"%s\")", path, nc_strerror(status));
#endif
}

void yac_nc_inq_dimid(int ncid, char const * name, int * dimidp) {

#ifndef YAC_NETCDF_ENABLED

  UNUSED(ncid);
  UNUSED(name);
  UNUSED(dimidp);
  die("ERROR(yac_nc_inq_dimid): YAC is built without the NetCDF support");
#else

  int status = nc_inq_dimid(ncid, name, dimidp);

  if (status == NC_EBADDIM) {
    // GCOVR_EXCL_START
    size_t pathlen;
    YAC_HANDLE_ERROR(nc_inq_path(ncid, &pathlen, NULL));
    char * path = xmalloc(pathlen * sizeof(*path));
    YAC_HANDLE_ERROR(nc_inq_path(ncid, NULL, path));
    YAC_ASSERT_F(
      0, "ERROR(yac_nc_inq_dimid): "
      "dimension \"%s\" could not be found in file \"%s\"", name, path);
    // GCOVR_EXCL_STOP
  } else YAC_HANDLE_ERROR(status);
#endif
}

void yac_nc_inq_varid(int ncid, char const * name, int * varidp) {

#ifndef YAC_NETCDF_ENABLED

  UNUSED(ncid);
  UNUSED(name);
  UNUSED(varidp);
  die("ERROR(yac_nc_inq_varid): YAC is built without the NetCDF support");
#else

  int status = nc_inq_varid(ncid, name, varidp);

  if (status == NC_ENOTVAR) {
    // GCOVR_EXCL_START
    size_t pathlen;
    YAC_HANDLE_ERROR(nc_inq_path(ncid, &pathlen, NULL));
    char * path = xmalloc(pathlen * sizeof(*path));
    YAC_HANDLE_ERROR(nc_inq_path(ncid, NULL, path));
    YAC_ASSERT_F(
      0, "ERROR(yac_nc_inq_varid): "
      "variable \"%s\" could not be found in file \"%s\"", name, path);
    // GCOVR_EXCL_STOP
  } else YAC_HANDLE_ERROR(status);
#endif
}
