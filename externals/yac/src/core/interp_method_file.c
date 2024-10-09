// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>

#include "interp_method_internal.h"
#include "interp_method_file.h"
#include "dist_grid.h"
#include "yac_mpi_internal.h"
#include "io_utils.h"

static size_t do_search_file(struct interp_method * method,
                              struct yac_interp_grid * interp_grid,
                              size_t * tgt_points, size_t count,
                              struct yac_interp_weights * weights);
static void delete_file(struct interp_method * method);

static struct interp_method_vtable
  interp_method_file_vtable = {
    .do_search = do_search_file,
    .delete = delete_file};

struct interp_method_file {

  struct interp_method_vtable * vtable;
  char const * weight_file_name;
};

struct link_data {
  union {
    yac_int global_id;
    size_t local_id;
  } src, tgt;
  double weight;
  size_t src_field_idx;
};

struct fixed_data {
  double value;
  struct {
    yac_int global_id;
  } * tgt;
  size_t num_tgt;
};

struct temp_fixed_data {
  double value;
  struct {
    yac_int global_id;
  } tgt;
};

struct temp_result {
  enum {
    LINK = 0,
    FIXED = 1,
    NO_INTERP = 2,
  } type;
  union {
    struct {
      struct link_data * links;
      size_t count;
    } link;
    struct {
      double value;
    } fixed;
  } data;
  struct {
    yac_int global_id;
    size_t local_id;
  } tgt;
  int owner;
  size_t reorder_idx;
};

static void read_weight_file(
  char const * weight_file_name, MPI_Comm comm,
  char const * src_grid_name, char const * tgt_grid_name,
  enum yac_location * src_locations, enum yac_location tgt_location,
  size_t num_src_fields,
  struct link_data ** links, size_t * num_links,
  struct fixed_data ** fixed, size_t * num_fixed) {

#ifndef YAC_NETCDF_ENABLED

  UNUSED(weight_file_name);
  UNUSED(comm);
  UNUSED(src_grid_name);
  UNUSED(tgt_grid_name);
  UNUSED(src_locations);
  UNUSED(tgt_location);
  UNUSED(num_src_fields);
  UNUSED(links);
  UNUSED(num_links);
  UNUSED(fixed);
  UNUSED(num_fixed);

  die("ERROR(read_weight_file): YAC is built without the NetCDF support");
#else

  *links = NULL;
  *num_links = 0;
  *fixed = NULL;
  *num_fixed = 0;

  int rank;
  yac_mpi_call(MPI_Comm_rank(comm, &rank), comm);

  // determine ranks that do input
  int local_is_io;
  int * io_ranks;
  int num_io_ranks;
  yac_get_io_ranks(comm, &local_is_io, &io_ranks, &num_io_ranks);

  int io_idx = INT_MAX;
  if (local_is_io)
    for (int i = 0; (i < num_io_ranks) && (io_idx == INT_MAX); ++i)
      if (io_ranks[i] == rank) io_idx = i;
  free(io_ranks);

  // if the local process does not have to read anything from the file
  if (!local_is_io) return;

  // open weight file on io processes
  int ncid;
  yac_nc_open(weight_file_name, NC_NOWRITE, &ncid);

  //\todo use parallel netcdf if available

  int dimid, var_id;

  // global attributes
  size_t version_len;
  char * version;
  YAC_HANDLE_ERROR(nc_inq_attlen(ncid, NC_GLOBAL, "version", &version_len));
  version = xmalloc(version_len + 1);
  version[version_len] = '\0';
  YAC_HANDLE_ERROR(nc_get_att_text(ncid, NC_GLOBAL, "version", version));
  YAC_ASSERT_F(
    (strlen(YAC_WEIGHT_FILE_VERSION_STRING) == version_len) &&
    !strncmp(version, YAC_WEIGHT_FILE_VERSION_STRING, version_len),
    "ERROR(read_weight_file): "
    "version string from weight file (\"%s\") does not match "
    "with YAC version \"%s\"", version, YAC_WEIGHT_FILE_VERSION_STRING)
  free(version);

  size_t grid_name_len;
  char * grid_name;
  YAC_HANDLE_ERROR(nc_inq_attlen(ncid, NC_GLOBAL, "src_grid_name", &grid_name_len));
  grid_name = xmalloc(grid_name_len + 1);
  grid_name[grid_name_len] = '\0';
  YAC_HANDLE_ERROR(nc_get_att_text(ncid, NC_GLOBAL, "src_grid_name", grid_name));
  YAC_ASSERT_F(
    (strlen(src_grid_name) == grid_name_len) &&
    !strncmp(src_grid_name, grid_name, grid_name_len),
    "ERROR(read_weight_file): source grid name from weight file (\"%s\") "
    "does not match with the one provided through the interface (\"%s\")",
    grid_name, src_grid_name)

  free(grid_name);

  YAC_HANDLE_ERROR(nc_inq_attlen(ncid, NC_GLOBAL, "dst_grid_name", &grid_name_len));
  grid_name = xmalloc(grid_name_len + 1);
  grid_name[grid_name_len] = '\0';
  YAC_HANDLE_ERROR(nc_get_att_text(ncid, NC_GLOBAL, "dst_grid_name", grid_name));
  YAC_ASSERT_F(
    (strlen(tgt_grid_name) == grid_name_len) &&
    !strncmp(tgt_grid_name, grid_name, grid_name_len),
    "ERROR(read_weight_file): target grid name from weight file (\"%s\") "
    "does not match with the one provided through the interface (\"%s\")",
    grid_name, tgt_grid_name)

  free(grid_name);

  //---------------
  // read link data
  //---------------

  size_t str_link_len;
  char * str_link;
  var_id = NC_GLOBAL;
  YAC_HANDLE_ERROR(nc_inq_attlen(ncid, var_id, "contains_links", &str_link_len));
  str_link = xmalloc(str_link_len + 1);
  str_link[str_link_len] = '\0';
  YAC_HANDLE_ERROR(nc_get_att_text(ncid, var_id, "contains_links", str_link));

  int contains_links = (strlen("TRUE") == str_link_len) &&
                       !strncmp("TRUE", str_link, str_link_len);
  YAC_ASSERT(
    contains_links ||
    ((strlen("FALSE") == str_link_len) &&
     !strncmp("FALSE", str_link, str_link_len)),
    "ERROR(read_weight_file): invalid global attribute contains_links")
  free(str_link);

  size_t num_wgts = 0;
  if (contains_links) {
    yac_nc_inq_dimid(ncid, "num_wgts", &dimid);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_wgts));
    YAC_ASSERT(
      num_wgts == 1, "ERROR(read_weight_file): YAC only supports num_wgts == 1")
  }

  // get number of links from file

  size_t num_links_in_file = 0;
  if (contains_links) {
    yac_nc_inq_dimid(ncid, "num_links", &dimid);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_links_in_file));
    YAC_ASSERT(
      num_links_in_file != 0, "ERROR(read_weight_file): no links defined")
  }

  // get number of source fields
  size_t tmp_num_src_fields = 0;
  if (contains_links) {
    yac_nc_inq_dimid(ncid, "num_src_fields", &dimid);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &tmp_num_src_fields));
    YAC_ASSERT(
      tmp_num_src_fields != 0,
      "ERROR(read_weight_file): no source fields in file")
    YAC_ASSERT_F(
      tmp_num_src_fields == num_src_fields,
      "ERROR(read_weight_file): number of source fields in file (%zu) does not "
      "match with the number provided through the interface (%zu)",
      tmp_num_src_fields, num_src_fields)
  }

  // get max location string length from file
  size_t max_loc_str_len;
  yac_nc_inq_dimid(ncid, "max_loc_str_len", &dimid);
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &max_loc_str_len));
  YAC_ASSERT(
    max_loc_str_len == YAC_MAX_LOC_STR_LEN,
    "ERROR(read_weight_file): wrong max location string length in weight file")

  // get source locations
  enum yac_location * tmp_src_locations =
    xmalloc(num_src_fields * sizeof(*tmp_src_locations));
  yac_nc_inq_varid(ncid, "src_locations", &var_id);

  for (size_t i = 0; i < num_src_fields; ++i) {

    char loc_str[YAC_MAX_LOC_STR_LEN];

    size_t str_start[2] = {i, 0};
    size_t str_count[2] = {1, YAC_MAX_LOC_STR_LEN};

    YAC_HANDLE_ERROR(
      nc_get_vara_text(ncid, var_id, str_start, str_count, loc_str));

    tmp_src_locations[i] = yac_str2loc(loc_str);

    YAC_ASSERT_F(
      tmp_src_locations[i] == src_locations[i],
      "ERROR(read_weight_file): source locations in file does not match with "
      "the locations provided through the interface\n"
      "location index:          %d\n"
      "location in weight file: %s\n"
      "location from interface: %s",
      (int)i, loc_str, yac_loc2str(src_locations[i]))
  }

  free(tmp_src_locations);

  // get target location
  yac_nc_inq_varid(ncid, "dst_location", &var_id);
  {
    char loc_str[YAC_MAX_LOC_STR_LEN];
    YAC_HANDLE_ERROR(nc_get_var_text(ncid, var_id, loc_str));
    YAC_ASSERT(
      tgt_location == yac_str2loc(loc_str),
      "ERROR(read_weight_file): target locations in file does not match with "
      "the locations provided through the interface")
  }

  if (contains_links) {

    // get number of links per source field
    unsigned * num_links_per_src_field =
      xmalloc(num_src_fields * sizeof(*num_links_per_src_field));
    yac_nc_inq_varid(ncid, "num_links_per_src_field", &var_id);
    YAC_HANDLE_ERROR(nc_get_var_uint(ncid, var_id, num_links_per_src_field));

    size_t offset =
      (size_t)(((long long)io_idx * (long long)num_links_in_file) /
               (long long)num_io_ranks);
    size_t count =
      (size_t)((((long long)io_idx + 1) * (long long)num_links_in_file) /
               (long long)num_io_ranks) - offset;

    *num_links = count;

    *links = xmalloc(count * sizeof(**links));

    if (count > 0) {
      size_t src_field_idx = 0, src_field_start_offset = 0;
      while ((src_field_idx < num_src_fields) &&
             (src_field_start_offset +
              (size_t)num_links_per_src_field[src_field_idx] <
              offset)) {

        src_field_start_offset += (size_t)num_links_per_src_field[src_field_idx];
        ++src_field_idx;
      };
      YAC_ASSERT(
        src_field_idx < num_src_fields,
        "ERROR(read_weight_file): problem in num_links_per_src_field")
      size_t src_field_offset = offset - src_field_start_offset;
      for (size_t k = 0; (src_field_idx < num_src_fields) && (k < count);
           src_field_offset = 0, ++src_field_idx) {

        for (; (src_field_offset <
                (size_t)num_links_per_src_field[src_field_idx]) &&
               (k < count); ++src_field_offset, ++k) {
          (*links)[k].src_field_idx = src_field_idx;
        }
      }
    }
    free(num_links_per_src_field);

    // get links
    int * address = xmalloc(count * sizeof(*address));
    yac_nc_inq_varid(ncid, "src_address", &var_id);
    YAC_HANDLE_ERROR(nc_get_vara_int(ncid, var_id, &offset, &count, address));
    for (size_t i = 0; i < count; ++i)
      (*links)[i].src.global_id = (yac_int)(address[i] - 1);

    yac_nc_inq_varid(ncid, "dst_address", &var_id);
    YAC_HANDLE_ERROR(nc_get_vara_int(ncid, var_id, &offset, &count, address));
    for (size_t i = 0; i < count; ++i)
      (*links)[i].tgt.global_id = (yac_int)(address[i] - 1);
    free(address);

    double * weights = xmalloc(count * sizeof(*weights));

    yac_nc_inq_varid(ncid, "remap_matrix", &var_id);
    size_t offsets[2] = {offset, 0};
    size_t counts[2] = {count, 1};
    YAC_HANDLE_ERROR(nc_get_vara_double(ncid, var_id, offsets, counts, weights));
    for (size_t i = 0; i < count; ++i)
      (*links)[i].weight = weights[i];

    free(weights);
  }

  //----------------
  // read fixed data
  //----------------

  // global attributes
  size_t str_fixed_len;
  char * str_fixed;
  var_id = NC_GLOBAL;
  YAC_HANDLE_ERROR(
    nc_inq_attlen(ncid, var_id, "contains_fixed_dst", &str_fixed_len));
  str_fixed = xmalloc(str_fixed_len + 1);
  str_fixed[str_fixed_len] = '\0';
  YAC_HANDLE_ERROR(nc_get_att_text(ncid, var_id, "contains_fixed_dst", str_fixed));

  int contains_fixed = (strlen("TRUE") == str_fixed_len) &&
                       !strncmp("TRUE", str_fixed, str_fixed_len);

  YAC_ASSERT(
    contains_fixed ||
    ((strlen("FALSE") == str_fixed_len) &&
     !strncmp("FALSE", str_fixed, str_fixed_len)),
    "ERROR(read_weight_file): invalid global attribute contains_fixed_dst")
  free(str_fixed);

  if (contains_fixed) {

    // get number of fixed values
    size_t num_fixed_values;
    yac_nc_inq_dimid(ncid, "num_fixed_values", &dimid);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_fixed_values));

    // get number of fixed target points
    size_t num_fixed_tgt;
    yac_nc_inq_dimid(ncid, "num_fixed_dst", &dimid);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_fixed_tgt));

    size_t offset =
      (size_t)(((long long)io_idx * (long long)num_fixed_tgt) /
               (long long)num_io_ranks);
    size_t count =
      (size_t)((((long long)io_idx + 1) * (long long)num_fixed_tgt) /
               (long long)num_io_ranks) - offset;

    *fixed = xmalloc(num_fixed_values * sizeof(**fixed));
    *num_fixed = num_fixed_values;

    double * fixed_values = xmalloc(num_fixed_values * sizeof(*fixed_values));
    int * num_tgt_indices_per_fixed_value =
      xmalloc(num_fixed_values * sizeof(*num_tgt_indices_per_fixed_value));
    int * global_tgt_indices = xmalloc(count * sizeof(*global_tgt_indices));

    // get number of fixed target points per fixed value
    yac_nc_inq_varid(ncid, "num_dst_per_fixed_value", &var_id);
    YAC_HANDLE_ERROR(nc_get_var_int(ncid, var_id, num_tgt_indices_per_fixed_value));

    // read data
    yac_nc_inq_varid(ncid, "fixed_values", &var_id);
    YAC_HANDLE_ERROR(nc_get_var_double(ncid, var_id, fixed_values));

    yac_nc_inq_varid(ncid, "dst_address_fixed", &var_id);
    YAC_HANDLE_ERROR(
      nc_get_vara_int(ncid, var_id, &offset, &count, global_tgt_indices));

    if (count > 0) {
      void * tgt_global_id_buffer = xmalloc(count * sizeof(*((*fixed)->tgt)));
      size_t fixed_value_idx = 0, fixed_value_start_offset = 0;
      while ((fixed_value_idx < num_fixed_values) &&
             (fixed_value_start_offset +
              (size_t)num_tgt_indices_per_fixed_value[fixed_value_idx] <
              offset)) {

        fixed_value_start_offset +=
          (size_t)num_tgt_indices_per_fixed_value[fixed_value_idx];
        (*fixed)[fixed_value_idx].value = fixed_values[fixed_value_idx];
        (*fixed)[fixed_value_idx].num_tgt = 0;
        (*fixed)[fixed_value_idx].tgt = tgt_global_id_buffer;
        ++fixed_value_idx;
      };
      YAC_ASSERT(
        fixed_value_idx < num_fixed_values,
        "ERROR(read_weight_file): problem in num_tgt_indices_per_fixed_value")
      size_t fixed_value_offset = offset - fixed_value_start_offset;
      for (size_t k = 0; (fixed_value_idx < num_fixed_values) && (k < count);
           fixed_value_offset = 0, ++fixed_value_idx) {

        size_t curr_tgt_count =
          MIN((size_t)num_tgt_indices_per_fixed_value[fixed_value_idx] -
              fixed_value_offset, count - k);

        (*fixed)[fixed_value_idx].value = fixed_values[fixed_value_idx];
        (*fixed)[fixed_value_idx].num_tgt = curr_tgt_count;
        (*fixed)[fixed_value_idx].tgt =
          (void*)(((unsigned char *)tgt_global_id_buffer) +
                  k * sizeof(*((*fixed)->tgt)));

        for (size_t i = 0; i < curr_tgt_count; ++fixed_value_offset, ++k, ++i)
          (*fixed)[fixed_value_idx].tgt[i].global_id =
            (yac_int)(global_tgt_indices[k] - 1);
      }
      for (; fixed_value_idx < num_fixed_values; ++fixed_value_idx) {
        (*fixed)[fixed_value_idx].value = fixed_values[fixed_value_idx];
        (*fixed)[fixed_value_idx].num_tgt = 0;
        (*fixed)[fixed_value_idx].tgt = NULL;
      }
    }
    free(fixed_values);
    free(global_tgt_indices);
    free(num_tgt_indices_per_fixed_value);
  }
#endif
}

static inline int
compute_bucket(yac_int value, int comm_size) {
  return (int)(value / 128) % comm_size;
}

static int get_tgt_pack_size(MPI_Comm comm) {

  int tgt_pack_size;
  yac_mpi_call(MPI_Pack_size(1, yac_int_dt, comm, &tgt_pack_size), comm);
  return tgt_pack_size;
}

static int get_link_pack_size(MPI_Comm comm) {

  int global_id_pack_size;
  int weight_pack_size;
  int src_field_idx_pack_size;

  yac_mpi_call(
    MPI_Pack_size(1, yac_int_dt, comm, &global_id_pack_size), comm);
  yac_mpi_call(
    MPI_Pack_size(1, MPI_DOUBLE, comm, &weight_pack_size), comm);
  yac_mpi_call(
    MPI_Pack_size(1, MPI_UINT64_T, comm, &src_field_idx_pack_size), comm);

  return 2 * global_id_pack_size + weight_pack_size +
         src_field_idx_pack_size;
}

static int get_fixed_pack_size(MPI_Comm comm) {

  int value_pack_size;
  int global_id_pack_size;
  int count_pack_size;

  yac_mpi_call(
    MPI_Pack_size(1, MPI_DOUBLE, comm, &value_pack_size), comm);
  yac_mpi_call(
    MPI_Pack_size(1, yac_int_dt, comm, &global_id_pack_size), comm);
  yac_mpi_call(
    MPI_Pack_size(1, MPI_UINT64_T, comm, &count_pack_size), comm);

  return value_pack_size + global_id_pack_size + count_pack_size;
}

static void pack_tgt(
  yac_int global_tgt_id, void * buffer, int buffer_size, MPI_Comm comm) {

  int position = 0;
  yac_mpi_call(
      MPI_Pack(&global_tgt_id, 1, yac_int_dt, buffer,
               buffer_size, &position, comm), comm);
}

static void pack_link(
  struct link_data * link, void * buffer, int buffer_size, MPI_Comm comm) {

  int position = 0;
  yac_mpi_call(
      MPI_Pack(&(link->src.global_id), 1, yac_int_dt, buffer,
               buffer_size, &position, comm), comm);
  yac_mpi_call(
      MPI_Pack(&(link->tgt.global_id), 1, yac_int_dt, buffer,
               buffer_size, &position, comm), comm);
  yac_mpi_call(
      MPI_Pack(&(link->weight), 1, MPI_DOUBLE, buffer,
               buffer_size, &position, comm), comm);
  uint64_t temp_src_field_idx = (uint64_t)(link->src_field_idx);
  yac_mpi_call(
      MPI_Pack(&temp_src_field_idx, 1, MPI_UINT64_T, buffer,
               buffer_size, &position, comm), comm);

}

static void pack_fixed(
  double fixed_value, yac_int global_tgt_id,
  void * buffer, int buffer_size, MPI_Comm comm) {

  int position = 0;
  yac_mpi_call(
      MPI_Pack(&fixed_value, 1, MPI_DOUBLE, buffer,
               buffer_size, &position, comm), comm);
  yac_mpi_call(
      MPI_Pack(&global_tgt_id, 1, yac_int_dt, buffer,
               buffer_size, &position, comm), comm);
}

static void unpack_tgt(
  void * buffer, int buffer_size, yac_int * global_tgt_id, MPI_Comm comm) {

  int position = 0;
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, &position, global_tgt_id, 1,
               yac_int_dt, comm), comm);
}

static void unpack_link(
  void * buffer, int buffer_size, struct link_data * link, MPI_Comm comm) {

  int position = 0;
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, &position, &(link->src.global_id), 1,
               yac_int_dt, comm), comm);
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, &position, &(link->tgt.global_id), 1,
               yac_int_dt, comm), comm);
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, &position, &(link->weight), 1,
               MPI_DOUBLE, comm), comm);
  uint64_t temp_src_field_idx;
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, &position, &temp_src_field_idx, 1,
               MPI_UINT64_T, comm), comm);
  link->src_field_idx = (size_t)temp_src_field_idx;
}

static void unpack_fixed(
  void * buffer, int buffer_size, struct temp_fixed_data * fixed,
  MPI_Comm comm) {

  int position = 0;
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, &position, &(fixed->value), 1,
               MPI_DOUBLE, comm), comm);
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, &position, &(fixed->tgt.global_id), 1,
               yac_int_dt, comm), comm);
}

static int compare_temp_result_tgt (void const * a, void const * b) {

  struct temp_result const * result_a = (struct temp_result const *)a;
  struct temp_result const * result_b = (struct temp_result const *)b;

  return (result_a->tgt.global_id > result_b->tgt.global_id) -
         (result_a->tgt.global_id < result_b->tgt.global_id);
}

static int compare_link_data_tgt (void const * a, void const * b) {

  struct link_data const * link_a = (struct link_data const *)a;
  struct link_data const * link_b = (struct link_data const *)b;

  return (link_a->tgt.global_id > link_b->tgt.global_id) -
         (link_a->tgt.global_id < link_b->tgt.global_id);
}

static int compare_temp_result_type (void const * a, void const * b) {

  struct temp_result const * result_a = (struct temp_result const *)a;
  struct temp_result const * result_b = (struct temp_result const *)b;

  return (result_a->type > result_b->type) - (result_a->type < result_b->type);
}

static int compare_temp_result_fixed (void const * a, void const * b) {

  struct temp_result const * result_a = (struct temp_result const *)a;
  struct temp_result const * result_b = (struct temp_result const *)b;

  return (result_a->data.fixed.value > result_b->data.fixed.value) -
         (result_a->data.fixed.value < result_b->data.fixed.value);
}

static int compare_temp_result_reorder_idx (void const * a, void const * b) {

  struct temp_result const * result_a = (struct temp_result const *)a;
  struct temp_result const * result_b = (struct temp_result const *)b;

  return (result_a->reorder_idx > result_b->reorder_idx) -
         (result_a->reorder_idx < result_b->reorder_idx);
}

static int compare_temp_fixed_data_tgt (void const * a, void const * b) {

  struct temp_fixed_data const * fixed_a = (struct temp_fixed_data const *)a;
  struct temp_fixed_data const * fixed_b = (struct temp_fixed_data const *)b;

  return (fixed_a->tgt.global_id > fixed_b->tgt.global_id) -
         (fixed_a->tgt.global_id < fixed_b->tgt.global_id);
}

static int get_temp_result_pack_size(
  struct temp_result * result, MPI_Comm comm) {

  int type_pack_size;
  int data_pack_size;
  int owner_pack_size;

  yac_mpi_call(
    MPI_Pack_size(1, MPI_INT, comm, &type_pack_size), comm);
  YAC_ASSERT(
    (result->type == LINK) ||
    (result->type == FIXED) ||
    (result->type == NO_INTERP),
    "ERROR(get_temp_result_pack_size): invalid result type")
  switch(result->type) {
    case(LINK): {
      int count_pack_size;
      int links_pack_size;
      yac_mpi_call(
        MPI_Pack_size(1, MPI_UINT64_T, comm, &count_pack_size), comm);
      links_pack_size =
        (int)(result->data.link.count) * get_link_pack_size(comm);
      data_pack_size = links_pack_size + count_pack_size;
      break;
    }
    case(FIXED): {
      int fixed_value_pack_size;
      yac_mpi_call(
        MPI_Pack_size(1, MPI_DOUBLE, comm, &fixed_value_pack_size), comm);
      data_pack_size = fixed_value_pack_size;
      break;
    }
    default:
    case(NO_INTERP): {
      data_pack_size = 0;
      break;
    }
  };
  int tgt_global_id_pack_size;
  yac_mpi_call(
    MPI_Pack_size(1, yac_int_dt, comm, &tgt_global_id_pack_size), comm);
  yac_mpi_call(
    MPI_Pack_size(1, MPI_INT, comm, &owner_pack_size), comm);

  return type_pack_size + data_pack_size + tgt_global_id_pack_size +
         owner_pack_size;
}

static void pack_temp_result(
  struct temp_result * result, void * buffer, int buffer_size, MPI_Comm comm) {

  int position = 0;
  int temp_type = (int)(result->type);
  yac_mpi_call(
      MPI_Pack(&temp_type, 1, MPI_INT, buffer,
               buffer_size, &position, comm), comm);
  YAC_ASSERT(
    (result->type == LINK) ||
    (result->type == FIXED) ||
    (result->type == NO_INTERP),
    "ERROR(get_temp_result_pack_size): invalid result type")
  switch(result->type) {
    case(LINK): {
      uint64_t temp_count = (uint64_t)(result->data.link.count);
      yac_mpi_call(
          MPI_Pack(&temp_count, 1, MPI_UINT64_T, buffer,
                   buffer_size, &position, comm), comm);
      for (uint64_t i = 0; i < temp_count; ++i) {
        yac_mpi_call(
            MPI_Pack(&(result->data.link.links[i].src.global_id), 1,
                     yac_int_dt, buffer, buffer_size, &position, comm), comm);
        yac_mpi_call(
            MPI_Pack(&(result->data.link.links[i].tgt.global_id), 1,
                     yac_int_dt, buffer, buffer_size, &position, comm), comm);
        yac_mpi_call(
            MPI_Pack(&(result->data.link.links[i].weight), 1,
                     MPI_DOUBLE, buffer, buffer_size, &position, comm), comm);
        uint64_t temp_src_field_idx =
          (uint64_t)result->data.link.links[i].src_field_idx;
        yac_mpi_call(
            MPI_Pack(&temp_src_field_idx, 1, MPI_UINT64_T,
                     buffer, buffer_size, &position, comm), comm);
      }
      break;
    }
    case(FIXED): {
      yac_mpi_call(
          MPI_Pack(&(result->data.fixed.value), 1,
                   MPI_DOUBLE, buffer, buffer_size, &position, comm), comm);
      break;
    }
    default:
    case(NO_INTERP): {
      break;
    }
  };
  yac_mpi_call(
      MPI_Pack(&(result->tgt.global_id), 1,
               yac_int_dt, buffer, buffer_size, &position, comm), comm);
  yac_mpi_call(
      MPI_Pack(&(result->owner), 1, MPI_INT,
               buffer, buffer_size, &position, comm), comm);
}

static void unpack_temp_result(
  void * buffer, int buffer_size, struct temp_result * result,
  MPI_Comm comm) {

  int position = 0;

  int temp_type;
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, &position, &temp_type, 1,
               MPI_INT, comm), comm);
  YAC_ASSERT(
    (temp_type == LINK) || (temp_type == FIXED) || (temp_type == NO_INTERP),
    "ERROR(unpack_temp_result): invalid result type")
  switch (temp_type) {
    case(LINK): {
      result->type = LINK;
      uint64_t temp_count;
      yac_mpi_call(
        MPI_Unpack(buffer, buffer_size, &position, &temp_count, 1,
                   MPI_UINT64_T, comm), comm);
      result->data.link.count = (size_t)temp_count;
      result->data.link.links =
        xmalloc((size_t)temp_count * sizeof(*(result->data.link.links)));
      for (uint64_t i = 0; i < temp_count; ++i) {

        yac_mpi_call(
          MPI_Unpack(buffer, buffer_size, &position,
                     &(result->data.link.links[i].src.global_id), 1,
                     yac_int_dt, comm), comm);
        yac_mpi_call(
          MPI_Unpack(buffer, buffer_size, &position,
                     &(result->data.link.links[i].tgt.global_id), 1,
                     yac_int_dt, comm), comm);
        yac_mpi_call(
          MPI_Unpack(buffer, buffer_size, &position,
                     &(result->data.link.links[i].weight), 1,
                     MPI_DOUBLE, comm), comm);
        uint64_t temp_src_field_idx;
        yac_mpi_call(
          MPI_Unpack(buffer, buffer_size, &position, &temp_src_field_idx, 1,
                     MPI_UINT64_T, comm), comm);
        result->data.link.links[i].src_field_idx = (size_t)temp_src_field_idx;
      }
      break;
    }
    case (FIXED): {
      result->type = FIXED;
      yac_mpi_call(
        MPI_Unpack(buffer, buffer_size, &position,
                   &(result->data.fixed.value), 1, MPI_DOUBLE, comm), comm);
      break;
    }
    default:
    case (NO_INTERP): {
      result->type = NO_INTERP;
      break;
    }
  };
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, &position,
               &(result->tgt.global_id), 1, yac_int_dt, comm), comm);
  yac_mpi_call(
    MPI_Unpack(buffer, buffer_size, &position,
               &(result->owner), 1, MPI_INT, comm), comm);
}

static void redist_weight_file_data(
  struct yac_interp_grid * interp_grid, size_t * tgt_points, size_t tgt_count,
  size_t * num_interpolated_points,
  struct link_data ** links, size_t * num_links,
  struct fixed_data ** fixed, size_t * num_fixed) {

  MPI_Comm comm = yac_interp_grid_get_MPI_Comm(interp_grid);

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // get global ids of target points that have to be interpolated
  yac_int * tgt_global_ids = xmalloc(tgt_count * sizeof(*tgt_global_ids));
  yac_interp_grid_get_tgt_global_ids(
    interp_grid, tgt_points, tgt_count, tgt_global_ids);

  yac_quicksort_index_yac_int_size_t(tgt_global_ids, tgt_count, tgt_points);

  int tgt_pack_size = get_tgt_pack_size(comm),
      link_pack_size = get_link_pack_size(comm),
      fixed_pack_size = get_fixed_pack_size(comm);

  size_t * size_t_buffer =
    xmalloc(4 * (size_t)comm_size * sizeof(*size_t_buffer));
  size_t * total_sendcounts = size_t_buffer + 0 * comm_size;
  size_t * total_recvcounts = size_t_buffer + 1 * comm_size;
  size_t * total_sdispls =    size_t_buffer + 2 * comm_size;
  size_t * total_rdispls =    size_t_buffer + 3 * comm_size;
  size_t * sendcounts, * recvcounts, * sdispls, *rdispls;
  yac_get_comm_buffers(
    3, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  for (size_t i = 0; i < tgt_count; ++i)
    sendcounts[3 * compute_bucket(tgt_global_ids[i], comm_size) + 0]++;
  for (size_t i = 0; i < *num_links; ++i)
    sendcounts[3 * compute_bucket((*links)[i].tgt.global_id, comm_size) + 1]++;
  for (size_t i = 0; i < *num_fixed; ++i)
    for (size_t j = 0; j < (*fixed)[i].num_tgt; ++j)
      sendcounts[
        3 * compute_bucket((*fixed)[i].tgt[j].global_id, comm_size) + 2]++;

  // exchange sendcounts
  yac_mpi_call(MPI_Alltoall(sendcounts, 3, YAC_MPI_SIZE_T,
                            recvcounts, 3, YAC_MPI_SIZE_T, comm), comm);

  size_t saccu = 0, raccu = 0;
  for (int i = 0; i < comm_size; ++i) {
    total_sdispls[i] = saccu;
    sdispls[3*i+0] = saccu;
    sdispls[3*i+1] =
      (saccu +
       (total_sendcounts[i] =
          sendcounts[3*i+0] * (size_t)tgt_pack_size));
    sdispls[3*i+2] =
      (saccu +
       (total_sendcounts[i] +=
          sendcounts[3*i+1] * (size_t)link_pack_size));
    total_sendcounts[i] +=
      sendcounts[3*i+2] * (size_t)fixed_pack_size;
    total_rdispls[i] = raccu;
    rdispls[3*i+0] = raccu;
    rdispls[3*i+1] =
      (raccu +
       (total_recvcounts[i] =
          recvcounts[3*i+0] * (size_t)tgt_pack_size));
    rdispls[3*i+2] =
      (raccu +
       (total_recvcounts[i] +=
          recvcounts[3*i+1] * (size_t)link_pack_size));
    total_recvcounts[i] +=
      recvcounts[3*i+2] * (size_t)fixed_pack_size;
    saccu += total_sendcounts[i];
    raccu += total_recvcounts[i];
  }

  size_t send_size = total_sendcounts[comm_size - 1] +
                     total_sdispls[comm_size - 1];
  size_t recv_size = total_recvcounts[comm_size - 1] +
                     total_rdispls[comm_size - 1];

  void * pack_buffer = xmalloc(send_size + recv_size);
  void * send_pack_buffer = pack_buffer;
  void * recv_pack_buffer = (void*)((unsigned char *)pack_buffer + send_size);

  // pack data
  for (size_t i = 0; i < tgt_count; ++i) {
    yac_int curr_global_id = tgt_global_ids[i];
    int rank = compute_bucket(curr_global_id, comm_size);
    pack_tgt(
      curr_global_id,
      (void*)((unsigned char*)send_pack_buffer + sdispls[3*rank+0]),
      tgt_pack_size, comm);
    sdispls[3*rank+0] += (size_t)tgt_pack_size;
  }
  for (size_t i = 0; i < *num_links; ++i) {
    struct link_data * curr_link = (*links) + i;
    int rank = compute_bucket(curr_link->tgt.global_id, comm_size);
    pack_link(
      curr_link,
      (void*)((unsigned char*)send_pack_buffer + sdispls[3*rank+1]),
      link_pack_size, comm);
    sdispls[3*rank+1] += (size_t)link_pack_size;
  }
  for (size_t i = 0; i < *num_fixed; ++i) {
    double curr_fixed_value = (*fixed)[i].value;
    size_t curr_num_tgt = (*fixed)[i].num_tgt;
    for (size_t j = 0; j < curr_num_tgt; ++j) {
      yac_int curr_global_id = (*fixed)[i].tgt[j].global_id;
      int rank = compute_bucket(curr_global_id, comm_size);
      pack_fixed(
        curr_fixed_value, curr_global_id,
        (void*)((unsigned char*)send_pack_buffer + sdispls[3*rank+2]),
        fixed_pack_size, comm);
      sdispls[3*rank+2] += (size_t)fixed_pack_size;
    }
  }

  // exchange data
  yac_alltoallv_packed_p2p(
    send_pack_buffer, total_sendcounts, total_sdispls,
    recv_pack_buffer, total_recvcounts, total_rdispls, comm);

  size_t recv_tgt_count = 0, recv_link_count = 0, recv_fixed_count = 0;
  for (int i = 0; i < comm_size; ++i) {
    recv_tgt_count += recvcounts[3*i+0];
    recv_link_count += recvcounts[3*i+1];
    recv_fixed_count += recvcounts[3*i+2];
  }
  struct link_data * recv_links =
    xrealloc(*links, recv_link_count * sizeof(*recv_links));
  struct temp_fixed_data * recv_fixed =
    xmalloc(recv_fixed_count * sizeof(*recv_fixed));
  struct temp_result * result_buffer =
    xmalloc((recv_tgt_count + tgt_count) * sizeof(*result_buffer));
  struct temp_result * send_results = result_buffer;
  struct temp_result * recv_results = result_buffer + recv_tgt_count;

  // unpack data
  recv_tgt_count = 0, recv_link_count = 0, recv_fixed_count = 0;
  size_t unpack_offset = 0;
  for (int i = 0; i < comm_size; ++i) {
    size_t curr_num_tgt = recvcounts[3 * i + 0];
    size_t curr_num_links = recvcounts[3 * i + 1];
    size_t curr_num_fixed = recvcounts[3 * i + 2];
    for (size_t j = 0; j < curr_num_tgt; ++j, ++recv_tgt_count) {
      unpack_tgt(
        (void*)((unsigned char *)recv_pack_buffer + unpack_offset),
        tgt_pack_size, &(send_results[recv_tgt_count].tgt.global_id),
        comm);
      unpack_offset += (size_t)tgt_pack_size;
      send_results[recv_tgt_count].owner = i;
      send_results[recv_tgt_count].reorder_idx = recv_tgt_count;
    }
    for (size_t j = 0; j < curr_num_links; ++j, ++recv_link_count) {
      unpack_link(
        (void*)((unsigned char *)recv_pack_buffer + unpack_offset),
        link_pack_size, recv_links + recv_link_count, comm);
      unpack_offset += (size_t)link_pack_size;
    }
    for (size_t j = 0; j < curr_num_fixed; ++j, ++recv_fixed_count) {
      unpack_fixed(
        (void*)((unsigned char *)recv_pack_buffer + unpack_offset),
        fixed_pack_size, recv_fixed + recv_fixed_count, comm);
      unpack_offset += (size_t)fixed_pack_size;
    }
  }

  // sort received target point inforation by global target id
  qsort(
    send_results, recv_tgt_count, sizeof(*send_results),
    compare_temp_result_tgt);

  // sort links by global target ids
  qsort(
    recv_links, recv_link_count, sizeof(*recv_links), compare_link_data_tgt);

  // sort fixed data by global target ids
  qsort(
    recv_fixed, recv_fixed_count, sizeof(*recv_fixed),
    compare_temp_fixed_data_tgt);

  // match target points with linked and fixed data from the weight file
  for (size_t tgt_idx = 0, link_idx = 0, fixed_idx = 0;
       tgt_idx < recv_tgt_count; ++tgt_idx) {

    yac_int curr_tgt = send_results[tgt_idx].tgt.global_id;

    while ((link_idx < recv_link_count) &&
           (recv_links[link_idx].tgt.global_id < curr_tgt)) ++link_idx;
    while ((fixed_idx < recv_fixed_count) &&
           (recv_fixed[fixed_idx].tgt.global_id < curr_tgt)) ++fixed_idx;

    YAC_ASSERT(
      (link_idx >= recv_link_count) ||
      (recv_links[link_idx].tgt.global_id != curr_tgt) ||
      (fixed_idx >= recv_fixed_count) ||
      (recv_fixed[fixed_idx].tgt.global_id != curr_tgt),
      "ERROR(redist_weight_file_data): inconsistent data in weight file;"
      "link and fixed data available for a target point")
    if ((link_idx < recv_link_count) &&
        (recv_links[link_idx].tgt.global_id == curr_tgt)) {

      size_t temp_link_idx = link_idx;
      while ((temp_link_idx < recv_link_count) &&
             (recv_links[temp_link_idx].tgt.global_id == curr_tgt))
        ++temp_link_idx;
      send_results[tgt_idx].type = LINK;
      send_results[tgt_idx].data.link.links = recv_links + link_idx;
      send_results[tgt_idx].data.link.count = temp_link_idx - link_idx;

    } else if ((fixed_idx < recv_fixed_count) &&
               (recv_fixed[fixed_idx].tgt.global_id == curr_tgt)) {

      send_results[tgt_idx].type = FIXED;
      send_results[tgt_idx].data.fixed.value = recv_fixed[fixed_idx].value;

    } else {
      send_results[tgt_idx].type = NO_INTERP;
    }
  }
  free(recv_fixed);

  // sort received target point inforatiom into origional receive order
  qsort(
    send_results, recv_tgt_count, sizeof(*send_results),
    compare_temp_result_reorder_idx);

  int * pack_size_buffer =
    xmalloc((recv_tgt_count + tgt_count) * sizeof(*pack_size_buffer));
  int * send_pack_sizes = pack_size_buffer;
  int * recv_pack_sizes = pack_size_buffer + recv_tgt_count;
  for (size_t i = 0; i < recv_tgt_count; ++i)
    send_pack_sizes[i] = get_temp_result_pack_size(send_results + i, comm);

  saccu = 0, raccu = 0;
  for (int i = 0; i < comm_size; ++i) {
    sdispls[i] = saccu;
    rdispls[i] = raccu;
    int sendcount = recvcounts[3*i];
    int recvcount = sendcounts[3*i];
    saccu += (sendcounts[i] = sendcount);
    raccu += (recvcounts[i] = recvcount);
  }

  yac_alltoallv_int_p2p(
    send_pack_sizes, sendcounts, sdispls,
    recv_pack_sizes, recvcounts, rdispls, comm);

  saccu = 0, raccu = 0;
  for (int i = 0, k = 0, l = 0; i < comm_size; ++i) {
    size_t sendcount = sendcounts[i];
    size_t recvcount = recvcounts[i];
    sendcounts[i] = 0;
    recvcounts[i] = 0;
    for (size_t j = 0; j < sendcount; ++j, ++k)
      sendcounts[i] += send_pack_sizes[k];
    for (size_t j = 0; j < recvcount; ++j, ++l)
      recvcounts[i] += recv_pack_sizes[l];
    sdispls[i] = saccu;
    rdispls[i] = raccu;
    saccu += sendcounts[i];
    raccu += recvcounts[i];
  }
  size_t total_send_pack_size =
    sdispls[comm_size-1] + sendcounts[comm_size-1];
  size_t total_recv_pack_size =
    rdispls[comm_size-1] + recvcounts[comm_size-1];

  pack_buffer =
    xrealloc(pack_buffer, total_send_pack_size + total_recv_pack_size);
  send_pack_buffer = pack_buffer;
  recv_pack_buffer =
    (void*)((unsigned char*)pack_buffer + total_send_pack_size);

  for (size_t i = 0, offset = 0; i < recv_tgt_count; ++i) {
    pack_temp_result(
      send_results + i, (void*)((unsigned char*)send_pack_buffer + offset),
      send_pack_sizes[i], comm);
    offset += (size_t)(send_pack_sizes[i]);
  }

  // redistribute link and fixed data to processes requiring this data for
  // their local target points
  yac_alltoallv_p2p(
    send_pack_buffer, sendcounts, sdispls,
    recv_pack_buffer, recvcounts, rdispls,
    1, MPI_PACKED, comm);

  // unpack results
  for (size_t i = 0, offset = 0; i < tgt_count; ++i) {
    unpack_temp_result(
      (void*)((unsigned char*)recv_pack_buffer + offset),
      recv_pack_sizes[i], recv_results + i, comm);
    offset += (size_t)(recv_pack_sizes[i]);
  }
  free(pack_size_buffer);

  // sort received results by global target id
  qsort(
    recv_results, tgt_count, sizeof(*recv_results),
    compare_temp_result_tgt);

  for (size_t i = 0; i < tgt_count; ++i)
    recv_results[i].tgt.local_id = tgt_points[i];

  // sort received results by type
  qsort(
    recv_results, tgt_count, sizeof(*recv_results),
    compare_temp_result_type);

  size_t new_num_links = 0;
  size_t new_total_num_links = 0;
  size_t new_num_fixed = 0;
  {
    size_t i = 0;
    while ((i < tgt_count) && (recv_results[i].type == LINK))
      new_total_num_links += recv_results[i++].data.link.count;
    new_num_links = i;
    while ((i < tgt_count) && (recv_results[i].type == FIXED)) ++i;
    new_num_fixed = i - new_num_links;
  }

  // extract links from results
  *links = xrealloc(recv_links, new_total_num_links * sizeof(**links));
  for (size_t i = 0, offset = 0; i < new_num_links;
       offset += recv_results[i++].data.link.count)
    memcpy(*links + offset, recv_results[i].data.link.links,
           recv_results[i].data.link.count * sizeof(**links));
  *num_links = new_total_num_links;

  // extract fixed data from results
  if (new_num_fixed > 0) {
    recv_results += new_num_links;
    qsort(
      recv_results, new_num_fixed, sizeof(*recv_results),
      compare_temp_result_fixed);
    double curr_fixed_value =
      (recv_results[0].data.fixed.value == 1337.0)?(-1337.0):(1337.0);
    size_t num_fixed_values = 0;
    void * tgt_global_id_buffer =
      xrealloc((*num_fixed > 0)?(*fixed)->tgt:NULL,
               new_num_fixed * sizeof(*((*fixed)->tgt)));
    struct fixed_data * curr_fixed = NULL;
    for (size_t i = 0, offset = 0; i < new_num_fixed; ++i) {
      if (recv_results[i].data.fixed.value != curr_fixed_value) {
        curr_fixed_value = recv_results[i].data.fixed.value;
        *fixed = xrealloc(*fixed, (num_fixed_values + 1) * sizeof(**fixed));
        curr_fixed = *fixed + num_fixed_values;
        ++num_fixed_values;
        curr_fixed->value = curr_fixed_value;
        curr_fixed->tgt = (void*)(((unsigned char *)tgt_global_id_buffer) +
                                  offset * sizeof(*((*fixed)->tgt)));
        curr_fixed->num_tgt = 0;
      }
      curr_fixed->tgt[curr_fixed->num_tgt++].global_id =
        recv_results[i].tgt.global_id;
    }
    *num_fixed = num_fixed_values;
    recv_results -= new_num_links;

  } else {
    if (*num_fixed > 0) free((*fixed)->tgt);
    free(*fixed);
    *fixed = NULL;
    *num_fixed = 0;
  }

  for (size_t i = 0; i < tgt_count; ++i)
    tgt_points[i] = recv_results[i].tgt.local_id;
  *num_interpolated_points = new_num_links + new_num_fixed;

  for (size_t i = 0; i < new_num_links; ++i)
    free(recv_results[i].data.link.links);
  free(result_buffer);
  free(pack_buffer);
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
  free(size_t_buffer);
  free(tgt_global_ids);
}

static int compare_link_data_field_idx_src(void const *a, void const *b) {

  struct link_data const * link_a = (struct link_data*)a;
  struct link_data const * link_b = (struct link_data*)b;

  int ret = (link_a->src_field_idx > link_b->src_field_idx) -
            (link_a->src_field_idx < link_b->src_field_idx);

  if (ret) return ret;

  else return (link_a->src.global_id > link_b->src.global_id) -
              (link_a->src.global_id < link_b->src.global_id);
}

// sort links by tgt local id and src ids (SIZE_MAX first; remaining by ascending src ids)
static int compare_link_data_tgt_src_mask(void const *a, void const *b) {

  struct link_data const * link_a = (struct link_data*)a;
  struct link_data const * link_b = (struct link_data*)b;

  int ret = (link_a->tgt.local_id > link_b->tgt.local_id) -
            (link_a->tgt.local_id < link_b->tgt.local_id);

  if (ret) return ret;

  if (link_a->src.local_id == SIZE_MAX) return -1;
  if (link_b->src.local_id == SIZE_MAX) return 1;

  return (link_a->src.local_id > link_b->src.local_id) -
         (link_a->src.local_id < link_b->src.local_id);
}

// sort links by tgt local id, source field index, and src local id
static int compare_link_data_tgt_src_field_src_id(void const *a, void const *b) {

  struct link_data const * link_a = (struct link_data*)a;
  struct link_data const * link_b = (struct link_data*)b;

  int ret = (link_a->tgt.local_id > link_b->tgt.local_id) -
            (link_a->tgt.local_id < link_b->tgt.local_id);

  if (ret) return ret;

  ret = (link_a->src_field_idx > link_b->src_field_idx) -
        (link_a->src_field_idx < link_b->src_field_idx);

  if (ret) return ret;

  if (link_a->src.local_id == SIZE_MAX) return -1;
  if (link_b->src.local_id == SIZE_MAX) return 1;

  return (link_a->src.local_id > link_b->src.local_id) -
         (link_a->src.local_id < link_b->src.local_id);
}


static size_t do_search_file(struct interp_method * method,
                             struct yac_interp_grid * interp_grid,
                             size_t * tgt_points, size_t count,
                             struct yac_interp_weights * weights) {

  struct interp_method_file * method_file =
    (struct interp_method_file*)(method);

  MPI_Comm comm = yac_interp_grid_get_MPI_Comm(interp_grid);

  size_t num_src_fields = yac_interp_grid_get_num_src_fields(interp_grid);
  enum yac_location * src_locations =
    xmalloc(num_src_fields * sizeof(*src_locations));
  for (size_t i = 0; i < num_src_fields; ++i)
    src_locations[i] = yac_interp_grid_get_src_field_location(interp_grid, i);
  enum yac_location tgt_location =
    yac_interp_grid_get_tgt_field_location(interp_grid);

  struct link_data * links = NULL;
  size_t num_links = 0;
  struct fixed_data * fixed = NULL;
  size_t num_fixed = 0;

  // read data from file
  read_weight_file(
    method_file->weight_file_name, comm,
    yac_interp_grid_get_src_grid_name(interp_grid),
    yac_interp_grid_get_tgt_grid_name(interp_grid),
    src_locations, tgt_location, num_src_fields,
    &links, &num_links, &fixed, &num_fixed);

  free(src_locations);

  // relocates links and fixed to the processes, which require them to
  // interpolate their local target points
  // (the tgt_points first contain the link-tgt_points, then the fixed-tgt_points,
  // and followed the non-interpolated tgt_points)
  size_t num_interpolated_points;
  redist_weight_file_data(
    interp_grid, tgt_points, count, &num_interpolated_points,
    &links, &num_links, &fixed, &num_fixed);

  size_t * tgt_points_reorder_idx =
    xmalloc(num_interpolated_points * sizeof(*tgt_points_reorder_idx));
  for (size_t i = 0; i < num_interpolated_points; ++i)
    tgt_points_reorder_idx[i] = i;

  size_t num_fixed_tgt = 0;
  if (num_fixed) {
    for (size_t i = 0; i < num_fixed; ++i) num_fixed_tgt += fixed[i].num_tgt;

    size_t tgt_points_offset = num_interpolated_points - num_fixed_tgt;

    for (size_t i = 0; i < num_fixed; ++i) {

      struct remote_points tgts = {
        .data =
          yac_interp_grid_get_tgt_remote_points(
            interp_grid, tgt_points + tgt_points_offset, fixed[i].num_tgt),
        .count = fixed[i].num_tgt};

      yac_interp_weights_add_fixed(weights, &tgts, fixed[i].value);
      free(tgts.data);
      tgt_points_offset += fixed[i].num_tgt;
    }
    free(fixed->tgt);
  }
  free(fixed);  

  size_t num_link_tgt = 0;
  if (num_links > 0) {

    // set local target ids for all links
    for (size_t link_idx = 0, tgt_idx = 0; link_idx < num_links; tgt_idx++) {

      yac_int curr_global_id = links[link_idx].tgt.global_id;
      size_t curr_local_id = tgt_points[tgt_idx];

      while ((link_idx < num_links) &&
             (links[link_idx].tgt.global_id == curr_global_id)) {
        links[link_idx].tgt.local_id = curr_local_id;
        link_idx++;
      }
    }

    // sort links by field idx and src id
    qsort(links, num_links, sizeof(*links), compare_link_data_field_idx_src);

    // determine number of source fields
    size_t num_src_fields = links[num_links-1].src_field_idx + 1;

    {
      int int_num_src_fields = (int)num_src_fields;
      MPI_Comm comm = yac_interp_grid_get_MPI_Comm(interp_grid);
      yac_mpi_call(
        MPI_Allreduce(MPI_IN_PLACE, &int_num_src_fields, 1, MPI_INT, MPI_MAX,
                      comm), comm);
      num_src_fields = (size_t)int_num_src_fields;
    }

    yac_int * src_global_ids = xmalloc(num_links * sizeof(*src_global_ids));
    size_t * src_local_ids = xmalloc(num_links * sizeof(*src_local_ids));
    size_t src_field_prefix_sum[num_src_fields];

    // check the mask of the source points and get the local ids for the
    // remaining ones
    for (size_t src_field_idx = 0, link_idx = 0, offset = 0;
         src_field_idx < num_src_fields; ++src_field_idx) {

      size_t curr_num_links = 0;
      while ((link_idx < num_links) &&
             (links[link_idx].src_field_idx == src_field_idx)) {

        // get global source ids
        src_global_ids[offset + curr_num_links] =
          links[link_idx].src.global_id;
        ++curr_num_links;
        ++link_idx;
      }

      // get local ids for all source points
      yac_interp_grid_src_global_to_local(
        interp_grid, src_field_idx, src_global_ids + offset, curr_num_links,
        src_local_ids + offset);

      // get the mask after the call to yac_interp_grid_src_global_to_local
      const_int_pointer src_mask =
        yac_interp_grid_get_src_field_mask(interp_grid, src_field_idx);

      // check source masks and mask src ids of marked links with SIZE_MAX
      link_idx -= curr_num_links;
      if (src_mask != NULL) {
        for (size_t i = 0; i < curr_num_links; ++i, ++link_idx)
          links[link_idx].src.local_id =
            (src_mask[src_local_ids[link_idx]])?
              src_local_ids[link_idx]:SIZE_MAX;
      } else {
        for (size_t i = 0; i < curr_num_links; ++i, ++link_idx)
          links[link_idx].src.local_id = src_local_ids[link_idx];
      }

      src_field_prefix_sum[src_field_idx] = offset;
      offset += curr_num_links;
    }

    // sort links by tgt local id and src ids (SIZE_MAX first; remaining by ascending src ids)
    qsort(links, num_links, sizeof(*links), compare_link_data_tgt_src_mask);

    // pack links (remove links with targets that have at least one masked source point)
    size_t new_num_links = 0;
    int contains_masked_tgt = 0;
    num_link_tgt = num_interpolated_points - num_fixed_tgt;
    size_t * tgt_local_ids = xmalloc(num_link_tgt * sizeof(*tgt_local_ids));
    size_t * num_src_per_field_per_tgt =
      xcalloc(num_link_tgt * num_src_fields, sizeof(*num_src_per_field_per_tgt));
    size_t num_links_per_src_field[num_src_fields];
    memset(num_links_per_src_field, 0, sizeof(num_links_per_src_field));
    num_link_tgt = 0;
    for (size_t link_idx = 0, tgt_idx = 0; link_idx < num_links; ++tgt_idx) {

      size_t curr_tgt_local_id = links[link_idx].tgt.local_id;
      int is_masked = links[link_idx].src.local_id == SIZE_MAX;
      tgt_points[tgt_idx] = curr_tgt_local_id;

      // if the current target point has a masked source point
      if (is_masked) {

        while ((link_idx < num_links) &&
               (links[link_idx].tgt.local_id == curr_tgt_local_id)) ++link_idx;
        contains_masked_tgt = 1;
        tgt_points_reorder_idx[tgt_idx] += num_interpolated_points;

      // if there already have been target points with masked source points
      } else {
        size_t * curr_num_src_per_field_per_tgt =
          num_src_per_field_per_tgt + num_link_tgt * num_src_fields;
        if (contains_masked_tgt) {
          while ((link_idx < num_links) &&
                 (links[link_idx].tgt.local_id == curr_tgt_local_id)) {
            size_t curr_src_field_idx = links[link_idx].src_field_idx;
            src_local_ids[
              src_field_prefix_sum[curr_src_field_idx] +
              num_links_per_src_field[curr_src_field_idx] +
              curr_num_src_per_field_per_tgt[curr_src_field_idx]++] =
                links[link_idx].src.local_id;
            links[new_num_links++] = links[link_idx++];
          }
        } else {

          while ((link_idx < num_links) &&
                 (links[link_idx].tgt.local_id == curr_tgt_local_id)) {
            size_t curr_src_field_idx = links[link_idx].src_field_idx;
            src_local_ids[
              src_field_prefix_sum[curr_src_field_idx] +
              num_links_per_src_field[curr_src_field_idx] +
              curr_num_src_per_field_per_tgt[curr_src_field_idx]++] =
                links[link_idx].src.local_id;
            ++link_idx, ++new_num_links;
          }
        }
        for (size_t i = 0; i < num_src_fields; ++i) {
          num_links_per_src_field[i] += curr_num_src_per_field_per_tgt[i];
        }
        tgt_local_ids[num_link_tgt++] = curr_tgt_local_id;
      }
    }
    num_links = new_num_links;
    tgt_local_ids =
      xrealloc(tgt_local_ids, num_link_tgt * sizeof(*tgt_local_ids));
    num_src_per_field_per_tgt =
      xrealloc(
        num_src_per_field_per_tgt,
        num_link_tgt * num_src_fields * sizeof(*num_src_per_field_per_tgt));

    // get the remote points for all required source points
    struct remote_point * srcs_per_field[num_src_fields];
    for (size_t src_field_idx = 0; src_field_idx < num_src_fields; ++src_field_idx)
      srcs_per_field[src_field_idx] =
        yac_interp_grid_get_src_remote_points(
          interp_grid, src_field_idx, src_local_ids + src_field_prefix_sum[src_field_idx],
          num_links_per_src_field[src_field_idx]);
    free(src_local_ids);

    // sort links by tgt local id, source field index, and src local id
    qsort(links, num_links, sizeof(*links), compare_link_data_tgt_src_field_src_id);

    // gather weights
    double * w = xmalloc(num_links * sizeof(*w));
    for (size_t link_idx = 0; link_idx < num_links; ++link_idx)
      w[link_idx] = links[link_idx].weight;

    struct remote_points tgts = {
      .data =
        yac_interp_grid_get_tgt_remote_points(
          interp_grid, tgt_local_ids, num_link_tgt),
      .count = num_link_tgt};
    free(tgt_local_ids);

    yac_interp_weights_add_wsum_mf(
      weights, &tgts, num_src_per_field_per_tgt, srcs_per_field, w, num_src_fields);

    free(tgts.data);
    for (size_t i = 0; i < num_src_fields; ++i) free(srcs_per_field[i]);
    free(num_src_per_field_per_tgt);
    free(w);
    free(src_global_ids);

  } else {
      int num_src_fields = 0;
      MPI_Comm comm = yac_interp_grid_get_MPI_Comm(interp_grid);
      yac_mpi_call(
        MPI_Allreduce(MPI_IN_PLACE, &num_src_fields, 1, MPI_INT, MPI_MAX,
                      comm), comm);
    // get local ids for all source points
    // (since this routine is collective for all processes, we have to call
    // this here if the local process received no links)
    for (int i = 0; i < num_src_fields; ++i)
      yac_interp_grid_src_global_to_local(interp_grid, 0, NULL, 0, NULL);
  }
  free(links);

  yac_quicksort_index_size_t_size_t(
    tgt_points_reorder_idx, num_interpolated_points, tgt_points);
  free(tgt_points_reorder_idx);

  return num_fixed_tgt + num_link_tgt;
}

struct interp_method * yac_interp_method_file_new(
  char const * weight_file_name) {

  struct interp_method_file * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_file_vtable;
  method->weight_file_name = strdup(weight_file_name);

  return (struct interp_method*)method;
}

static void delete_file(struct interp_method * method) {

  struct interp_method_file * method_file =
    (struct interp_method_file *)method;
  free((void*)(method_file->weight_file_name));
  free(method_file);
}
