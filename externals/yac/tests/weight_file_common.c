/**
 * @file weight_file_common.c
 *
 * @copyright Copyright  (C)  2016 DKRZ, MPI-M
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
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
 
#include <string.h>

#include "tests.h"
#include "weight_file_common.h"
#include "interp_method_file.h"
#include "utils.h"
#include "io_utils.h"

#include <netcdf.h>

void write_weight_file(char const * file_name, int const * src_id,
                       int const * tgt_id, double const * weights,
                       unsigned num_links,
                       enum yac_location const * src_locations,
                       unsigned num_src_fields,
                       int const * num_links_per_src_field,
                       int * tgt_id_fixed, unsigned num_fixed_tgt,
                       double * fixed_values,
                       int * num_tgt_per_fixed_value,
                       unsigned num_fixed_values,
                       enum yac_location tgt_location,
                       char const * src_grid_name,
                       char const * tgt_grid_name) {

  int ncid;

  // create file
  yac_nc_create(file_name, NC_CLOBBER, &ncid);

  int dim_weight_id[6];

  // define dimensions
  if (num_links > 0) {
    HANDLE_ERROR(nc_def_dim(ncid, "num_links", num_links, &dim_weight_id[0]));
    HANDLE_ERROR(nc_def_dim(ncid, "num_wgts", 1, &dim_weight_id[1]));
  }
  HANDLE_ERROR(
    nc_def_dim(ncid, "num_src_fields", num_src_fields, &dim_weight_id[2]));
  HANDLE_ERROR(
    nc_def_dim(ncid, "max_loc_str_len", YAC_MAX_LOC_STR_LEN, &dim_weight_id[3]));
  if (num_fixed_values > 0) {
    HANDLE_ERROR(
      nc_def_dim(ncid, "num_fixed_values", num_fixed_values, &dim_weight_id[4]));
    HANDLE_ERROR(
      nc_def_dim(ncid, "num_fixed_dst", num_fixed_tgt, &dim_weight_id[5]));
  }

  int var_src_add_id, var_dst_add_id, var_weight_id, var_num_links_id,
      src_var_locs_id, tgt_var_loc_id, var_fixed_values_id,
      var_num_dst_per_fixed_value_id, var_dst_add_fixed_id;

  // define variables
  if (num_links > 0) {
    HANDLE_ERROR(
      nc_def_var(
        ncid, "src_address", NC_INT, 1, dim_weight_id, &var_src_add_id));
    HANDLE_ERROR(
      nc_def_var(
        ncid, "dst_address", NC_INT, 1, dim_weight_id, &var_dst_add_id));
    HANDLE_ERROR(
      nc_def_var(
        ncid, "remap_matrix", NC_DOUBLE, 2, dim_weight_id, &var_weight_id));
    HANDLE_ERROR(
      nc_def_var(
        ncid, "num_links_per_src_field", NC_INT, 1, dim_weight_id + 2,
        &var_num_links_id));
  }
  HANDLE_ERROR(
    nc_def_var(
      ncid, "src_locations", NC_CHAR, 2, &dim_weight_id[2], &src_var_locs_id));
  HANDLE_ERROR(
    nc_def_var(
      ncid, "dst_location", NC_CHAR, 1, &dim_weight_id[3], &tgt_var_loc_id));
  if (num_fixed_values > 0) {
    HANDLE_ERROR(
      nc_def_var(
        ncid, "fixed_values", NC_DOUBLE, 1, &dim_weight_id[4],
        &var_fixed_values_id));
    HANDLE_ERROR(
      nc_def_var(
        ncid, "num_dst_per_fixed_value", NC_INT, 1, &dim_weight_id[4],
        &var_num_dst_per_fixed_value_id));
    HANDLE_ERROR(
      nc_def_var(
        ncid, "dst_address_fixed", NC_INT, 1, &dim_weight_id[5],
        &var_dst_add_fixed_id));
  }

  // put attributes
  HANDLE_ERROR(nc_put_att_text(ncid, NC_GLOBAL, "version",
                               strlen(YAC_WEIGHT_FILE_VERSION_STRING),
                               YAC_WEIGHT_FILE_VERSION_STRING));
  HANDLE_ERROR(nc_put_att_text(ncid, NC_GLOBAL, "src_grid_name",
                               strlen(src_grid_name), src_grid_name));
  HANDLE_ERROR(nc_put_att_text(ncid, NC_GLOBAL, "dst_grid_name",
                               strlen(tgt_grid_name), tgt_grid_name));
  {
    char const * str_logical[2] = {"FALSE", "TRUE"};
    HANDLE_ERROR(nc_put_att_text(ncid, NC_GLOBAL, "contains_fixed_dst",
                                 strlen(str_logical[num_fixed_values > 0]),
                                 str_logical[num_fixed_values > 0]));
    HANDLE_ERROR(nc_put_att_text(ncid, NC_GLOBAL, "contains_links",
                                 strlen(str_logical[num_links > 0]),
                                 str_logical[num_links > 0]));
  }

  // end definition
  HANDLE_ERROR(nc_enddef(ncid));

  // write mapping data

  if (num_links > 0) {
    int * src_address = xmalloc(num_links * sizeof(*src_address));
    int * tgt_address = xmalloc(num_links * sizeof(*tgt_address));

    for (unsigned i = 0; i < num_links; i++) {
      src_address[i] = src_id[i] + 1;
      tgt_address[i] = tgt_id[i] + 1;
    }

    HANDLE_ERROR(nc_put_var_int(ncid, var_src_add_id, src_address));
    HANDLE_ERROR(nc_put_var_int(ncid, var_dst_add_id, tgt_address));

    free(src_address);
    free(tgt_address);

    HANDLE_ERROR(nc_put_var_double(ncid, var_weight_id, weights));
    HANDLE_ERROR(
      nc_put_var_int(ncid, var_num_links_id, num_links_per_src_field));
  }

  for (unsigned i = 0; i < num_src_fields; ++i) {
    char const * loc_str = yac_loc2str(src_locations[i]);
    size_t str_start[2] = {i, 0};
    size_t str_count[2] = {1, strlen(loc_str)};
    HANDLE_ERROR(
      nc_put_vara_text(ncid, src_var_locs_id, str_start, str_count, loc_str));
  }

  {
    char const * loc_str = yac_loc2str(tgt_location);
    size_t str_start[1] = {0};
    size_t str_count[1] = {strlen(loc_str)};
    HANDLE_ERROR(
      nc_put_vara_text(ncid, tgt_var_loc_id, str_start, str_count, loc_str));
  }

  if (num_fixed_values > 0) {

    int * tgt_address_fixed =
      xmalloc(num_fixed_tgt * sizeof(*tgt_address_fixed));
    for (unsigned i = 0; i < num_fixed_tgt; i++)
      tgt_address_fixed[i] = tgt_id_fixed[i] + 1;

    HANDLE_ERROR(nc_put_var_int(ncid, var_dst_add_fixed_id, tgt_address_fixed));
    HANDLE_ERROR(nc_put_var_double(ncid, var_fixed_values_id, fixed_values));
    HANDLE_ERROR(nc_put_var_int(ncid, var_num_dst_per_fixed_value_id,
                                num_tgt_per_fixed_value));

    free(tgt_address_fixed);
  }

  HANDLE_ERROR(nc_close(ncid));
}

void check_weight_file(char const * file_name, int const * ref_src_address,
                       int const * ref_tgt_address, double const * ref_weights,
                       unsigned ref_num_links,
                       enum yac_location const * ref_src_locations,
                       unsigned ref_num_src_fields,
                       int const * ref_num_links_per_src_field,
                       int const * ref_tgt_address_fixed,
                       double const * ref_fixed_values,
                       int const * ref_num_tgt_per_fixed_value,
                       unsigned ref_num_fixed_values,
                       enum yac_location ref_tgt_location,
                       char const * ref_src_grid_name,
                       char const * ref_tgt_grid_name) {

  int ncid;

  // open file
  yac_nc_open(file_name, NC_NOWRITE, &ncid);

  int dimid, var_id;
  size_t num_wgts;
  HANDLE_ERROR(nc_inq_dimid(ncid, "num_wgts", &dimid));
  HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_wgts));
  YAC_ASSERT(
    num_wgts == 1, "ERROR(check_weight_file): test only supports num_wgts == 1")

  // global attributes
  size_t str_version_len;
  char * str_version;
  var_id = NC_GLOBAL;
  HANDLE_ERROR(nc_inq_attlen(ncid, var_id, "version", &str_version_len));
  str_version = xmalloc(str_version_len + 1);
  str_version[str_version_len] = '\0';
  HANDLE_ERROR(nc_get_att_text(ncid, var_id, "version", str_version));

  size_t str_fixed_len;
  char * str_fixed;
  HANDLE_ERROR(
    nc_inq_attlen(ncid, var_id, "contains_fixed_dst", &str_fixed_len));
  str_fixed = xmalloc(str_fixed_len + 1);
  str_fixed[str_fixed_len] = '\0';
  HANDLE_ERROR(nc_get_att_text(ncid, var_id, "contains_fixed_dst", str_fixed));

  unsigned contains_fixed_dst =
    ((strlen("TRUE") == str_fixed_len) &&
     !strncmp("TRUE", str_fixed, str_fixed_len));
  YAC_ASSERT(
    contains_fixed_dst ||
    ((strlen("FALSE") == str_fixed_len) &&
     !strncmp("FALSE", str_fixed, str_fixed_len)),
    "ERROR(check_weight_file): invalid global attribute contains_fixed_dst")

  char * src_grid_name, * tgt_grid_name;
  size_t src_grid_name_len = 0, tgt_grid_name_len = 0;
  var_id = NC_GLOBAL;
  HANDLE_ERROR(nc_inq_attlen(ncid, var_id, "src_grid_name", &src_grid_name_len));
  src_grid_name = xmalloc(src_grid_name_len + 1);
  src_grid_name[src_grid_name_len] = '\0';
  HANDLE_ERROR(nc_get_att_text(ncid, var_id, "src_grid_name", src_grid_name));

  HANDLE_ERROR(nc_inq_attlen(ncid, var_id, "dst_grid_name", &tgt_grid_name_len));
  tgt_grid_name = xmalloc(tgt_grid_name_len + 1);
  tgt_grid_name[tgt_grid_name_len] = '\0';
  HANDLE_ERROR(nc_get_att_text(ncid, var_id, "dst_grid_name", tgt_grid_name));

  size_t str_link_len;
  char * str_link;
  var_id = NC_GLOBAL;
  HANDLE_ERROR(nc_inq_attlen(ncid, var_id, "contains_links", &str_link_len));
  str_link = xmalloc(str_link_len + 1);
  str_link[str_link_len] = '\0';
  HANDLE_ERROR(nc_get_att_text(ncid, var_id, "contains_links", str_link));

  unsigned contains_links =
    (strlen("TRUE") == str_link_len) &&
    !strncmp("TRUE", str_link, str_link_len);
  YAC_ASSERT(
    contains_links ||
    ((strlen("FALSE") == str_link_len) &&
     !strncmp("FALSE", str_link, str_link_len)),
    "ERROR(check_weight_file): invalid global attribute contains_links")
  free(str_link);

  // get number of links from file
  size_t num_links = 0;
  if (contains_links) {
    HANDLE_ERROR(nc_inq_dimid(ncid, "num_links", &dimid));
    HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_links));
    YAC_ASSERT(num_links != 0, "ERROR(check_weight_file): no links defined")
  }
  if (ref_num_links != num_links)
    PUT_ERR("wrong number of links in weight file\n");

  // get number of source fields from file
  size_t num_src_fields = 0;
  if (contains_links) {
    HANDLE_ERROR(nc_inq_dimid(ncid, "num_src_fields", &dimid));
    HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_src_fields));
    YAC_ASSERT(
      num_src_fields != 0, "ERROR(check_weight_file): no source fields")
    if (ref_num_src_fields != num_src_fields)
      PUT_ERR("wrong number of source fields in weight file\n");
  }

  // get max location string length from file
  size_t max_loc_str_len;
  HANDLE_ERROR(nc_inq_dimid(ncid, "max_loc_str_len", &dimid));
  HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &max_loc_str_len));
  YAC_ASSERT(
    max_loc_str_len == YAC_MAX_LOC_STR_LEN,
    "ERROR(check_weight_file): wrong max location string length")

  size_t num_fixed_values = 0;
  size_t num_fixed_tgt = 0;
  if (contains_fixed_dst) {

    // get number of fixed values
    HANDLE_ERROR(nc_inq_dimid(ncid, "num_fixed_values", &dimid));
    HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_fixed_values));

    // get number for fixed target points
    HANDLE_ERROR(nc_inq_dimid(ncid, "num_fixed_dst", &dimid));
    HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, &num_fixed_tgt));
  }

  int * src_address = xmalloc(num_links * sizeof(*src_address));
  int * tgt_address = xmalloc(num_links * sizeof(*tgt_address));
  double * weights = xmalloc(num_links * sizeof(*weights));
  int * num_links_per_src_field =
    xmalloc(num_src_fields * sizeof(*num_links_per_src_field));
  enum yac_location * src_locations =
    xmalloc(num_src_fields * sizeof(*src_locations));
  enum yac_location tgt_location;
  double * fixed_values = xmalloc(num_fixed_values * sizeof(*fixed_values));
  int * num_tgt_per_fixed_value =
    xmalloc(num_fixed_values * sizeof(*num_tgt_per_fixed_value));
  int * tgt_address_fixed = xmalloc(num_fixed_tgt * sizeof(*tgt_address_fixed));

  // read data
  if (contains_links) {
    yac_nc_inq_varid(ncid, "src_address", &var_id);
    HANDLE_ERROR(nc_get_var_int(ncid, var_id, src_address));

    yac_nc_inq_varid(ncid, "dst_address", &var_id);
    HANDLE_ERROR(nc_get_var_int(ncid, var_id, tgt_address));

    yac_nc_inq_varid(ncid, "remap_matrix", &var_id);
    HANDLE_ERROR(nc_get_var_double(ncid, var_id, weights));

    yac_nc_inq_varid(ncid, "num_links_per_src_field", &var_id);
    HANDLE_ERROR(nc_get_var_int(ncid, var_id, num_links_per_src_field));
  }

  yac_nc_inq_varid(ncid, "src_locations", &var_id);
  for (unsigned i = 0; i < num_src_fields; ++i) {
    char loc_str[YAC_MAX_LOC_STR_LEN];
    size_t str_start[2] = {i, 0};
    size_t str_count[2] = {1, YAC_MAX_LOC_STR_LEN};
    HANDLE_ERROR(nc_get_vara_text(ncid, var_id, str_start, str_count, loc_str));
    src_locations[i] = yac_str2loc(loc_str);
  }

  yac_nc_inq_varid(ncid, "dst_location", &var_id);
  {
    char loc_str[YAC_MAX_LOC_STR_LEN];
    HANDLE_ERROR(nc_get_var_text(ncid, var_id, loc_str));
    tgt_location = yac_str2loc(loc_str);
  }

  if (contains_fixed_dst) {

    yac_nc_inq_varid(ncid, "fixed_values", &var_id);
    HANDLE_ERROR(nc_get_var_double(ncid, var_id, fixed_values));

    yac_nc_inq_varid(ncid, "num_dst_per_fixed_value", &var_id);
    HANDLE_ERROR(nc_get_var_int(ncid, var_id, num_tgt_per_fixed_value));

    yac_nc_inq_varid(ncid, "dst_address_fixed", &var_id);
    HANDLE_ERROR(nc_get_var_int(ncid, var_id, tgt_address_fixed));
  }

  // check data
  if ((strlen(YAC_WEIGHT_FILE_VERSION_STRING) != str_version_len) ||
      strncmp(YAC_WEIGHT_FILE_VERSION_STRING, str_version, str_version_len))
    PUT_ERR("wrong version string\n");
  if ((strlen(ref_src_grid_name) != src_grid_name_len) ||
      strncmp(ref_src_grid_name, src_grid_name, src_grid_name_len))
    PUT_ERR("wrong src_grid_name\n");
  if ((strlen(ref_tgt_grid_name) != tgt_grid_name_len) ||
      strncmp(ref_tgt_grid_name, tgt_grid_name, tgt_grid_name_len))
    PUT_ERR("wrong tgt_grid_name\n");

  if (contains_links != (ref_num_links > 0)) {
    PUT_ERR("file contains links, but reference data does not\n");
  } else if (contains_links) {
    for (unsigned i = 0; i < MIN(num_links, ref_num_links); ++i)
      if (ref_src_address[i] != src_address[i] - 1)
        PUT_ERR("wrong src_address\n");
    if ((strlen(ref_src_grid_name) != src_grid_name_len) ||
        strncmp(ref_src_grid_name, src_grid_name, src_grid_name_len))
      PUT_ERR("wrong src_grid_name\n");
    for (unsigned i = 0; i < MIN(num_links, ref_num_links); ++i)
      if (ref_tgt_address[i] != tgt_address[i] - 1)
        PUT_ERR("wrong tgt_address\n");
    if ((strlen(ref_tgt_grid_name) != tgt_grid_name_len) ||
        strncmp(ref_tgt_grid_name, tgt_grid_name, tgt_grid_name_len))
      PUT_ERR("wrong src_grid_name\n");
    for (unsigned i = 0; i < MIN(num_links, ref_num_links); ++i)
      if (fabs(ref_weights[i] - weights[i]) > 1e-10)
        PUT_ERR("wrong weight\n");
    for (unsigned i = 0; i < MIN(num_src_fields, ref_num_src_fields); ++i)
      if (ref_num_links_per_src_field[i] != num_links_per_src_field[i])
        PUT_ERR("wrong number of links per source field\n");
  }
  for (unsigned i = 0; i < MIN(num_src_fields, ref_num_src_fields); ++i)
    if (ref_src_locations[i] != src_locations[i])
      PUT_ERR("wrong source location\n");
  if (ref_tgt_location != tgt_location) PUT_ERR("wrong target location\n");

  if (contains_fixed_dst != (ref_num_fixed_values > 0)) {
    PUT_ERR("file contains fixed target data, but reference data does not\n");
  } else if (contains_fixed_dst) {
    if (ref_num_fixed_values != num_fixed_values) {
      PUT_ERR("wrong number of fixed values\n");
    } else {
      unsigned * used_flag = xcalloc(num_fixed_values, sizeof(*used_flag));
      unsigned match_count = 0;
      for (unsigned i = 0, ref_offset = 0; i < ref_num_fixed_values;
           ref_offset += ref_num_tgt_per_fixed_value[i++]) {
        for (unsigned j = 0, offset = 0; j < num_fixed_values;
             offset += num_tgt_per_fixed_value[j++]) {
          if ((!used_flag[j]) &&
              !(fabs(ref_fixed_values[i] - fixed_values[j]) > 0.0) &&
              (ref_num_tgt_per_fixed_value[i] == num_tgt_per_fixed_value[j])) {

            for (int k = 0; k < num_tgt_per_fixed_value[j]; ++k) {
              if (ref_tgt_address_fixed[ref_offset + k] !=
                  tgt_address_fixed[offset + k] - 1)
                PUT_ERR("wrong fixed target address\n");
            }

            used_flag[j] = 1;
            match_count++;
            break;
          }
        }
      }
      if (match_count != ref_num_fixed_values)
        PUT_ERR("wrong fixed values data");
      free(used_flag);
    }
  }

  // clean up
  free(str_version);
  free(fixed_values);
  free(num_tgt_per_fixed_value);
  free(tgt_address_fixed);
  free(src_grid_name);
  free(tgt_grid_name);
  free(src_address);
  free(tgt_address);
  free(weights);
  free(num_links_per_src_field);
  free(src_locations);
  free(str_fixed);

  // close file
  HANDLE_ERROR(nc_close(ncid));
}
