// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>
#include <math.h>
#include <netcdf.h>

#include "tests.h"
#include "grid_file_common.h"
#include "io_utils.h"
#include "utils_common.h"

static int compare_cell_coords(double * a, double * b, int count) {

  for (int order = -1; order <= 1; order += 2) {
    for (int start = 0; start < count; ++start) {

      int differences = 0;

      for (int i = 0; i < count; ++i) {

        int j = (count + start + order * i) % count;

        if (fabs(a[i] - b[j]) > 1e-6) ++differences;
      }

      if (!differences) return 0;
    }
  }
  return 1;
}

static int compare_ints(int * a, int * b, int count) {

  for (int order = -1; order <= 1; order += 2) {
    for (int start = 0; start < count; ++start) {

      int differences = 0;

      for (int i = 0; i < count; ++i) {

        int j = (count + start + order * i) % count;

        if (a[i] != b[j]) ++differences;
      }

      if (!differences) return 0;
    }
  }
  return 1;
}

void check_grid_file(
  char const * filename, char const * grid_name,
  size_t ref_num_cells, size_t ref_num_corners_per_cell,
  double * ref_cla, double * ref_clo, double * ref_lat, double * ref_lon,
  int * ref_cell_global_ids,int * ref_core_cell_mask,
  int * ref_vertex_global_ids,int * ref_core_vertex_mask,
  int * ref_edge_global_ids,int * ref_core_edge_mask) {

  if (!yac_file_exists(filename)) {
    PUT_ERR("error file does not exist");
    return;
  }

  size_t grid_name_len = strlen(grid_name) + 1;

  char nv_dim_name[3 + grid_name_len];
  char nc_dim_name[3 + grid_name_len];

  char cla_var_name[4 + grid_name_len];
  char clo_var_name[4 + grid_name_len];
  char lat_var_name[4 + grid_name_len];
  char lon_var_name[4 + grid_name_len];
  char gid_var_name[4 + grid_name_len];
  char cmk_var_name[4 + grid_name_len];
  char vgid_var_name[5 + grid_name_len];
  char vcmk_var_name[5 + grid_name_len];
  char egid_var_name[5 + grid_name_len];
  char ecmk_var_name[5 + grid_name_len];

  snprintf(nv_dim_name, 3 + grid_name_len, "nv_%s", grid_name);
  snprintf(nc_dim_name, 3 + grid_name_len, "nc_%s", grid_name);

  snprintf(cla_var_name, 4 + grid_name_len, "%s.cla", grid_name);
  snprintf(clo_var_name, 4 + grid_name_len, "%s.clo", grid_name);
  snprintf(lat_var_name, 4 + grid_name_len, "%s.lat", grid_name);
  snprintf(lon_var_name, 4 + grid_name_len, "%s.lon", grid_name);
  snprintf(gid_var_name, 4 + grid_name_len, "%s.gid", grid_name);
  snprintf(cmk_var_name, 4 + grid_name_len, "%s.cmk", grid_name);
  snprintf(vgid_var_name, 5 + grid_name_len, "%s.vgid", grid_name);
  snprintf(vcmk_var_name, 5 + grid_name_len, "%s.vcmk", grid_name);
  snprintf(egid_var_name, 5 + grid_name_len, "%s.egid", grid_name);
  snprintf(ecmk_var_name, 5 + grid_name_len, "%s.ecmk", grid_name);

  int ncid;
  yac_nc_open(filename, NC_NOWRITE, &ncid);

  int nv_dim_id, nc_dim_id;
  yac_nc_inq_dimid(ncid, nv_dim_name, &nv_dim_id);
  yac_nc_inq_dimid(ncid, nc_dim_name, &nc_dim_id);

  int cla_var_id, clo_var_id, lat_var_id, lon_var_id,
      gid_var_id, cmk_var_id, vgid_var_id, vcmk_var_id,
      egid_var_id, ecmk_var_id;
  yac_nc_inq_varid(ncid, cla_var_name, &cla_var_id);
  yac_nc_inq_varid(ncid, clo_var_name, &clo_var_id);
  if (ref_lat)
    yac_nc_inq_varid(ncid, lat_var_name, &lat_var_id);
  if (ref_lon)
    yac_nc_inq_varid(ncid, lon_var_name, &lon_var_id);
  if (ref_cell_global_ids)
    yac_nc_inq_varid(ncid, gid_var_name, &gid_var_id);
  if (ref_core_cell_mask)
    yac_nc_inq_varid(ncid, cmk_var_name, &cmk_var_id);
  if (ref_vertex_global_ids)
    yac_nc_inq_varid(ncid, vgid_var_name, &vgid_var_id);
  if (ref_core_vertex_mask)
    yac_nc_inq_varid(ncid, vcmk_var_name, &vcmk_var_id);
  if (ref_edge_global_ids)
    yac_nc_inq_varid(ncid, egid_var_name, &egid_var_id);
  if (ref_core_edge_mask)
    yac_nc_inq_varid(ncid, ecmk_var_name, &ecmk_var_id);

  size_t nv_dim_len, nc_dim_len;
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, nv_dim_id, &nv_dim_len));
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, nc_dim_id, &nc_dim_len));

  if (nv_dim_len != ref_num_corners_per_cell)
    PUT_ERR("wrong nv dim size");
  if (nc_dim_len < ref_num_cells) PUT_ERR("wrong nc dim size");

  nc_type type;
  int ndims;
  int dimids[NC_MAX_VAR_DIMS];

  YAC_HANDLE_ERROR(
    nc_inq_var(ncid, cla_var_id, NULL, &type, &ndims, dimids, NULL));
  if (type != NC_DOUBLE) PUT_ERR("wrong type for cla");
  if (ndims != 2) PUT_ERR("wrong ndims for cla");
  if ((dimids[0] != nc_dim_id) ||
      (dimids[1] != nv_dim_id)) PUT_ERR("wrong dimensions for cla");

  YAC_HANDLE_ERROR(
    nc_inq_var(ncid, clo_var_id, NULL, &type, &ndims, dimids, NULL));
  if (type != NC_DOUBLE) PUT_ERR("wrong type for clo");
  if (ndims != 2) PUT_ERR("wrong ndims for clo");
  if ((dimids[0] != nc_dim_id) ||
      (dimids[1] != nv_dim_id)) PUT_ERR("wrong dimensions for clo");

  if (ref_lat) {
    YAC_HANDLE_ERROR(
      nc_inq_var(ncid, lat_var_id, NULL, &type, &ndims, dimids, NULL));
    if (type != NC_DOUBLE) PUT_ERR("wrong type for lat");
    if (ndims != 1) PUT_ERR("wrong ndims for lat");
    if (dimids[0] != nc_dim_id) PUT_ERR("wrong dimensions for lat");
  }

  if (ref_lon) {
    YAC_HANDLE_ERROR(
      nc_inq_var(ncid, lon_var_id, NULL, &type, &ndims, dimids, NULL));
    if (type != NC_DOUBLE) PUT_ERR("wrong type for lon");
    if (ndims != 1) PUT_ERR("wrong ndims for lon");
    if (dimids[0] != nc_dim_id) PUT_ERR("wrong dimensions for lon");
  }

  if (ref_cell_global_ids) {
    YAC_HANDLE_ERROR(
      nc_inq_var(ncid, gid_var_id, NULL, &type, &ndims, dimids, NULL));
    if (type != NC_INT) PUT_ERR("wrong type for gid");
    if (ndims != 1) PUT_ERR("wrong ndims for gid");
    if (dimids[0] != nc_dim_id) PUT_ERR("wrong dimensions for gid");
  }

  if (ref_core_cell_mask) {
    YAC_HANDLE_ERROR(
      nc_inq_var(ncid, cmk_var_id, NULL, &type, &ndims, dimids, NULL));
    if (type != NC_INT) PUT_ERR("wrong type for cmk");
    if (ndims != 1) PUT_ERR("wrong ndims for cmk");
    if (dimids[0] != nc_dim_id) PUT_ERR("wrong dimensions for cmk");
  }

  if (ref_vertex_global_ids) {
    YAC_HANDLE_ERROR(
      nc_inq_var(ncid, vgid_var_id, NULL, &type, &ndims, dimids, NULL));
    if (type != NC_INT) PUT_ERR("wrong type for vgid");
  if (ndims != 2) PUT_ERR("wrong ndims for vgid");
  if ((dimids[0] != nc_dim_id) ||
      (dimids[1] != nv_dim_id)) PUT_ERR("wrong dimensions for vgid");
  }

  if (ref_core_vertex_mask) {
    YAC_HANDLE_ERROR(
      nc_inq_var(ncid, vcmk_var_id, NULL, &type, &ndims, dimids, NULL));
    if (type != NC_INT) PUT_ERR("wrong type for vcmk");
  if (ndims != 2) PUT_ERR("wrong ndims for vcmk");
  if ((dimids[0] != nc_dim_id) ||
      (dimids[1] != nv_dim_id)) PUT_ERR("wrong dimensions for vcmk");
  }

  if (ref_edge_global_ids) {
    YAC_HANDLE_ERROR(
      nc_inq_var(ncid, egid_var_id, NULL, &type, &ndims, dimids, NULL));
    if (type != NC_INT) PUT_ERR("wrong type for egid");
  if (ndims != 2) PUT_ERR("wrong ndims for egid");
  if ((dimids[0] != nc_dim_id) ||
      (dimids[1] != nv_dim_id)) PUT_ERR("wrong dimensions for egid");
  }

  if (ref_core_edge_mask) {
    YAC_HANDLE_ERROR(
      nc_inq_var(ncid, ecmk_var_id, NULL, &type, &ndims, dimids, NULL));
    if (type != NC_INT) PUT_ERR("wrong type for ecmk");
  if (ndims != 2) PUT_ERR("wrong ndims for ecmk");
  if ((dimids[0] != nc_dim_id) ||
      (dimids[1] != nv_dim_id)) PUT_ERR("wrong dimensions for ecmk");
  }

  double * clo = malloc(nc_dim_len * nv_dim_len * sizeof(*clo));
  double * cla = malloc(nc_dim_len * nv_dim_len * sizeof(*cla));
  double * lon = ref_lon?malloc(nc_dim_len * sizeof(*lon)):NULL;
  double * lat = ref_lat?malloc(nc_dim_len * sizeof(*lat)):NULL;
  int * gid = ref_cell_global_ids?malloc(nc_dim_len * sizeof(*gid)):NULL;
  int * cmk = ref_core_cell_mask?malloc(nc_dim_len * sizeof(*cmk)):NULL;
  int * vgid = ref_vertex_global_ids?malloc(nc_dim_len * nv_dim_len * sizeof(*vgid)):NULL;
  int * vcmk = ref_core_vertex_mask?malloc(nc_dim_len * nv_dim_len * sizeof(*vcmk)):NULL;
  int * egid = ref_edge_global_ids?malloc(nc_dim_len * nv_dim_len * sizeof(*egid)):NULL;
  int * ecmk = ref_core_edge_mask?malloc(nc_dim_len * nv_dim_len * sizeof(*ecmk)):NULL;

  YAC_HANDLE_ERROR(nc_get_var_double(ncid, cla_var_id, cla));
  YAC_HANDLE_ERROR(nc_get_var_double(ncid, clo_var_id, clo));
  if (ref_lat)
    YAC_HANDLE_ERROR(nc_get_var_double(ncid, lat_var_id, lat));
  if (ref_lon)
    YAC_HANDLE_ERROR(nc_get_var_double(ncid, lon_var_id, lon));
  if (ref_cell_global_ids)
    YAC_HANDLE_ERROR(nc_get_var_int(ncid, gid_var_id, gid));
  if (ref_core_cell_mask)
    YAC_HANDLE_ERROR(nc_get_var_int(ncid, cmk_var_id, cmk));
  if (ref_vertex_global_ids)
    YAC_HANDLE_ERROR(nc_get_var_int(ncid, vgid_var_id, vgid));
  if (ref_core_vertex_mask)
    YAC_HANDLE_ERROR(nc_get_var_int(ncid, vcmk_var_id, vcmk));
  if (ref_edge_global_ids)
    YAC_HANDLE_ERROR(nc_get_var_int(ncid, egid_var_id, egid));
  if (ref_core_edge_mask)
    YAC_HANDLE_ERROR(nc_get_var_int(ncid, ecmk_var_id, ecmk));

  int * ref_cell_flag = calloc(ref_num_cells, sizeof(*ref_cell_flag));

  // for all cells in the grid file
  for (size_t i = 0; i < nc_dim_len; ++i) {

    // get the coordinates of the current cell
    double * curr_cla = cla + i * nv_dim_len;
    double * curr_clo = clo + i * nv_dim_len;

    // match the current cell with the reference cells
    size_t match_idx = SIZE_MAX;
    {
      if (ref_cell_global_ids) {
        for (size_t j = 0; (j < ref_num_cells) && (match_idx == SIZE_MAX); ++j)
          if (gid[i] == ref_cell_global_ids[j])
            match_idx = j;
      } else if (ref_lat && ref_lon) {
        for (size_t j = 0; (j < ref_num_cells) && (match_idx == SIZE_MAX); ++j)
          if ((fabs(lat[i] - ref_lat[j]) < 1e-3) &&
              (fabs(lon[i] - ref_lon[j]) < 1e-3))
            match_idx = j;
      } else {
        for (size_t j = 0; (j < ref_num_cells) && (match_idx == SIZE_MAX); ++j)
          if (!compare_cell_coords(
                curr_cla, ref_cla + j * ref_num_corners_per_cell,
                ref_num_corners_per_cell) &&
              !compare_cell_coords(
                curr_clo, ref_clo + j * ref_num_corners_per_cell,
                ref_num_corners_per_cell))
            match_idx = j;
      }
    }

    if (match_idx == SIZE_MAX) {

      PUT_ERR("error no matching cell");

    } else {

      ref_cell_flag[match_idx] = 1;

      if (compare_cell_coords(
            curr_cla, ref_cla + match_idx * ref_num_corners_per_cell,
            ref_num_corners_per_cell))
        PUT_ERR("wrong cla");
      if (compare_cell_coords(
            curr_clo, ref_clo + match_idx * ref_num_corners_per_cell,
            ref_num_corners_per_cell))
        PUT_ERR("wrong clo");
      if (ref_lat)
        if (fabs(lat[i] - ref_lat[match_idx]) > 1e-3) PUT_ERR("wrong lat");
      if (ref_lon)
        if (fabs(lon[i] - ref_lon[match_idx]) > 1e-3) PUT_ERR("wrong lon");
      if (ref_cell_global_ids)
        if (gid[i] != (int)(ref_cell_global_ids[match_idx])) PUT_ERR("wrong gid");
      if (ref_core_cell_mask)
        if (cmk[i] != (ref_core_cell_mask[match_idx]))
          PUT_ERR("wrong cmk");
      if (ref_vertex_global_ids)
        if (compare_ints(
              vgid + i * nv_dim_len,
              ref_vertex_global_ids + match_idx * ref_num_corners_per_cell,
              ref_num_corners_per_cell))
          PUT_ERR("wrong vgid");
      if (ref_core_vertex_mask)
        if (compare_ints(
              vcmk + i * nv_dim_len,
              ref_core_vertex_mask + match_idx * ref_num_corners_per_cell,
              ref_num_corners_per_cell))
          PUT_ERR("wrong vcmk");
      if (ref_edge_global_ids)
        if (compare_ints(
              egid + i * nv_dim_len,
              ref_edge_global_ids + match_idx * ref_num_corners_per_cell,
              ref_num_corners_per_cell))
          PUT_ERR("wrong egid");
      if (ref_core_edge_mask)
        if (compare_ints(
              ecmk + i * nv_dim_len,
              ref_core_edge_mask + match_idx * ref_num_corners_per_cell,
              ref_num_corners_per_cell))
          PUT_ERR("wrong ecmk");
    }
  }

  size_t match_count = 0;
  for (size_t i = 0; i < ref_num_cells; ++i)
    if (ref_cell_flag[i]) match_count++;

  if (match_count != ref_num_cells) PUT_ERR("missing cells");

  free(ref_cell_flag);
  free(ecmk);
  free(egid);
  free(vcmk);
  free(vgid);
  free(cmk);
  free(gid);
  free(lat);
  free(lon);
  free(clo);
  free(cla);

  YAC_HANDLE_ERROR(nc_close(ncid));
}
