// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef YAC_NETCDF_ENABLED
#include <netcdf.h>
#endif

#include <string.h>

#include "geometry.h"
#include "utils_common.h"
#include "io_utils.h"

#ifdef YAC_NETCDF_ENABLED
static size_t check_dimension(int ncid, int varids[2]) {

  for (int i = 0; i < 2; ++i) {
    int ndims;
    YAC_HANDLE_ERROR(nc_inq_varndims(ncid, varids[i], &ndims));
    YAC_ASSERT(
      ndims == 1,
      "ERROR(check_dimension): coordinate array has more than one dimension")
  }

  int dimids[2];
  for (int i = 0; i < 2; ++i)
    YAC_HANDLE_ERROR(nc_inq_vardimid(ncid, varids[i], &(dimids[i])));

  YAC_ASSERT(
    dimids[0] == dimids[1],
    "ERROR(check_dimension): "
    "lon lat coordinate arrays have differing dimensions")

  size_t dimlen;
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimids[0], &dimlen));
  return dimlen;
}
#endif

int yac_check_coord_units(int ncid, int varid) {

#ifndef YAC_NETCDF_ENABLED

  UNUSED(ncid);
  UNUSED(varid);
  die("ERROR(yac_check_coord_units): YAC is built without the NetCDF support");

  return -1;
#else

  int is_degree = 0;
  nc_type att_type;
  size_t att_len;
  int status = nc_inq_att(ncid, varid, "units", &att_type, &att_len);
  // if the status is not "attribute not found"
  if (status != NC_ENOTATT) {
    YAC_HANDLE_ERROR(status);
    // if the attribute is not a string or too long
    YAC_ASSERT(
      (att_type == NC_CHAR) && (att_len <= 8),
      "ERROR(yac_check_coord_units): invalid units type or len")
    char units[8];
    memset(units, 0, 8 * sizeof(units[0]));
    YAC_HANDLE_ERROR(nc_get_att_text(ncid, varid, "units", units));
    is_degree = !strcmp(units, "degree");
    YAC_ASSERT(
      is_degree || !strcmp(units, "radian"),
      "ERROR(yac_check_coord_units): unsupported units type")
  }
  return is_degree;
#endif
}

#ifdef YAC_NETCDF_ENABLED
static double * read_coord(int ncid, int varid, size_t varlen) {

  int is_degree = yac_check_coord_units(ncid, varid);

  double * coord = xmalloc(varlen * sizeof(*coord));
  YAC_HANDLE_ERROR(nc_get_var_double (ncid, varid, coord));

  // convert to radiant if necessary
  if (is_degree) for (size_t i = 0; i < varlen; ++i) coord[i] *= YAC_RAD;

  return coord;
}
#endif

void yac_read_coords(
  int ncid, char const * lon_name, char const * lat_name,
  double ** lon, double ** lat, size_t * len) {

#ifndef YAC_NETCDF_ENABLED

  UNUSED(ncid);
  UNUSED(lon_name);
  UNUSED(lat_name);
  UNUSED(lon);
  UNUSED(lat);
  UNUSED(len);
  die("ERROR(yac_read_coords): YAC is built without the NetCDF support");
#else

  int vlonid, vlatid;
  yac_nc_inq_varid(ncid, lon_name, &vlonid);
  yac_nc_inq_varid(ncid, lat_name, &vlatid);

  size_t varlen = (*len = check_dimension(ncid, (int[]){vlonid, vlatid}));
  *lon = read_coord(ncid, vlonid, varlen);
  *lat = read_coord(ncid, vlatid, varlen);
#endif
}
